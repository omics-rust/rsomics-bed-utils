//! Compute the disjoint-interval union of N sorted BED files and annotate each
//! sub-interval with which files cover it.
//!
//! Matches `bedtools multiinter`: output columns are
//!   chrom  start  end  count  list  [0|1 per file …]
//! where `list` is a comma-separated list of 1-based file indices (or names
//! when `-names` is supplied) and the per-file indicator columns follow.
//! Input files must be sorted (chrom, then start) — no sorting is performed.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Debug)]
struct Interval {
    start: u64,
    end: u64,
}

/// Load one BED file into a per-chrom interval list and merge overlapping
/// intervals within each chrom into their union. This ensures a single file
/// does not contribute internal boundary points to the sweep-line endpoint
/// set, which would cause spurious splits not produced by `bedtools multiinter`.
fn load_file(path: &Path) -> Result<BTreeMap<String, Vec<Interval>>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut by_chrom: BTreeMap<String, Vec<Interval>> = BTreeMap::new();

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let mut fields = line.splitn(4, '\t');
        let chrom = match fields.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        let start: u64 = match fields.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let end: u64 = match fields.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        by_chrom
            .entry(chrom)
            .or_default()
            .push(Interval { start, end });
    }

    // Merge overlapping/adjacent intervals within each chrom so that
    // intra-file boundaries do not produce extra sweep-line endpoints.
    for ivs in by_chrom.values_mut() {
        ivs.sort_unstable_by_key(|iv| (iv.start, iv.end));
        let mut merged: Vec<Interval> = Vec::with_capacity(ivs.len());
        for iv in ivs.drain(..) {
            if let Some(last) = merged.last_mut()
                && iv.start <= last.end
            {
                last.end = last.end.max(iv.end);
                continue;
            }
            merged.push(iv);
        }
        *ivs = merged;
    }

    Ok(by_chrom)
}

/// Event kind: +1 when a file's merged interval opens, -1 when it closes.
#[derive(Clone, Copy)]
struct Event {
    pos: u64,
    file_idx: usize,
    /// +1 = interval starts, -1 = interval ends.
    delta: i8,
}

/// O(E log E) event sweep-line over a single chromosome across all files.
///
/// For each file we push two events per merged interval (start/end). We sort by
/// position, then walk events in groups of equal position. Between each pair of
/// consecutive distinct positions we maintain per-file active counts; when the
/// active *set* is non-empty we emit one output segment. Counts are updated
/// incrementally — no per-segment interval rescan.
fn sweep_chrom(
    chrom: &str,
    file_intervals: &[&[Interval]],
    names: &[String],
    out: &mut impl Write,
) -> Result<()> {
    let n = file_intervals.len();

    // Build the event list: one open + one close per merged interval per file.
    let total_events: usize = file_intervals.iter().map(|ivs| ivs.len() * 2).sum();
    let mut events: Vec<Event> = Vec::with_capacity(total_events);
    for (fi, ivs) in file_intervals.iter().enumerate() {
        for iv in *ivs {
            events.push(Event {
                pos: iv.start,
                file_idx: fi,
                delta: 1,
            });
            events.push(Event {
                pos: iv.end,
                file_idx: fi,
                delta: -1,
            });
        }
    }
    if events.is_empty() {
        return Ok(());
    }

    // Sort by position; tie-break: closes (-1) before opens (+1) so that a
    // close at position P and an open at P are handled correctly (the gap
    // [prev..P] is emitted with the old set, then P starts the new set).
    events.sort_unstable_by_key(|e| (e.pos, e.delta));

    // Per-file active depth; a file is "active" when its depth > 0.
    let mut depth: Vec<u32> = vec![0; n];
    // How many files are currently active.
    let mut n_active: usize = 0;

    let mut prev_pos: u64 = 0;
    let mut started = false;

    // Reusable scratch for building the output line.
    let mut list_buf = String::new();
    let mut indicators: Vec<u8> = vec![b'0'; n];

    let mut i = 0;
    while i < events.len() {
        let cur_pos = events[i].pos;

        // Emit the segment [prev_pos, cur_pos) if any files are active.
        if started && n_active > 0 && prev_pos < cur_pos {
            list_buf.clear();
            let mut first = true;
            let mut count = 0usize;
            for fi in 0..n {
                if depth[fi] > 0 {
                    if !first {
                        list_buf.push(',');
                    }
                    list_buf.push_str(names[fi].as_str());
                    indicators[fi] = b'1';
                    first = false;
                    count += 1;
                } else {
                    indicators[fi] = b'0';
                }
            }
            write!(out, "{chrom}\t{prev_pos}\t{cur_pos}\t{count}\t{list_buf}")
                .map_err(RsomicsError::Io)?;
            for b in &indicators {
                write!(out, "\t{}", *b as char).map_err(RsomicsError::Io)?;
            }
            writeln!(out).map_err(RsomicsError::Io)?;
        }

        // Process all events at cur_pos.
        while i < events.len() && events[i].pos == cur_pos {
            let e = events[i];
            let was_active = depth[e.file_idx] > 0;
            if e.delta > 0 {
                depth[e.file_idx] += 1;
                if !was_active {
                    n_active += 1;
                }
            } else {
                depth[e.file_idx] -= 1;
                if was_active && depth[e.file_idx] == 0 {
                    n_active -= 1;
                }
            }
            i += 1;
        }

        prev_pos = cur_pos;
        started = true;
    }

    Ok(())
}

pub fn multiinter(
    inputs: &[&Path],
    names: Option<&[String]>,
    output: &mut dyn Write,
) -> Result<()> {
    let n = inputs.len();

    // Build label list: either user-supplied names or 1-based indices.
    let default_names: Vec<String> = (1..=n).map(|i| i.to_string()).collect();
    let labels: &[String] = names.unwrap_or(&default_names);

    if labels.len() != n {
        return Err(RsomicsError::InvalidInput(
            "-names must supply one name per input file".to_string(),
        ));
    }

    // Load all files.
    let loaded: Vec<BTreeMap<String, Vec<Interval>>> =
        inputs.iter().map(|p| load_file(p)).collect::<Result<_>>()?;

    // Collect all chroms across all files, preserving first-seen order via a
    // BTreeMap (sorted lexicographically — matches bedtools).
    let mut all_chroms: BTreeMap<String, ()> = BTreeMap::new();
    for m in &loaded {
        for k in m.keys() {
            all_chroms.insert(k.clone(), ());
        }
    }

    let mut out = BufWriter::with_capacity(64 * 1024, output);

    for chrom in all_chroms.keys() {
        let per_file: Vec<&[Interval]> = loaded
            .iter()
            .map(|m| m.get(chrom).map(|v| v.as_slice()).unwrap_or(&[]))
            .collect();

        sweep_chrom(chrom, &per_file, labels, &mut out)?;
    }

    out.flush().map_err(RsomicsError::Io)
}
