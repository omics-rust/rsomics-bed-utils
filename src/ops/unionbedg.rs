//! Combine N bedGraph files into a matrix of disjoint maximal segments.
//!
//! Matches `bedtools unionbedg -i f1 f2 ... [-names n1 n2 ...] [-header]`.
//!
//! Algorithm — O(E log E) event-sweep (same shape as multiinter):
//! 1. For each file i, parse all (chrom, start, end, value) records. The value
//!    is the raw string from column 4 (bedtools keeps it as-is); default = "0".
//! 2. Build one open-event and one close-event per record per file. Sort all
//!    events by (chrom, pos, delta) where closes (-1) sort before opens (+1)
//!    at the same position — same tie-break as bedtools so that a close at P
//!    and an open at P emit the [prev..P] segment with the old values, then
//!    reset.
//! 3. Walk the events. Between consecutive distinct positions, if any file is
//!    active, emit one line: chrom start end val[0] val[1] … where inactive
//!    files emit "0".
//! 4. Each file's "current value" is the value of the most recently opened
//!    record. Because bedGraph is non-overlapping within a file, at most one
//!    record per file is active at any moment.
//! 5. Chromosomes are processed in the order they first appear in the first input
//!    file, with any additional chroms (from subsequent files) appended in their
//!    first-seen order — matching bedtools' output order.

use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

type BgMap = BTreeMap<String, Vec<Record>>;

struct Record {
    start: u64,
    end: u64,
    value: String,
}

/// Returns (chrom_order, by_chrom) where chrom_order lists chroms in first-seen order.
fn load_bedgraph(path: &Path) -> Result<(Vec<String>, BgMap)> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut by_chrom: BgMap = BTreeMap::new();
    let mut chrom_order: Vec<String> = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let mut f = line.splitn(5, '\t');
        let chrom = match f.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        let start: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let end: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let value = f.next().unwrap_or("0").trim_end().to_string();
        if !by_chrom.contains_key(&chrom) {
            chrom_order.push(chrom.clone());
        }
        by_chrom
            .entry(chrom)
            .or_default()
            .push(Record { start, end, value });
    }
    Ok((chrom_order, by_chrom))
}

// ---------------------------------------------------------------------------
// Event sweep
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Event<'a> {
    pos: u64,
    /// +1 = interval opens, -1 = interval closes.
    delta: i8,
    file_idx: usize,
    /// Value string (only meaningful for open events; ignored on close).
    value: &'a str,
}

fn sweep_chrom(
    chrom: &str,
    file_records: &[&[Record]],
    n_files: usize,
    out: &mut impl Write,
) -> Result<()> {
    // Build event list.
    let total = file_records.iter().map(|rs| rs.len() * 2).sum();
    let mut events: Vec<Event<'_>> = Vec::with_capacity(total);
    for (fi, records) in file_records.iter().enumerate() {
        for rec in *records {
            events.push(Event {
                pos: rec.start,
                delta: 1,
                file_idx: fi,
                value: &rec.value,
            });
            events.push(Event {
                pos: rec.end,
                delta: -1,
                file_idx: fi,
                value: "",
            });
        }
    }
    if events.is_empty() {
        return Ok(());
    }

    // Sort: by position, then closes (-1) before opens (+1).
    events.sort_unstable_by_key(|e| (e.pos, e.delta));

    // Per-file current value ("0" = not active).
    let mut cur_values: Vec<&str> = vec!["0"; n_files];
    let mut active: Vec<bool> = vec![false; n_files];
    let mut n_active: usize = 0;

    let mut prev_pos: u64 = 0;
    let mut started = false;

    let mut i = 0;
    while i < events.len() {
        let cur_pos = events[i].pos;

        // Emit segment [prev_pos, cur_pos) if any file is active.
        if started && n_active > 0 && prev_pos < cur_pos {
            write!(out, "{chrom}\t{prev_pos}\t{cur_pos}").map_err(RsomicsError::Io)?;
            for fi in 0..n_files {
                let v = if active[fi] { cur_values[fi] } else { "0" };
                write!(out, "\t{v}").map_err(RsomicsError::Io)?;
            }
            writeln!(out).map_err(RsomicsError::Io)?;
        }

        // Process all events at cur_pos.
        while i < events.len() && events[i].pos == cur_pos {
            let e = events[i];
            if e.delta > 0 {
                // Open: update value and mark active.
                cur_values[e.file_idx] = e.value;
                if !active[e.file_idx] {
                    active[e.file_idx] = true;
                    n_active += 1;
                }
            } else {
                // Close.
                if active[e.file_idx] {
                    active[e.file_idx] = false;
                    cur_values[e.file_idx] = "0";
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

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn unionbedg(
    inputs: &[&Path],
    names: Option<&[String]>,
    header: bool,
    out: &mut dyn Write,
) -> Result<()> {
    let n = inputs.len();

    if let Some(ns) = names
        && ns.len() != n
    {
        return Err(RsomicsError::InvalidInput(
            "-names must supply one name per input file".to_string(),
        ));
    }

    // Load all files; each returns (chrom_order, by_chrom).
    let loaded: Vec<(Vec<String>, BgMap)> = inputs
        .iter()
        .map(|p| load_bedgraph(p))
        .collect::<Result<_>>()?;

    // Collect chroms in first-seen order across all files, preserving per-file order.
    let mut seen: HashSet<String> = HashSet::new();
    let mut all_chroms: Vec<String> = Vec::new();
    for (order, _) in &loaded {
        for chrom in order {
            if seen.insert(chrom.clone()) {
                all_chroms.push(chrom.clone());
            }
        }
    }

    let mut writer = BufWriter::with_capacity(64 * 1024, out);

    // Optional header line.
    if header {
        write!(writer, "chrom\tstart\tend").map_err(RsomicsError::Io)?;
        if let Some(ns) = names {
            for name in ns {
                write!(writer, "\t{name}").map_err(RsomicsError::Io)?;
            }
        } else {
            for i in 1..=n {
                write!(writer, "\t{i}").map_err(RsomicsError::Io)?;
            }
        }
        writeln!(writer).map_err(RsomicsError::Io)?;
    }

    for chrom in &all_chroms {
        let per_file: Vec<&[Record]> = loaded
            .iter()
            .map(|(_, m)| m.get(chrom).map(|v| v.as_slice()).unwrap_or(&[]))
            .collect();
        sweep_chrom(chrom, &per_file, n, &mut writer)?;
    }

    writer.flush().map_err(RsomicsError::Io)
}
