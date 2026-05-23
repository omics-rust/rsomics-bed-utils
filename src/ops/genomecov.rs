#![allow(clippy::cast_precision_loss)]

//! Per-base genome coverage histogram (`bedtools genomecov` default mode).
//!
//! For each chromosome (in input order) and then for the whole genome, emits:
//!   chrom  depth  count  chrom_len  fraction
//!
//! where `count` is the number of bases at that depth.  Rows with depth > 0
//! come from a linear sweep over the sorted depth signal; depth-0 row fills the
//! remainder.  Fractions use C `%.6g` (6 significant digits, trailing zeros
//! stripped).  The genome-level summary adds depth 0/1/2/… lines pooled across
//! all chromosomes, using the total genome length as denominator.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// Format a fraction with C `%.6g` semantics (6 significant digits, trailing
/// zeros stripped, decimal point removed when not needed).
fn fmt_g6(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    if v == 1.0 {
        return "1".to_string();
    }
    let mag = v.abs().log10().floor() as i32;
    if (-4..6).contains(&mag) {
        let decimals = (5 - mag).max(0) as usize;
        let s = format!("{:.prec$}", v, prec = decimals);
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    } else {
        let raw = format!("{:.5e}", v);
        normalize_exp_2digits(raw)
    }
}

fn normalize_exp_2digits(s: String) -> String {
    if let Some(e_pos) = s.find('e') {
        let (mantissa, exp_part) = s.split_at(e_pos);
        let after_e = &exp_part[1..];
        let (sign, digits) = if after_e.starts_with(['+', '-']) {
            (&after_e[..1], &after_e[1..])
        } else {
            ("+", after_e)
        };
        if digits.len() < 2 {
            format!("{mantissa}e{sign}{digits:0>2}")
        } else {
            s
        }
    } else {
        s
    }
}

pub fn genomecov(bed_path: &Path, genome_path: &Path, output: &mut dyn Write) -> Result<()> {
    // Load genome in declaration order (order matters for output).
    let genome_order = load_genome_ordered(genome_path)?;
    let genome_map: HashMap<&str, u64> =
        genome_order.iter().map(|(c, l)| (c.as_str(), *l)).collect();

    // Load all BED intervals grouped by chrom.
    let mut ivs_by_chrom: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    let file = File::open(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let chrom = f[0].to_string();
        let start: u64 = f[1]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("start: {e}")))?;
        let end: u64 = f[2]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("end: {e}")))?;
        if start < end {
            ivs_by_chrom.entry(chrom).or_default().push((start, end));
        }
    }

    let mut out = BufWriter::with_capacity(64 * 1024, output);

    // Depth histogram pooled across the whole genome.
    let mut genome_hist: HashMap<u32, u64> = HashMap::new();
    let total_genome: u64 = genome_order.iter().map(|(_, l)| l).sum();

    for (chrom, chrom_len) in &genome_order {
        let chrom_len = *chrom_len;
        let hist = depth_histogram(
            ivs_by_chrom.get(chrom.as_str()).map_or(&[], |v| v),
            chrom_len,
            &genome_map,
            chrom,
        );

        for (depth, count) in &hist {
            let frac = *count as f64 / chrom_len as f64;
            writeln!(
                out,
                "{chrom}\t{depth}\t{count}\t{chrom_len}\t{}",
                fmt_g6(frac)
            )
            .map_err(RsomicsError::Io)?;
            *genome_hist.entry(*depth).or_insert(0) += count;
        }
    }

    // Genome-level summary.
    let mut depths: Vec<u32> = genome_hist.keys().copied().collect();
    depths.sort_unstable();
    for depth in depths {
        let count = genome_hist[&depth];
        let frac = count as f64 / total_genome as f64;
        writeln!(
            out,
            "genome\t{depth}\t{count}\t{total_genome}\t{}",
            fmt_g6(frac)
        )
        .map_err(RsomicsError::Io)?;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

/// Compute a depth histogram for one chromosome using a coordinate-compression
/// sweep.  Returns `(depth → base_count)` sorted by depth ascending.
///
/// Algorithm: sort intervals by start, then run an event sweep to compute a
/// depth array compactly (via coordinate-compressed "events" not per-base
/// allocation).  O(N log N) time, O(N) space.
fn depth_histogram(
    ivs: &[(u64, u64)],
    chrom_len: u64,
    genome_map: &HashMap<&str, u64>,
    chrom: &str,
) -> Vec<(u32, u64)> {
    // Clamp intervals to chromosome bounds.
    let chrom_len = genome_map.get(chrom).copied().unwrap_or(chrom_len);

    // Build start/end events.
    let mut events: Vec<(u64, i32)> = Vec::with_capacity(ivs.len() * 2);
    for &(s, e) in ivs {
        let s = s.min(chrom_len);
        let e = e.min(chrom_len);
        if s < e {
            events.push((s, 1));
            events.push((e, -1));
        }
    }
    // Sort by position; at ties, closes (-1) before opens (+1) so [a,b) semantics hold.
    events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    // Sweep to compute (start, depth) segments.
    let mut hist: HashMap<u32, u64> = HashMap::new();
    let mut prev_pos: u64 = 0;
    let mut depth: i32 = 0;
    let mut i = 0;
    while i < events.len() {
        let pos = events[i].0;
        // Consume all events at the same position.
        if pos > prev_pos && depth >= 0 {
            let bases = pos - prev_pos;
            *hist.entry(depth as u32).or_insert(0) += bases;
        }
        while i < events.len() && events[i].0 == pos {
            depth += events[i].1;
            i += 1;
        }
        prev_pos = pos;
    }
    // Tail: from last event to end of chromosome (depth should be 0).
    if prev_pos < chrom_len && depth >= 0 {
        let bases = chrom_len - prev_pos;
        *hist.entry(depth as u32).or_insert(0) += bases;
    }

    let mut v: Vec<(u32, u64)> = hist.into_iter().collect();
    v.sort_unstable_by_key(|(d, _)| *d);
    v
}

fn load_genome_ordered(path: &Path) -> Result<Vec<(String, u64)>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut out = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 2 {
            let len: u64 = f[1]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad length: {e}")))?;
            out.push((f[0].to_string(), len));
        }
    }
    Ok(out)
}
