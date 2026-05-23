//! Relative distance distribution between two BED files.
//!
//! Matches `bedtools reldist -a A -b B` (summary mode, the default).
//!
//! Algorithm (from bedtools RelDist.cpp, GPL — implemented here from the
//! published method and black-box behaviour):
//!
//! For each interval in A:
//!   1. Compute its midpoint: midA = (start + end) / 2 (integer division).
//!   2. Binary-search the sorted B-midpoint array for the same chromosome to
//!      find the adjacent pair [left, right] that straddles midA.
//!      - low_idx = lower_bound(midA) - 1 (clamped to 0).
//!      - high_idx = low_idx + 1.
//!      - Skip if low_idx == last index (no right neighbour).
//!      - Skip if left > midA (midA is before all B midpoints on this chrom).
//!   3. reldist = min(|midA - left|, |midA - right|) / (right - left)
//!      Special case: if the minimum distance is 0, reldist = 0.0.
//!   4. Accumulate: floor(reldist * 100) / 100 (two-decimal bin).
//!
//! Output: `reldist count total fraction` with `printf("%.2f\t%lu\t%lu\t%.3lf\n", ...)`.
//! The header line is always printed. If no A interval qualifies, only the
//! header is output (matching bedtools).

#![allow(clippy::cast_precision_loss)]

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

fn load_bed_midpoints(path: &Path) -> Result<BTreeMap<String, Vec<i64>>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut by_chrom: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut f = line.splitn(4, '\t');
        let chrom = match f.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        let start: i64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let end: i64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let mid = (start + end) / 2;
        by_chrom.entry(chrom).or_default().push(mid);
    }
    for mids in by_chrom.values_mut() {
        mids.sort_unstable();
    }
    Ok(by_chrom)
}

pub fn reldist(a_path: &Path, b_path: &Path, out: &mut dyn Write) -> Result<()> {
    // Load B midpoints sorted per chromosome.
    let db_mids = load_bed_midpoints(b_path)?;

    // Accumulate relative-distance bins: bin_key = floor(reldist * 100) / 100.
    let mut reldists: BTreeMap<u64, usize> = BTreeMap::new();
    let mut tot_queries: usize = 0;

    let file = File::open(a_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", a_path.display())))?;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut f = line.splitn(4, '\t');
        let chrom = match f.next() {
            Some(c) => c,
            None => continue,
        };
        let start: i64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let end: i64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };

        let chrom_mids = match db_mids.get(chrom) {
            Some(m) => m,
            None => continue,
        };

        let mid_a = (start + end) / 2;

        // Binary search: lower_bound gives the first index where chrom_mids[i] >= mid_a.
        let lb = chrom_mids.partition_point(|&m| m < mid_a);

        // low_idx = lb - 1 (or lb if lb==0, matching bedtools' clamp).
        let low_idx = if lb == 0 { 0 } else { lb - 1 };
        let high_idx = low_idx + 1;

        // Skip if there is no right neighbour.
        if low_idx == chrom_mids.len() - 1 {
            continue;
        }

        let left = chrom_mids[low_idx];
        let right = chrom_mids[high_idx];

        // Skip if mid_a is before the left B midpoint (bedtools skips this too).
        if left > mid_a {
            continue;
        }

        let left_dist = (mid_a - left).unsigned_abs();
        let right_dist = (mid_a - right).unsigned_abs();
        let min_dist = left_dist.min(right_dist);
        let span = (right - left) as u64;

        let rel = if min_dist == 0 || span == 0 {
            0.0f64
        } else {
            min_dist as f64 / span as f64
        };

        // Floor to two decimal places (matching bedtools: floor(rel*100)/100).
        let bin = (rel * 100.0).floor() as u64;
        *reldists.entry(bin).or_insert(0) += 1;
        tot_queries += 1;
    }

    // Output header + bins.
    writeln!(out, "reldist\tcount\ttotal\tfraction").map_err(RsomicsError::Io)?;
    for (bin_key, count) in &reldists {
        let bin_val = *bin_key as f64 / 100.0;
        let fraction = *count as f64 / tot_queries as f64;
        writeln!(out, "{bin_val:.2}\t{count}\t{tot_queries}\t{fraction:.3}")
            .map_err(RsomicsError::Io)?;
    }

    Ok(())
}
