//! Assign a cluster ID to each sorted BED interval such that overlapping or
//! bookended intervals (within `-d` bp) on the same chrom share one ID.
//!
//! Matches `bedtools cluster` semantics: IDs are globally incrementing integers
//! starting at 1; the cluster ID is appended as a trailing tab-separated column.
//! Input must be sorted (chrom, then start). No sorting is performed here.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// Append a cluster-ID column to each BED record.
///
/// `max_dist` = maximum gap between intervals that still counts as the same
/// cluster (bedtools `-d`; 0 means overlap/bookend only).
pub fn cluster(input: &Path, output: &mut dyn Write, max_dist: i64) -> Result<()> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);

    let mut prev_chrom = String::new();
    let mut cluster_end: i64 = 0;
    let mut cluster_id: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let chrom = fields[0];
        let start: i64 = fields[1]
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("bad start: {:?}", fields[1])))?;
        let end: i64 = fields[2]
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("bad end: {:?}", fields[2])))?;

        // New cluster when chrom changes or gap exceeds max_dist.
        // bedtools cluster: a gap of exactly max_dist still clusters (<=).
        if chrom != prev_chrom || start > cluster_end + max_dist {
            cluster_id += 1;
            cluster_end = end;
            prev_chrom = chrom.to_string();
        } else {
            // Extend the current cluster's reach.
            if end > cluster_end {
                cluster_end = end;
            }
        }

        writeln!(out, "{line}\t{cluster_id}").map_err(RsomicsError::Io)?;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
