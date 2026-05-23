use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone)]
struct Interval {
    chrom: String,
    start: u64,
    end: u64,
    line: String,
}

fn load_sorted_bed(path: &Path) -> Result<Vec<Interval>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
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
        intervals.push(Interval {
            chrom,
            start,
            end,
            line: line.clone(),
        });
    }
    Ok(intervals)
}

/// Compute the non-overlapping gap between A and B (0 if adjacent, positive if
/// there is a gap).  Returns `None` if A and B strictly overlap (share bases).
fn gap(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> Option<u64> {
    if a_end <= b_start {
        Some(b_start - a_end)
    } else if b_end <= a_start {
        Some(a_start - b_end)
    } else {
        // Strict overlap (shared bases).
        None
    }
}

/// Find closest B feature(s) for each A interval, matching `bedtools closest`.
///
/// Priority rule (mirrors bedtools):
///   1. If any B strictly overlaps A, emit ALL strictly-overlapping B.
///   2. Otherwise emit all B at the minimum non-overlapping gap distance.
///
/// No distance column is appended (bedtools default omits it).
/// When no B exists on the same chromosome: emits `A\t.\t-1\t-1`.
pub fn closest(a_path: &Path, b_path: &Path, output: &mut dyn Write) -> Result<u64> {
    let a_intervals = load_sorted_bed(a_path)?;
    let b_intervals = load_sorted_bed(b_path)?;
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    // Group B by chromosome.
    let mut b_by_chrom: std::collections::HashMap<&str, Vec<&Interval>> =
        std::collections::HashMap::new();
    for b in &b_intervals {
        b_by_chrom.entry(b.chrom.as_str()).or_default().push(b);
    }

    for a in &a_intervals {
        let a_fields = a.line.as_str();
        match b_by_chrom.get(a.chrom.as_str()) {
            None => {
                writeln!(out, "{a_fields}\t.\t-1\t-1").map_err(RsomicsError::Io)?;
                count += 1;
            }
            Some(bs) => {
                // Collect all strictly-overlapping B intervals.
                let overlapping: Vec<&&Interval> = bs
                    .iter()
                    .filter(|b| gap(a.start, a.end, b.start, b.end).is_none())
                    .collect();

                if !overlapping.is_empty() {
                    // Emit all overlapping B.
                    for b in overlapping {
                        writeln!(out, "{a_fields}\t{}", b.line).map_err(RsomicsError::Io)?;
                        count += 1;
                    }
                } else {
                    // No overlap — find minimum gap distance, emit all at that distance.
                    let min_gap = bs
                        .iter()
                        .filter_map(|b| gap(a.start, a.end, b.start, b.end))
                        .min();
                    if let Some(min_d) = min_gap {
                        for b in bs.iter() {
                            if gap(a.start, a.end, b.start, b.end) == Some(min_d) {
                                writeln!(out, "{a_fields}\t{}", b.line)
                                    .map_err(RsomicsError::Io)?;
                                count += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
