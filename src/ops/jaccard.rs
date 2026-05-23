#![allow(clippy::cast_precision_loss)]

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Clone)]
struct Iv {
    chrom: String,
    start: u64,
    end: u64,
}

fn load_bed(path: &Path) -> Result<Vec<Iv>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut out: Vec<Iv> = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim_end_matches('\r');
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let mut cols = line.splitn(4, '\t');
        let chrom = cols
            .next()
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!("line {}: missing chrom", lineno + 1))
            })?
            .to_string();
        let start_s = cols.next().ok_or_else(|| {
            RsomicsError::InvalidInput(format!("line {}: missing start", lineno + 1))
        })?;
        let end_s = cols.next().ok_or_else(|| {
            RsomicsError::InvalidInput(format!("line {}: missing end", lineno + 1))
        })?;
        let start: u64 = start_s.parse().map_err(|e| {
            RsomicsError::InvalidInput(format!("line {}: bad start {start_s:?}: {e}", lineno + 1))
        })?;
        let end: u64 = end_s.parse().map_err(|e| {
            RsomicsError::InvalidInput(format!("line {}: bad end {end_s:?}: {e}", lineno + 1))
        })?;
        out.push(Iv { chrom, start, end });
    }
    out.sort_unstable_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    Ok(out)
}

/// Merge sorted intervals into non-overlapping segments per chromosome.
///
/// Adjacent or overlapping intervals on the same chromosome are fused.
fn merge_sorted(ivs: &[Iv]) -> BTreeMap<String, Vec<(u64, u64)>> {
    let mut map: BTreeMap<String, Vec<(u64, u64)>> = BTreeMap::new();
    for iv in ivs {
        let segs = map.entry(iv.chrom.clone()).or_default();
        if let Some(last) = segs.last_mut()
            && iv.start <= last.1
        {
            last.1 = last.1.max(iv.end);
            continue;
        }
        segs.push((iv.start, iv.end));
    }
    map
}

/// Format a float with C `%.6g` semantics (6 significant digits, trailing
/// zeros stripped, decimal point removed when not needed).
fn format_6g(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    if v == 1.0 {
        return "1".to_string();
    }
    let s = format!("{:.6}", v);
    let s = s.trim_end_matches('0');
    let s = s.trim_end_matches('.');
    s.to_string()
}

#[derive(Debug)]
pub struct JaccardResult {
    pub intersection: u64,
    pub union: u64,
    pub jaccard: f64,
    pub n_intersections: u64,
}

/// Compute Jaccard similarity between two BED files.
///
/// Mirrors `bedtools jaccard`: each file's intervals are merged independently
/// before computing intersection and union base-pair counts.  The ratio
/// `intersection / union` is the Jaccard index.
pub fn jaccard(a_path: &Path, b_path: &Path) -> Result<JaccardResult> {
    let a_raw = load_bed(a_path)?;
    let b_raw = load_bed(b_path)?;

    let a_merged = merge_sorted(&a_raw);
    let b_merged = merge_sorted(&b_raw);

    let sum_a: u64 = a_merged
        .values()
        .flat_map(|v| v.iter())
        .map(|(s, e)| e - s)
        .sum();
    let sum_b: u64 = b_merged
        .values()
        .flat_map(|v| v.iter())
        .map(|(s, e)| e - s)
        .sum();

    let mut intersection: u64 = 0;
    let mut n_intersections: u64 = 0;

    for (chrom, a_segs) in &a_merged {
        let Some(b_segs) = b_merged.get(chrom) else {
            continue;
        };
        let mut j = 0usize;
        for &(as_, ae) in a_segs {
            while j < b_segs.len() && b_segs[j].1 <= as_ {
                j += 1;
            }
            let mut k = j;
            while k < b_segs.len() && b_segs[k].0 < ae {
                let (bs, be) = b_segs[k];
                let ov_s = as_.max(bs);
                let ov_e = ae.min(be);
                if ov_s < ov_e {
                    intersection += ov_e - ov_s;
                    n_intersections += 1;
                }
                k += 1;
            }
        }
    }

    let union = sum_a + sum_b - intersection;
    let jaccard_val = if union == 0 {
        0.0f64
    } else {
        intersection as f64 / union as f64
    };

    Ok(JaccardResult {
        intersection,
        union,
        jaccard: jaccard_val,
        n_intersections,
    })
}

/// Write jaccard result to output, matching bedtools' `header\nrow` format.
pub fn write_result(result: &JaccardResult, out: &mut dyn std::io::Write) -> Result<()> {
    writeln!(out, "intersection\tunion\tjaccard\tn_intersections").map_err(RsomicsError::Io)?;
    writeln!(
        out,
        "{}\t{}\t{}\t{}",
        result.intersection,
        result.union,
        format_6g(result.jaccard),
        result.n_intersections
    )
    .map_err(RsomicsError::Io)
}
