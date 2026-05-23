//! Randomly reposition BED intervals within chromosome bounds.
//!
//! Matches `bedtools shuffle -i in.bed -g genome.txt` semantics:
//! - Each interval is placed at a random start position such that the entire
//!   interval fits within the chromosome. The interval length is preserved.
//! - `-seed N` makes placement reproducible.
//! - `-excl excl.bed` rejects positions that overlap any excluded region.
//! - `-chrom` keeps each interval on its original chromosome; without it the
//!   destination chromosome is chosen proportionally to chrom length minus
//!   interval length (same weighting as bedtools).
//!
//! The RNG output intentionally differs from bedtools (different algorithm),
//! so compat tests verify INVARIANTS rather than byte identity:
//! - Same number of output intervals as input.
//! - Each output interval has the same length as its input.
//! - Each output interval lands within its chrom's bounds.
//! - With `-excl`, no output interval overlaps an excluded region.
//! - With `-chrom`, the output chrom equals the input chrom.
//!
//! Uses `rand::rngs::StdRng` (ChaCha-based; already a dep in Cargo.toml).
//! Maximum retry limit per interval: 10 000 attempts before error.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rsomics_common::{Result, RsomicsError};

const MAX_TRIES: u64 = 10_000;

/// Parsed BED record (first 3 columns + trailing raw fields for passthrough).
struct Rec {
    chrom: String,
    start: u64,
    end: u64,
    /// Columns 4+ verbatim (empty string if absent).
    rest: String,
}

pub fn shuffle(
    input_path: &Path,
    genome_path: &Path,
    seed: Option<u64>,
    excl_path: Option<&Path>,
    same_chrom: bool,
    output: &mut dyn Write,
) -> Result<()> {
    let genome = load_genome(genome_path)?;
    let excl = match excl_path {
        Some(p) => load_intervals(p)?,
        None => HashMap::new(),
    };
    let records = load_bed(input_path)?;

    // Chroms as an ordered vec for weighted selection (by max-placeable length).
    let chrom_list: Vec<(&String, u64)> = genome.iter().map(|(k, &v)| (k, v)).collect();

    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    let mut out = BufWriter::with_capacity(64 * 1024, output);

    for rec in &records {
        let len = rec.end.saturating_sub(rec.start);

        // Build the candidate chrom set: chroms large enough to hold `len`.
        // With -chrom: only the record's own chrom (if large enough).
        let candidates: Vec<(&str, u64)> = if same_chrom {
            let chrom_len = genome.get(&rec.chrom).copied().unwrap_or(0);
            if chrom_len < len {
                return Err(RsomicsError::InvalidInput(format!(
                    "interval {}:{}-{} (len={len}) longer than chrom {chrom_len}",
                    rec.chrom, rec.start, rec.end
                )));
            }
            vec![(rec.chrom.as_str(), chrom_len)]
        } else {
            chrom_list
                .iter()
                .filter_map(|&(name, clen)| {
                    if clen >= len {
                        Some((name.as_str(), clen))
                    } else {
                        None
                    }
                })
                .collect()
        };

        if candidates.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "no chromosome long enough to place interval of length {len}"
            )));
        }

        // Weighted chrom selection: weight = chrom_len - len + 1 (placeable positions).
        let total_weight: u64 = candidates.iter().map(|(_, clen)| clen - len + 1).sum();

        let mut placed = false;
        for _ in 0..MAX_TRIES {
            // Choose a chrom proportional to placeable positions.
            let chosen_chrom = if candidates.len() == 1 {
                candidates[0].0
            } else {
                let pick: u64 = rng.gen_range(0..total_weight);
                let mut acc = 0u64;
                let mut chosen = candidates[0].0;
                for &(name, clen) in &candidates {
                    acc += clen - len + 1;
                    if pick < acc {
                        chosen = name;
                        break;
                    }
                }
                chosen
            };

            let chrom_len = genome[chosen_chrom];
            let new_start: u64 = rng.gen_range(0..=(chrom_len - len));
            let new_end = new_start + len;

            if let Some(excl_ivs) = excl.get(chosen_chrom)
                && overlaps_any(new_start, new_end, excl_ivs)
            {
                continue;
            }

            // Emit the interval.
            if rec.rest.is_empty() {
                writeln!(out, "{chosen_chrom}\t{new_start}\t{new_end}")
            } else {
                writeln!(out, "{chosen_chrom}\t{new_start}\t{new_end}\t{}", rec.rest)
            }
            .map_err(RsomicsError::Io)?;

            placed = true;
            break;
        }

        if !placed {
            return Err(RsomicsError::InvalidInput(format!(
                "could not place interval {}:{}-{} after {MAX_TRIES} tries (too many excl regions?)",
                rec.chrom, rec.start, rec.end
            )));
        }
    }

    out.flush().map_err(RsomicsError::Io)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn load_genome(path: &Path) -> Result<HashMap<String, u64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut map = HashMap::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut f = line.splitn(3, '\t');
        let chrom = match f.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        let len: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        map.insert(chrom, len);
    }
    Ok(map)
}

type Iv = (u64, u64);

fn load_intervals(path: &Path) -> Result<HashMap<String, Vec<Iv>>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut by_chrom: HashMap<String, Vec<Iv>> = HashMap::new();

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
        let start: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let end: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        by_chrom.entry(chrom).or_default().push((start, end));
    }

    // Sort each chrom's intervals for binary-search overlap check.
    for ivs in by_chrom.values_mut() {
        ivs.sort_unstable_by_key(|&(s, e)| (s, e));
    }

    Ok(by_chrom)
}

fn load_bed(path: &Path) -> Result<Vec<Rec>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut recs = Vec::new();
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
        let start: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let end: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let rest = f.next().unwrap_or("").trim_end_matches('\n').to_string();
        recs.push(Rec {
            chrom,
            start,
            end,
            rest,
        });
    }
    Ok(recs)
}

/// Return true if [start, end) overlaps any interval in the sorted list.
fn overlaps_any(start: u64, end: u64, ivs: &[Iv]) -> bool {
    // Binary search for the first interval whose end > start.
    let idx = ivs.partition_point(|&(_, e)| e <= start);
    if idx < ivs.len() && ivs[idx].0 < end {
        return true;
    }
    false
}
