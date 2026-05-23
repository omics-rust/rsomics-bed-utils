//! Mask FASTA bases overlapping BED intervals.
//!
//! Matches `bedtools maskfasta -fi ref.fa -bed regions.bed -fo out.fa`:
//! - Default: replace masked bases with `N`.
//! - `-soft`: replace masked bases with their lowercase equivalent.
//! - `-mc C`: replace masked bases with the given character.
//!
//! Output FASTA preserves the original per-sequence line width (taken from
//! the first data line of each sequence in the input), matching bedtools
//! exactly — output is byte-identical to `bedtools maskfasta`.
//!
//! Algorithm: O(N + M log M) where N = FASTA length, M = BED intervals.
//! 1. Parse BED into a per-chrom interval list; merge overlapping intervals.
//! 2. Stream the FASTA line by line. For each sequence, track the flat byte
//!    offset and emit bases, masking those that fall inside the merged
//!    intervals for that chrom. Re-wrap using the input line width.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// Masking mode.
pub enum MaskMode {
    /// Replace with `N` (or the given character).
    Hard(u8),
    /// Replace with the lowercase equivalent.
    Soft,
}

pub fn maskfasta(
    fasta_path: &Path,
    bed_path: &Path,
    mode: &MaskMode,
    output: &mut dyn Write,
) -> Result<()> {
    let intervals = load_intervals(bed_path)?;
    let out = &mut BufWriter::with_capacity(256 * 1024, output);
    stream_and_mask(fasta_path, &intervals, mode, out)
}

// ---------------------------------------------------------------------------
// BED loading — per-chrom merged intervals
// ---------------------------------------------------------------------------

/// Sorted, merged half-open intervals [start, end) for one chrom.
type Iv = (u64, u64);

fn load_intervals(bed_path: &Path) -> Result<HashMap<String, Vec<Iv>>> {
    let file = File::open(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
    let mut by_chrom: HashMap<String, Vec<Iv>> = HashMap::new();

    for line in BufReader::new(file).lines() {
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
        by_chrom.entry(chrom).or_default().push((start, end));
    }

    // Sort and merge each chrom's intervals.
    for ivs in by_chrom.values_mut() {
        ivs.sort_unstable_by_key(|&(s, e)| (s, e));
        let mut merged: Vec<Iv> = Vec::with_capacity(ivs.len());
        for &(s, e) in ivs.iter() {
            if let Some(last) = merged.last_mut()
                && s < last.1
            {
                last.1 = last.1.max(e);
                continue;
            }
            merged.push((s, e));
        }
        *ivs = merged;
    }

    Ok(by_chrom)
}

// ---------------------------------------------------------------------------
// FASTA streaming with in-place masking
// ---------------------------------------------------------------------------

fn stream_and_mask(
    fasta_path: &Path,
    intervals: &HashMap<String, Vec<Iv>>,
    mode: &MaskMode,
    out: &mut impl Write,
) -> Result<()> {
    let file = File::open(fasta_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fasta_path.display())))?;

    // State for the current sequence being accumulated.
    let mut cur_name: Option<String> = None;
    // The sequence bytes, flat, for the current record.
    let mut cur_seq: Vec<u8> = Vec::with_capacity(1 << 20);
    // Line width (from first data line) of the current sequence.
    let mut cur_line_width: usize = 0;

    for line_res in BufReader::new(file).lines() {
        let line = line_res.map_err(RsomicsError::Io)?;
        if let Some(rest) = line.strip_prefix('>') {
            // Flush previous sequence.
            if let Some(name) = cur_name.take() {
                flush_seq(out, &name, &cur_seq, cur_line_width, intervals, mode)?;
                cur_seq.clear();
                cur_line_width = 0;
            }
            // Extract name (everything before first space).
            let name = rest.split_whitespace().next().unwrap_or("").to_string();
            cur_name = Some(name);
        } else if cur_name.is_some() {
            let data = line.as_bytes();
            if cur_line_width == 0 && !data.is_empty() {
                cur_line_width = data.len();
            }
            cur_seq.extend_from_slice(data);
        }
    }

    // Flush final sequence.
    if let Some(name) = cur_name {
        flush_seq(out, &name, &cur_seq, cur_line_width, intervals, mode)?;
    }

    out.flush().map_err(RsomicsError::Io)
}

fn flush_seq(
    out: &mut impl Write,
    name: &str,
    seq: &[u8],
    line_width: usize,
    intervals: &HashMap<String, Vec<Iv>>,
    mode: &MaskMode,
) -> Result<()> {
    writeln!(out, ">{name}").map_err(RsomicsError::Io)?;

    if seq.is_empty() {
        return Ok(());
    }

    let ivs = intervals.get(name).map(|v| v.as_slice()).unwrap_or(&[]);
    let wrap = if line_width == 0 {
        seq.len()
    } else {
        line_width
    };

    // Walk through seq in wrap-sized chunks, masking as we go.
    // `flat_pos` tracks the absolute position in the sequence for interval lookup.
    let mut flat_pos: u64 = 0;
    // iv_idx: pointer into the sorted merged interval list — O(N+M) total.
    let mut iv_idx: usize = 0;

    for chunk in seq.chunks(wrap) {
        // Emit each byte in the chunk, masking those inside an interval.
        for &b in chunk {
            // Advance interval pointer past intervals that end before flat_pos.
            while iv_idx < ivs.len() && ivs[iv_idx].1 <= flat_pos {
                iv_idx += 1;
            }
            let masked = iv_idx < ivs.len() && ivs[iv_idx].0 <= flat_pos;
            let out_byte = if masked {
                match mode {
                    MaskMode::Hard(ch) => *ch,
                    MaskMode::Soft => b.to_ascii_lowercase(),
                }
            } else {
                b
            };
            out.write_all(&[out_byte]).map_err(RsomicsError::Io)?;
            flat_pos += 1;
        }
        writeln!(out).map_err(RsomicsError::Io)?;
    }

    Ok(())
}
