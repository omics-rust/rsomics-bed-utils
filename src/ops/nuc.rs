#![allow(clippy::cast_precision_loss)]

//! Nucleotide content per BED interval (`bedtools nuc`).
//!
//! Output columns match bedtools exactly:
//!   - Dynamic header: one `#N_usercol` / `N_usercol` column per BED column
//!     (the first gets the `#` prefix), then pct_at, pct_gc, num_A/C/G/T/N/oth,
//!     seq_len.  Column numbers are 1-based and continue from the BED column count.
//!   - pct_at and pct_gc use `%.6f` (always 6 decimal places).
//!   - count columns are integers.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

pub fn bed_nuc(bed_path: &Path, fasta_path: &Path, output: &mut dyn Write) -> Result<u64> {
    let seqs = load_fasta(fasta_path)?;

    // Two-pass: first count max BED columns to build correct header, then emit.
    let raw_lines: Vec<String> = {
        let file = File::open(bed_path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
        BufReader::new(file)
            .lines()
            .map(|l| l.map_err(RsomicsError::Io))
            .collect::<Result<Vec<_>>>()?
    };

    let bed_lines: Vec<&str> = raw_lines
        .iter()
        .map(String::as_str)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();

    let max_cols = bed_lines
        .iter()
        .map(|l| l.split('\t').count())
        .max()
        .unwrap_or(3)
        .max(3);

    let mut out = BufWriter::with_capacity(64 * 1024, output);

    // Build and write header.
    let stats_start = max_cols + 1;
    let mut hdr = String::new();
    for i in 1..=max_cols {
        if i == 1 {
            hdr.push_str(&format!("#{i}_usercol"));
        } else {
            hdr.push_str(&format!("\t{i}_usercol"));
        }
    }
    let n = stats_start;
    hdr.push_str(&format!(
        "\t{}_pct_at\t{}_pct_gc\t{}_num_A\t{}_num_C\t{}_num_G\t{}_num_T\t{}_num_N\t{}_num_oth\t{}_seq_len",
        n, n+1, n+2, n+3, n+4, n+5, n+6, n+7, n+8
    ));
    writeln!(out, "{hdr}").map_err(RsomicsError::Io)?;

    let mut count: u64 = 0;

    for line in &bed_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let chrom = fields[0];
        let start: usize = fields[1].parse().unwrap_or(0);
        let end: usize = fields[2].parse().unwrap_or(0);

        let bc = if let Some(seq) = seqs.get(chrom) {
            let s = start.min(seq.len());
            let e = end.min(seq.len());
            count_bases(&seq[s..e])
        } else {
            BaseCounts::default()
        };

        let len = bc.adenine + bc.cytosine + bc.guanine + bc.thymine + bc.ambiguous + bc.other;
        let pct_at = if len > 0 {
            (bc.adenine + bc.thymine) as f64 / len as f64
        } else {
            0.0
        };
        let pct_gc = if len > 0 {
            (bc.guanine + bc.cytosine) as f64 / len as f64
        } else {
            0.0
        };

        // All BED columns (exactly as in input).
        for (i, &f) in fields.iter().enumerate() {
            if i > 0 {
                write!(out, "\t").map_err(RsomicsError::Io)?;
            }
            write!(out, "{f}").map_err(RsomicsError::Io)?;
        }
        // Pad to max_cols if this record has fewer columns.
        for _ in fields.len()..max_cols {
            write!(out, "\t").map_err(RsomicsError::Io)?;
        }
        // Stats columns.
        writeln!(
            out,
            "\t{pct_at:.6}\t{pct_gc:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{len}",
            bc.adenine, bc.cytosine, bc.guanine, bc.thymine, bc.ambiguous, bc.other
        )
        .map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}

#[derive(Default)]
struct BaseCounts {
    adenine: u64,
    cytosine: u64,
    guanine: u64,
    thymine: u64,
    ambiguous: u64,
    other: u64,
}

fn count_bases(seq: &[u8]) -> BaseCounts {
    let mut bc = BaseCounts::default();
    for &base in seq {
        match base.to_ascii_uppercase() {
            b'A' => bc.adenine += 1,
            b'C' => bc.cytosine += 1,
            b'G' => bc.guanine += 1,
            b'T' => bc.thymine += 1,
            b'N' => bc.ambiguous += 1,
            _ => bc.other += 1,
        }
    }
    bc
}

fn load_fasta(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    if std::fs::metadata(path).is_ok_and(|m| m.len() == 0) {
        return Err(RsomicsError::InvalidInput("empty FASTA".into()));
    }

    let mut reader = parse_fastx_file(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;

    let mut seqs = HashMap::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(|e| RsomicsError::InvalidInput(format!("reading: {e}")))?;
        let id = std::str::from_utf8(record.id())
            .unwrap_or("unknown")
            .to_string();
        seqs.insert(id, record.seq().to_vec());
    }

    Ok(seqs)
}
