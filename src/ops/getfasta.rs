use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_fasta_index::fetch_seq;

pub fn getfasta(bed_path: &Path, fasta_path: &Path, output: &mut dyn Write) -> Result<u64> {
    let fai_path = fasta_path.with_extension(format!(
        "{}.fai",
        fasta_path.extension().unwrap_or_default().to_string_lossy()
    ));

    let file = File::open(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let chrom = fields[0];
        let start: usize = fields[1]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("bad start: {e}")))?;
        let end: usize = fields[2]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("bad end: {e}")))?;

        // Use fetch_seq with explicit 0-based coordinates (BED is 0-based).
        let seq = fetch_seq(fasta_path, &fai_path, chrom, start, end)?;
        writeln!(out, ">{chrom}:{start}-{end}").map_err(RsomicsError::Io)?;
        out.write_all(&seq).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
