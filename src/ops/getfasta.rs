use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

struct SeqGeom {
    length: usize,
    /// Byte offset of the sequence's first base in the FASTA buffer.
    offset: usize,
    line_bases: usize,
    line_width: usize,
}

/// Extract a region's bases out of the in-memory FASTA buffer.
///
/// Returns `None` for features bedtools getfasta skips entirely (emitting no
/// FASTA record, only a stderr warning): zero-length (`start >= end`) and those
/// running past the sequence end (`end > length`). Matching that skip set
/// keeps stdout byte-identical to bedtools.
///
/// The `.fai`-style geometry (offset of the sequence's first base, line layout)
/// lets us compute byte ranges with zero per-region file IO — the whole
/// extraction is in-RAM slicing. Wrapped sequences (`line_width > line_bases`)
/// carry an EOL separator every `line_bases` characters; we copy base runs
/// span-by-span into `scratch`, skipping separators. Single-line sequences hit
/// the contiguous fast path and borrow the buffer directly with no copy.
fn extract<'a>(
    bytes: &'a [u8],
    geom: &SeqGeom,
    start: usize,
    end: usize,
    scratch: &'a mut Vec<u8>,
) -> Option<&'a [u8]> {
    if start >= end || end > geom.length {
        return None;
    }

    if geom.line_width == geom.line_bases {
        // No newlines inside this sequence: the requested bases are a single
        // contiguous slice of the file buffer — borrow it directly, no copy.
        return Some(&bytes[geom.offset + start..geom.offset + end]);
    }

    scratch.clear();
    scratch.reserve(end - start);
    let mut pos = start;
    while pos < end {
        let line_idx = pos / geom.line_bases;
        let col = pos % geom.line_bases;
        let take = (geom.line_bases - col).min(end - pos);
        let from = geom.offset + line_idx * geom.line_width + col;
        scratch.extend_from_slice(&bytes[from..from + take]);
        pos += take;
    }
    Some(&scratch[..])
}

/// Scan FASTA geometry (`.fai`-equivalent layout) directly from the in-memory
/// buffer.
///
/// Reuses the buffer we already slurped instead of re-reading the file — one
/// FASTA-sized allocation total, not two. Per-sequence layout only needs the
/// first line's width, so the scan jumps to the next header once geometry is
/// known.
fn scan_geometry(data: &[u8]) -> HashMap<&str, SeqGeom> {
    let n = data.len();
    let find_nl = |from: usize| -> usize { memchr_nl(&data[from..]).map_or(n, |p| from + p) };
    let mut geom = HashMap::new();
    let mut i = 0;
    while i < n {
        if data[i] != b'>' {
            let nl = find_nl(i);
            i = if nl < n { nl + 1 } else { n };
            continue;
        }
        let hdr_nl = find_nl(i);
        let header = &data[i + 1..hdr_nl];
        let name_end = header
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(header.len());
        let name = std::str::from_utf8(&header[..name_end]).unwrap_or("");

        let mut p = if hdr_nl < n { hdr_nl + 1 } else { n };
        let offset = p;
        let (mut length, mut line_bases, mut line_width) = (0usize, 0usize, 0usize);
        let mut first = true;
        while p < n && data[p] != b'>' {
            let nl = find_nl(p);
            let raw_len = if nl < n { nl - p + 1 } else { nl - p };
            let mut content_end = nl;
            if content_end > p && data[content_end - 1] == b'\r' {
                content_end -= 1;
            }
            length += content_end - p;
            if first {
                line_bases = content_end - p;
                line_width = raw_len;
                first = false;
            }
            p = if nl < n { nl + 1 } else { n };
        }

        geom.insert(
            name,
            SeqGeom {
                length,
                offset,
                line_bases,
                line_width,
            },
        );
        i = p;
    }
    geom
}

pub fn getfasta(bed_path: &Path, fasta_path: &Path, output: &mut dyn Write) -> Result<u64> {
    // Slurp the whole FASTA once. bedtools getfasta loads the reference into
    // memory and slices regions from RAM; 100k BED rows against a per-region
    // fseek+read lose to that in-memory slicing, so we mirror it.
    let mut fasta = Vec::new();
    File::open(fasta_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fasta_path.display())))?
        .read_to_end(&mut fasta)
        .map_err(RsomicsError::Io)?;

    let geom = scan_geometry(&fasta);

    let bed = std::fs::read(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;

    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut scratch: Vec<u8> = Vec::new();
    let mut count: u64 = 0;

    for line in bed.split(|&b| b == b'\n') {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let mut fields = line.splitn(3, |&b| b == b'\t');
        let chrom = match fields.next() {
            Some(c) if !c.is_empty() => c,
            _ => continue,
        };
        let Some(start_b) = fields.next() else {
            continue;
        };
        let Some(end_b) = fields.next() else {
            continue;
        };
        // The end field may carry trailing BED columns; take up to the next tab.
        let end_b = end_b.split(|&b| b == b'\t').next().unwrap_or(end_b);

        let chrom_str = std::str::from_utf8(chrom)
            .map_err(|e| RsomicsError::InvalidInput(format!("bad chrom: {e}")))?;
        let start = parse_usize(start_b, "start")?;
        let end = parse_usize(end_b, "end")?;

        let g = geom.get(chrom_str).ok_or_else(|| {
            RsomicsError::InvalidInput(format!("sequence '{chrom_str}' not found in index"))
        })?;

        // bedtools skips (warns on stderr, emits nothing) zero-length features
        // and features running past the sequence end; mirror that so stdout
        // stays byte-identical.
        let Some(seq) = extract(&fasta, g, start, end, &mut scratch) else {
            if start >= end {
                eprintln!("Feature ({chrom_str}:{start}-{end}) has length = 0, Skipping.");
            } else {
                eprintln!(
                    "Feature ({chrom_str}:{start}-{end}) beyond the length of {chrom_str} size ({} bp).  Skipping.",
                    g.length
                );
            }
            continue;
        };
        writeln!(out, ">{chrom_str}:{start}-{end}").map_err(RsomicsError::Io)?;
        out.write_all(seq).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}

fn parse_usize(bytes: &[u8], what: &str) -> Result<usize> {
    let s = std::str::from_utf8(bytes)
        .map_err(|e| RsomicsError::InvalidInput(format!("bad {what}: {e}")))?
        .trim_end_matches('\r');
    s.parse()
        .map_err(|e| RsomicsError::InvalidInput(format!("bad {what}: {e}")))
}

fn memchr_nl(haystack: &[u8]) -> Option<usize> {
    haystack.iter().position(|&b| b == b'\n')
}
