use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::{IntervalIndex, IntervalSet, bed};

use crate::ops::byteparse::{is_skippable, parse_coord, rest_after_end};

/// Report the intersected region of every A interval with each overlapping B
/// interval (`bedtools intersect -a A -b B` default). B is read once and put in
/// a coitrees index; A is streamed as byte slices and queried against it, so
/// the large side never builds per-record `Interval`/`IntervalSet` or clones
/// chrom strings. A's columns past end are preserved per overlap. Output is
/// A-file order; per-A overlaps are emitted in coordinate order.
pub fn intersect(a_path: &Path, b_path: &Path, output: &mut dyn Write) -> Result<()> {
    let b_ivs = bed::read(File::open(b_path).map_err(RsomicsError::Io)?)?;
    let b_set: IntervalSet = b_ivs.into_iter().collect();
    let index = IntervalIndex::build(&b_set);

    let mut data = Vec::new();
    File::open(a_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", a_path.display())))?
        .read_to_end(&mut data)
        .map_err(RsomicsError::Io)?;

    let mut out = BufWriter::new(output);
    let mut hits: Vec<(u64, u64)> = Vec::new();
    let mut lineno = 0usize;
    for raw in data.split(|&b| b == b'\n') {
        let line = match raw.last() {
            Some(b'\r') => &raw[..raw.len() - 1],
            _ => raw,
        };
        if is_skippable(line) {
            continue;
        }
        lineno += 1;
        let mut fields = line.split(|&c| c == b'\t');
        let chrom = fields.next().unwrap_or(b"");
        let start = parse_coord(fields.next(), lineno, "start")?;
        let end = parse_coord(fields.next(), lineno, "end")?;
        let chrom = std::str::from_utf8(chrom).map_err(|e| {
            RsomicsError::InvalidInput(format!("BED line {lineno}: non-UTF8 chrom: {e}"))
        })?;
        let rest = rest_after_end(line);

        hits.clear();
        index.for_each_overlap(chrom, start, end, |bi| {
            let lo = start.max(bi.start);
            let hi = end.min(bi.end);
            if hi > lo {
                hits.push((lo, hi));
            }
        });
        hits.sort_unstable();
        for &(lo, hi) in &hits {
            write!(out, "{chrom}\t{lo}\t{hi}").map_err(RsomicsError::Io)?;
            out.write_all(rest).map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
        }
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
