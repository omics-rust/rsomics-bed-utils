use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::{IntervalIndex, IntervalSet, bed};

/// Report the intersected region of every A interval with each overlapping B
/// interval (`bedtools intersect -a A -b B` default). B is read once and put in
/// a coitrees index; A is streamed as byte slices and queried against it, so
/// the large side never builds per-record `Interval`/`IntervalSet` or clones
/// chrom strings. Output is A-file order, and per-A overlaps are emitted in
/// coordinate order — byte-identical to bedtools.
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
        if line.is_empty()
            || line[0] == b'#'
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }
        lineno += 1;
        let mut fields = line.split(|&c| c == b'\t');
        let chrom = fields.next().unwrap_or(b"");
        let start = parse_field(fields.next(), lineno, "start")?;
        let end = parse_field(fields.next(), lineno, "end")?;
        let chrom = std::str::from_utf8(chrom).map_err(|e| {
            RsomicsError::InvalidInput(format!("BED line {lineno}: non-UTF8 chrom: {e}"))
        })?;
        // A's columns past end (BED4+: name/score/strand/...), kept verbatim and
        // re-emitted per overlap with the coordinates replaced — matches bedtools.
        let rest = nth_tab(line, 3).map_or(&b""[..], |i| &line[i..]);

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

/// Byte index of the `n`-th tab (1-based) in `line`, or None if there are fewer.
fn nth_tab(line: &[u8], n: usize) -> Option<usize> {
    line.iter()
        .enumerate()
        .filter(|&(_, &b)| b == b'\t')
        .map(|(i, _)| i)
        .nth(n - 1)
}

fn parse_field(f: Option<&[u8]>, lineno: usize, what: &str) -> Result<u64> {
    let bytes =
        f.ok_or_else(|| RsomicsError::InvalidInput(format!("BED line {lineno}: missing {what}")))?;
    std::str::from_utf8(bytes)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "BED line {lineno}: bad {what} {:?}",
                String::from_utf8_lossy(bytes)
            ))
        })
}
