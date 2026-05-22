//! Byte-slice helpers shared by the streaming BED ops (intersect, subtract):
//! parse coordinate fields and locate A's trailing columns without allocating.

use rsomics_common::{Result, RsomicsError};

/// Byte index of the `n`-th tab (1-based) in `line`, or None if there are fewer.
pub fn nth_tab(line: &[u8], n: usize) -> Option<usize> {
    line.iter()
        .enumerate()
        .filter(|&(_, &b)| b == b'\t')
        .map(|(i, _)| i)
        .nth(n - 1)
}

/// A line's columns past `end` (BED4+: name/score/strand/...), kept verbatim
/// (leading tab included) so a streaming op can re-emit them with replaced
/// coordinates — matches bedtools. Empty for BED3.
pub fn rest_after_end(line: &[u8]) -> &[u8] {
    nth_tab(line, 3).map_or(&b""[..], |i| &line[i..])
}

/// Parse a required unsigned coordinate field, failing loud with line context.
pub fn parse_coord(f: Option<&[u8]>, lineno: usize, what: &str) -> Result<u64> {
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

/// True for lines a BED reader skips: blank, comment, or track/browser headers.
pub fn is_skippable(line: &[u8]) -> bool {
    line.is_empty() || line[0] == b'#' || line.starts_with(b"track") || line.starts_with(b"browser")
}
