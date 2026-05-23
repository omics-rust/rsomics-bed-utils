//! Fisher's exact test on whether two BED files overlap more than expected.
//!
//! Matches `bedtools fisher -a A -b B -g genome.txt` exactly.
//!
//! Algorithm (from bedtools Fisher.cpp / kfunc.cpp, MIT licensed):
//!
//! 1. Count A intervals (queryCounts) and B intervals (dbCounts).
//! 2. Count overlapping pairs: overlapCounts (each A × each B it hits).
//! 3. Sum raw lengths: queryUnion (sum of A lengths) and dbUnion (sum of B lengths).
//! 4. Compute n_possible via the bedtools heuristic:
//!    - dMean = 1 + dbUnion / dbCounts
//!    - qMean = 1 + queryUnion / queryCounts
//!    - bMean = qMean + dMean
//!    - n22_full = max(n11+n12+n21, floor(genomeSize / bMean))
//! 5. Build the 2×2 contingency table:
//!    - n11 = overlapCounts, n12 = max(0, queryCounts - overlapCounts)
//!    - n21 = max(0, dbCounts - overlapCounts), n22 = n22_full - n11 - n12 - n21
//! 6. Run kt_fisher_exact (hypergeometric-based Fisher's exact test) for left /
//!    right / two-tail p-values and the odds ratio.
//! 7. Print the bedtools-identical multi-line report.

#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::ops::byteparse::is_skippable;

// ---------------------------------------------------------------------------
// BED loading helpers
// ---------------------------------------------------------------------------

/// One BED set parsed into struct-of-arrays form. Chrom names are interned to a
/// dense `u32` id shared across both A and B (`interner`), so the overlap sweep
/// works on integer ids instead of re-comparing strings. Coordinates are `i32`
/// (genomic positions fit comfortably below `i32::MAX` ≈ 2.1 Gb), halving the
/// per-interval footprint. `union` is Σ(end−start) and `count` the interval
/// count, both accumulated during the single parse pass so fisher needs no
/// second walk over the records.
struct BedSet {
    chrom_id: Vec<u32>,
    start: Vec<i32>,
    end: Vec<i32>,
    union: i64,
    count: i64,
}

/// Parse one coordinate field from raw bytes. BED coordinates are non-negative
/// integers below `i32::MAX`; anything else is a malformed file and fails loud.
fn parse_coord(bytes: &[u8], lineno: usize, what: &str) -> Result<i32> {
    if bytes.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "BED line {lineno}: empty {what}"
        )));
    }
    let mut v: i64 = 0;
    for &c in bytes {
        if !c.is_ascii_digit() {
            return Err(RsomicsError::InvalidInput(format!(
                "BED line {lineno}: bad {what} {:?}",
                String::from_utf8_lossy(bytes)
            )));
        }
        v = v * 10 + i64::from(c - b'0');
        if v > i64::from(i32::MAX) {
            return Err(RsomicsError::InvalidInput(format!(
                "BED line {lineno}: {what} exceeds i32::MAX {:?}",
                String::from_utf8_lossy(bytes)
            )));
        }
    }
    Ok(v as i32)
}

/// Load a BED file in one `read_to_end` + byte-slice pass: no per-line `String`
/// allocation, no `to_string()` per chrom, integer fields parsed directly from
/// bytes. Chrom names are interned through the shared `interner`. Sorted input
/// is grouped by chrom contiguously, so the previous chrom's id is reused
/// without a hash lookup on the common case.
fn load_bed(path: &Path, interner: &mut HashMap<Vec<u8>, u32>) -> Result<BedSet> {
    let mut data = Vec::new();
    File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?
        .read_to_end(&mut data)
        .map_err(RsomicsError::Io)?;

    let mut chrom_id = Vec::new();
    let mut start = Vec::new();
    let mut end = Vec::new();
    let mut union: i64 = 0;
    let mut count: i64 = 0;

    let mut last_chrom: Vec<u8> = Vec::new();
    let mut last_id: u32 = u32::MAX;
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

        let t1 = line.iter().position(|&c| c == b'\t').ok_or_else(|| {
            RsomicsError::InvalidInput(format!("BED line {lineno}: missing start"))
        })?;
        let chrom = &line[..t1];
        let rest = &line[t1 + 1..];
        let t2 = rest
            .iter()
            .position(|&c| c == b'\t')
            .ok_or_else(|| RsomicsError::InvalidInput(format!("BED line {lineno}: missing end")))?;
        let s = parse_coord(&rest[..t2], lineno, "start")?;
        let rest2 = &rest[t2 + 1..];
        let t3 = rest2
            .iter()
            .position(|&c| c == b'\t')
            .unwrap_or(rest2.len());
        let e = parse_coord(&rest2[..t3], lineno, "end")?;

        let id = if chrom == last_chrom.as_slice() {
            last_id
        } else {
            let next = interner.len() as u32;
            let id = *interner.entry(chrom.to_vec()).or_insert(next);
            last_chrom.clear();
            last_chrom.extend_from_slice(chrom);
            last_id = id;
            id
        };

        chrom_id.push(id);
        start.push(s);
        end.push(e);
        union += i64::from(e - s);
        count += 1;
    }

    Ok(BedSet {
        chrom_id,
        start,
        end,
        union,
        count,
    })
}

/// Count every (A, B) pair that overlaps under half-open `[start, end)`
/// semantics, summed across all chroms shared by both sets.
///
/// Per chrom, A and B are each sorted by start, then swept by merging the two
/// start streams. Two min-heaps hold the end coordinates of the currently-open
/// A and B intervals; an interval expires (is popped) once the sweep coordinate
/// reaches its end. When an A interval opens it overlaps every still-open B
/// (and vice versa). Coincident starts are processed A-first then B-first so a
/// pair sharing a start coordinate is counted exactly once — matching the
/// close-before-open / A-before-B event ordering bedtools' sort implies.
///
/// This avoids the 2·(|A|+|B|) endpoint-event vector and its sort: it allocates
/// only the small per-side active heaps (bounded by the depth of overlap, not
/// the interval count).
fn count_overlaps(a: &BedSet, b: &BedSet, nchrom: usize) -> i64 {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;

    let mut a_groups: Vec<Vec<u32>> = vec![Vec::new(); nchrom];
    let mut b_groups: Vec<Vec<u32>> = vec![Vec::new(); nchrom];
    for (i, &c) in a.chrom_id.iter().enumerate() {
        a_groups[c as usize].push(i as u32);
    }
    for (i, &c) in b.chrom_id.iter().enumerate() {
        b_groups[c as usize].push(i as u32);
    }

    let mut overlap: i64 = 0;
    let mut a_ends: BinaryHeap<Reverse<i32>> = BinaryHeap::new();
    let mut b_ends: BinaryHeap<Reverse<i32>> = BinaryHeap::new();

    for chrom in 0..nchrom {
        let ag = &mut a_groups[chrom];
        let bg = &mut b_groups[chrom];
        if ag.is_empty() || bg.is_empty() {
            continue;
        }
        ag.sort_unstable_by_key(|&i| a.start[i as usize]);
        bg.sort_unstable_by_key(|&i| b.start[i as usize]);

        a_ends.clear();
        b_ends.clear();
        let (na, nb) = (ag.len(), bg.len());
        let mut ia = 0usize;
        let mut ib = 0usize;
        loop {
            let a_next = ag.get(ia).map(|&i| a.start[i as usize]);
            let b_next = bg.get(ib).map(|&i| b.start[i as usize]);
            let coord = match (a_next, b_next) {
                (None, None) => break,
                (Some(x), None) => x,
                (None, Some(y)) => y,
                (Some(x), Some(y)) => x.min(y),
            };
            while a_ends.peek().is_some_and(|&Reverse(e)| e <= coord) {
                a_ends.pop();
            }
            while b_ends.peek().is_some_and(|&Reverse(e)| e <= coord) {
                b_ends.pop();
            }
            while ia < na && a.start[ag[ia] as usize] == coord {
                overlap += b_ends.len() as i64;
                a_ends.push(Reverse(a.end[ag[ia] as usize]));
                ia += 1;
            }
            while ib < nb && b.start[bg[ib] as usize] == coord {
                overlap += a_ends.len() as i64;
                b_ends.push(Reverse(b.end[bg[ib] as usize]));
                ib += 1;
            }
        }
    }
    overlap
}

fn load_genome(path: &Path) -> Result<i64> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut total: i64 = 0;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut f = line.splitn(2, '\t');
        f.next(); // chrom
        if let Some(sz) = f.next().and_then(|s| s.trim().parse::<i64>().ok()) {
            total += sz;
        }
    }
    Ok(total)
}

// ---------------------------------------------------------------------------
// Fisher's exact test (kt_fisher_exact) — ported from bedtools kfunc.cpp
// (MIT license, Heng Li).  The hypergeometric recursion is kept as-is.
// ---------------------------------------------------------------------------

fn lbinom(n: i64, k: i64) -> f64 {
    if k == 0 || n == k {
        return 0.0;
    }
    lgamma((n + 1) as f64) - lgamma((k + 1) as f64) - lgamma((n - k + 1) as f64)
}

fn lgamma(x: f64) -> f64 {
    // Stirling via the standard `libm` lgamma.  In Rust we can call f64's own.
    statrs_lgamma(x)
}

// Use Lanczos approximation matching the kf_lgamma from bedtools exactly.
fn statrs_lgamma(z: f64) -> f64 {
    let mut x = 0.0f64;
    x += 0.165_947_018_740_846_2e-6 / (z + 7.0);
    x += 0.993_493_711_393_074_8e-5 / (z + 6.0);
    x -= 0.138_571_033_129_652_6 / (z + 5.0);
    x += 12.507_343_240_090_56 / (z + 4.0);
    x -= 176.615_029_149_838_6 / (z + 3.0);
    x += 771.323_428_775_767_4 / (z + 2.0);
    x -= 1_259.139_216_722_289 / (z + 1.0);
    x += 676.520_368_121_883_5 / z;
    x += 0.999_999_999_999_518_3;
    x.ln() - 5.581_061_466_795_328 - z + (z - 0.5) * (z + 6.5).ln()
}

fn hypergeo(n11: i64, n1_: i64, n_1: i64, n: i64) -> f64 {
    (lbinom(n1_, n11) + lbinom(n - n1_, n_1 - n11) - lbinom(n, n_1)).exp()
}

struct HgAcc {
    n11: i64,
    n1_: i64,
    n_1: i64,
    n: i64,
    p: f64,
}

fn hypergeo_acc(n11: i64, n1_: i64, n_1: i64, n: i64, aux: &mut HgAcc) -> f64 {
    if n1_ != 0 || n_1 != 0 || n != 0 {
        aux.n11 = n11;
        aux.n1_ = n1_;
        aux.n_1 = n_1;
        aux.n = n;
    } else if n11 % 11 != 0 && n11 + aux.n - aux.n1_ - aux.n_1 != 0 {
        if n11 == aux.n11 + 1 {
            aux.p *= (aux.n1_ - aux.n11) as f64 / n11 as f64 * (aux.n_1 - aux.n11) as f64
                / (n11 + aux.n - aux.n1_ - aux.n_1) as f64;
            aux.n11 = n11;
            return aux.p;
        }
        if n11 == aux.n11 - 1 {
            aux.p *= aux.n11 as f64 / (aux.n1_ - n11) as f64
                * (aux.n11 + aux.n - aux.n1_ - aux.n_1) as f64
                / (aux.n_1 - n11) as f64;
            aux.n11 = n11;
            return aux.p;
        }
        aux.n11 = n11;
    } else {
        aux.n11 = n11;
    }
    aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n);
    aux.p
}

fn kt_fisher_exact(n11: i64, n12: i64, n21: i64, n22: i64) -> (f64, f64, f64) {
    let n1_ = n11 + n12;
    let n_1 = n11 + n21;
    let n = n11 + n12 + n21 + n22;
    let max = n_1.min(n1_);
    let min = (n1_ + n_1 - n).max(0);

    if min == max {
        return (1.0, 1.0, 1.0);
    }

    let mut aux = HgAcc {
        n11: 0,
        n1_: 0,
        n_1: 0,
        n: 0,
        p: 0.0,
    };
    let q = hypergeo_acc(n11, n1_, n_1, n, &mut aux);

    if q == 0.0 {
        // n11 left or right of mode
        if n11 * (n + 2) < (n_1 + 1) * (n1_ + 1) {
            return (0.0, 1.0, 0.0);
        } else {
            return (1.0, 0.0, 0.0);
        }
    }

    // left tail
    let mut p = hypergeo_acc(min, 0, 0, 0, &mut aux);
    let mut i = min + 1;
    let mut left = 0.0f64;
    while p < 0.999_999_99 * q && i <= max {
        left += p;
        p = hypergeo_acc(i, 0, 0, 0, &mut aux);
        i += 1;
    }
    i -= 1;
    if p < 1.000_000_01 * q {
        left += p;
    } else {
        i -= 1;
    }

    // right tail
    p = hypergeo_acc(max, 0, 0, 0, &mut aux);
    let mut j = max - 1;
    let mut right = 0.0f64;
    while p < 0.999_999_99 * q && j >= 0 {
        right += p;
        p = hypergeo_acc(j, 0, 0, 0, &mut aux);
        j -= 1;
    }
    j += 1;
    if p < 1.000_000_01 * q {
        right += p;
    } else {
        j += 1;
    }

    // two-tail is computed before left/right adjustment (matches bedtools).
    let two = (left + right).min(1.0);

    // Adjust the non-dominant tail so the reported left/right are consistent.
    if (i - n11).abs() < (j - n11).abs() {
        right = 1.0 - left + q;
    } else {
        left = 1.0 - right + q;
    }

    (left, right, two)
}

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

/// Format a float using C's `%.5g` semantics: 5 significant figures, no
/// trailing zeros, no trailing decimal point.  Special cases: exactly 0 → "0",
/// exactly 1 → "1".
///
/// Scientific notation uses at least two exponent digits (e.g. `1.3582e-06`),
/// matching C's printf which always writes `e+XX` / `e-XX` with a 2-digit minimum.
fn fmt_g5(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    if v == 1.0 {
        return "1".to_string();
    }
    // Determine magnitude for %g-style formatting.
    let mag = v.abs().log10().floor() as i32;
    // %g uses fixed notation when -4 <= exp < precision (5).
    if (-4..5).contains(&mag) {
        // fixed: number of decimals = 4 - mag (so total sig figs = 5)
        let decimals = (4 - mag).max(0) as usize;
        let s = format!("{:.prec$}", v, prec = decimals);
        // Strip trailing zeros after decimal point, then bare decimal point.
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    } else {
        // Scientific notation: 4 decimal places = 5 sig figs.
        // Rust's {:.4e} writes e.g. "1.3582e-6"; C printf writes "1.3582e-06".
        // Normalise exponent to minimum 2 digits to match C behaviour.
        let raw = format!("{:.4e}", v);
        normalize_exp_2digits(raw)
    }
}

/// Ensure the exponent in a scientific-notation string has at least 2 digits,
/// matching C's `printf` style (`e-6` → `e-06`, `e+6` → `e+06`).
fn normalize_exp_2digits(s: String) -> String {
    if let Some(e_pos) = s.find('e') {
        let (mantissa, exp_part) = s.split_at(e_pos);
        // exp_part starts with 'e', followed by sign (+/-) and digits.
        let after_e = &exp_part[1..]; // skip 'e'
        let (sign, digits) = if after_e.starts_with(['+', '-']) {
            (&after_e[..1], &after_e[1..])
        } else {
            ("+", after_e)
        };
        if digits.len() < 2 {
            format!("{mantissa}e{sign}{digits:0>2}")
        } else {
            s
        }
    } else {
        s
    }
}

// ---------------------------------------------------------------------------
// Main fisher function
// ---------------------------------------------------------------------------

pub fn fisher(a_path: &Path, b_path: &Path, genome_path: &Path, out: &mut dyn Write) -> Result<()> {
    let mut interner: HashMap<Vec<u8>, u32> = HashMap::new();
    let a = load_bed(a_path, &mut interner)?;
    let b = load_bed(b_path, &mut interner)?;
    let genome_size = load_genome(genome_path)?;

    let query_counts = a.count;
    let db_counts = b.count;
    let query_union = a.union;
    let db_union = b.union;

    let overlap_counts = count_overlaps(&a, &b, interner.len());

    // Contingency table entries.
    let n11 = overlap_counts;
    let n12 = (query_counts - overlap_counts).max(0);
    let n21 = (db_counts - overlap_counts).max(0);

    // Bedtools heuristic for n_possible (n22_full).
    let d_mean = 1.0 + db_union as f64 / db_counts as f64;
    let q_mean = 1.0 + query_union as f64 / query_counts as f64;
    let b_mean = q_mean + d_mean;
    let n22_full = (n21 + n12 + n11).max((genome_size as f64 / b_mean) as i64);
    let n22 = (n22_full - n12 - n21 - n11).max(0);

    // Fisher's exact test.
    let (left, right, two) = kt_fisher_exact(n11, n12, n21, n22);

    // Odds ratio (can be inf when n12=0 or n21=0).
    let ratio = if n12 == 0 || n21 == 0 {
        f64::INFINITY
    } else {
        (n11 as f64 / n12 as f64) / (n21 as f64 / n22 as f64)
    };

    // Print report matching bedtools format exactly.
    writeln!(out, "# Number of query intervals: {query_counts}").map_err(RsomicsError::Io)?;
    writeln!(out, "# Number of db intervals: {db_counts}").map_err(RsomicsError::Io)?;
    writeln!(out, "# Number of overlaps: {overlap_counts}").map_err(RsomicsError::Io)?;
    writeln!(
        out,
        "# Number of possible intervals (estimated): {n22_full}"
    )
    .map_err(RsomicsError::Io)?;
    writeln!(
        out,
        "# phyper({n11} - 1, {query_counts}, {n22_full} - {query_counts}, {db_counts}, lower.tail=F)"
    )
    .map_err(RsomicsError::Io)?;
    writeln!(out, "# Contingency Table Of Counts").map_err(RsomicsError::Io)?;
    writeln!(out, "#_________________________________________").map_err(RsomicsError::Io)?;
    writeln!(
        out,
        "#           | {:<12} | {:<12} |",
        " in -b", "not in -b"
    )
    .map_err(RsomicsError::Io)?;
    writeln!(out, "#     in -a | {n11:<12} | {n12:<12} |").map_err(RsomicsError::Io)?;
    writeln!(out, "# not in -a | {n21:<12} | {n22:<12} |").map_err(RsomicsError::Io)?;
    writeln!(out, "#_________________________________________").map_err(RsomicsError::Io)?;
    writeln!(out, "# p-values for fisher's exact test").map_err(RsomicsError::Io)?;
    writeln!(out, "left\tright\ttwo-tail\tratio").map_err(RsomicsError::Io)?;

    if ratio.is_nan() || ratio.is_infinite() {
        writeln!(
            out,
            "{}\t{}\t{}\tinf",
            fmt_g5(left),
            fmt_g5(right),
            fmt_g5(two)
        )
        .map_err(RsomicsError::Io)?;
    } else {
        writeln!(
            out,
            "{}\t{}\t{}\t{:.3}",
            fmt_g5(left),
            fmt_g5(right),
            fmt_g5(two),
            ratio
        )
        .map_err(RsomicsError::Io)?;
    }

    Ok(())
}
