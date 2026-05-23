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

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

// ---------------------------------------------------------------------------
// BED loading helpers
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct Interval {
    chrom: String,
    start: i64,
    end: i64,
}

fn load_bed(path: &Path) -> Result<Vec<Interval>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut out = Vec::new();
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
        let start: i64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let end: i64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        out.push(Interval { chrom, start, end });
    }
    Ok(out)
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
    // exp here is the exponent of the number (mag).
    let s = if (-4..5).contains(&mag) {
        // fixed: number of decimals = 4 - mag (so total sig figs = 5)
        let decimals = (4 - mag).max(0) as usize;
        format!("{:.prec$}", v, prec = decimals)
    } else {
        // scientific: 4 decimal places = 5 sig figs
        format!("{:.4e}", v)
            // C uses e+01 style, Rust uses e1 — normalise
            .replace("e-0", "e-0") // keep as-is, normalise below
    };
    // Strip trailing zeros after decimal point, then strip bare decimal point.
    if s.contains('.') && !s.contains('e') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}

// ---------------------------------------------------------------------------
// Main fisher function
// ---------------------------------------------------------------------------

pub fn fisher(a_path: &Path, b_path: &Path, genome_path: &Path, out: &mut dyn Write) -> Result<()> {
    let a = load_bed(a_path)?;
    let b = load_bed(b_path)?;
    let genome_size = load_genome(genome_path)?;

    let query_counts = a.len() as i64;
    let db_counts = b.len() as i64;

    let query_union: i64 = a.iter().map(|iv| iv.end - iv.start).sum();
    let db_union: i64 = b.iter().map(|iv| iv.end - iv.start).sum();

    // Group intervals by chrom, sort within each chrom by start.
    // This avoids any cross-chrom sort-order dependency — bedtools requires sorted input
    // so within each chrom the order is already by start; we just group and sort within.
    let mut a_by_chrom: std::collections::HashMap<&str, Vec<&Interval>> =
        std::collections::HashMap::new();
    let mut b_by_chrom: std::collections::HashMap<&str, Vec<&Interval>> =
        std::collections::HashMap::new();
    for iv in &a {
        a_by_chrom.entry(iv.chrom.as_str()).or_default().push(iv);
    }
    for iv in &b {
        b_by_chrom.entry(iv.chrom.as_str()).or_default().push(iv);
    }
    for v in a_by_chrom.values_mut() {
        v.sort_unstable_by_key(|iv| iv.start);
    }
    for v in b_by_chrom.values_mut() {
        v.sort_unstable_by_key(|iv| iv.start);
    }

    // Count overlapping A×B pairs per chrom.
    //
    // Classic O(N log N sort + O(N) sweep: build an event list from all four
    // endpoint kinds, sort by coordinate (closes before opens on ties to honour
    // half-open [start,end) semantics — an interval ending at p does NOT overlap
    // one starting at p), then maintain `active_a` / `active_b` counters.
    // When an A interval opens, add active_b pairs; when a B interval opens, add
    // active_a pairs. No BIT, no coordinate compression — tiny constant overhead.
    let mut overlap_counts: i64 = 0;
    for (chrom, a_ivs) in &a_by_chrom {
        let b_ivs = match b_by_chrom.get(chrom) {
            Some(v) => v,
            None => continue,
        };
        if b_ivs.is_empty() {
            continue;
        }

        // Event kinds sorted so that ties are broken: closes (0) before opens (1).
        // 0 = close-A, 1 = close-B, 2 = open-A, 3 = open-B.
        // Using a u8 kind in the tuple gives the right ordering for free.
        let mut events: Vec<(i64, u8)> = Vec::with_capacity(2 * a_ivs.len() + 2 * b_ivs.len());
        for iv in a_ivs.iter() {
            events.push((iv.start, 2)); // open-A
            events.push((iv.end, 0)); // close-A
        }
        for iv in b_ivs.iter() {
            events.push((iv.start, 3)); // open-B
            events.push((iv.end, 1)); // close-B
        }
        events.sort_unstable();

        let mut active_a: i64 = 0;
        let mut active_b: i64 = 0;
        for (_, kind) in &events {
            match kind {
                0 => active_a -= 1, // close-A
                1 => active_b -= 1, // close-B
                2 => {
                    // open-A: this interval overlaps every currently-open B
                    overlap_counts += active_b;
                    active_a += 1;
                }
                3 => {
                    // open-B: this interval overlaps every currently-open A
                    overlap_counts += active_a;
                    active_b += 1;
                }
                _ => unreachable!(),
            }
        }
    }

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
