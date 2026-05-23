use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone)]
struct Interval {
    chrom: String,
    start: u64,
    end: u64,
    line: String,
}

fn load_sorted_bed(path: &Path) -> Result<Vec<Interval>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    for line in reader.lines() {
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
        intervals.push(Interval {
            chrom,
            start,
            end,
            line: line.clone(),
        });
    }
    Ok(intervals)
}

/// Per-chromosome sorted B list plus a prefix-max-end array.
///
/// `prefix_max_end[i]` = max(B[0].end, …, B[i].end).
/// Non-decreasing; enables O(1) pruning of left scans.
struct ChromB<'a> {
    bs: Vec<&'a Interval>,
    /// `prefix_max_end[i]` = max of `bs[0].end .. bs[i].end`. Same length as bs.
    prefix_max_end: Vec<u64>,
}

impl<'a> ChromB<'a> {
    fn new(mut bs: Vec<&'a Interval>) -> Self {
        // Stable sort by start: preserves input-file order within same-start ties,
        // matching the output order bedtools produces.
        bs.sort_by_key(|b| b.start);
        let mut prefix_max_end = Vec::with_capacity(bs.len());
        let mut running_max = 0u64;
        for b in &bs {
            running_max = running_max.max(b.end);
            prefix_max_end.push(running_max);
        }
        Self { bs, prefix_max_end }
    }
}

/// Find closest B feature(s) for each A interval, matching `bedtools closest`.
///
/// Priority rule (mirrors bedtools):
///   1. If any B strictly overlaps A, emit ALL strictly-overlapping B.
///   2. Otherwise emit all B at the minimum non-overlapping gap distance.
///
/// No distance column is appended (bedtools default omits it).
/// When no B exists on the same chromosome: emits `A\t.\t-1\t-1`.
///
/// Algorithm: O(log m + k) per query using binary search + prefix-max-end pruning.
pub fn closest(a_path: &Path, b_path: &Path, output: &mut dyn Write) -> Result<u64> {
    let a_intervals = load_sorted_bed(a_path)?;
    let b_intervals = load_sorted_bed(b_path)?;
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    // Group B by chromosome; build ChromB (sorted + prefix_max_end) per chrom.
    let mut b_raw: HashMap<&str, Vec<&Interval>> = HashMap::new();
    for b in &b_intervals {
        b_raw.entry(b.chrom.as_str()).or_default().push(b);
    }
    let b_by_chrom: HashMap<&str, ChromB<'_>> = b_raw
        .into_iter()
        .map(|(c, v)| (c, ChromB::new(v)))
        .collect();

    // Candidate indices in bs; reused across queries to avoid per-query allocation.
    let mut cand_idx: Vec<usize> = Vec::new();

    for a in &a_intervals {
        let a_fields = a.line.as_str();
        let Some(cb) = b_by_chrom.get(a.chrom.as_str()) else {
            writeln!(out, "{a_fields}\t.\t-1\t-1").map_err(RsomicsError::Io)?;
            count += 1;
            continue;
        };

        let bs = &cb.bs;
        let pmax = &cb.prefix_max_end;

        cand_idx.clear();

        // lb = first index where B.start >= a.start.
        let lb = bs.partition_point(|b| b.start < a.start);

        // ── Overlap search ─────────────────────────────────────────────────
        //
        // RIGHT of lb: B[lb..ri_end) with start in [a.start, a.end).
        // These automatically have end > start >= a.start → all overlap A.
        let ri_end = bs.partition_point(|b| b.start < a.end);
        for idx in lb..ri_end {
            cand_idx.push(idx);
        }

        // LEFT of lb: B[..lb] with start < a.start; overlap iff end > a.start.
        // Scan leftward from lb-1; prune when pmax[li] <= a.start
        // (prefix_max_end is non-decreasing, so no further-left B can have end > a.start).
        if lb > 0 {
            let mut li = lb;
            loop {
                li -= 1;
                if pmax[li] <= a.start {
                    break;
                }
                if bs[li].end > a.start {
                    cand_idx.push(li);
                }
                if li == 0 {
                    break;
                }
            }
        }

        if !cand_idx.is_empty() {
            // Emit in ascending index order (= bs sort order = start-sorted,
            // original-input-order tiebreaking within same start).
            cand_idx.sort_unstable();
            for &idx in &cand_idx {
                writeln!(out, "{a_fields}\t{}", bs[idx].line).map_err(RsomicsError::Io)?;
                count += 1;
            }
            continue;
        }

        // ── Non-overlap: find minimum gap, collect all ties ────────────────
        //
        // RIGHT side [ri_end..]: starts >= a.end; gap = start - a.end (non-decreasing).
        // Scan right until start > a.end + min_gap (gap increases monotonically).
        //
        // LEFT side [..lb]: starts < a.start, ends <= a.start; gap = a.start - end.
        // Scan leftward from lb-1; prune via prefix_max_end:
        //   once pmax[li] < a.start - min_gap, all gaps in [0..=li] exceed min_gap.
        // (pmax[li] < a.start - min_gap means max end in prefix < a.start - min_gap,
        //  so a.start - any_end > min_gap.)

        let mut min_gap: Option<u64> = None;

        // Right scan.
        let mut ri = ri_end;
        while ri < bs.len() {
            let g = bs[ri].start - a.end;
            match min_gap {
                None => {
                    min_gap = Some(g);
                    cand_idx.push(ri);
                }
                Some(mg) if g == mg => {
                    cand_idx.push(ri);
                }
                Some(_) => break, // starts non-decreasing → farther right, done
            }
            ri += 1;
        }

        // Left scan.
        if lb > 0 {
            let mut li = lb;
            loop {
                li -= 1;
                // Prune: if pmax[li] < a.start - min_gap, all left-side gaps > min_gap.
                if let Some(mg) = min_gap
                    && pmax[li] < a.start.saturating_sub(mg)
                {
                    break;
                }
                debug_assert!(bs[li].end <= a.start, "left-gap B overlaps A");
                let g = a.start - bs[li].end;
                match min_gap {
                    None => {
                        min_gap = Some(g);
                        cand_idx.push(li);
                    }
                    Some(mg) if g < mg => {
                        // Strictly better: discard all previously collected candidates.
                        cand_idx.clear();
                        min_gap = Some(g);
                        cand_idx.push(li);
                    }
                    Some(mg) if g == mg => {
                        cand_idx.push(li);
                    }
                    Some(_) => {} // g > mg; keep scanning (end not monotone, pmax prune handles it)
                }
                if li == 0 {
                    break;
                }
            }
        }

        // Emit in ascending index order.
        if !cand_idx.is_empty() {
            cand_idx.sort_unstable();
            for &idx in &cand_idx {
                writeln!(out, "{a_fields}\t{}", bs[idx].line).map_err(RsomicsError::Io)?;
                count += 1;
            }
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
