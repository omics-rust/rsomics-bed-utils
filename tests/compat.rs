use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bed-utils"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bedtools_available() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn run(mut cmd: Command) -> Vec<u8> {
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn sorted_lines(bytes: &[u8]) -> Vec<&[u8]> {
    let mut v: Vec<&[u8]> = bytes
        .split(|&b| b == b'\n')
        .filter(|l| !l.is_empty())
        .collect();
    v.sort_unstable();
    v
}

// `sort` must be byte-identical to `bedtools sort` (stable tie-break, all columns).
#[test]
fn sort_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let small = golden("small.bed");
    let ours = run({
        let mut c = bin();
        c.arg("sort").arg("-i").arg(&small);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("sort").arg("-i").arg(&small);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs)
    );
}

// `merge` (on sorted input) must be byte-identical to `bedtools merge`.
#[test]
fn merge_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let small = golden("small.bed");
    let sorted = run({
        let mut c = Command::new("bedtools");
        c.arg("sort").arg("-i").arg(&small);
        c
    });
    let dir = std::env::temp_dir().join("rsomics-bed-merge-compat");
    std::fs::create_dir_all(&dir).unwrap();
    let sorted_path = dir.join("sorted.bed");
    std::fs::write(&sorted_path, &sorted).unwrap();

    let ours = run({
        let mut c = bin();
        c.arg("merge").arg("-i").arg(&sorted_path);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("merge").arg("-i").arg(&sorted_path);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs)
    );
}

// `intersect` must produce the same overlap SET as `bedtools intersect` (A's
// columns preserved). Per-A multi-overlap order is implementation-defined, so
// compare as sorted sets, not byte-for-byte.
#[test]
fn intersect_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("small.bed");
    let b = golden("b.bed");
    let ours = run({
        let mut c = bin();
        c.arg("intersect").arg("-a").arg(&a).arg("-b").arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("intersect").arg("-a").arg(&a).arg("-b").arg(&b);
        c
    });
    assert_eq!(sorted_lines(&ours), sorted_lines(&theirs));
}

// `cluster` must be byte-identical to `bedtools cluster` (default d=0 and d=20).
#[test]
fn cluster_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let input = golden("cluster.bed");

    // default distance (0)
    let ours = run({
        let mut c = bin();
        c.arg("cluster").arg("-i").arg(&input);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("cluster").arg("-i").arg(&input);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "cluster d=0 mismatch"
    );

    // with -d 20: geneC (80-100) is 30 bp from geneB (20-50), not within 20 → own cluster;
    // geneD and geneE on chr2 are 20 bp apart (60-40=20) → same cluster.
    let ours_d20 = run({
        let mut c = bin();
        c.arg("cluster").arg("-i").arg(&input).arg("-d").arg("20");
        c
    });
    let theirs_d20 = run({
        let mut c = Command::new("bedtools");
        c.arg("cluster").arg("-i").arg(&input).arg("-d").arg("20");
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours_d20),
        String::from_utf8_lossy(&theirs_d20),
        "cluster d=20 mismatch"
    );
}

// `groupby` must be byte-identical to `bedtools groupby` for the core ops.
// bedtools formats integer-valued results without a decimal and non-integer
// results with 10 significant figures (%.10g).
#[test]
fn groupby_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let input = golden("groupby.bed");

    let ops = [
        "sum", "mean", "min", "max", "count", "collapse", "distinct", "first", "last",
    ];
    for op in ops {
        let ours = run({
            let mut c = bin();
            c.arg("groupby")
                .arg("-i")
                .arg(&input)
                .arg("-g")
                .arg("1")
                .arg("-c")
                .arg("5")
                .arg("--op")
                .arg(op);
            c
        });
        let theirs = run({
            let mut c = Command::new("bedtools");
            c.arg("groupby")
                .arg("-i")
                .arg(&input)
                .arg("-g")
                .arg("1")
                .arg("-c")
                .arg("5")
                .arg("-o")
                .arg(op);
            c
        });
        assert_eq!(
            String::from_utf8_lossy(&ours),
            String::from_utf8_lossy(&theirs),
            "groupby op={op} mismatch"
        );
    }
}

// `multiinter` must be byte-identical to `bedtools multiinter`.
// Output columns: chrom start end count list [0|1 per file].
#[test]
fn multiinter_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("multiinter_a.bed");
    let b = golden("multiinter_b.bed");

    let ours = run({
        let mut c = bin();
        c.arg("multiinter").arg("-i").arg(&a).arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("multiinter").arg("-i").arg(&a).arg(&b);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "multiinter mismatch"
    );
}

// `multiinter` with intra-file overlapping intervals must produce the same
// maximal segments as `bedtools multiinter` (2-file overlapping case).
#[test]
fn multiinter_overlap_2file_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("multiinter_overlap_a.bed");
    let b = golden("multiinter_overlap_b.bed");

    let ours = run({
        let mut c = bin();
        c.arg("multiinter").arg("-i").arg(&a).arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("multiinter").arg("-i").arg(&a).arg(&b);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "multiinter overlapping 2-file mismatch"
    );
}

// `multiinter` with intra-file overlapping intervals, 3-file case.
#[test]
fn multiinter_overlap_3file_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("multiinter_overlap_a.bed");
    let b = golden("multiinter_overlap_b.bed");
    let c_file = golden("multiinter_overlap_c.bed");

    let ours = run({
        let mut c = bin();
        c.arg("multiinter").arg("-i").arg(&a).arg(&b).arg(&c_file);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("multiinter").arg("-i").arg(&a).arg(&b).arg(&c_file);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "multiinter overlapping 3-file mismatch"
    );
}

// `subtract` (gaps of A not covered by B) must be byte-identical to
// `bedtools subtract` — A-file order, gaps in coordinate order, columns kept.
#[test]
fn subtract_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("small.bed");
    let b = golden("b.bed");
    let ours = run({
        let mut c = bin();
        c.arg("subtract").arg("-a").arg(&a).arg("-b").arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("subtract").arg("-a").arg(&a).arg("-b").arg(&b);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs)
    );
}

// `fisher` must be byte-identical to `bedtools fisher`.
// Bedtools version at time of writing: v2.31.1.
// Number format: p-values use C's %.5g (5 sig figs, trailing zeros stripped);
// ratio uses %.3f; "inf" printed literally when either marginal is zero.
#[test]
fn fisher_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("fisher_a.bed");
    let b = golden("fisher_b.bed");
    let g = golden("fisher_genome.txt");

    let ours = run({
        let mut c = bin();
        c.arg("fisher")
            .arg("-a")
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .arg("-g")
            .arg(&g);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("fisher")
            .arg("-a")
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .arg("-g")
            .arg(&g);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "fisher mismatch"
    );
}

// `fisher` on inputs where all A intervals overlap B (high-overlap edge case).
#[test]
fn fisher_high_overlap_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("small.bed");
    let b = golden("b.bed");
    let g = golden("genome.txt");

    let ours = run({
        let mut c = bin();
        c.arg("fisher")
            .arg("-a")
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .arg("-g")
            .arg(&g);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("fisher")
            .arg("-a")
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .arg("-g")
            .arg(&g);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "fisher high-overlap mismatch"
    );
}

// `reldist` must be byte-identical to `bedtools reldist`.
// Output columns: reldist count total fraction.
// Bins use floor(reldist*100)/100 (two decimal places).
// A intervals without two adjacent B midpoints are silently skipped.
#[test]
fn reldist_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("reldist_a.bed");
    let b = golden("reldist_b.bed");

    let ours = run({
        let mut c = bin();
        c.arg("reldist").arg("-a").arg(&a).arg("-b").arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("reldist").arg("-a").arg(&a).arg("-b").arg(&b);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "reldist mismatch"
    );
}

// `unionbedg` 2-file case must be byte-identical to `bedtools unionbedg`.
#[test]
fn unionbedg_2file_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("unionbedg_a.bg");
    let b = golden("unionbedg_b.bg");

    let ours = run({
        let mut c = bin();
        c.arg("unionbedg").arg("-i").arg(&a).arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("unionbedg").arg("-i").arg(&a).arg(&b);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "unionbedg 2-file mismatch"
    );
}

// `unionbedg` 3-file case with overlapping intervals across files.
#[test]
fn unionbedg_3file_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("unionbedg_a.bg");
    let b = golden("unionbedg_b.bg");
    let c_file = golden("unionbedg_c.bg");

    let ours = run({
        let mut c = bin();
        c.arg("unionbedg").arg("-i").arg(&a).arg(&b).arg(&c_file);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("unionbedg").arg("-i").arg(&a).arg(&b).arg(&c_file);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "unionbedg 3-file mismatch"
    );
}

// `unionbedg` with -names and -header.
#[test]
fn unionbedg_header_names_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("unionbedg_a.bg");
    let b = golden("unionbedg_b.bg");

    let ours = run({
        let mut c = bin();
        c.arg("unionbedg")
            .arg("--header")
            .arg("--names")
            .arg("FileA")
            .arg("FileB")
            .arg("-i")
            .arg(&a)
            .arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("unionbedg")
            .arg("-header")
            .arg("-names")
            .arg("FileA")
            .arg("FileB")
            .arg("-i")
            .arg(&a)
            .arg(&b);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "unionbedg header+names mismatch"
    );
}

// `maskfasta` (default N-masking) must be BYTE-IDENTICAL to `bedtools maskfasta`.
// Bedtools version at time of writing: v2.31.1.
#[test]
fn maskfasta_n_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let fasta = golden("mask_ref.fa");
    let bed = golden("mask_regions.bed");
    let tmpdir = std::env::temp_dir().join("rsomics-maskfasta-compat");
    std::fs::create_dir_all(&tmpdir).unwrap();
    let theirs_out = tmpdir.join("bedtools_masked.fa");

    let bt_status = Command::new("bedtools")
        .arg("maskfasta")
        .arg("-fi")
        .arg(&fasta)
        .arg("-bed")
        .arg(&bed)
        .arg("-fo")
        .arg(&theirs_out)
        .status()
        .unwrap();
    assert!(bt_status.success(), "bedtools maskfasta failed");
    let theirs = std::fs::read(&theirs_out).unwrap();

    let ours = run({
        let mut c = bin();
        c.arg("maskfasta")
            .arg("--fasta")
            .arg(&fasta)
            .arg("--bed")
            .arg(&bed);
        c
    });

    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "maskfasta N-mask mismatch"
    );
}

// `maskfasta --soft` must be BYTE-IDENTICAL to `bedtools maskfasta -soft`.
#[test]
fn maskfasta_soft_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let fasta = golden("mask_ref.fa");
    let bed = golden("mask_regions.bed");
    let tmpdir = std::env::temp_dir().join("rsomics-maskfasta-soft-compat");
    std::fs::create_dir_all(&tmpdir).unwrap();
    let theirs_out = tmpdir.join("bedtools_soft.fa");

    let bt_status = Command::new("bedtools")
        .arg("maskfasta")
        .arg("-fi")
        .arg(&fasta)
        .arg("-bed")
        .arg(&bed)
        .arg("-fo")
        .arg(&theirs_out)
        .arg("-soft")
        .status()
        .unwrap();
    assert!(bt_status.success(), "bedtools maskfasta -soft failed");
    let theirs = std::fs::read(&theirs_out).unwrap();

    let ours = run({
        let mut c = bin();
        c.arg("maskfasta")
            .arg("--fasta")
            .arg(&fasta)
            .arg("--bed")
            .arg(&bed)
            .arg("--soft");
        c
    });

    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "maskfasta soft-mask mismatch"
    );
}

// `shuffle` invariant tests (RNG differs from bedtools, so no byte-diff):
// - Same number of output intervals as input.
// - Each output interval has the same length as its input.
// - Each output interval lands within its chrom's bounds.
#[test]
fn shuffle_invariants() {
    let input = golden("shuffle_in.bed");
    let genome = golden("shuffle_genome.txt");

    if !input.exists() || !genome.exists() {
        eprintln!("skipping: shuffle fixtures not found");
        return;
    }

    let ours = run({
        let mut c = bin();
        c.arg("shuffle")
            .arg("-i")
            .arg(&input)
            .arg("-g")
            .arg(&genome)
            .arg("--seed")
            .arg("1");
        c
    });

    // Parse input intervals.
    let input_bytes = std::fs::read(&input).unwrap();
    let input_ivs: Vec<(String, u64, u64)> = String::from_utf8_lossy(&input_bytes)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            let mut f = l.splitn(4, '\t');
            let c = f.next().unwrap().to_string();
            let s: u64 = f.next().unwrap().parse().unwrap();
            let e: u64 = f.next().unwrap().parse().unwrap();
            (c, s, e)
        })
        .collect();

    // Parse genome.
    let genome_bytes = std::fs::read(&genome).unwrap();
    let genome_map: std::collections::HashMap<String, u64> = String::from_utf8_lossy(&genome_bytes)
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let mut f = l.splitn(3, '\t');
            let c = f.next().unwrap().to_string();
            let len: u64 = f.next().unwrap().parse().unwrap();
            (c, len)
        })
        .collect();

    // Parse output intervals.
    let out_ivs: Vec<(String, u64, u64)> = String::from_utf8_lossy(&ours)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            let mut f = l.splitn(4, '\t');
            let c = f.next().unwrap().to_string();
            let s: u64 = f.next().unwrap().parse().unwrap();
            let e: u64 = f.next().unwrap().parse().unwrap();
            (c, s, e)
        })
        .collect();

    assert_eq!(
        out_ivs.len(),
        input_ivs.len(),
        "shuffle: output interval count mismatch"
    );

    for (i, ((ic, is, ie), (oc, os, oe))) in input_ivs.iter().zip(out_ivs.iter()).enumerate() {
        let in_len = ie - is;
        let out_len = oe - os;
        assert_eq!(
            in_len, out_len,
            "interval {i}: length changed {in_len} -> {out_len}"
        );

        let chrom_len = *genome_map
            .get(oc)
            .unwrap_or_else(|| panic!("interval {i}: output chrom '{oc}' not in genome"));
        assert!(
            *oe <= chrom_len,
            "interval {i}: output {oc}:{os}-{oe} exceeds chrom len {chrom_len}"
        );
        let _ = ic; // input chrom not checked here (no -chrom flag)
    }
}

// `shuffle --excl`: output intervals must not overlap excluded regions.
#[test]
fn shuffle_excl_invariant() {
    let input = golden("shuffle_in.bed");
    let genome = golden("shuffle_genome.txt");
    let excl = golden("shuffle_excl.bed");

    if !input.exists() || !genome.exists() || !excl.exists() {
        eprintln!("skipping: shuffle_excl fixtures not found");
        return;
    }

    let ours = run({
        let mut c = bin();
        c.arg("shuffle")
            .arg("-i")
            .arg(&input)
            .arg("-g")
            .arg(&genome)
            .arg("--excl")
            .arg(&excl)
            .arg("--seed")
            .arg("42");
        c
    });

    // Parse excluded intervals.
    let excl_bytes = std::fs::read(&excl).unwrap();
    let excl_map: std::collections::HashMap<String, Vec<(u64, u64)>> = {
        let mut m: std::collections::HashMap<String, Vec<(u64, u64)>> =
            std::collections::HashMap::new();
        for line in String::from_utf8_lossy(&excl_bytes).lines() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut f = line.splitn(4, '\t');
            let c = f.next().unwrap().to_string();
            let s: u64 = f.next().unwrap().parse().unwrap();
            let e: u64 = f.next().unwrap().parse().unwrap();
            m.entry(c).or_default().push((s, e));
        }
        m
    };

    for line in String::from_utf8_lossy(&ours).lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut f = line.splitn(4, '\t');
        let chrom = f.next().unwrap().to_string();
        let start: u64 = f.next().unwrap().parse().unwrap();
        let end: u64 = f.next().unwrap().parse().unwrap();

        if let Some(ivs) = excl_map.get(&chrom) {
            for &(es, ee) in ivs {
                assert!(
                    end <= es || start >= ee,
                    "shuffle --excl: output {chrom}:{start}-{end} overlaps excluded {chrom}:{es}-{ee}"
                );
            }
        }
    }
}

// `shuffle --chrom`: each output interval must stay on the same chromosome.
#[test]
fn shuffle_chrom_invariant() {
    let input = golden("shuffle_in.bed");
    let genome = golden("shuffle_genome.txt");

    if !input.exists() || !genome.exists() {
        eprintln!("skipping: shuffle fixtures not found");
        return;
    }

    let ours = run({
        let mut c = bin();
        c.arg("shuffle")
            .arg("-i")
            .arg(&input)
            .arg("-g")
            .arg(&genome)
            .arg("--chrom")
            .arg("--seed")
            .arg("7");
        c
    });

    let input_bytes = std::fs::read(&input).unwrap();
    let input_chroms: Vec<String> = String::from_utf8_lossy(&input_bytes)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();

    let out_chroms: Vec<String> = String::from_utf8_lossy(&ours)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();

    assert_eq!(
        input_chroms.len(),
        out_chroms.len(),
        "shuffle --chrom: count mismatch"
    );

    for (i, (ic, oc)) in input_chroms.iter().zip(out_chroms.iter()).enumerate() {
        assert_eq!(
            ic, oc,
            "interval {i}: --chrom violated: input chrom '{ic}' != output chrom '{oc}'"
        );
    }
}

// `jaccard` must be byte-identical to `bedtools jaccard` on a fixture with
// self-overlapping intervals (which must be merged before computing statistics).
#[test]
fn jaccard_overlapping_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("jaccard_a.bed");
    let b = golden("jaccard_b.bed");

    let ours = run({
        let mut c = bin();
        c.arg("jaccard").arg(&a).arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("jaccard").arg("-a").arg(&a).arg("-b").arg(&b);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "jaccard overlapping mismatch"
    );
}

// `closest` must be byte-identical to `bedtools closest` including:
//   - Ties at the same distance (emit all B at minimum distance).
//   - Overlapping-B priority over adjacent-B (both at numeric distance 0).
//   - No-chrom-match emits ".\t-1\t-1" without a distance column.
#[test]
fn closest_overlapping_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("closest_a.bed");
    let b = golden("closest_b.bed");

    let ours = run({
        let mut c = bin();
        c.arg("closest").arg(&a).arg(&b);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("closest").arg("-a").arg(&a).arg("-b").arg(&b);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "closest overlapping mismatch"
    );
}

// `map` number formatting must match bedtools' `%.10g` convention:
// integer-valued floats print without a decimal point; others strip trailing zeros.
#[test]
fn map_format_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let a = golden("map_a.bed");
    let b = golden("map_b.bed");

    for op in ["mean", "sum", "min", "max"] {
        let ours = run({
            let mut c = bin();
            c.arg("map")
                .arg(&a)
                .arg(&b)
                .arg("--operation")
                .arg(op)
                .arg("-c")
                .arg("4");
            c
        });
        let theirs = run({
            let mut c = Command::new("bedtools");
            c.arg("map")
                .arg("-a")
                .arg(&a)
                .arg("-b")
                .arg(&b)
                .arg("-o")
                .arg(op)
                .arg("-c")
                .arg("4");
            c
        });
        assert_eq!(
            String::from_utf8_lossy(&ours),
            String::from_utf8_lossy(&theirs),
            "map op={op} format mismatch"
        );
    }
}

// `genomecov` must be byte-identical to `bedtools genomecov` on an input with
// self-overlapping intervals (requiring per-base depth accounting).
#[test]
fn genomecov_overlapping_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    // small.bed has chr1:10-30 and chr1:20-50 (overlapping), so depth 2 exists.
    let input = golden("small.bed");
    let genome = golden("genome.txt");

    let ours = run({
        let mut c = bin();
        c.arg("genomecov").arg(&input).arg("-g").arg(&genome);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("genomecov")
            .arg("-i")
            .arg(&input)
            .arg("-g")
            .arg(&genome);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "genomecov overlapping mismatch"
    );
}

// `nuc` header and data columns must be byte-identical to `bedtools nuc`.
#[test]
fn nuc_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }

    // Build .fai for the reference (bedtools nuc reads the FASTA directly).
    let fasta = golden("nuc_ref.fa");
    let bed = golden("nuc_regions.bed");

    let ours = run({
        let mut c = bin();
        c.arg("nuc").arg(&bed).arg("-f").arg(&fasta);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("nuc").arg("-bed").arg(&bed).arg("-fi").arg(&fasta);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "nuc mismatch"
    );
}

// `getfasta` must extract correct sequences (0-based BED coords).
#[test]
fn getfasta_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skipping: bedtools not found");
        return;
    }
    let fasta = golden("nuc_ref.fa");
    let bed = golden("nuc_regions.bed");

    // bedtools getfasta needs an .fai; build it if absent.
    let fai = fasta.with_extension("fa.fai");
    if !fai.exists() {
        Command::new("samtools")
            .arg("faidx")
            .arg(&fasta)
            .status()
            .ok();
    }

    let ours = run({
        let mut c = bin();
        c.arg("getfasta").arg(&bed).arg("-f").arg(&fasta);
        c
    });
    let theirs = run({
        let mut c = Command::new("bedtools");
        c.arg("getfasta")
            .arg("-fi")
            .arg(&fasta)
            .arg("-bed")
            .arg(&bed);
        c
    });
    assert_eq!(
        String::from_utf8_lossy(&ours),
        String::from_utf8_lossy(&theirs),
        "getfasta mismatch"
    );
}
