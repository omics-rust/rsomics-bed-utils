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
