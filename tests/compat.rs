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
