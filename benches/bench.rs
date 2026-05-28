use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bed_utils(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-utils");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bed = manifest.join("tests/golden/b.bed");
    c.bench_function("rsomics-bed-utils sort golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["sort", bed.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bed_utils);
criterion_main!(benches);
