use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_RECORDS: usize = 100_000;
const SEED: u64 = 0x00CAFE_0042;

fn xorshift(x: &mut u64) -> u64 {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    *x
}

fn synth_vcf(path: &PathBuf) {
    let f = File::create(path).expect("create vcf");
    let mut w = BufWriter::new(f);
    writeln!(w, "##fileformat=VCFv4.2").unwrap();
    writeln!(w, "##contig=<ID=chr1,length=248956422>").unwrap();
    writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    let mut rng = SEED;
    let mut pos = 1u64;
    let bases = b"ACGT";
    for i in 0..N_RECORDS {
        pos += 1 + (xorshift(&mut rng) % 1000);
        let r = bases[(xorshift(&mut rng) % 4) as usize] as char;
        let a = bases[(xorshift(&mut rng) % 4) as usize] as char;
        writeln!(w, "chr1\t{pos}\trs{i}\t{r}\t{a}\t.\tPASS\t.").unwrap();
    }
}

fn ensure_fixture() -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("rsomics-vcf-sample-bench-{N_RECORDS}.vcf"));
    if !p.exists() {
        synth_vcf(&p);
    }
    p
}

fn bench(c: &mut Criterion) {
    let fixture = ensure_fixture();
    let ours = env!("CARGO_BIN_EXE_rsomics-vcf-sample");
    let mut group = c.benchmark_group(format!("vcf_sample/{N_RECORDS}"));
    group.sample_size(10);

    // Fraction mode: keep 50%.
    group.bench_function("rsomics-vcf-sample-frac50", |bm| {
        bm.iter(|| {
            let out = Command::new(ours)
                .arg(&fixture)
                .args(["-p", "0.5", "--seed", "42"])
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    // Exact mode: reservoir sample 10k.
    group.bench_function("rsomics-vcf-sample-exact10k", |bm| {
        bm.iter(|| {
            let out = Command::new(ours)
                .arg(&fixture)
                .args(["-n", "10000", "--seed", "42"])
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    // bcftools view -s for comparison (subset columns, not records — closest bcftools analog is view with prob).
    if Command::new("bcftools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
    {
        group.bench_function("bcftools-view-pass", |bm| {
            bm.iter(|| {
                // bcftools view -f PASS — selects by filter, not by fraction; used as baseline cost.
                let out = Command::new("bcftools")
                    .args(["view", "-f", "PASS"])
                    .arg(&fixture)
                    .output()
                    .expect("bcftools run");
                assert!(
                    out.status.success(),
                    "{}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
    } else {
        eprintln!("bcftools not on PATH — skipping upstream comparison");
    }

    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
