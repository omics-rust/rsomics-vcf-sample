use std::io::Write;
use std::path::Path;

use rsomics_vcf_sample::{SampleMode, sample_vcf};

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn data_lines(s: &str) -> Vec<&str> {
    s.lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect()
}

fn header_lines(s: &str) -> Vec<&str> {
    s.lines().filter(|l| l.starts_with('#')).collect()
}

#[test]
fn header_always_emitted_fraction() {
    let input = golden("input.vcf");
    let mut out = Vec::new();
    sample_vcf(&input, &mut out, SampleMode::Fraction(1.0), 42).unwrap();
    let result = String::from_utf8(out).unwrap();
    let headers = header_lines(&result);
    assert_eq!(headers.len(), 5, "expected 5 header lines: {result}");
    assert!(
        headers.last().unwrap().starts_with("#CHROM"),
        "last header must be #CHROM"
    );
}

#[test]
fn fraction_one_keeps_all() {
    let input = golden("input.vcf");
    let mut out = Vec::new();
    sample_vcf(&input, &mut out, SampleMode::Fraction(1.0), 42).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines = data_lines(&result);
    assert_eq!(
        lines.len(),
        10,
        "p=1.0 should keep all 10 records: {result}"
    );
}

#[test]
fn fraction_zero_keeps_none() {
    let input = golden("input.vcf");
    let mut out = Vec::new();
    sample_vcf(&input, &mut out, SampleMode::Fraction(0.0), 42).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines = data_lines(&result);
    assert_eq!(lines.len(), 0, "p=0.0 should keep no records: {result}");
}

#[test]
fn exact_n_fewer_than_total() {
    let input = golden("input.vcf");
    let mut out = Vec::new();
    sample_vcf(&input, &mut out, SampleMode::Exact(5), 42).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines = data_lines(&result);
    assert_eq!(
        lines.len(),
        5,
        "reservoir(5) should yield exactly 5 records: {result}"
    );
}

#[test]
fn exact_n_more_than_total_clamps() {
    let input = golden("input.vcf");
    let mut out = Vec::new();
    sample_vcf(&input, &mut out, SampleMode::Exact(100), 42).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines = data_lines(&result);
    assert_eq!(
        lines.len(),
        10,
        "reservoir(100) should yield all 10 records: {result}"
    );
}

#[test]
fn exact_preserves_input_order() {
    let input = golden("input.vcf");
    let mut out = Vec::new();
    sample_vcf(&input, &mut out, SampleMode::Exact(5), 42).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines = data_lines(&result);
    // Output must appear in the same relative order as the input records.
    // We verify this by checking that the selected records are a subsequence
    // of the input, preserving the original record order (not just POS order,
    // since multiple chromosomes have independent coordinate spaces).
    let input_lines: Vec<String> = std::fs::read_to_string(&input)
        .unwrap()
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_owned)
        .collect();
    let mut last_pos = 0usize;
    for line in &lines {
        let pos = input_lines.iter().position(|l| l == line).unwrap();
        assert!(
            pos >= last_pos,
            "output record at input index {pos} precedes previous index {last_pos}: {line}"
        );
        last_pos = pos + 1;
    }
}

#[test]
fn reproducible_with_same_seed() {
    let input = golden("input.vcf");
    let mut out1 = Vec::new();
    let mut out2 = Vec::new();
    sample_vcf(&input, &mut out1, SampleMode::Fraction(0.5), 99).unwrap();
    sample_vcf(&input, &mut out2, SampleMode::Fraction(0.5), 99).unwrap();
    assert_eq!(out1, out2, "same seed must produce identical output");
}

#[test]
fn empty_input_produces_only_header() {
    use tempfile::NamedTempFile;
    let mut f = NamedTempFile::new().unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    let mut out = Vec::new();
    sample_vcf(f.path(), &mut out, SampleMode::Exact(10), 42).unwrap();
    let result = String::from_utf8(out).unwrap();
    let data = data_lines(&result);
    assert_eq!(data.len(), 0, "empty input → no data records: {result}");
    let headers = header_lines(&result);
    assert_eq!(headers.len(), 2, "header must still emit: {result}");
}
