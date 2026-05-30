//! VCF variant subsampler — fraction-based or reservoir-sampled exact count.
//!
//! Equivalent to `bcftools view --min-ac 0 | shuf` for fraction sampling, but
//! deterministic (seeded) and single-pass (Vitter's Algorithm R for exact-N).
//! Header lines (##, #CHROM) are always emitted verbatim. Data records are
//! selected using one of two modes:
//!
//! - Fraction mode (`-p <f>`): each record is kept independently with
//!   probability `p`. Output size is approximately `p × total_records`.
//! - Exact mode (`-n <n>`): reservoir sampling over a single pass; exactly
//!   `min(n, total_records)` records are emitted, in input order (the
//!   reservoir is sorted by original index before output).
//!
//! Both modes accept a `--seed` for reproducibility. Without `--seed`, a
//! fixed default seed (42) is used so the tool is reproducible by default;
//! pass `--seed 0` to get a time-seeded run.

use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rand::Rng;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone, Copy)]
pub enum SampleMode {
    /// Keep each record with probability `p`.
    Fraction(f64),
    /// Reservoir sample exactly `n` records.
    Exact(usize),
}

/// Subsample VCF records from `input`, writing to `output`.
///
/// Header lines (starting with `#`) are always passed through. Data records
/// are selected according to `mode`. The `seed` controls reproducibility.
pub fn sample_vcf(
    input: &Path,
    output: &mut dyn Write,
    mode: SampleMode,
    seed: u64,
) -> Result<SampleStats> {
    let file = std::fs::File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    sample_vcf_from_reader(BufReader::new(file), output, mode, seed)
}

/// Same as [`sample_vcf`] but reads from stdin.
pub fn sample_vcf_stdin(
    output: &mut dyn Write,
    mode: SampleMode,
    seed: u64,
) -> Result<SampleStats> {
    sample_vcf_from_reader(BufReader::new(io::stdin()), output, mode, seed)
}

pub struct SampleStats {
    pub total: u64,
    pub kept: u64,
}

fn sample_vcf_from_reader<R: io::Read>(
    reader: BufReader<R>,
    output: &mut dyn Write,
    mode: SampleMode,
    seed: u64,
) -> Result<SampleStats> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut out = BufWriter::new(output);

    match mode {
        SampleMode::Fraction(p) => fraction_pass(reader, &mut out, &mut rng, p),
        SampleMode::Exact(n) => exact_pass(reader, &mut out, &mut rng, n),
    }
}

/// Single-pass fraction sampling: keep each record with probability `p`.
fn fraction_pass<R: io::Read>(
    reader: BufReader<R>,
    out: &mut BufWriter<impl Write>,
    rng: &mut SmallRng,
    p: f64,
) -> Result<SampleStats> {
    let mut total = 0u64;
    let mut kept = 0u64;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let bytes = line.as_bytes();
        if bytes.first() == Some(&b'#') {
            out.write_all(bytes).map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
        } else if !bytes.is_empty() {
            total += 1;
            if rng.random::<f64>() < p {
                kept += 1;
                out.write_all(bytes).map_err(RsomicsError::Io)?;
                out.write_all(b"\n").map_err(RsomicsError::Io)?;
            }
        }
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(SampleStats { total, kept })
}

/// Single-pass reservoir sampling (Vitter's Algorithm R): exactly `n` records.
/// Header lines are streamed immediately; data records are buffered and
/// emitted in input order after the full pass.
fn exact_pass<R: io::Read>(
    reader: BufReader<R>,
    out: &mut BufWriter<impl Write>,
    rng: &mut SmallRng,
    n: usize,
) -> Result<SampleStats> {
    let mut header: Vec<Vec<u8>> = Vec::new();
    let mut reservoir: Vec<(usize, Vec<u8>)> = Vec::with_capacity(n);
    let mut total = 0usize;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let bytes = line.into_bytes();
        if bytes.first() == Some(&b'#') {
            header.push(bytes);
        } else if !bytes.is_empty() {
            let idx = total;
            total += 1;
            if reservoir.len() < n {
                reservoir.push((idx, bytes));
            } else {
                // Replace a random entry with probability n / total.
                let j = rng.random_range(0..total);
                if j < n {
                    reservoir[j] = (idx, bytes);
                }
            }
        }
    }

    for h in &header {
        out.write_all(h).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }

    reservoir.sort_unstable_by_key(|(i, _)| *i);
    let kept = reservoir.len() as u64;
    for (_, bytes) in reservoir {
        out.write_all(&bytes).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(SampleStats {
        total: total as u64,
        kept,
    })
}
