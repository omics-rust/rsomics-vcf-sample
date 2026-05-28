use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_vcf_sample::{SampleMode, sample_vcf, sample_vcf_stdin};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-sample", disable_help_flag = true)]
pub struct Cli {
    /// Input VCF (default: stdin)
    input: Option<PathBuf>,
    /// Output VCF (default: stdout)
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
    /// Sampling fraction (0 < p ≤ 1); mutually exclusive with -n
    #[arg(short = 'p', long, conflicts_with = "count")]
    fraction: Option<f64>,
    /// Exact record count via reservoir sampling; mutually exclusive with -p
    #[arg(short = 'n', long, value_name = "INT")]
    count: Option<usize>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mode = match (self.fraction, self.count) {
            (Some(p), None) => {
                if !(0.0..=1.0).contains(&p) {
                    return Err(RsomicsError::InvalidInput(
                        "-p/--fraction must be in range (0, 1]".into(),
                    ));
                }
                SampleMode::Fraction(p)
            }
            (None, Some(n)) => SampleMode::Exact(n),
            (None, None) => {
                return Err(RsomicsError::InvalidInput(
                    "one of -p/--fraction or -n/--count is required".into(),
                ));
            }
            (Some(_), Some(_)) => unreachable!("clap conflicts_with enforces mutual exclusion"),
        };

        let seed = self.common.seed_rng();

        let mut stdout_lock;
        let mut file_out;
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = std::fs::File::create(p).map_err(RsomicsError::Io)?;
            &mut file_out
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };

        match self.input {
            Some(ref p) => sample_vcf(p.as_path(), out, mode, seed)?,
            None => sample_vcf_stdin(out, mode, seed)?,
        };
        Ok(())
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Random subsample VCF variants by fraction or exact count.",
    origin: Some(Origin {
        upstream: "bcftools view",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/gigascience/giab008"),
    }),
    usage_lines: &["[OPTIONS] [INPUT]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('p'),
                long: "fraction",
                aliases: &[],
                value: Some("<FLOAT>"),
                type_hint: Some("f64"),
                required: false,
                default: None,
                description: "Keep each record with probability p (0 < p ≤ 1)",
                why_default: None,
            },
            FlagSpec {
                short: Some('n'),
                long: "count",
                aliases: &[],
                value: Some("<INT>"),
                type_hint: Some("usize"),
                required: false,
                default: None,
                description: "Reservoir-sample exactly N records",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output VCF path",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Keep 10% of variants",
            command: "rsomics-vcf-sample -p 0.1 input.vcf",
        },
        Example {
            description: "Reservoir-sample exactly 1000 records",
            command: "rsomics-vcf-sample -n 1000 input.vcf",
        },
        Example {
            description: "Reproducible fraction sample with custom seed",
            command: "rsomics-vcf-sample -p 0.5 --seed 123 input.vcf",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
