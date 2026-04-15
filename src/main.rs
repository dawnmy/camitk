mod cami;
mod commands;
mod expression;
mod processing;
mod taxonomy;

use anyhow::{Result, bail};
use clap::{Parser, Subcommand};
use std::io;
use std::path::PathBuf;

use commands::{
    benchmark::{self as benchmark_cmd, BenchmarkConfig as BenchmarkRunConfig},
    convert::{self as convert_cmd, ConvertConfig as ConvertRunConfig},
    fillup::{self as fillup_cmd, FillupConfig},
    filter::{self as filter_cmd, FilterConfig},
    list::{self as list_cmd, ListConfig},
    preview::{self as preview_cmd, PreviewConfig},
    renorm::{self as renorm_cmd, RenormConfig},
    sort::{self as sort_cmd, SortConfig, SortMode, TaxPathField},
};

#[derive(Parser)]
#[command(
    author,
    version,
    about = "Explore and post-process CAMI profiling tables",
    long_about = "camitk reads CAMI-format abundance tables and provides commands to inspect, filter, and reshape them.",
    disable_help_subcommand = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

fn split_labels(input: &str) -> Vec<String> {
    input
        .split(',')
        .map(|part| part.trim().to_string())
        .filter(|part| !part.is_empty())
        .collect()
}

fn split_ranks(input: &str) -> Vec<String> {
    input
        .split(|c: char| c == ',' || c.is_whitespace())
        .map(|part| part.trim().to_string())
        .filter(|part| !part.is_empty())
        .collect()
}

#[derive(Subcommand)]
enum Commands {
    #[command(
        about = "Benchmark predicted profiles against a ground truth",
        long_about = "Compare predicted CAMI abundance tables to a ground truth profile per sample and rank. Reports TP/FP/FN counts, precision (purity), recall (completeness), F1-score, Jaccard index, L1 error, Bray-Curtis distance, Shannon diversity, equitability, Pearson and Spearman correlations, weighted and unweighted UniFrac differences, the Abundance Rank Error (ARE), and the mass-weighted Abundance Rank Error (mARE)."
    )]
    Benchmark {
        #[arg(
            short = 'g',
            long = "ground-truth",
            value_name = "FILE",
            help = "Ground truth CAMI profile to evaluate against."
        )]
        ground_truth: PathBuf,
        #[arg(
            value_name = "PREDICTED",
            num_args = 1..,
            help = "Predicted CAMI profiles to benchmark."
        )]
        predictions: Vec<PathBuf>,
        #[arg(
            short = 'u',
            long = "update-taxonomy",
            help = "Refresh taxonomy fields from the local NCBI taxdump when benchmarking."
        )]
        update_taxonomy: bool,
        #[arg(
            short = 'l',
            long = "labels",
            value_name = "LABELS",
            help = "Comma-separated labels for the predicted profiles."
        )]
        labels: Option<String>,
        #[arg(
            long = "af",
            value_name = "EXPR",
            help = "Filter expression applied to both ground truth and predicted profiles before scoring.",
            long_help = "Expression syntax matches `camitk filter`, combining rank (r), sample (s), abundance (a), taxonomy (t/tax), and cumulative-sum (c) predicates with & (and), | (or), and parentheses."
        )]
        all_filter: Option<String>,
        #[arg(
            long = "gf",
            value_name = "EXPR",
            help = "Filter expression applied to the ground truth profile before scoring.",
            long_help = "Expression syntax matches `camitk filter`, combining rank (r), sample (s), abundance (a), taxonomy (t/tax), and cumulative-sum (c) predicates with & (and), | (or), and parentheses."
        )]
        ground_filter: Option<String>,
        #[arg(
            long = "pf",
            value_name = "EXPR",
            help = "Filter expression applied to each predicted profile before scoring.",
            long_help = "Expression syntax matches `camitk filter`, combining rank (r), sample (s), abundance (a), taxonomy (t/tax), and cumulative-sum (c) predicates with & (and), | (or), and parentheses."
        )]
        pred_filter: Option<String>,
        #[arg(
            short = 'n',
            long = "normalize",
            help = "Normalize abundances within each sample/rank to 100 before scoring."
        )]
        normalize: bool,
        #[arg(
            long = "by-domain",
            help = "Also write domain-specific reports for Bacteria, Archaea, Eukarya, and Viruses."
        )]
        by_domain: bool,
        #[arg(
            long = "group-realms",
            action = clap::ArgAction::SetTrue,
            default_value_t = true,
            help = "Group viral realms under the Viruses superkingdom when scoring and filtering (default: enabled)."
        )]
        group_realms: bool,
        #[arg(
            long = "no-group-realms",
            action = clap::ArgAction::SetTrue,
            conflicts_with = "group_realms",
            help = "Disable viral realm grouping so realms are not folded into Viruses."
        )]
        no_group_realms: bool,
        #[arg(
            short = 'o',
            long = "output",
            value_name = "DIR",
            help = "Directory where benchmark TSV reports will be written."
        )]
        output: PathBuf,
        #[arg(
            short = 'r',
            long = "ranks",
            value_name = "RANKS",
            help = "Comma-separated list of ranks to evaluate (mix short and full names)."
        )]
        ranks: Option<String>,
        #[arg(
            long = "dmp-dir",
            value_name = "DIR",
            help = "Directory containing nodes.dmp and names.dmp (defaults to ~/.camitk)."
        )]
        dmp_dir: Option<PathBuf>,
    },
    #[command(
        about = "Filter CAMI profiling data with logical expressions",
        long_about = "Apply boolean filter expressions to CAMI profiling tables. Combine rank, sample, abundance, taxonomy, and cumulative-sum predicates with & (and), | (or), and parentheses. When --fill-up is provided the command completes missing lineages using the NCBI taxdump before writing a filtered CAMI table."
    )]
    Filter {
        #[arg(
            value_name = "EXPR",
            help = "Filter expression combining rank (r), sample (s), abundance (a), taxonomy (t/tax), and cumsum (c) tests.",
            long_help = concat!(
                "Filter expression combining rank (r), sample (s), abundance (a), taxonomy (t/tax), and cumulative-sum tests. Use & (and), | (or), and parentheses.\n\n",
                "Rank selectors: r==rank, r!=rank, r<=rank (current or more specific), r>=rank (current or more general), and range comparisons follow the order declared in @Ranks.\n",
                "Sample selectors: s==id accepts IDs, 1-based indices, comma-separated lists, or inclusive ranges (e.g. s==1:3); use '.' or ':' for all samples; s!= negates; s~'regex' matches sample IDs with a regular expression.\n",
                "Abundance selectors: a>=value, a<=value, a>value, a<value, a==value, and a!=value compare percentages.\n",
                "Taxonomy selectors: t==taxid, t!=taxid, t<=taxid (ancestor or self), t<taxid (strict ancestor); prefix with ! to invert.\n",
                "Cumulative sums: c<=threshold keeps the least-abundant taxa per rank whose cumulative percentage is within the threshold; prefix with ! to discard them instead."
            )
        )]
        expression: String,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
        output: Option<PathBuf>,
        #[arg(
            long,
            help = "Fill in missing higher ranks using the NCBI taxdump (downloaded to ~/.camitk if absent)."
        )]
        fill_up: bool,
        #[arg(
            long = "from",
            help = "Rank to aggregate from when filling (defaults to species when available)."
        )]
        from_rank: Option<String>,
        #[arg(
            long = "to",
            alias = "to-rank",
            default_value = "phylum",
            help = "Highest rank to build during fill-up (inclusive)."
        )]
        to_rank: String,
        #[arg(
            long,
            help = "Renormalize each rank to 100% after filtering and filling."
        )]
        renorm: bool,
        #[arg(
            value_name = "INPUT",
            index = 2,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
        #[arg(
            long = "dmp-dir",
            value_name = "DIR",
            help = "Directory containing nodes.dmp and names.dmp (defaults to ~/.camitk)."
        )]
        dmp_dir: Option<PathBuf>,
    },
    #[command(
        about = "List samples and per-rank summaries",
        long_about = "Report each sample's declared ranks, total taxa, and the summed abundance at every rank so you can spot coverage gaps before further processing."
    )]
    List {
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
    },
    #[command(
        about = "Preview the first entries per sample",
        long_about = "Show the first N taxonomic entries for each sample, preserving the CAMI header format so you can quickly verify parsing and taxonomy strings."
    )]
    Preview {
        #[arg(
            short = 'n',
            long,
            default_value_t = 5,
            help = "Number of entries per sample to display.",
            value_name = "N"
        )]
        n: usize,
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
    },
    #[command(
        about = "Renormalize abundances per rank",
        long_about = "Scale positive abundances within each rank of every sample so they sum to 100 without altering their precision."
    )]
    Renorm {
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
        output: Option<PathBuf>,
    },
    #[command(
        about = "Fill missing ranks using taxonomy",
        long_about = "Complete partial lineages in each sample by consulting the NCBI taxdump while preserving the numeric precision of existing abundances."
    )]
    Fillup {
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
        output: Option<PathBuf>,
        #[arg(
            long = "to",
            alias = "to-rank",
            help = "Highest rank to fill to (inclusive); defaults to the first declared rank."
        )]
        to_rank: Option<String>,
        #[arg(
            long = "from",
            help = "Rank to aggregate from when filling (defaults to species when available)."
        )]
        from_rank: Option<String>,
        #[arg(
            long = "dmp-dir",
            value_name = "DIR",
            help = "Directory containing nodes.dmp and names.dmp (defaults to ~/.camitk)."
        )]
        dmp_dir: Option<PathBuf>,
    },
    #[command(
        about = "Convert TSV profiling results into CAMI format",
        long_about = "Read taxonomic abundances from a TSV file and produce a CAMI-formatted profile with ranks completed from the NCBI taxdump."
    )]
    Convert {
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input TSV file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
        output: Option<PathBuf>,
        #[arg(
            short = 'i',
            long = "taxid-column",
            value_name = "INDEX",
            default_value_t = 1,
            help = "1-based column index containing NCBI taxids."
        )]
        taxid_column: usize,
        #[arg(
            short = 'a',
            long = "abundance-column",
            value_name = "INDEX",
            default_value_t = 2,
            help = "1-based column index containing abundances."
        )]
        abundance_column: usize,
        #[arg(
            long = "input-percent",
            help = "Treat abundances as percentages instead of fractions (0-1)."
        )]
        input_percent: bool,
        #[arg(
            long = "norm",
            help = "Normalize abundances to total 100 before converting."
        )]
        normalize: bool,
        #[arg(
            short = 's',
            long = "sample-id",
            value_name = "ID",
            default_value = "sample",
            help = "Sample ID to use in the generated CAMI profile."
        )]
        sample_id: String,
        #[arg(
            short = 'T',
            long = "taxonomy-tag",
            value_name = "TAG",
            help = "@TaxonomyID value describing the NCBI taxonomy snapshot (e.g. 2025-06-19)."
        )]
        taxonomy_tag: Option<String>,
        #[arg(
            long = "dmp-dir",
            value_name = "DIR",
            help = "Directory containing nodes.dmp and names.dmp (defaults to ~/.camitk)."
        )]
        dmp_dir: Option<PathBuf>,
    },
    #[command(
        about = "Sort taxa within each rank",
        long_about = "Reorder taxa in every sample either by decreasing abundance (dropping zero-abundance entries) or by their taxonomy paths so related lineages stay adjacent."
    )]
    Sort {
        #[arg(
            short = 'a',
            long,
            conflicts_with = "taxpath",
            help = "Sort taxa by descending abundance and drop entries with zero abundance."
        )]
        abundance: bool,
        #[arg(
            short = 't',
            value_enum,
            num_args = 0..=1,
            default_missing_value = "taxpath",
            help = "Sort taxa by their taxonomy path (TAXPATH or TAXPATHSN)."
        )]
        taxpath: Option<TaxPathField>,
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
        output: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let result = match &cli.command {
        Commands::Benchmark {
            ground_truth,
            predictions,
            update_taxonomy,
            labels,
            all_filter,
            ground_filter,
            pred_filter,
            normalize,
            by_domain,
            group_realms,
            no_group_realms,
            output,
            ranks,
            dmp_dir,
        } => {
            let label_vec = labels.as_ref().map(|s| split_labels(s)).unwrap_or_default();
            let rank_vec = ranks.as_ref().map(|s| split_ranks(s));
            let cfg = BenchmarkRunConfig {
                ground_truth: ground_truth.clone(),
                predictions: predictions.clone(),
                update_taxonomy: *update_taxonomy,
                labels: label_vec,
                all_filter: all_filter.clone(),
                ground_filter: ground_filter.clone(),
                pred_filter: pred_filter.clone(),
                normalize: *normalize,
                by_domain: *by_domain,
                group_realms: *group_realms && !*no_group_realms,
                output: output.clone(),
                ranks: rank_vec,
                dmp_dir: dmp_dir.clone(),
            };
            benchmark_cmd::run(&cfg)
        }
        Commands::Filter {
            expression,
            output,
            fill_up,
            from_rank,
            to_rank,
            renorm,
            input,
            dmp_dir,
        } => {
            let cfg = FilterConfig {
                expression,
                output: output.as_ref(),
                fill_up: *fill_up,
                from_rank: from_rank.as_deref(),
                to_rank,
                renorm: *renorm,
                input: input.as_ref(),
                dmp_dir: dmp_dir.as_ref(),
            };
            filter_cmd::run(&cfg)
        }
        Commands::List { input } => {
            let cfg = ListConfig {
                input: input.as_ref(),
            };
            list_cmd::run(&cfg)
        }
        Commands::Preview { n, input } => {
            let cfg = PreviewConfig {
                n: *n,
                input: input.as_ref(),
            };
            preview_cmd::run(&cfg)
        }
        Commands::Renorm { input, output } => {
            let cfg = RenormConfig {
                input: input.as_ref(),
                output: output.as_ref(),
            };
            renorm_cmd::run(&cfg)
        }
        Commands::Fillup {
            input,
            output,
            to_rank,
            from_rank,
            dmp_dir,
        } => {
            let cfg = FillupConfig {
                input: input.as_ref(),
                output: output.as_ref(),
                to_rank: to_rank.as_deref(),
                from_rank: from_rank.as_deref(),
                dmp_dir: dmp_dir.as_ref(),
            };
            fillup_cmd::run(&cfg)
        }
        Commands::Convert {
            input,
            output,
            taxid_column,
            abundance_column,
            input_percent,
            normalize,
            sample_id,
            taxonomy_tag,
            dmp_dir,
        } => {
            let cfg = ConvertRunConfig {
                input: input.as_ref(),
                output: output.as_ref(),
                taxid_column: *taxid_column,
                abundance_column: *abundance_column,
                input_is_percent: *input_percent,
                normalize: *normalize,
                sample_id,
                taxonomy_tag: taxonomy_tag.as_deref(),
                dmp_dir: dmp_dir.as_ref(),
            };
            convert_cmd::run(&cfg)
        }
        Commands::Sort {
            abundance,
            taxpath,
            input,
            output,
        } => {
            let mode = if *abundance {
                SortMode::Abundance
            } else if let Some(field) = taxpath {
                SortMode::TaxPath(*field)
            } else {
                bail!("either -a/--abundance or -t must be provided");
            };
            let cfg = SortConfig {
                input: input.as_ref(),
                output: output.as_ref(),
                mode,
            };
            sort_cmd::run(&cfg)
        }
    };

    match result {
        Ok(()) => Ok(()),
        Err(err) if is_broken_pipe(&err) => Ok(()),
        Err(err) => Err(err),
    }
}

fn is_broken_pipe(err: &anyhow::Error) -> bool {
    err.chain()
        .any(|cause| match cause.downcast_ref::<io::Error>() {
            Some(io_err) => io_err.kind() == io::ErrorKind::BrokenPipe,
            None => false,
        })
}
