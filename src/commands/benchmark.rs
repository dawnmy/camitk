use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail, ensure};

use crate::cami;
use crate::cami::{Entry, Sample};
use crate::expression::{apply_filter, expr_needs_taxdump, parse_expression};
use crate::taxonomy::{Taxonomy, default_taxdump_dir, ensure_taxdump};

#[derive(Clone)]
pub struct BenchmarkConfig {
    pub ground_truth: PathBuf,
    pub predictions: Vec<PathBuf>,
    pub update_taxonomy: bool,
    pub labels: Vec<String>,
    pub all_filter: Option<String>,
    pub ground_filter: Option<String>,
    pub pred_filter: Option<String>,
    pub normalize: bool,
    pub by_domain: bool,
    pub group_realms: bool,
    pub output: PathBuf,
    pub ranks: Option<Vec<String>>,
    pub dmp_dir: Option<PathBuf>,
}

#[derive(Clone)]
struct ProfileEntry {
    taxid: String,
    percentage: f64,
    lineage: Arc<[Option<String>]>,
}

pub fn run(cfg: &BenchmarkConfig) -> Result<()> {
    if cfg.predictions.is_empty() {
        bail!("at least one predicted profile must be provided");
    }

    if !cfg.labels.is_empty() && cfg.labels.len() != cfg.predictions.len() {
        bail!("labels count must match number of predicted profiles");
    }

    let all_expr = cfg
        .all_filter
        .as_ref()
        .map(|s| parse_expression(s))
        .transpose()
        .context("parsing global filter expression")?;
    let ground_expr = cfg
        .ground_filter
        .as_ref()
        .map(|s| parse_expression(s))
        .transpose()
        .context("parsing ground-truth filter expression")?;
    let pred_expr = cfg
        .pred_filter
        .as_ref()
        .map(|s| parse_expression(s))
        .transpose()
        .context("parsing predicted-profile filter expression")?;

    let needs_taxdump = cfg.update_taxonomy
        || cfg.by_domain
        || cfg.group_realms
        || all_expr
            .as_ref()
            .is_some_and(|expr| expr_needs_taxdump(expr))
        || ground_expr
            .as_ref()
            .is_some_and(|expr| expr_needs_taxdump(expr))
        || pred_expr
            .as_ref()
            .is_some_and(|expr| expr_needs_taxdump(expr));

    let tax_dir = cfg.dmp_dir.clone().unwrap_or_else(default_taxdump_dir);

    let mut taxonomy = if needs_taxdump {
        ensure_taxdump(&tax_dir)
            .with_context(|| format!("ensuring taxdump in {}", tax_dir.display()))?;
        Some(Taxonomy::load(&tax_dir)?)
    } else {
        None
    };

    let mut gt_samples = cami::parse_cami(&cfg.ground_truth)?;
    if gt_samples.is_empty() {
        bail!("ground truth profile has no samples");
    }

    let mut taxonomy_cache: HashMap<u32, LineageInfo> = HashMap::new();

    if cfg.update_taxonomy || samples_need_superkingdom(&gt_samples) {
        let taxonomy_ref = ensure_taxonomy_loaded(&mut taxonomy, &tax_dir)?;
        update_samples_taxonomy(
            &mut gt_samples,
            taxonomy_ref,
            cfg.update_taxonomy,
            &mut taxonomy_cache,
            cfg.group_realms,
        );
    }

    if let Some(expr) = all_expr.as_ref() {
        gt_samples = apply_filter(&gt_samples, expr, taxonomy.as_ref());
    }

    if let Some(expr) = ground_expr.as_ref() {
        gt_samples = apply_filter(&gt_samples, expr, taxonomy.as_ref());
        if gt_samples.is_empty() {
            bail!("ground truth profile has no samples after applying filters");
        }
    } else if gt_samples.is_empty() {
        bail!("ground truth profile has no samples after applying filters");
    }

    fs::create_dir_all(&cfg.output)
        .with_context(|| format!("creating output directory {}", cfg.output.display()))?;

    let ranks = cfg.ranks.as_ref().map(|v| canonical_ranks(v)).transpose()?;

    let labels: Vec<String> = if cfg.labels.is_empty() {
        cfg.predictions
            .iter()
            .map(|p| {
                p.file_stem()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_string())
                    .unwrap_or_else(|| p.to_string_lossy().into_owned())
            })
            .collect()
    } else {
        cfg.labels.clone()
    };

    let domains = if cfg.by_domain {
        vec![
            None,
            Some("Bacteria".to_string()),
            Some("Archaea".to_string()),
            Some("Eukarya".to_string()),
            Some("Viruses".to_string()),
        ]
    } else {
        vec![None]
    };

    for domain in domains {
        let suffix = domain
            .as_ref()
            .map(|d| format!("_{}", d.to_lowercase()))
            .unwrap_or_else(|| "".to_string());
        let path = cfg.output.join(format!("benchmark{}.tsv", suffix));
        let mut writer = File::create(&path)
            .with_context(|| format!("creating benchmark report {}", path.display()))?;
        writeln!(
            writer,
            "profile\tsample\trank\ttp\tfp\tfn\tprecision\trecall\tf1\tjaccard\tl1_error\tbray_curtis\tshannon_pred\tshannon_truth\tevenness_pred\tevenness_truth\tpearson\tspearman\tweighted_unifrac\tunweighted_unifrac\tabundance_rank_error\tmass_weighted_abundance_rank_error"
        )?;

        let gt_map = build_profile_map(
            &gt_samples,
            ranks.as_ref(),
            cfg.normalize,
            domain.as_deref(),
            taxonomy.as_ref(),
            cfg.group_realms,
        );

        let mut sample_ids: Vec<_> = gt_map.keys().cloned().collect();
        sample_ids.sort();

        for (pred_path, label) in cfg.predictions.iter().zip(labels.iter()) {
            let mut pred_samples = cami::parse_cami(pred_path)?;

            if cfg.update_taxonomy || samples_need_superkingdom(&pred_samples) {
                let taxonomy_ref = ensure_taxonomy_loaded(&mut taxonomy, &tax_dir)?;
                update_samples_taxonomy(
                    &mut pred_samples,
                    taxonomy_ref,
                    cfg.update_taxonomy,
                    &mut taxonomy_cache,
                    cfg.group_realms,
                );
            }

            if let Some(expr) = all_expr.as_ref() {
                pred_samples = apply_filter(&pred_samples, expr, taxonomy.as_ref());
            }
            if let Some(expr) = pred_expr.as_ref() {
                pred_samples = apply_filter(&pred_samples, expr, taxonomy.as_ref());
            }
            let pred_map = build_profile_map(
                &pred_samples,
                ranks.as_ref(),
                cfg.normalize,
                domain.as_deref(),
                taxonomy.as_ref(),
                cfg.group_realms,
            );

            for sample_id in &sample_ids {
                let Some(gt_ranks) = gt_map.get(sample_id) else {
                    continue;
                };
                let mut rank_names: Vec<_> = gt_ranks.keys().cloned().collect();
                rank_names.sort();
                let pred_ranks = pred_map.get(sample_id);
                for rank in rank_names {
                    if rank_is_above_phylum(&rank) {
                        continue;
                    }
                    let gt_entries = gt_ranks.get(&rank).unwrap();
                    let pred_entries = pred_ranks
                        .and_then(|map| map.get(&rank))
                        .map(|v| v.as_slice())
                        .unwrap_or(&[]);
                    let metrics = compute_metrics(&rank, gt_entries, pred_entries)?;
                    write_metrics(&mut writer, label, sample_id, &rank, &metrics)?;
                }
            }
        }
    }

    Ok(())
}

fn canonical_ranks(ranks: &[String]) -> Result<Vec<String>> {
    ranks.iter().map(|r| canonical_rank(r)).collect()
}

fn rank_is_above_phylum(rank: &str) -> bool {
    matches!(
        rank.trim().to_lowercase().as_str(),
        "superkingdom" | "domain" | "kingdom" | "realm" | "cellular root" | "acellular root"
    )
}

fn canonical_rank(rank: &str) -> Result<String> {
    let trimmed = rank.trim();
    if trimmed.is_empty() {
        bail!("rank names must not be empty");
    }
    let lower = trimmed.to_lowercase();
    let canonical = match lower.as_str() {
        "d" | "domain" | "superkingdom" | "realm" | "cellular root" | "acellular root" => {
            "superkingdom"
        }
        "k" | "kingdom" => "kingdom",
        "p" | "phylum" => "phylum",
        "c" | "class" => "class",
        "o" | "order" => "order",
        "f" | "family" => "family",
        "g" | "genus" => "genus",
        "s" | "species" => "species",
        "t" | "strain" => "strain",
        other => other,
    };
    Ok(canonical.to_string())
}

fn build_profile_map(
    samples: &[Sample],
    ranks: Option<&Vec<String>>,
    normalize: bool,
    domain: Option<&str>,
    taxonomy: Option<&Taxonomy>,
    group_realms: bool,
) -> HashMap<String, HashMap<String, Vec<ProfileEntry>>> {
    let mut map: HashMap<String, HashMap<String, Vec<ProfileEntry>>> = HashMap::new();
    let rank_filter: Option<BTreeSet<String>> = ranks.map(|r| r.iter().cloned().collect());
    let domain_lower = domain.map(|d| d.to_lowercase());
    let mut lineage_cache: HashMap<String, Arc<[Option<String>]>> = HashMap::new();

    for sample in samples {
        let mut sample_map: HashMap<String, Vec<ProfileEntry>> = HashMap::new();
        for entry in &sample.entries {
            if entry.percentage <= 0.0 {
                continue;
            }
            if let Some(filter) = &rank_filter {
                if !filter.contains(&entry.rank) {
                    continue;
                }
            }
            if let Some(domain) = &domain_lower {
                if !entry_belongs_to_domain(entry, domain, taxonomy, group_realms) {
                    continue;
                }
            }
            let rank_name = canonical_rank(&entry.rank).unwrap_or_else(|_| entry.rank.clone());
            let cache_key = lineage_cache_key(entry, &rank_name);
            let lineage = lineage_cache
                .entry(cache_key)
                .or_insert_with(|| {
                    Arc::from(
                        compute_lineage(entry, &rank_name, taxonomy, group_realms)
                            .into_boxed_slice(),
                    )
                })
                .clone();
            sample_map
                .entry(rank_name.clone())
                .or_default()
                .push(ProfileEntry {
                    taxid: entry.taxid.clone(),
                    percentage: entry.percentage,
                    lineage,
                });
        }

        if normalize {
            for entries in sample_map.values_mut() {
                let sum: f64 = entries.iter().map(|e| e.percentage).sum();
                if sum > 0.0 {
                    for entry in entries.iter_mut() {
                        entry.percentage = entry.percentage / sum * 100.0;
                    }
                }
            }
        }

        if !sample_map.is_empty() {
            map.insert(sample.id.clone(), sample_map);
        }
    }

    map
}

fn entry_belongs_to_domain(
    entry: &Entry,
    domain: &str,
    taxonomy: Option<&Taxonomy>,
    group_realms: bool,
) -> bool {
    let Some(target) = canonical_domain_name(domain, group_realms) else {
        return false;
    };

    entry_domain_candidates(entry, taxonomy, group_realms)
        .into_iter()
        .any(|candidate| candidate.eq_ignore_ascii_case(&target))
}

fn entry_domain_candidates(
    entry: &Entry,
    taxonomy: Option<&Taxonomy>,
    group_realms: bool,
) -> Vec<String> {
    let mut candidates: Vec<String> = Vec::new();

    let mut push_candidate = |value: String| {
        if !candidates
            .iter()
            .any(|existing| existing.eq_ignore_ascii_case(&value))
        {
            candidates.push(value);
        }
    };

    for component in entry.taxpathsn.split('|') {
        if let Some(name) = canonical_domain_name(component, group_realms) {
            push_candidate(name);
        }
    }

    if group_realms {
        if entry.taxid.trim() == "10239" {
            push_candidate("Viruses".to_string());
        }
        for component in split_taxpath(&entry.taxpath) {
            if component.trim() == "10239" {
                push_candidate("Viruses".to_string());
            }
        }
    }

    if let Some(tax) = taxonomy {
        let mut visit_taxid = |tid: u32| {
            if let Some(name) = tax.domain_of(tid) {
                if let Some(normalized) = canonical_domain_name(&name, group_realms) {
                    push_candidate(normalized);
                }
            }
            if group_realms {
                if let Some(realm) = tax.realm_of(tid) {
                    if component_is_virus_realm(&realm) {
                        push_candidate("Viruses".to_string());
                    }
                }
            }
        };

        if let Some(tid) = highest_taxid_in_taxpath(entry, Some(tax)) {
            visit_taxid(tid);
        }
        if let Some(tid) = tax.resolve_taxid_str(&entry.taxid) {
            visit_taxid(tid);
        }
    }

    candidates
}

fn canonical_domain_name(raw: &str, group_realms: bool) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let lower = trimmed.to_ascii_lowercase();
    if lower.contains("bacteria") {
        Some("Bacteria".to_string())
    } else if lower.contains("archaea") || lower.contains("archaeon") || lower.contains("archaeota")
    {
        Some("Archaea".to_string())
    } else if lower.contains("eukary") {
        Some("Eukarya".to_string())
    } else if lower.contains("viruses") {
        Some("Viruses".to_string())
    } else if group_realms && component_is_virus_realm(trimmed) {
        Some("Viruses".to_string())
    } else {
        None
    }
}

fn component_is_virus_realm(name: &str) -> bool {
    let lower = name.trim().to_ascii_lowercase();
    lower.contains("viria") || lower.contains("virae") || lower.contains("viruses")
}

fn highest_taxid_in_taxpath(entry: &Entry, taxonomy: Option<&Taxonomy>) -> Option<u32> {
    let raw = entry
        .taxpath
        .split('|')
        .map(|tid| tid.trim())
        .find_map(|tid| {
            if tid.is_empty() {
                None
            } else {
                tid.parse::<u32>().ok().filter(|id| *id > 1)
            }
        })
        .or_else(|| entry.taxid.trim().parse::<u32>().ok());

    match (taxonomy, raw) {
        (Some(tax), Some(id)) => tax.resolve_taxid(id),
        (_, value) => value,
    }
}

#[derive(Default)]
struct Metrics {
    tp: usize,
    fp: usize,
    fn_: usize,
    precision: Option<f64>,
    recall: Option<f64>,
    f1: Option<f64>,
    jaccard: Option<f64>,
    l1_error: f64,
    bray_curtis: Option<f64>,
    shannon_pred: f64,
    shannon_truth: f64,
    evenness_pred: Option<f64>,
    evenness_truth: Option<f64>,
    pearson: Option<f64>,
    spearman: Option<f64>,
    weighted_unifrac: Option<f64>,
    unweighted_unifrac: Option<f64>,
    abundance_rank_error: Option<f64>,
    mass_weighted_abundance_rank_error: Option<f64>,
}

fn compute_metrics(
    rank: &str,
    gt_entries: &[ProfileEntry],
    pred_entries: &[ProfileEntry],
) -> Result<Metrics> {
    let mut metrics = Metrics::default();
    let mut gt_map: HashMap<String, f64> = HashMap::new();
    for entry in gt_entries {
        ensure!(
            !entry.percentage.is_nan(),
            "ground truth abundance for taxid {} is NaN",
            entry.taxid
        );
        ensure!(
            entry.percentage >= 0.0,
            "ground truth abundance for taxid {} is negative",
            entry.taxid
        );
        if entry.percentage > 0.0 {
            gt_map.insert(entry.taxid.clone(), entry.percentage);
        }
    }
    let mut pred_map: HashMap<String, f64> = HashMap::new();
    for entry in pred_entries {
        ensure!(
            !entry.percentage.is_nan(),
            "predicted abundance for taxid {} is NaN",
            entry.taxid
        );
        ensure!(
            entry.percentage >= 0.0,
            "predicted abundance for taxid {} is negative",
            entry.taxid
        );
        if entry.percentage > 0.0 {
            pred_map.insert(entry.taxid.clone(), entry.percentage);
        }
    }

    let gt_total: f64 = gt_entries.iter().map(|e| e.percentage).sum();
    let pred_total: f64 = pred_entries.iter().map(|e| e.percentage).sum();

    let mut union: BTreeSet<String> = BTreeSet::new();
    for taxid in gt_map.keys() {
        union.insert(taxid.clone());
    }
    for taxid in pred_map.keys() {
        union.insert(taxid.clone());
    }

    let mut gt_values: Vec<f64> = Vec::new();
    let mut pred_values: Vec<f64> = Vec::new();

    let mut sum_abs = 0.0;
    let mut sum_tot = 0.0;

    for taxid in &union {
        let g = gt_map.get(taxid).copied().unwrap_or(0.0);
        let p = pred_map.get(taxid).copied().unwrap_or(0.0);
        if g > 0.0 && p > 0.0 {
            metrics.tp += 1;
        } else if p > 0.0 {
            metrics.fp += 1;
        } else if g > 0.0 {
            metrics.fn_ += 1;
        }

        let g_norm = if gt_total > 0.0 { g / gt_total } else { 0.0 };
        let p_norm = if pred_total > 0.0 {
            p / pred_total
        } else {
            0.0
        };

        sum_abs += (g_norm - p_norm).abs();
        sum_tot += g_norm + p_norm;
        gt_values.push(g_norm);
        pred_values.push(p_norm);
    }

    metrics.l1_error = sum_abs;
    metrics.precision = ratio(metrics.tp as f64, (metrics.tp + metrics.fp) as f64);
    metrics.recall = ratio(metrics.tp as f64, (metrics.tp + metrics.fn_) as f64);
    metrics.f1 = match (metrics.precision, metrics.recall) {
        (Some(p), Some(r)) if p + r > 0.0 => Some(2.0 * p * r / (p + r)),
        _ => None,
    };
    metrics.jaccard = ratio(
        metrics.tp as f64,
        (metrics.tp + metrics.fp + metrics.fn_) as f64,
    );
    metrics.bray_curtis = if sum_tot > 0.0 {
        Some(sum_abs / sum_tot)
    } else {
        None
    };

    metrics.shannon_truth = shannon(&gt_values);
    metrics.shannon_pred = shannon(&pred_values);
    metrics.evenness_truth = evenness(&gt_values, metrics.shannon_truth);
    metrics.evenness_pred = evenness(&pred_values, metrics.shannon_pred);
    metrics.pearson = pearson(&gt_values, &pred_values);
    metrics.spearman = spearman(&gt_values, &pred_values);

    let unifrac = unifrac(rank, gt_entries, pred_entries);
    metrics.weighted_unifrac = unifrac.weighted;
    metrics.unweighted_unifrac = unifrac.unweighted;
    metrics.abundance_rank_error = Some(abundance_rank_error(&gt_map, &pred_map)?);
    metrics.mass_weighted_abundance_rank_error =
        Some(mass_weighted_abundance_rank_error(&gt_map, &pred_map, 1.0)?);

    Ok(metrics)
}

fn ratio(num: f64, denom: f64) -> Option<f64> {
    if denom > 0.0 { Some(num / denom) } else { None }
}

fn shannon(values: &[f64]) -> f64 {
    let sum: f64 = values.iter().copied().sum();
    if sum <= 0.0 {
        return 0.0;
    }
    let mut h = 0.0;
    for v in values {
        if *v <= 0.0 {
            continue;
        }
        let p = v / sum;
        h -= p * p.ln();
    }
    h
}

fn evenness(values: &[f64], shannon: f64) -> Option<f64> {
    let count = values.iter().filter(|v| **v > 0.0).count();
    if count <= 1 || shannon <= 0.0 {
        return None;
    }
    let denom = (count as f64).ln();
    if denom <= 0.0 {
        None
    } else {
        Some(shannon / denom)
    }
}

fn pearson(x: &[f64], y: &[f64]) -> Option<f64> {
    if x.len() != y.len() || x.len() < 2 {
        return None;
    }
    let mean_x = x.iter().copied().sum::<f64>() / x.len() as f64;
    let mean_y = y.iter().copied().sum::<f64>() / y.len() as f64;
    let mut num = 0.0;
    let mut denom_x = 0.0;
    let mut denom_y = 0.0;
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let dx = xi - mean_x;
        let dy = yi - mean_y;
        num += dx * dy;
        denom_x += dx * dx;
        denom_y += dy * dy;
    }
    if denom_x <= 0.0 || denom_y <= 0.0 {
        None
    } else {
        Some(num / (denom_x.sqrt() * denom_y.sqrt()))
    }
}

fn spearman(x: &[f64], y: &[f64]) -> Option<f64> {
    if x.len() != y.len() || x.len() < 2 {
        return None;
    }
    let rx = rank_values(x);
    let ry = rank_values(y);
    pearson(&rx, &ry)
}

fn rank_values(values: &[f64]) -> Vec<f64> {
    let mut pairs: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let mut ranks = vec![0.0; values.len()];
    let mut i = 0;
    while i < pairs.len() {
        let value = pairs[i].1;
        let mut j = i + 1;
        while j < pairs.len() && (pairs[j].1 - value).abs() < f64::EPSILON {
            j += 1;
        }
        let rank = (i + j - 1) as f64 / 2.0 + 1.0;
        for k in i..j {
            ranks[pairs[k].0] = rank;
        }
        i = j;
    }
    ranks
}

#[derive(Clone, Copy, Debug, Default)]
pub struct UniFracResult {
    pub weighted: Option<f64>,
    pub unweighted: Option<f64>,
}

struct UnifracComponents {
    weighted_raw: f64,
    unweighted_raw: f64,
    unweighted_max: f64,
    max_weighted: f64,
}

const CANONICAL_RANKS: [&str; 8] = [
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
];

const BRANCH_LENGTHS: [f64; CANONICAL_RANKS.len()] = [1.0; CANONICAL_RANKS.len()];
const MASS_EPSILON: f64 = 1e-12;

fn branch_length_for_level(level: usize) -> f64 {
    BRANCH_LENGTHS
        .get(level)
        .copied()
        .unwrap_or_else(|| BRANCH_LENGTHS.last().copied().unwrap_or(1.0))
}

fn unifrac(
    rank: &str,
    gt_entries: &[ProfileEntry],
    pred_entries: &[ProfileEntry],
) -> UniFracResult {
    let Some(rank_index) = canonical_rank_index(rank) else {
        return UniFracResult::default();
    };

    match unifrac_components(rank_index, gt_entries, pred_entries) {
        Some(components) => {
            let weighted = if components.max_weighted > 0.0 {
                Some((components.weighted_raw / components.max_weighted).clamp(0.0, 1.0))
            } else {
                None
            };
            let unweighted = if components.unweighted_max > MASS_EPSILON {
                Some((components.unweighted_raw / components.unweighted_max).clamp(0.0, 1.0))
            } else {
                None
            };
            UniFracResult {
                weighted,
                unweighted,
            }
        }
        None => UniFracResult::default(),
    }
}

fn unifrac_components(
    rank_index: usize,
    gt_entries: &[ProfileEntry],
    pred_entries: &[ProfileEntry],
) -> Option<UnifracComponents> {
    let gt_total: f64 = gt_entries
        .iter()
        .filter(|e| e.percentage > 0.0)
        .map(|e| e.percentage)
        .sum();
    let pred_total: f64 = pred_entries
        .iter()
        .filter(|e| e.percentage > 0.0)
        .map(|e| e.percentage)
        .sum();

    if gt_total <= 0.0 && pred_total <= 0.0 {
        return None;
    }

    let mut tree = TaxTree::new();

    if gt_total > 0.0 {
        for entry in gt_entries.iter().filter(|e| e.percentage > 0.0) {
            let mass = entry.percentage / gt_total;
            if !mass.is_finite() || mass <= 0.0 {
                continue;
            }
            tree.add_mass(entry.lineage.as_ref(), rank_index, mass, true);
        }
    }

    if pred_total > 0.0 {
        for entry in pred_entries.iter().filter(|e| e.percentage > 0.0) {
            let mass = entry.percentage / pred_total;
            if !mass.is_finite() || mass <= 0.0 {
                continue;
            }
            tree.add_mass(entry.lineage.as_ref(), rank_index, mass, false);
        }
    }

    let (weighted_raw, unweighted_raw, unweighted_max) = tree.compute_edge_flows();
    let path_length = BRANCH_LENGTHS
        .iter()
        .take(rank_index + 1)
        .copied()
        .sum::<f64>();
    let max_weighted = 2.0 * path_length;
    let clamped_weighted = if max_weighted > 0.0 {
        weighted_raw.max(0.0).min(max_weighted)
    } else {
        weighted_raw.max(0.0)
    };

    Some(UnifracComponents {
        weighted_raw: clamped_weighted,
        unweighted_raw: unweighted_raw.max(0.0),
        unweighted_max,
        max_weighted,
    })
}

fn canonical_rank_index(rank: &str) -> Option<usize> {
    let canonical = canonical_rank(rank).ok()?;
    match canonical.as_str() {
        "superkingdom" => Some(0),
        "phylum" => Some(1),
        "class" => Some(2),
        "order" => Some(3),
        "family" => Some(4),
        "genus" => Some(5),
        "species" => Some(6),
        "strain" => Some(7),
        _ => None,
    }
}

fn lineage_cache_key(entry: &Entry, canonical_rank: &str) -> String {
    let base = if let Ok(tid) = entry.taxid.trim().parse::<u32>() {
        format!("taxid:{}", tid)
    } else if !entry.taxpath.trim().is_empty() {
        entry.taxpath.trim().to_string()
    } else if !entry.taxpathsn.trim().is_empty() {
        entry.taxpathsn.trim().to_string()
    } else {
        entry.taxid.trim().to_string()
    };

    format!("{}|{}", canonical_rank, base)
}

fn compute_lineage(
    entry: &Entry,
    canonical_rank: &str,
    taxonomy: Option<&Taxonomy>,
    group_realms: bool,
) -> Vec<Option<String>> {
    let entry_rank_index = canonical_rank_index(canonical_rank)
        .unwrap_or_else(|| CANONICAL_RANKS.len().saturating_sub(1));
    let mut result = vec![None; CANONICAL_RANKS.len()];

    if let Some(taxonomy) = taxonomy {
        if let Some(taxid) = taxonomy
            .resolve_taxid_str(&entry.taxid)
            .or_else(|| highest_taxid_in_taxpath(entry, Some(taxonomy)))
        {
            fill_lineage_from_taxonomy(&mut result, entry_rank_index, taxonomy, taxid);
        }
    }

    let mut normalized = taxid_lineage_from_path(&entry.taxpath, &entry.taxid);

    if normalized.is_empty() {
        normalized = split_taxpath(&entry.taxpathsn)
            .into_iter()
            .filter_map(|part| normalize_taxpath_component(&part, group_realms))
            .collect();
    }

    if normalized.is_empty() {
        normalized = split_taxpath(&entry.taxpath)
            .into_iter()
            .filter_map(|part| normalize_taxpath_component(&part, group_realms))
            .collect();
    }

    if normalized.is_empty() {
        if let Some(value) = normalize_taxid_component(&entry.taxid) {
            normalized.push(value);
        } else if let Some(value) = normalize_taxpath_component(&entry.taxid, group_realms) {
            normalized.push(value);
        }
    }

    if normalized.len() > entry_rank_index + 1 {
        normalized = normalized.split_off(normalized.len() - (entry_rank_index + 1));
    }

    let offset = entry_rank_index + 1 - normalized.len();
    for (idx, name) in normalized.into_iter().enumerate() {
        let target = offset + idx;
        if result[target].is_none() {
            result[target] = Some(name);
        }
    }

    if result[0].is_none() {
        if let Some(superkingdom) = infer_superkingdom(entry, taxonomy, group_realms) {
            result[0] = Some(superkingdom);
        }
    }

    result
}

fn fill_lineage_from_taxonomy(
    result: &mut [Option<String>],
    entry_rank_index: usize,
    taxonomy: &Taxonomy,
    taxid: u32,
) {
    let lineage = taxonomy.lineage(taxid);
    for (tid, rank, _) in lineage.iter() {
        if let Some(idx) = canonical_rank_index(rank) {
            if idx <= entry_rank_index {
                result[idx] = Some(format!("taxid:{}", tid));
            }
        }
    }
}

fn taxid_lineage_from_path(taxpath: &str, taxid: &str) -> Vec<String> {
    let mut lineage: Vec<String> = split_taxpath(taxpath)
        .into_iter()
        .filter_map(|component| normalize_taxid_component(&component))
        .collect();

    if let Some(taxid_component) = normalize_taxid_component(taxid) {
        if lineage
            .last()
            .map(|last| last != &taxid_component)
            .unwrap_or(true)
        {
            lineage.push(taxid_component);
        }
    }

    lineage
}

fn split_taxpath(path: &str) -> Vec<String> {
    path.split('|')
        .map(|part| part.trim())
        .filter(|part| !part.is_empty())
        .map(|part| part.to_string())
        .collect()
}

fn ensure_taxonomy_loaded<'a>(
    taxonomy: &'a mut Option<Taxonomy>,
    dir: &Path,
) -> Result<&'a Taxonomy> {
    if taxonomy.is_none() {
        ensure_taxdump(dir).with_context(|| format!("ensuring taxdump in {}", dir.display()))?;
        *taxonomy = Some(Taxonomy::load(dir)?);
    }
    Ok(taxonomy.as_ref().unwrap())
}

fn samples_need_superkingdom(samples: &[Sample]) -> bool {
    samples
        .iter()
        .flat_map(|sample| sample.entries.iter())
        .any(|entry| entry_missing_superkingdom(entry))
}

fn entry_missing_superkingdom(entry: &Entry) -> bool {
    let mut id_parts = entry
        .taxpath
        .split('|')
        .map(|part| part.trim())
        .filter(|part| !part.is_empty());
    let mut name_parts = entry
        .taxpathsn
        .split('|')
        .map(|part| part.trim())
        .filter(|part| !part.is_empty());

    let id_present = id_parts
        .next()
        .and_then(|value| value.parse::<u32>().ok())
        .map(|tid| matches!(tid, 2 | 2157 | 2759 | 10239))
        .unwrap_or(false);

    let name_present = name_parts
        .next()
        .map(|value| {
            let lower = value.to_lowercase();
            lower.contains("bacteria")
                || lower.contains("archaea")
                || lower.contains("eukary")
                || lower.contains("viruses")
        })
        .unwrap_or(false);

    !(id_present || name_present)
}

fn update_samples_taxonomy(
    samples: &mut [Sample],
    taxonomy: &Taxonomy,
    update_paths: bool,
    cache: &mut HashMap<u32, LineageInfo>,
    group_realms: bool,
) {
    for sample in samples.iter_mut() {
        for entry in sample.entries.iter_mut() {
            let Some(taxid) = taxonomy.resolve_taxid_str(&entry.taxid) else {
                continue;
            };
            let info = cache
                .entry(taxid)
                .or_insert_with(|| canonical_lineage(taxonomy, taxid, group_realms));
            if update_paths {
                apply_full_update(entry, info);
            } else {
                ensure_superkingdom_only(entry, info);
            }
        }
    }
}

fn apply_full_update(entry: &mut Entry, info: &LineageInfo) {
    entry.taxid.clear();
    entry.taxid.push_str(info.resolved_taxid());
    if let Some(rank_name) = info.rank_name.as_ref() {
        entry.rank = rank_name.clone();
    } else if let Ok(canonical) = canonical_rank(&entry.rank) {
        entry.rank = canonical;
    }

    let target_idx = info
        .rank_index
        .or_else(|| canonical_rank_index(&entry.rank))
        .unwrap_or(info.deepest_index);

    let (taxpath, taxpathsn) = info.path_up_to(target_idx);
    if !taxpath.is_empty() {
        entry.taxpath = taxpath;
    }
    if !taxpathsn.is_empty() {
        entry.taxpathsn = taxpathsn;
    }
}

fn ensure_superkingdom_only(entry: &mut Entry, info: &LineageInfo) {
    if let Some((taxid, name)) = info.superkingdom() {
        let mut ids = split_taxpath(&entry.taxpath);
        if !ids.first().map(|first| first == &taxid).unwrap_or(false) {
            ids.insert(0, taxid.clone());
        }

        let mut names = split_taxpath(&entry.taxpathsn);
        if !names
            .first()
            .map(|first| first.eq_ignore_ascii_case(&name))
            .unwrap_or(false)
        {
            names.insert(0, name.clone());
        }

        entry.taxpath = ids.join("|");
        entry.taxpathsn = names.join("|");
    }
}

#[derive(Clone)]
struct LineageInfo {
    canonical: Arc<[Option<(String, String)>]>,
    rank_index: Option<usize>,
    rank_name: Option<String>,
    deepest_index: usize,
    resolved_taxid: String,
}

impl LineageInfo {
    fn path_up_to(&self, upto: usize) -> (String, String) {
        if self.canonical.is_empty() {
            return (String::new(), String::new());
        }
        let limit = upto.min(self.canonical.len().saturating_sub(1));
        let mut ids = Vec::new();
        let mut names = Vec::new();
        for idx in 0..=limit {
            if let Some((tid, name)) = self.canonical[idx].as_ref() {
                ids.push(tid.clone());
                names.push(name.clone());
            }
        }
        (ids.join("|"), names.join("|"))
    }

    fn superkingdom(&self) -> Option<(String, String)> {
        self.canonical.get(0).and_then(|entry| {
            entry
                .as_ref()
                .map(|(tid, name)| (tid.clone(), name.clone()))
        })
    }

    fn resolved_taxid(&self) -> &str {
        &self.resolved_taxid
    }
}

fn canonical_lineage(taxonomy: &Taxonomy, taxid: u32, group_realms: bool) -> LineageInfo {
    let lineage = taxonomy.lineage(taxid);
    let mut canonical: Vec<Option<(String, String)>> = vec![None; CANONICAL_RANKS.len()];
    let mut rank_index: Option<usize> = None;
    let mut rank_name: Option<String> = None;

    for (tid, rank, name) in lineage.iter() {
        let canonical_rank_name = canonical_rank(rank).unwrap_or_else(|_| rank.clone());
        if let Some(idx) = canonical_rank_index(&canonical_rank_name) {
            canonical[idx] = Some((tid.to_string(), name.clone()));
            if *tid == taxid {
                rank_index = Some(idx);
                rank_name = Some(CANONICAL_RANKS[idx].to_string());
            }
        } else if *tid == taxid {
            rank_name = Some(canonical_rank_name.clone());
        }
    }

    let taxid_str = taxid.to_string();
    let tip_name = taxonomy.name_of(taxid).unwrap_or_else(|| taxid_str.clone());

    if let Some(idx) = rank_index {
        let needs_update = canonical[idx]
            .as_ref()
            .map(|(tid, _)| tid != &taxid_str)
            .unwrap_or(true);
        if needs_update {
            canonical[idx] = Some((taxid_str.clone(), tip_name.clone()));
        }
    } else {
        let fallback_idx = canonical
            .iter()
            .rposition(|v| v.is_some())
            .map(|idx| (idx + 1).min(CANONICAL_RANKS.len() - 1))
            .unwrap_or(CANONICAL_RANKS.len() - 1);
        canonical[fallback_idx] = Some((taxid_str.clone(), tip_name.clone()));
        rank_index = Some(fallback_idx);
        rank_name = Some(CANONICAL_RANKS[fallback_idx].to_string());
    }

    let mut domain_name = taxonomy.domain_of(taxid);
    let belongs_to_virus = group_realms
        && lineage.iter().any(|(tid, rank, name)| {
            *tid == 10239
                || rank.eq_ignore_ascii_case("realm") && component_is_virus_realm(name)
                || name.eq_ignore_ascii_case("viruses")
        });

    if belongs_to_virus {
        canonical[0] = Some(("10239".to_string(), "Viruses".to_string()));
        domain_name = Some("Viruses".to_string());
    }

    if canonical[0].is_none() {
        if let Some(domain) = domain_name.clone() {
            if let Some((tid, _, name)) = lineage
                .iter()
                .find(|(_, _, n)| n.eq_ignore_ascii_case(&domain))
            {
                canonical[0] = Some((tid.to_string(), name.clone()))
            } else if domain.eq_ignore_ascii_case("Viruses") {
                canonical[0] = Some(("10239".to_string(), domain));
            } else {
                canonical[0] = Some((taxid_str.clone(), domain));
            }
        }
    } else if let Some((tid, name)) = canonical[0].as_mut() {
        if name.eq_ignore_ascii_case("acellular root")
            || name.eq_ignore_ascii_case("acellular organisms")
            || name.eq_ignore_ascii_case("cellular organisms")
        {
            if let Some(domain) = domain_name.clone() {
                if !domain.eq_ignore_ascii_case(name) {
                    if let Some((domain_tid, _, _)) = lineage
                        .iter()
                        .find(|(_, _, n)| n.eq_ignore_ascii_case(&domain))
                    {
                        *tid = domain_tid.to_string();
                        *name = domain;
                    } else if domain.eq_ignore_ascii_case("Viruses") {
                        *tid = "10239".to_string();
                        *name = domain;
                    } else {
                        *name = domain;
                    }
                }
            }
        }
    }

    let deepest_index = canonical
        .iter()
        .rposition(|v| v.is_some())
        .unwrap_or(rank_index.unwrap_or(0));

    LineageInfo {
        canonical: Arc::from(canonical.into_boxed_slice()),
        rank_index,
        rank_name,
        deepest_index,
        resolved_taxid: taxid_str,
    }
}

fn normalize_taxpath_component(raw: &str, group_realms: bool) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let lower = trimmed.to_lowercase();
    if matches!(
        lower.as_str(),
        "root" | "cellular organisms" | "cellular root" | "acellular organisms" | "acellular root"
    ) {
        return None;
    }
    if lower.contains("viruses") && !lower.eq("viruses") {
        return Some("Viruses".to_string());
    }
    if group_realms && component_is_virus_realm(trimmed) {
        return Some("Viruses".to_string());
    }
    Some(trimmed.to_string())
}

fn normalize_taxid_component(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    trimmed
        .parse::<u32>()
        .ok()
        .filter(|id| *id != 0)
        .map(|id| format!("taxid:{}", id))
}

fn infer_superkingdom(
    entry: &Entry,
    taxonomy: Option<&Taxonomy>,
    group_realms: bool,
) -> Option<String> {
    if let Some(taxonomy) = taxonomy {
        if let Some(tid) = taxonomy
            .resolve_taxid_str(&entry.taxid)
            .or_else(|| highest_taxid_in_taxpath(entry, Some(taxonomy)))
        {
            if group_realms {
                if let Some(realm) = taxonomy.realm_of(tid) {
                    if component_is_virus_realm(&realm) {
                        return Some("Viruses".to_string());
                    }
                }
            }
            if let Some(name) = taxonomy.domain_of(tid) {
                if let Some(normalized) = canonical_domain_name(&name, group_realms) {
                    return Some(normalized);
                }
                return Some(name);
            }
        }
    }

    for component in split_taxpath(&entry.taxpathsn) {
        if let Some(label) = normalize_taxpath_component(&component, group_realms)
            .and_then(|value| detect_superkingdom(&value, group_realms))
        {
            return Some(label);
        }
        if let Some(label) = detect_superkingdom(&component, group_realms) {
            return Some(label);
        }
    }

    for component in split_taxpath(&entry.taxpath) {
        if let Some(label) = normalize_taxpath_component(&component, group_realms)
            .and_then(|value| detect_superkingdom(&value, group_realms))
        {
            return Some(label);
        }
        if let Some(label) = detect_superkingdom(&component, group_realms) {
            return Some(label);
        }
        if let Ok(taxid) = component.parse::<u32>() {
            if let Some(name) = superkingdom_from_taxid(taxid) {
                return Some(name.to_string());
            }
        }
    }

    detect_superkingdom(&entry.taxid, group_realms)
}

fn detect_superkingdom(raw: &str, group_realms: bool) -> Option<String> {
    let trimmed = raw
        .trim()
        .split("__")
        .last()
        .unwrap_or(raw)
        .split(':')
        .last()
        .unwrap_or(raw)
        .trim();
    canonical_domain_name(trimmed, group_realms)
}

fn superkingdom_from_taxid(taxid: u32) -> Option<&'static str> {
    match taxid {
        2 => Some("Bacteria"),
        2157 => Some("Archaea"),
        2759 => Some("Eukarya"),
        10239 => Some("Viruses"),
        _ => None,
    }
}

struct TaxTree {
    nodes: Vec<TaxNode>,
    index: HashMap<(usize, String), usize>,
}

impl TaxTree {
    fn new() -> Self {
        let root = TaxNode {
            name: "root".to_string(),
            parent: None,
            children: Vec::new(),
            gt_mass: 0.0,
            pred_mass: 0.0,
            branch_length: 0.0,
        };
        Self {
            nodes: vec![root],
            index: HashMap::new(),
        }
    }

    fn add_mass(&mut self, path: &[Option<String>], rank_index: usize, mass: f64, is_gt: bool) {
        if mass <= 0.0 {
            return;
        }

        let mut parent = 0usize;
        let mut parent_name = self.nodes[parent].name.clone();
        let mut target_node = parent;

        for level in 0..CANONICAL_RANKS.len() {
            let node_name = match path.get(level).and_then(|p| p.clone()) {
                Some(name) => name,
                None => format!("__missing_{}_under_{}", CANONICAL_RANKS[level], parent_name),
            };
            let child = self.get_or_create_child(parent, level, node_name.clone());
            parent = child;
            if level == rank_index {
                target_node = child;
            }
            parent_name = node_name;
        }

        if is_gt {
            self.nodes[target_node].gt_mass += mass;
        } else {
            self.nodes[target_node].pred_mass += mass;
        }
    }

    fn get_or_create_child(&mut self, parent: usize, rank_index: usize, name: String) -> usize {
        let key = (parent, name.clone());
        if let Some(id) = self.index.get(&key) {
            return *id;
        }

        let id = self.nodes.len();
        self.index.insert(key, id);
        self.nodes[parent].children.push(id);
        self.nodes.push(TaxNode {
            name,
            parent: Some(parent),
            children: Vec::new(),
            gt_mass: 0.0,
            pred_mass: 0.0,
            branch_length: branch_length_for_level(rank_index),
        });
        id
    }

    fn compute_edge_flows(&self) -> (f64, f64, f64) {
        let mut gt_subtree: Vec<f64> = self.nodes.iter().map(|node| node.gt_mass).collect();
        let mut pred_subtree: Vec<f64> = self.nodes.iter().map(|node| node.pred_mass).collect();

        for node_id in (1..self.nodes.len()).rev() {
            if let Some(parent) = self.nodes[node_id].parent {
                gt_subtree[parent] += gt_subtree[node_id];
                pred_subtree[parent] += pred_subtree[node_id];
            }
        }

        let mut weighted = 0.0;
        let mut unweighted = 0.0;
        let mut union_unweighted = 0.0;

        for node_id in 1..self.nodes.len() {
            let branch_length = self.nodes[node_id].branch_length;
            if branch_length <= 0.0 {
                continue;
            }

            let gt_mass = gt_subtree[node_id];
            let pred_mass = pred_subtree[node_id];
            let diff = (gt_mass - pred_mass).abs();
            if diff > MASS_EPSILON {
                weighted += branch_length * diff;
            }

            let gt_present = gt_mass > MASS_EPSILON;
            let pred_present = pred_mass > MASS_EPSILON;
            let presence_diff = match (gt_present, pred_present) {
                (true, false) | (false, true) => 1.0,
                _ => 0.0,
            };
            if presence_diff > 0.0 {
                unweighted += branch_length * presence_diff;
            }
            if gt_present || pred_present {
                union_unweighted += branch_length;
            }
        }

        (weighted, unweighted, union_unweighted)
    }
}

struct TaxNode {
    name: String,
    parent: Option<usize>,
    children: Vec<usize>,
    gt_mass: f64,
    pred_mass: f64,
    branch_length: f64,
}

fn abundance_rank_error(
    gt_map: &HashMap<String, f64>,
    pred_map: &HashMap<String, f64>,
) -> Result<f64> {
    let num_gt = gt_map.len();
    let num_pred = pred_map.len();

    if num_gt == 0 && num_pred == 0 {
        return Ok(0.0);
    }

    let mut gt_sorted: Vec<(String, f64)> = gt_map
        .iter()
        .map(|(taxon, value)| (taxon.clone(), *value))
        .collect();
    gt_sorted.sort_by(
        |a, b| match b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal) {
            Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        },
    );

    let mut pred_sorted: Vec<(String, f64)> = pred_map
        .iter()
        .map(|(taxon, value)| (taxon.clone(), *value))
        .collect();
    pred_sorted.sort_by(
        |a, b| match b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal) {
            Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        },
    );

    let mut gt_ranks: HashMap<String, usize> = HashMap::new();
    for (idx, (taxon, _)) in gt_sorted.iter().enumerate() {
        gt_ranks.insert(taxon.clone(), idx + 1);
    }

    let mut pred_ranks: HashMap<String, usize> = HashMap::new();
    for (idx, (taxon, _)) in pred_sorted.iter().enumerate() {
        pred_ranks.insert(taxon.clone(), idx + 1);
    }

    let mut total_error = 0.0;

    if num_gt > 0 {
        let num_gt_f = num_gt as f64;
        for (taxon, &gt_rank) in &gt_ranks {
            if let Some(&pred_rank) = pred_ranks.get(taxon) {
                let weight = (num_gt - gt_rank + 1) as f64 / num_gt_f;
                let rank_diff = if gt_rank > pred_rank {
                    (gt_rank - pred_rank) as f64
                } else {
                    (pred_rank - gt_rank) as f64
                };
                total_error += rank_diff * weight;
            }
        }

        for (taxon, &gt_rank) in &gt_ranks {
            if !pred_ranks.contains_key(taxon) {
                let offset = (num_gt - gt_rank + 1) as f64;
                let weight = offset / num_gt_f;
                total_error += offset * weight;
            }
        }
    }

    if num_pred > 0 {
        let num_pred_f = num_pred as f64;
        for (taxon, &pred_rank) in &pred_ranks {
            if !gt_ranks.contains_key(taxon) {
                let offset = (num_pred - pred_rank + 1) as f64;
                let weight = offset / num_pred_f;
                total_error += offset * weight;
            }
        }
    }

    let mut normalization = 0.0;
    if num_gt > 0 {
        let num_gt_f = num_gt as f64;
        for j in 1..=num_gt {
            let term = (num_gt + 1 - j) as f64;
            normalization += (term * term) / num_gt_f;
        }
    }
    if num_pred > 0 {
        let num_pred_f = num_pred as f64;
        for k in 1..=num_pred {
            let term = (num_pred + 1 - k) as f64;
            normalization += (term * term) / num_pred_f;
        }
    }

    if normalization <= 0.0 {
        return Ok(0.0);
    }

    let ratio = total_error / normalization;
    if ratio < 0.0 {
        Ok(0.0)
    } else if ratio > 1.0 {
        Ok(1.0)
    } else {
        Ok(ratio)
    }
}

fn mass_weighted_abundance_rank_error(
    gt_map: &HashMap<String, f64>,
    pred_map: &HashMap<String, f64>,
    p: f64,
) -> Result<f64> {
    ensure!(p > 0.0, "mARE exponent must be positive");

    let (gt_weights, gt_ranks, gt_count, gt_mass) = prepare_profile(gt_map);
    let (pred_weights, pred_ranks, pred_count, pred_mass) = prepare_profile(pred_map);

    let denominator = gt_mass + pred_mass;
    if denominator <= 0.0 {
        return Ok(0.0);
    }

    let mut union: BTreeSet<String> = BTreeSet::new();
    for taxon in gt_weights.keys() {
        union.insert(taxon.clone());
    }
    for taxon in pred_weights.keys() {
        union.insert(taxon.clone());
    }

    if union.is_empty() {
        return Ok(0.0);
    }

    let max_rank_diff = (gt_count.max(pred_count).saturating_sub(1)) as f64;

    let mut numerator = 0.0;
    for taxon in union {
        match (gt_weights.get(&taxon), pred_weights.get(&taxon)) {
            (Some(&p_gt), Some(&p_pred)) => {
                let rank_penalty = if max_rank_diff > 0.0 {
                    let gt_rank = gt_ranks.get(&taxon).copied().unwrap_or(1);
                    let pred_rank = pred_ranks.get(&taxon).copied().unwrap_or(1);
                    let diff = (gt_rank as i64 - pred_rank as i64).abs() as f64;
                    (diff / max_rank_diff).powf(p)
                } else {
                    0.0
                };
                numerator += rank_penalty * (p_gt + p_pred);
            }
            (Some(&p_gt), None) => {
                numerator += p_gt;
            }
            (None, Some(&p_pred)) => {
                numerator += p_pred;
            }
            (None, None) => {}
        }
    }

    let ratio = numerator / denominator;
    Ok(ratio.max(0.0).min(1.0))
}

fn prepare_profile(
    abundances: &HashMap<String, f64>,
) -> (HashMap<String, f64>, HashMap<String, usize>, usize, f64) {
    let mut entries: Vec<(String, f64)> = abundances
        .iter()
        .map(|(taxon, value)| (taxon.clone(), *value))
        .collect();

    let total: f64 = entries.iter().map(|(_, v)| *v).sum();
    if !total.is_finite() || total <= 0.0 {
        return (HashMap::new(), HashMap::new(), 0, 0.0);
    }

    for (_, value) in entries.iter_mut() {
        *value /= total;
    }

    entries.sort_by(
        |a, b| match b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal) {
            Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        },
    );

    let mut ranks = HashMap::new();
    let mut weights = HashMap::new();
    for (idx, (taxon, value)) in entries.into_iter().enumerate() {
        ranks.insert(taxon.clone(), idx + 1);
        weights.insert(taxon, value);
    }

    let sum_norm = weights.values().copied().sum();
    let count = ranks.len();
    (weights, ranks, count, sum_norm)
}

fn write_metrics<W: Write>(
    writer: &mut W,
    label: &str,
    sample: &str,
    rank: &str,
    metrics: &Metrics,
) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        label,
        sample,
        rank,
        metrics.tp,
        metrics.fp,
        metrics.fn_,
        format_opt(metrics.precision),
        format_opt(metrics.recall),
        format_opt(metrics.f1),
        format_opt(metrics.jaccard),
        format_float(metrics.l1_error),
        format_opt(metrics.bray_curtis),
        format_float(metrics.shannon_pred),
        format_float(metrics.shannon_truth),
        format_opt(metrics.evenness_pred),
        format_opt(metrics.evenness_truth),
        format_opt(metrics.pearson),
        format_opt(metrics.spearman),
        format_opt(metrics.weighted_unifrac),
        format_opt(metrics.unweighted_unifrac),
        format_opt(metrics.abundance_rank_error),
        format_opt(metrics.mass_weighted_abundance_rank_error)
    )?;
    Ok(())
}

fn format_opt(value: Option<f64>) -> String {
    value
        .map(|v| format_float(v))
        .unwrap_or_else(|| "NA".to_string())
}

fn format_float(value: f64) -> String {
    format!("{:.5}", value)
}

#[cfg(test)]
mod tests {
    use super::{
        CANONICAL_RANKS, LineageInfo, ProfileEntry, abundance_rank_error, canonical_rank_index,
        compute_lineage, ensure_superkingdom_only, entry_missing_superkingdom,
        highest_taxid_in_taxpath, mass_weighted_abundance_rank_error, samples_need_superkingdom,
        unifrac, unifrac_components,
    };
    use crate::cami::{Entry, Sample};
    use std::collections::HashMap;
    use std::sync::Arc;

    fn map(entries: &[(&str, f64)]) -> HashMap<String, f64> {
        entries
            .iter()
            .map(|(name, value)| ((*name).to_string(), *value))
            .collect()
    }

    #[test]
    fn entry_missing_superkingdom_detects_empty_paths() {
        let missing = Entry {
            taxid: "1".to_string(),
            rank: "superkingdom".to_string(),
            taxpath: "".to_string(),
            taxpathsn: "".to_string(),
            percentage: 100.0,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        };
        assert!(entry_missing_superkingdom(&missing));

        let present = Entry {
            taxid: "2".to_string(),
            rank: "superkingdom".to_string(),
            taxpath: "2".to_string(),
            taxpathsn: "Bacteria".to_string(),
            percentage: 100.0,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        };
        assert!(!entry_missing_superkingdom(&present));
    }

    #[test]
    fn ensure_superkingdom_only_sets_prefix() {
        let mut entry = Entry {
            taxid: "123".to_string(),
            rank: "phylum".to_string(),
            taxpath: "123".to_string(),
            taxpathsn: "Firmicutes".to_string(),
            percentage: 10.0,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        };
        let mut canonical = vec![None; CANONICAL_RANKS.len()];
        canonical[0] = Some(("2".to_string(), "Bacteria".to_string()));
        canonical[1] = Some(("123".to_string(), "Firmicutes".to_string()));
        let info = LineageInfo {
            canonical: Arc::from(canonical.into_boxed_slice()),
            rank_index: Some(1),
            rank_name: Some("phylum".to_string()),
            deepest_index: 1,
            resolved_taxid: "123".to_string(),
        };

        ensure_superkingdom_only(&mut entry, &info);
        assert_eq!(entry.taxpath, "2|123");
        assert_eq!(entry.taxpathsn, "Bacteria|Firmicutes");
    }

    #[test]
    fn samples_need_superkingdom_flags_entries() {
        let sample_missing = Sample {
            id: "s1".to_string(),
            version: None,
            taxonomy_tag: None,
            ranks: Vec::new(),
            rank_groups: Vec::new(),
            rank_aliases: HashMap::new(),
            entries: vec![Entry {
                taxid: "1".to_string(),
                rank: "superkingdom".to_string(),
                taxpath: "".to_string(),
                taxpathsn: "".to_string(),
                percentage: 100.0,
                cami_genome_id: None,
                cami_otu: None,
                hosts: None,
            }],
        };
        let sample_ok = Sample {
            id: "s2".to_string(),
            version: None,
            taxonomy_tag: None,
            ranks: Vec::new(),
            rank_groups: Vec::new(),
            rank_aliases: HashMap::new(),
            entries: vec![Entry {
                taxid: "2".to_string(),
                rank: "superkingdom".to_string(),
                taxpath: "2".to_string(),
                taxpathsn: "Bacteria".to_string(),
                percentage: 100.0,
                cami_genome_id: None,
                cami_otu: None,
                hosts: None,
            }],
        };

        assert!(samples_need_superkingdom(&[sample_missing.clone()]));
        assert!(!samples_need_superkingdom(&[sample_ok.clone()]));
        assert!(samples_need_superkingdom(&[sample_ok, sample_missing]));
    }

    fn profile_entry_for_test(taxpath: &str, percentage: f64) -> ProfileEntry {
        let parts: Vec<&str> = taxpath.split('|').filter(|s| !s.is_empty()).collect();
        let taxid = parts.last().copied().unwrap_or(taxpath);
        let rank_index = parts
            .len()
            .checked_sub(1)
            .unwrap_or(0)
            .min(CANONICAL_RANKS.len() - 1);
        let rank = CANONICAL_RANKS[rank_index].to_string();
        let entry = Entry {
            taxid: taxid.to_string(),
            rank: rank.clone(),
            taxpath: taxpath.to_string(),
            taxpathsn: taxpath.to_string(),
            percentage,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        };
        let lineage = Arc::from(compute_lineage(&entry, &rank, None, false).into_boxed_slice());
        ProfileEntry {
            taxid: entry.taxid.clone(),
            percentage,
            lineage,
        }
    }

    #[test]
    fn mass_weighted_abundance_rank_error_perfect_match() {
        let gt = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let pred = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 0.0).abs() < 1e-9);
    }

    #[test]
    fn mass_weighted_abundance_rank_error_only_missed_taxa() {
        let gt = map(&[("A", 0.6), ("B", 0.4)]);
        let pred = map(&[]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn unifrac_identical_profiles_are_zero() {
        let rank = "species";
        let gt = vec![
            profile_entry_for_test("sk|p|c|o|f|g|s", 60.0),
            profile_entry_for_test("sk|p|c|o|f|h|i", 40.0),
        ];
        let pred = vec![
            profile_entry_for_test("sk|p|c|o|f|g|s", 60.0),
            profile_entry_for_test("sk|p|c|o|f|h|i", 40.0),
        ];
        let result = unifrac(rank, &gt, &pred);
        assert!(result.weighted.unwrap_or(0.0) < 1e-12);
        assert!(result.unweighted.unwrap_or(0.0) < 1e-12);
    }

    #[test]
    fn unifrac_uses_taxid_lineage_even_when_names_differ() {
        let rank = "species";
        let gt_entry = Entry {
            taxid: "123".to_string(),
            rank: rank.to_string(),
            taxpath: "1|2|123".to_string(),
            taxpathsn: "Root|Clade|Species alpha".to_string(),
            percentage: 100.0,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        };
        let pred_entry = Entry {
            taxid: "123".to_string(),
            rank: rank.to_string(),
            taxpath: "1|2|123".to_string(),
            taxpathsn: "Root|Synonym|Species beta".to_string(),
            percentage: 100.0,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        };
        let gt = vec![ProfileEntry {
            taxid: gt_entry.taxid.clone(),
            percentage: gt_entry.percentage,
            lineage: Arc::from(compute_lineage(&gt_entry, rank, None, false).into_boxed_slice()),
        }];
        let pred = vec![ProfileEntry {
            taxid: pred_entry.taxid.clone(),
            percentage: pred_entry.percentage,
            lineage: Arc::from(compute_lineage(&pred_entry, rank, None, false).into_boxed_slice()),
        }];
        let result = unifrac(rank, &gt, &pred);
        assert!(result.weighted.unwrap_or(1.0) < 1e-12);
        assert!(result.unweighted.unwrap_or(1.0) < 1e-12);
    }

    #[test]
    fn unifrac_strain_singletons_on_opposite_branches_reach_max() {
        let rank = "strain";
        let gt = vec![profile_entry_for_test("a|b|c|d|e|f|g|t1", 100.0)];
        let pred = vec![profile_entry_for_test("z|y|x|w|v|u|s|t2", 100.0)];
        let result = unifrac(rank, &gt, &pred);
        assert!((result.weighted.unwrap() - 1.0).abs() < 1e-9);
        assert!((result.unweighted.unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn unweighted_unifrac_max_follows_opal_formula() {
        let rank = "genus";
        let gt = vec![
            profile_entry_for_test("sk|p|c|o|f|g1", 60.0),
            profile_entry_for_test("sk|p|c|o|f|g2", 40.0),
        ];
        let pred = vec![
            profile_entry_for_test("sk|p|c|o|f|g1", 60.0),
            profile_entry_for_test("sk|p|c|o|f|g3", 40.0),
        ];
        let result = unifrac(rank, &gt, &pred);
        let rank_index = canonical_rank_index(rank).unwrap();
        let components = unifrac_components(rank_index, &gt, &pred).unwrap();
        let gt_baseline = unifrac_components(rank_index, &gt, &[]).unwrap();
        let pred_baseline = unifrac_components(rank_index, &[], &pred).unwrap();
        let shared = vec![profile_entry_for_test("sk|p|c|o|f|g1", 60.0)];
        let shared_components = unifrac_components(rank_index, &shared, &shared).unwrap();
        let union_length = gt_baseline.unweighted_max + pred_baseline.unweighted_max
            - shared_components.unweighted_max;
        let expected_normalized = components.unweighted_raw / union_length;
        assert!((result.unweighted.unwrap() - expected_normalized).abs() < 1e-9);
    }

    #[test]
    fn weighted_unifrac_is_monotonic_with_deeper_mismatch() {
        let rank = "strain";
        let gt = vec![profile_entry_for_test("sk|p|c|o|f|g1|s1|t1", 100.0)];
        let pred_genus = vec![profile_entry_for_test("sk|p|c|o|f|g2", 100.0)];
        let pred_strain = vec![profile_entry_for_test("sk|p|c|o|f|g2|s2|t2", 100.0)];

        let rank_index = canonical_rank_index(rank).unwrap();
        let genus_components = unifrac_components(rank_index, &gt, &pred_genus).unwrap();
        let strain_components = unifrac_components(rank_index, &gt, &pred_strain).unwrap();

        assert!(genus_components.weighted_raw > 0.0);
        assert!(strain_components.weighted_raw >= genus_components.weighted_raw - 1e-9);
        assert!(strain_components.weighted_raw <= genus_components.max_weighted + 1e-9);
    }

    #[test]
    fn mass_weighted_abundance_rank_error_only_false_positives() {
        let gt = map(&[]);
        let pred = map(&[("X", 0.7), ("Y", 0.3)]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn mass_weighted_abundance_rank_error_mixed_case() {
        let gt = map(&[("A", 0.6), ("B", 0.3), ("C", 0.1)]);
        let pred = map(&[("A", 0.55), ("C", 0.35), ("D", 0.10)]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 0.3125).abs() < 1e-9);
    }

    #[test]
    fn mass_weighted_abundance_rank_error_rank_reversal() {
        let gt = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let pred = map(&[("A", 0.2), ("B", 0.3), ("C", 0.5)]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 0.7).abs() < 1e-9);
    }

    #[test]
    fn abundance_rank_error_perfect_match() {
        let gt = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let pred = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!((score - 0.0).abs() < 1e-9);
    }

    #[test]
    fn abundance_rank_error_rank_swap() {
        let gt = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let pred = map(&[("B", 0.5), ("A", 0.3), ("C", 0.2)]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!(score > 0.0 && score < 1.0);
    }

    #[test]
    fn abundance_rank_error_only_missed_taxa() {
        let gt = map(&[("A", 0.6), ("B", 0.4)]);
        let pred = map(&[]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn abundance_rank_error_only_false_positives() {
        let gt = map(&[]);
        let pred = map(&[("X", 0.7), ("Y", 0.3)]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn abundance_rank_error_mixed_case() {
        let gt = map(&[("A", 0.6), ("B", 0.3), ("C", 0.1)]);
        let pred = map(&[("A", 0.55), ("C", 0.35), ("D", 0.10)]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!(score > 0.0 && score < 1.0);
    }

    #[test]
    fn highest_taxid_prefers_first_taxpath_entry() {
        let entry = Entry {
            taxid: "2956277".to_string(),
            rank: "species".to_string(),
            taxpath: "2731618|2731619|3424659|1198980|2956277".to_string(),
            taxpathsn: "".to_string(),
            percentage: 1.0,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        };

        assert_eq!(highest_taxid_in_taxpath(&entry, None), Some(2_731_618));
    }

    #[test]
    fn highest_taxid_falls_back_to_entry_taxid() {
        let entry = Entry {
            taxid: "12345".to_string(),
            rank: "species".to_string(),
            taxpath: "||abc||".to_string(),
            taxpathsn: "".to_string(),
            percentage: 1.0,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        };

        assert_eq!(highest_taxid_in_taxpath(&entry, None), Some(12_345));
    }
}
