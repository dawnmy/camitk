use crate::cami::{Entry, Sample};
use crate::taxonomy::Taxonomy;
use anyhow::{Result, anyhow, bail};
use regex::Regex;
use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap};

#[derive(Debug, Clone)]
pub enum Expr {
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Atom(String),
}

pub fn parse_expression(s: &str) -> Result<Expr> {
    let normalized = s.replace("&&", "&").replace("||", "|");
    let chars: Vec<char> = normalized.chars().collect();
    let mut pos = 0;

    fn skip_ws(chars: &[char], pos: &mut usize) {
        while *pos < chars.len() && chars[*pos].is_whitespace() {
            *pos += 1;
        }
    }

    fn parse_primary(chars: &[char], pos: &mut usize) -> Result<Expr> {
        skip_ws(chars, pos);
        if *pos >= chars.len() {
            return Err(anyhow!("unexpected end of expression"));
        }
        if chars[*pos] == '(' {
            *pos += 1;
            let expr = parse_expr(chars, pos)?;
            skip_ws(chars, pos);
            if *pos >= chars.len() || chars[*pos] != ')' {
                return Err(anyhow!("missing closing parenthesis"));
            }
            *pos += 1;
            Ok(expr)
        } else {
            let start = *pos;
            while *pos < chars.len()
                && chars[*pos] != '&'
                && chars[*pos] != '|'
                && chars[*pos] != ')'
            {
                *pos += 1;
            }
            let atom: String = chars[start..*pos].iter().collect();
            Ok(Expr::Atom(atom.trim().to_string()))
        }
    }

    fn parse_term(chars: &[char], pos: &mut usize) -> Result<Expr> {
        let mut left = parse_primary(chars, pos)?;
        loop {
            skip_ws(chars, pos);
            if *pos + 1 < chars.len() && chars[*pos] == '&' && chars[*pos + 1] == '&' {
                *pos += 1;
            }
            if *pos >= chars.len() || chars[*pos] != '&' {
                break;
            }
            *pos += 1;
            let right = parse_primary(chars, pos)?;
            left = Expr::And(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn parse_expr(chars: &[char], pos: &mut usize) -> Result<Expr> {
        let mut left = parse_term(chars, pos)?;
        loop {
            skip_ws(chars, pos);
            if *pos + 1 < chars.len() && chars[*pos] == '|' && chars[*pos + 1] == '|' {
                *pos += 1;
            }
            if *pos >= chars.len() || chars[*pos] != '|' {
                break;
            }
            *pos += 1;
            let right = parse_term(chars, pos)?;
            left = Expr::Or(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    let expr = parse_expr(&chars, &mut pos)?;
    skip_ws(&chars, &mut pos);
    if pos != chars.len() {
        return Err(anyhow!("unexpected characters after expression"));
    }
    Ok(expr)
}

pub fn expr_needs_taxdump(expr: &Expr) -> bool {
    match expr {
        Expr::And(a, b) | Expr::Or(a, b) => expr_needs_taxdump(a) || expr_needs_taxdump(b),
        Expr::Atom(s) => {
            let trimmed = s.trim().trim_start_matches('!').trim_start();
            trimmed.starts_with("tax") || trimmed.starts_with('t')
        }
    }
}

pub fn apply_filter(samples: &[Sample], expr: &Expr, taxonomy: Option<&Taxonomy>) -> Vec<Sample> {
    let cumsums = compute_cumsums(samples);
    let mut filtered = Vec::new();
    for (sample_index, sample) in samples.iter().enumerate() {
        let mut ns = sample.clone();
        ns.entries.clear();
        for (entry_index, entry) in sample.entries.iter().enumerate() {
            if eval_expr(
                expr,
                samples,
                sample,
                entry,
                sample_index,
                entry_index,
                taxonomy,
                &cumsums,
            ) {
                ns.entries.push(entry.clone());
            }
        }
        if !ns.entries.is_empty() {
            filtered.push(ns);
        }
    }
    filtered
}

fn compute_cumsums(samples: &[Sample]) -> Vec<Vec<f64>> {
    samples
        .iter()
        .map(|sample| {
            let mut by_rank: HashMap<&str, Vec<(usize, f64)>> = HashMap::new();
            for (idx, entry) in sample.entries.iter().enumerate() {
                by_rank
                    .entry(entry.rank.as_str())
                    .or_default()
                    .push((idx, entry.percentage));
            }

            let mut cumsums = vec![0.0; sample.entries.len()];
            for entries in by_rank.values_mut() {
                entries.sort_by(|a, b| {
                    a.1.partial_cmp(&b.1)
                        .unwrap_or(Ordering::Equal)
                        .then_with(|| a.0.cmp(&b.0))
                });
                let mut running = 0.0;
                for (idx, value) in entries.iter() {
                    running += value;
                    cumsums[*idx] = running;
                }
            }
            cumsums
        })
        .collect()
}

fn eval_expr(
    e: &Expr,
    samples: &[Sample],
    sample: &Sample,
    entry: &Entry,
    sample_index: usize,
    entry_index: usize,
    taxonomy: Option<&Taxonomy>,
    cumsums: &[Vec<f64>],
) -> bool {
    match e {
        Expr::And(a, b) => {
            eval_expr(
                a,
                samples,
                sample,
                entry,
                sample_index,
                entry_index,
                taxonomy,
                cumsums,
            ) && eval_expr(
                b,
                samples,
                sample,
                entry,
                sample_index,
                entry_index,
                taxonomy,
                cumsums,
            )
        }
        Expr::Or(a, b) => {
            eval_expr(
                a,
                samples,
                sample,
                entry,
                sample_index,
                entry_index,
                taxonomy,
                cumsums,
            ) || eval_expr(
                b,
                samples,
                sample,
                entry,
                sample_index,
                entry_index,
                taxonomy,
                cumsums,
            )
        }
        Expr::Atom(s) => eval_atom(
            s,
            samples,
            sample,
            entry,
            sample_index,
            entry_index,
            taxonomy,
            cumsums,
        ),
    }
}

fn eval_atom(
    atom: &str,
    samples: &[Sample],
    sample: &Sample,
    entry: &Entry,
    sample_index: usize,
    entry_index: usize,
    taxonomy: Option<&Taxonomy>,
    cumsums: &[Vec<f64>],
) -> bool {
    let atom = atom.trim();

    if let Some(res) = eval_rank(atom, sample, entry) {
        return res;
    }
    if let Some(res) = eval_sample(atom, samples, sample, sample_index) {
        return res;
    }
    if let Some(res) = eval_abundance(atom, entry) {
        return res;
    }
    if let Some(res) = eval_tax(atom, entry, taxonomy) {
        return res;
    }
    if let Some(res) = eval_cumsum(atom, sample_index, entry_index, cumsums) {
        return res;
    }

    false
}

fn eval_rank(atom: &str, sample: &Sample, entry: &Entry) -> Option<bool> {
    let rest = if let Some(r) = atom.strip_prefix("rank") {
        r
    } else if let Some(r) = atom.strip_prefix('r') {
        r
    } else {
        return None;
    };
    let rest = rest.trim_start();
    if let Some(v) = rest.strip_prefix("!=") {
        let values = parse_rank_list(v);
        if values.is_empty() {
            let value = strip_quotes(v.trim());
            return Some(entry.rank != value);
        }
        return Some(values.iter().all(|value| entry.rank != *value));
    }
    if let Some(v) = rest.strip_prefix("==") {
        let values = parse_rank_list(v);
        if values.is_empty() {
            let value = strip_quotes(v.trim());
            return Some(entry.rank == value);
        }
        return Some(values.iter().any(|value| entry.rank == *value));
    }
    if let Some(v) = rest.strip_prefix("<=") {
        let value = strip_quotes(v.trim());
        return Some(rank_compare(sample, &entry.rank, value.as_str(), |a, b| {
            a >= b
        }));
    }
    if let Some(v) = rest.strip_prefix("<") {
        let value = strip_quotes(v.trim());
        return Some(rank_compare(sample, &entry.rank, value.as_str(), |a, b| {
            a > b
        }));
    }
    if let Some(v) = rest.strip_prefix(">=") {
        let value = strip_quotes(v.trim());
        return Some(rank_compare(sample, &entry.rank, value.as_str(), |a, b| {
            a <= b
        }));
    }
    if let Some(v) = rest.strip_prefix('>') {
        let value = strip_quotes(v.trim());
        return Some(rank_compare(sample, &entry.rank, value.as_str(), |a, b| {
            a < b
        }));
    }
    None
}

fn rank_compare<F>(sample: &Sample, entry_rank: &str, other: &str, cmp: F) -> bool
where
    F: Fn(usize, usize) -> bool,
{
    let Some(entry_idx) = sample.rank_index(entry_rank) else {
        return false;
    };
    let Some(other_idx) = sample.rank_index(other) else {
        return false;
    };
    cmp(entry_idx, other_idx)
}

fn eval_sample(
    atom: &str,
    samples: &[Sample],
    sample: &Sample,
    sample_index: usize,
) -> Option<bool> {
    let rest = if let Some(r) = atom.strip_prefix("sample") {
        r
    } else if let Some(r) = atom.strip_prefix('s') {
        r
    } else {
        return None;
    };
    let rest = rest.trim_start();
    if let Some(v) = rest.strip_prefix("==") {
        return Some(match_sample(v.trim(), samples, sample, sample_index));
    }
    if let Some(v) = rest.strip_prefix("!=") {
        return Some(!match_sample(v.trim(), samples, sample, sample_index));
    }
    if let Some(v) = rest.strip_prefix('~') {
        return Some(match_sample_regex(v.trim(), sample));
    }
    None
}

fn match_sample(selector: &str, samples: &[Sample], sample: &Sample, sample_index: usize) -> bool {
    if selector.is_empty() || selector == ":" {
        return true;
    }
    let mut matched = false;
    for part in split_unquoted(selector, ',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if part == ":" || part == "." {
            matched = true;
            break;
        }
        let colon_parts = split_unquoted(part, ':');
        if colon_parts.len() == 2 {
            if match_sample_range_parts(&colon_parts[0], &colon_parts[1], samples, sample_index) {
                matched = true;
                break;
            }
        } else if colon_parts.len() == 1 {
            if match_sample_value(&colon_parts[0], samples, sample, sample_index) {
                matched = true;
                break;
            }
        } else if match_sample_value(part, samples, sample, sample_index) {
            matched = true;
            break;
        }
    }
    matched
}

fn match_sample_regex(pattern_raw: &str, sample: &Sample) -> bool {
    let Some((pattern, _quoted)) = parse_selector_value(pattern_raw) else {
        return false;
    };
    if pattern.is_empty() {
        return false;
    }
    Regex::new(&pattern)
        .map(|re| re.is_match(&sample.id))
        .unwrap_or(false)
}

fn match_sample_value(
    value: &str,
    samples: &[Sample],
    sample: &Sample,
    sample_index: usize,
) -> bool {
    let Some((parsed, quoted)) = parse_selector_value(value) else {
        return true;
    };
    if !quoted {
        if let Ok(idx) = parsed.parse::<usize>() {
            return idx == sample_index + 1;
        }
    }
    if let Some(idx) = samples.iter().position(|s| s.id == parsed) {
        return idx == sample_index;
    }
    parsed == sample.id
}

fn match_sample_range_parts(
    start_raw: &str,
    end_raw: &str,
    samples: &[Sample],
    sample_index: usize,
) -> bool {
    let start_val = parse_selector_value(start_raw);
    let end_val = parse_selector_value(end_raw);

    if start_val.is_none() && end_val.is_none() {
        return true;
    }

    let start_idx = match start_val {
        Some(ref value) => selector_value_to_index(value, samples),
        None => Some(0),
    };
    if start_val.is_some() && start_idx.is_none() {
        return false;
    }

    let end_idx = match end_val {
        Some(ref value) => selector_value_to_index(value, samples),
        None => samples.len().checked_sub(1),
    };
    if end_val.is_some() && end_idx.is_none() {
        return false;
    }

    let (Some(start), Some(end)) = (start_idx, end_idx) else {
        return false;
    };

    let (low, high) = if start <= end {
        (start, end)
    } else {
        (end, start)
    };
    (low..=high).contains(&sample_index)
}

fn eval_abundance(atom: &str, entry: &Entry) -> Option<bool> {
    let rest = if let Some(r) = atom.strip_prefix("abundance") {
        r
    } else if let Some(r) = atom.strip_prefix('a') {
        r
    } else {
        return None;
    };
    let rest = rest.trim_start();
    if let Some(v) = rest.strip_prefix("<=") {
        return Some(entry.percentage <= parse_f64(v.trim()));
    }
    if let Some(v) = rest.strip_prefix(">=") {
        return Some(entry.percentage >= parse_f64(v.trim()));
    }
    if let Some(v) = rest.strip_prefix("==") {
        let target = parse_f64(v.trim());
        return Some((entry.percentage - target).abs() < 1e-9);
    }
    if let Some(v) = rest.strip_prefix("!=") {
        let target = parse_f64(v.trim());
        return Some((entry.percentage - target).abs() >= 1e-9);
    }
    if let Some(v) = rest.strip_prefix('>') {
        return Some(entry.percentage > parse_f64(v.trim()));
    }
    if let Some(v) = rest.strip_prefix('<') {
        return Some(entry.percentage < parse_f64(v.trim()));
    }
    None
}

fn eval_cumsum(
    atom: &str,
    sample_index: usize,
    entry_index: usize,
    cumsums: &[Vec<f64>],
) -> Option<bool> {
    let mut working = atom.trim();
    let mut negate = false;
    while let Some(rest) = working.strip_prefix('!') {
        negate = !negate;
        working = rest.trim_start();
    }

    let rest = if let Some(r) = working.strip_prefix("cumsum") {
        r
    } else if let Some(r) = working.strip_prefix('c') {
        r
    } else {
        return None;
    };

    let rest = rest.trim_start();
    if let Some(v) = rest.strip_prefix("<=") {
        let threshold = parse_f64(v.trim());
        if let Some(value) = cumsums
            .get(sample_index)
            .and_then(|v| v.get(entry_index))
            .copied()
        {
            let result = value <= threshold + 1e-12;
            return Some(if negate { !result } else { result });
        }
        return Some(false);
    }

    Some(false)
}

fn split_unquoted(input: &str, delimiter: char) -> Vec<String> {
    let mut parts = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;
    let mut quote_char = '\0';
    for ch in input.chars() {
        match ch {
            '"' | '\'' => {
                if in_quotes {
                    if ch == quote_char {
                        in_quotes = false;
                    }
                } else {
                    in_quotes = true;
                    quote_char = ch;
                }
                current.push(ch);
            }
            _ if ch == delimiter && !in_quotes => {
                parts.push(current.trim().to_string());
                current.clear();
            }
            _ => current.push(ch),
        }
    }
    if !current.is_empty() {
        parts.push(current.trim().to_string());
    } else if input.ends_with(delimiter) {
        parts.push(String::new());
    }
    if parts.is_empty() {
        parts.push(String::new());
    }
    parts
}

fn strip_quotes(value: &str) -> String {
    let trimmed = value.trim();
    if trimmed.len() >= 2
        && ((trimmed.starts_with('"') && trimmed.ends_with('"'))
            || (trimmed.starts_with('\'') && trimmed.ends_with('\'')))
    {
        trimmed[1..trimmed.len() - 1].to_string()
    } else {
        trimmed.to_string()
    }
}

fn parse_rank_list(raw: &str) -> Vec<String> {
    split_unquoted(raw, ',')
        .into_iter()
        .map(|s| strip_quotes(&s))
        .filter(|s| !s.is_empty())
        .collect()
}

pub fn validate_rank_selectors(expr: &Expr, samples: &[Sample]) -> Result<()> {
    let mut ranks = BTreeSet::new();
    for sample in samples {
        for group in &sample.rank_groups {
            for rank in group {
                ranks.insert(rank.clone());
            }
        }
    }

    if ranks.is_empty() {
        return Ok(());
    }

    validate_rank_expr(expr, &ranks)
}

fn validate_rank_expr(expr: &Expr, ranks: &BTreeSet<String>) -> Result<()> {
    match expr {
        Expr::And(a, b) | Expr::Or(a, b) => {
            validate_rank_expr(a, ranks)?;
            validate_rank_expr(b, ranks)?;
            Ok(())
        }
        Expr::Atom(atom) => validate_rank_atom(atom, ranks),
    }
}

fn validate_rank_atom(atom: &str, ranks: &BTreeSet<String>) -> Result<()> {
    if let Some(values) = extract_rank_values(atom) {
        for value in values {
            if !value.is_empty() && !ranks.contains(&value) {
                let available = ranks.iter().cloned().collect::<Vec<_>>().join(", ");
                bail!(
                    "rank selector references '{}' which is not present in the @Ranks header; available ranks: {}",
                    value,
                    available
                );
            }
        }
    }
    Ok(())
}

fn extract_rank_values(atom: &str) -> Option<Vec<String>> {
    let rest = if let Some(r) = atom.strip_prefix("rank") {
        r
    } else if let Some(r) = atom.strip_prefix('r') {
        r
    } else {
        return None;
    };

    let rest = rest.trim_start();
    if let Some(v) = rest.strip_prefix("==") {
        let values = parse_rank_list(v);
        if values.is_empty() {
            let single = strip_quotes(v.trim());
            if single.is_empty() {
                return Some(Vec::new());
            }
            return Some(vec![single]);
        }
        return Some(values);
    }
    if let Some(v) = rest.strip_prefix("!=") {
        let values = parse_rank_list(v);
        if values.is_empty() {
            let single = strip_quotes(v.trim());
            if single.is_empty() {
                return Some(Vec::new());
            }
            return Some(vec![single]);
        }
        return Some(values);
    }
    for prefix in ["<=", "<", ">=", ">"] {
        if let Some(v) = rest.strip_prefix(prefix) {
            let value = strip_quotes(v.trim());
            if value.is_empty() {
                return Some(Vec::new());
            }
            return Some(vec![value]);
        }
    }

    None
}

fn parse_selector_value(raw: &str) -> Option<(String, bool)> {
    let trimmed = raw.trim();
    if trimmed.is_empty() || trimmed == "." {
        return None;
    }
    let first = trimmed.chars().next();
    let last = trimmed.chars().last();
    if trimmed.len() >= 2 {
        if let (Some(f), Some(l)) = (first, last) {
            if (f == '"' && l == '"') || (f == '\'' && l == '\'') {
                return Some((trimmed[1..trimmed.len() - 1].to_string(), true));
            }
        }
    }
    Some((trimmed.to_string(), false))
}

fn selector_value_to_index(value: &(String, bool), samples: &[Sample]) -> Option<usize> {
    if !value.1 {
        if let Ok(idx) = value.0.parse::<usize>() {
            return Some(idx.saturating_sub(1));
        }
    }
    samples.iter().position(|s| s.id == value.0)
}

fn parse_f64(s: &str) -> f64 {
    s.parse().unwrap_or(0.0)
}

fn eval_tax(atom: &str, entry: &Entry, taxonomy: Option<&Taxonomy>) -> Option<bool> {
    let mut working = atom.trim();
    let mut negate = false;
    if let Some(rest) = working.strip_prefix('!') {
        negate = true;
        working = rest.trim_start();
    }

    let rest = if let Some(r) = working.strip_prefix("tax") {
        r
    } else if let Some(r) = working.strip_prefix('t') {
        r
    } else {
        return None;
    };

    let rest = rest.trim_start();
    let mut invert_op = false;
    let (op, value) = if let Some(v) = rest.strip_prefix("!=") {
        invert_op = true;
        ("==", v.trim())
    } else if let Some(v) = rest.strip_prefix("==") {
        ("==", v.trim())
    } else if let Some(v) = rest.strip_prefix("<=") {
        ("<=", v.trim())
    } else if let Some(v) = rest.strip_prefix('<') {
        ("<", v.trim())
    } else {
        return Some(false);
    };

    let result = if let Some(taxonomy) = taxonomy {
        eval_taxonomy(entry, taxonomy, op, value)
    } else {
        eval_taxpath(entry, op, value)
    };
    let result = if invert_op { !result } else { result };
    Some(if negate { !result } else { result })
}

fn eval_taxonomy(entry: &Entry, taxonomy: &Taxonomy, op: &str, value: &str) -> bool {
    let entry_taxid = taxonomy.resolve_taxid_str(&entry.taxid);
    let target_taxid = taxonomy.resolve_taxid_str(value);

    match op {
        "==" => match (entry_taxid, target_taxid) {
            (Some(e), Some(t)) => e == t,
            _ => entry.taxid == value,
        },
        "<=" => {
            if let (Some(e), Some(t)) = (entry_taxid, target_taxid) {
                if e == t {
                    return true;
                }
                let ancestors = taxonomy.ancestors_of(e);
                ancestors.iter().any(|ancestor| *ancestor == t)
            } else {
                entry.taxpath.split('|').any(|tid| tid == value)
            }
        }
        "<" => {
            if let (Some(e), Some(t)) = (entry_taxid, target_taxid) {
                if e == t {
                    return false;
                }
                let ancestors = taxonomy.ancestors_of(e);
                ancestors.iter().any(|ancestor| *ancestor == t)
            } else {
                entry
                    .taxpath
                    .split('|')
                    .any(|tid| tid == value && tid != entry.taxid)
            }
        }
        _ => false,
    }
}

fn eval_taxpath(entry: &Entry, op: &str, value: &str) -> bool {
    match op {
        "==" => entry
            .taxpath
            .split('|')
            .last()
            .map(|tid| tid == value)
            .unwrap_or(false),
        "<=" => entry.taxpath.split('|').any(|tid| tid == value),
        "<" => {
            let mut parts: Vec<&str> = entry.taxpath.split('|').collect();
            if let Some(last) = parts.pop() {
                parts.contains(&value) && last != value
            } else {
                false
            }
        }
        _ => false,
    }
}
