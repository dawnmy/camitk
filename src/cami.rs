use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct Sample {
    pub id: String,
    pub version: Option<String>,
    pub taxonomy_tag: Option<String>,
    pub ranks: Vec<String>,
    pub rank_groups: Vec<Vec<String>>,
    pub rank_aliases: HashMap<String, usize>,
    pub entries: Vec<Entry>,
}

#[derive(Debug, Clone)]
pub struct Entry {
    pub taxid: String,
    pub rank: String,
    pub taxpath: String,
    pub taxpathsn: String,
    pub percentage: f64,
    pub cami_genome_id: Option<String>,
    pub cami_otu: Option<String>,
    pub hosts: Option<String>,
}

impl Sample {
    pub fn rank_index(&self, rank: &str) -> Option<usize> {
        self.rank_aliases.get(rank).copied()
    }

    pub fn header_rank_tokens(&self) -> Vec<String> {
        self.rank_groups
            .iter()
            .map(|group| group.join(","))
            .collect()
    }

    pub fn is_modern_format(&self) -> bool {
        self.rank_groups.iter().any(|g| g.len() > 1)
    }

    pub fn set_rank_groups(&mut self, groups: Vec<Vec<String>>) {
        self.rank_groups = groups
            .into_iter()
            .map(|group| {
                group
                    .into_iter()
                    .map(|name| name.trim().to_string())
                    .filter(|name| !name.is_empty())
                    .collect::<Vec<_>>()
            })
            .filter(|group| !group.is_empty())
            .collect();
        self.ranks = self
            .rank_groups
            .iter()
            .map(|group| group.first().cloned().unwrap_or_default())
            .collect();
        self.rank_aliases.clear();
        for (idx, group) in self.rank_groups.iter().enumerate() {
            for name in group {
                self.rank_aliases.insert(name.clone(), idx);
            }
        }
    }
}

fn parse_rank_groups(spec: &str) -> Vec<Vec<String>> {
    spec.split('|')
        .map(|part| {
            let trimmed = part.trim();
            if let Some(stripped) = trimmed.strip_prefix('{') {
                if let Some(inner) = stripped.strip_suffix('}') {
                    return inner
                        .split(',')
                        .map(|name| name.trim().to_string())
                        .filter(|name| !name.is_empty())
                        .collect();
                }
            }
            let names: Vec<String> = trimmed
                .split(',')
                .map(|name| name.trim().to_string())
                .filter(|name| !name.is_empty())
                .collect();
            if names.is_empty() {
                Vec::new()
            } else if names.len() == 1 {
                vec![trimmed.to_string()]
            } else {
                names
            }
        })
        .collect()
}

pub fn parse_cami_reader<R: BufRead>(reader: R) -> Result<Vec<Sample>> {
    let mut samples = Vec::new();
    let mut current: Option<Sample> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with("@SampleID:") {
            if let Some(s) = current.take() {
                samples.push(s);
            }
            let id = line[10..].trim().to_string();
            current = Some(Sample {
                id,
                version: None,
                taxonomy_tag: None,
                ranks: Vec::new(),
                rank_groups: Vec::new(),
                rank_aliases: HashMap::new(),
                entries: Vec::new(),
            });
            continue;
        }
        if line.starts_with("@Version:") {
            if let Some(s) = current.as_mut() {
                s.version = Some(line[9..].trim().to_string());
            }
            continue;
        }
        if line.starts_with("@TaxonomyID:") {
            if let Some(s) = current.as_mut() {
                s.taxonomy_tag = Some(line[12..].trim().to_string());
            }
            continue;
        }
        if line.starts_with("@Ranks:") {
            if let Some(s) = current.as_mut() {
                let groups = parse_rank_groups(line[7..].trim());
                s.set_rank_groups(groups);
            }
            continue;
        }
        if line.starts_with("@@TAXID") {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let (taxid, rank, taxpath, taxpathsn, percentage, cami_genome_id, cami_otu, hosts) =
            if fields.len() >= 5 {
                let genome_id = fields.get(5).map(|v| v.trim()).filter(|v| !v.is_empty());
                let otu = fields.get(6).map(|v| v.trim()).filter(|v| !v.is_empty());
                let hosts = fields.get(7).map(|v| v.trim()).filter(|v| !v.is_empty());
                (
                    fields[0], fields[1], fields[2], fields[3], fields[4], genome_id, otu, hosts,
                )
            } else {
                let fields_ws: Vec<&str> = line.split_whitespace().collect();
                if fields_ws.len() < 5 {
                    continue;
                }
                (
                    fields_ws[0],
                    fields_ws[1],
                    fields_ws[2],
                    fields_ws[3],
                    fields_ws[4],
                    fields_ws.get(5).map(|v| v.trim()).filter(|v| !v.is_empty()),
                    fields_ws.get(6).map(|v| v.trim()).filter(|v| !v.is_empty()),
                    fields_ws.get(7).map(|v| v.trim()).filter(|v| !v.is_empty()),
                )
            };

        if let Some(s) = current.as_mut() {
            let entry = Entry {
                taxid: taxid.to_string(),
                rank: rank.to_string(),
                taxpath: taxpath.to_string(),
                taxpathsn: taxpathsn.to_string(),
                percentage: percentage.parse().unwrap_or(0.0),
                cami_genome_id: cami_genome_id.map(|v| v.to_string()),
                cami_otu: cami_otu.map(|v| v.to_string()),
                hosts: hosts.map(|v| v.to_string()),
            };
            s.entries.push(entry);
        }
    }

    if let Some(s) = current.take() {
        samples.push(s);
    }
    Ok(samples)
}

pub fn parse_cami(path: &Path) -> Result<Vec<Sample>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    parse_cami_reader(reader)
}

pub fn load_samples(input: Option<&PathBuf>) -> Result<Vec<Sample>> {
    match input {
        Some(path) if path != Path::new("-") => parse_cami(path),
        _ => {
            let stdin = io::stdin();
            let handle = stdin.lock();
            parse_cami_reader(handle)
        }
    }
}

pub fn write_cami(samples: &[Sample], out: &mut dyn Write) -> Result<()> {
    for s in samples {
        writeln!(out, "@SampleID:{}", s.id)?;
        if let Some(v) = &s.version {
            writeln!(out, "@Version:{}", v)?;
        }
        let rank_tokens = s.header_rank_tokens();
        if !rank_tokens.is_empty() {
            writeln!(out, "@Ranks:{}", rank_tokens.join("|"))?;
        }
        if let Some(tag) = &s.taxonomy_tag {
            writeln!(out, "@TaxonomyID:{}", tag)?;
        }
        let extended = s.is_modern_format()
            || s.entries
                .iter()
                .any(|e| e.cami_genome_id.is_some() || e.cami_otu.is_some() || e.hosts.is_some());
        if extended {
            writeln!(
                out,
                "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\tHOSTS"
            )?;
        } else {
            writeln!(out, "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE")?;
        }
        for e in &s.entries {
            if extended {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    e.taxid,
                    e.rank,
                    e.taxpath,
                    e.taxpathsn,
                    e.percentage,
                    e.cami_genome_id.as_deref().unwrap_or(""),
                    e.cami_otu.as_deref().unwrap_or(""),
                    e.hosts.as_deref().unwrap_or("")
                )?;
            } else {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    e.taxid, e.rank, e.taxpath, e.taxpathsn, e.percentage
                )?;
            }
        }
    }
    Ok(())
}

pub fn open_output(path: Option<&PathBuf>) -> Result<Box<dyn Write>> {
    let writer: Box<dyn Write> = match path {
        Some(p) if p != Path::new("-") => Box::new(File::create(p)?),
        _ => Box::new(io::stdout()),
    };
    Ok(writer)
}
