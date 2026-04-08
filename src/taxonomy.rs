use anyhow::{Context, Result, anyhow};
use flate2::read::GzDecoder;
use reqwest::blocking::get;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};
use std::thread;
use tar::Archive;

#[derive(Debug, Clone)]
pub struct TaxNode {
    pub parent: u32,
    pub rank: String,
}

type LineageEntry = (u32, String, String);

#[derive(Debug)]
pub struct Taxonomy {
    nodes: HashMap<u32, TaxNode>,
    names: HashMap<u32, String>,
    merged: HashMap<u32, u32>,
    deleted: HashSet<u32>,
    ancestors: RwLock<HashMap<u32, Arc<[u32]>>>,
    lineages: RwLock<HashMap<u32, Arc<[LineageEntry]>>>,
    domains: RwLock<HashMap<u32, Option<String>>>,
    realms: RwLock<HashMap<u32, Option<String>>>,
    is_modern: bool,
}

impl Taxonomy {
    pub fn load(dir: &Path) -> Result<Self> {
        let nodes_path = dir.join("nodes.dmp");
        let names_path = dir.join("names.dmp");
        let merged_path = dir.join("merged.dmp");
        let delnodes_path = dir.join("delnodes.dmp");

        let nodes_handle = {
            let path = nodes_path.clone();
            thread::spawn(move || parse_nodes(&path))
        };
        let names_handle = {
            let path = names_path.clone();
            thread::spawn(move || parse_names(&path))
        };

        let merged_handle = {
            let path = merged_path.clone();
            thread::spawn(move || parse_merged(&path))
        };

        let deleted_handle = {
            let path = delnodes_path.clone();
            thread::spawn(move || parse_delnodes(&path))
        };

        let nodes = nodes_handle
            .join()
            .map_err(|_| anyhow!("failed to parse nodes.dmp"))??;
        let names = names_handle
            .join()
            .map_err(|_| anyhow!("failed to parse names.dmp"))??;
        let merged = merged_handle
            .join()
            .map_err(|_| anyhow!("failed to parse merged.dmp"))??;
        let deleted = deleted_handle
            .join()
            .map_err(|_| anyhow!("failed to parse delnodes.dmp"))??;

        let is_modern = nodes.values().any(|node| {
            matches!(
                node.rank.to_ascii_lowercase().as_str(),
                "acellular root" | "cellular root" | "realm"
            )
        });

        Ok(Self {
            nodes,
            names,
            merged,
            deleted,
            ancestors: RwLock::new(HashMap::new()),
            lineages: RwLock::new(HashMap::new()),
            domains: RwLock::new(HashMap::new()),
            realms: RwLock::new(HashMap::new()),
            is_modern,
        })
    }

    pub fn ancestors_of(&self, taxid: u32) -> Arc<[u32]> {
        if let Some(cached) = self.ancestors.read().unwrap().get(&taxid) {
            return cached.clone();
        }

        let mut lineage = Vec::new();
        let mut current = taxid;
        let mut seen = HashSet::new();
        while let Some(node) = self.nodes.get(&current) {
            if node.parent == current || seen.contains(&current) {
                break;
            }
            seen.insert(current);
            lineage.push(node.parent);
            current = node.parent;
        }

        let cached: Arc<[u32]> = lineage.into_boxed_slice().into();
        self.ancestors
            .write()
            .unwrap()
            .insert(taxid, cached.clone());
        cached
    }

    pub fn lineage(&self, taxid: u32) -> Arc<[LineageEntry]> {
        if let Some(cached) = self.lineages.read().unwrap().get(&taxid) {
            return cached.clone();
        }
        let mut stack = Vec::new();
        let mut current = Some(taxid);
        let mut visited = HashSet::new();
        while let Some(tid) = current {
            if visited.contains(&tid) {
                break;
            }
            visited.insert(tid);
            let rank = self
                .nodes
                .get(&tid)
                .map(|n| n.rank.clone())
                .unwrap_or_else(|| "no_rank".to_string());
            let name = self
                .names
                .get(&tid)
                .cloned()
                .unwrap_or_else(|| tid.to_string());
            stack.push((tid, rank, name));
            current = self.nodes.get(&tid).and_then(|n| {
                if n.parent == tid {
                    None
                } else {
                    Some(n.parent)
                }
            });
        }
        stack.reverse();
        let cached: Arc<[LineageEntry]> = stack.into_boxed_slice().into();
        self.lineages.write().unwrap().insert(taxid, cached.clone());
        cached
    }

    pub fn name_of(&self, taxid: u32) -> Option<String> {
        self.names.get(&taxid).cloned()
    }

    pub fn rank_of(&self, taxid: u32) -> Option<String> {
        self.nodes.get(&taxid).map(|n| n.rank.clone())
    }

    pub fn domain_of(&self, taxid: u32) -> Option<String> {
        if let Some(cached) = self.domains.read().unwrap().get(&taxid) {
            return cached.clone();
        }

        let lineage = self.lineage(taxid);
        let result = lineage.iter().find_map(|(_, rank, name)| {
            if rank.eq_ignore_ascii_case("superkingdom")
                || rank.eq_ignore_ascii_case("domain")
                || rank.eq_ignore_ascii_case("acellular root")
            {
                Some(name.clone())
            } else {
                None
            }
        });

        let final_result = result.or_else(|| self.name_of(taxid));
        self.domains
            .write()
            .unwrap()
            .insert(taxid, final_result.clone());
        final_result
    }

    pub fn realm_of(&self, taxid: u32) -> Option<String> {
        if let Some(cached) = self.realms.read().unwrap().get(&taxid) {
            return cached.clone();
        }

        let lineage = self.lineage(taxid);
        let result = lineage.iter().find_map(|(_, rank, name)| {
            if rank.eq_ignore_ascii_case("realm") {
                Some(name.clone())
            } else {
                None
            }
        });

        self.realms.write().unwrap().insert(taxid, result.clone());

        result
    }

    pub fn uses_modern_ranks(&self) -> bool {
        self.is_modern
    }

    pub fn resolve_taxid(&self, taxid: u32) -> Option<u32> {
        if self.nodes.contains_key(&taxid) {
            return Some(taxid);
        }

        if let Some(&merged) = self.merged.get(&taxid) {
            let mut visited = HashSet::new();
            let mut current = merged;
            visited.insert(taxid);

            while !self.nodes.contains_key(&current) {
                if !visited.insert(current) {
                    eprintln!(
                        "warning: taxid {taxid} could not be resolved because merged targets form a cycle"
                    );
                    return None;
                }
                if let Some(&next) = self.merged.get(&current) {
                    current = next;
                    continue;
                }
                if self.deleted.contains(&current) {
                    eprintln!("warning: taxid {taxid} was merged into deleted taxid {current}");
                    return None;
                }
                eprintln!("warning: taxid {taxid} was merged into unknown taxid {current}");
                return None;
            }

            let name = self
                .names
                .get(&current)
                .cloned()
                .unwrap_or_else(|| current.to_string());
            eprintln!("warning: taxid {taxid} has been merged into {current} ({name})");
            return Some(current);
        }

        if self.deleted.contains(&taxid) {
            eprintln!("warning: taxid {taxid} has been deleted from the NCBI taxonomy");
            return None;
        }

        eprintln!("warning: taxid {taxid} is not present in the NCBI taxonomy");
        None
    }

    pub fn resolve_taxid_str(&self, value: &str) -> Option<u32> {
        let Ok(raw) = value.trim().parse::<u32>() else {
            return None;
        };
        self.resolve_taxid(raw)
    }
}

pub fn ensure_taxdump(dir: &Path) -> Result<()> {
    fs::create_dir_all(dir)?;
    let nodes = dir.join("nodes.dmp");
    let names = dir.join("names.dmp");
    if nodes.exists() && names.exists() {
        return Ok(());
    }
    let url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
    let resp = get(url)?;
    let bytes = resp.bytes()?;
    let gz = GzDecoder::new(&bytes[..]);
    let mut ar = Archive::new(gz);
    ar.unpack(dir)?;
    Ok(())
}

pub fn default_taxdump_dir() -> PathBuf {
    dirs::home_dir()
        .map(|p| p.join(".camitk"))
        .unwrap_or_else(|| PathBuf::from(".camitk"))
}

fn parse_nodes(path: &Path) -> Result<HashMap<u32, TaxNode>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut nodes = HashMap::new();
    for line in reader.lines() {
        let l = line?;
        let parts: Vec<&str> = l.split('|').collect();
        if parts.len() < 3 {
            continue;
        }
        let taxid: u32 = parts[0]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(0);
        let parent: u32 = parts[1]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(taxid);
        let rank = parts[2].trim_matches(|c: char| c.is_whitespace());
        nodes.insert(
            taxid,
            TaxNode {
                parent,
                rank: canonicalize_rank(rank),
            },
        );
    }
    Ok(nodes)
}

fn parse_names(path: &Path) -> Result<HashMap<u32, String>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut names = HashMap::new();
    for line in reader.lines() {
        let l = line?;
        let parts: Vec<&str> = l.split('|').collect();
        if parts.len() < 4 {
            continue;
        }
        let taxid: u32 = parts[0]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(0);
        let name = parts[1]
            .trim_matches(|c: char| c.is_whitespace())
            .to_string();
        let class = parts[3].trim_matches(|c: char| c.is_whitespace());
        if class != "scientific name" {
            continue;
        }
        names.entry(taxid).or_insert(name);
    }
    Ok(names)
}

fn canonicalize_rank(rank: &str) -> String {
    rank.trim().to_string()
}

fn parse_merged(path: &Path) -> Result<HashMap<u32, u32>> {
    if !path.exists() {
        return Ok(HashMap::new());
    }
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut merged = HashMap::new();
    for line in reader.lines() {
        let l = line?;
        let parts: Vec<&str> = l.split('|').collect();
        if parts.len() < 2 {
            continue;
        }
        let from: u32 = parts[0]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(0);
        let to: u32 = parts[1]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(0);
        if from > 0 && to > 0 {
            merged.insert(from, to);
        }
    }
    Ok(merged)
}

fn parse_delnodes(path: &Path) -> Result<HashSet<u32>> {
    if !path.exists() {
        return Ok(HashSet::new());
    }
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut deleted = HashSet::new();
    for line in reader.lines() {
        let l = line?;
        let parts: Vec<&str> = l.split('|').collect();
        if parts.is_empty() {
            continue;
        }
        let taxid: u32 = parts[0]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(0);
        if taxid > 0 {
            deleted.insert(taxid);
        }
    }
    Ok(deleted)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn domain_of_includes_acellular_root() {
        let mut nodes = HashMap::new();
        nodes.insert(
            1,
            TaxNode {
                parent: 1,
                rank: "no rank".to_string(),
            },
        );
        nodes.insert(
            10239,
            TaxNode {
                parent: 1,
                rank: "acellular root".to_string(),
            },
        );

        let mut names = HashMap::new();
        names.insert(1, "root".to_string());
        names.insert(10239, "Viruses".to_string());

        let taxonomy = Taxonomy {
            nodes,
            names,
            merged: HashMap::new(),
            deleted: HashSet::new(),
            ancestors: RwLock::new(HashMap::new()),
            lineages: RwLock::new(HashMap::new()),
            domains: RwLock::new(HashMap::new()),
            realms: RwLock::new(HashMap::new()),
            is_modern: true,
        };

        assert_eq!(taxonomy.domain_of(10239), Some("Viruses".to_string()));
    }
}
