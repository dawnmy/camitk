use crate::camitk::{Entry, Sample};
use crate::taxonomy::Taxonomy;
use std::collections::{BTreeSet, HashMap};

#[derive(Clone)]
struct RankDetail {
    rank: String,
    taxid: String,
    name: String,
}

type RankMap = HashMap<usize, RankDetail>;

pub fn renormalize(samples: &mut [Sample]) {
    for sample in samples.iter_mut() {
        let mut by_rank: HashMap<String, Vec<usize>> = HashMap::new();
        for (idx, entry) in sample.entries.iter().enumerate() {
            by_rank.entry(entry.rank.clone()).or_default().push(idx);
        }
        for idxs in by_rank.values() {
            let sum: f64 = idxs
                .iter()
                .map(|&i| {
                    let v = sample.entries[i].percentage;
                    if v > 0.0 { v } else { 0.0 }
                })
                .sum();
            if sum <= 0.0 {
                continue;
            }
            for &i in idxs {
                if sample.entries[i].percentage > 0.0 {
                    sample.entries[i].percentage = sample.entries[i].percentage / sum * 100.0;
                }
            }
        }
    }
}

pub fn round_percentages(samples: &mut [Sample]) {
    for sample in samples.iter_mut() {
        for entry in &mut sample.entries {
            entry.percentage = (entry.percentage * 100000.0).round() / 100000.0;
        }
    }
}

pub fn fill_up_to(
    samples: &mut [Sample],
    from_rank: Option<&str>,
    to_rank: &str,
    taxonomy: &Taxonomy,
) {
    for sample in samples.iter_mut() {
        if sample.entries.is_empty() {
            continue;
        }
        let Some(target_idx) = sample.rank_index(to_rank) else {
            continue;
        };
        let Some(base_idx) = select_base_rank(sample, from_rank) else {
            continue;
        };
        if base_idx < target_idx {
            continue;
        }

        let start_idx = target_idx.min(base_idx);
        let end_idx = target_idx.max(base_idx);

        let mut existing_by_key: HashMap<(usize, String), Entry> = HashMap::new();
        let mut existing_by_rank: HashMap<usize, Vec<String>> = HashMap::new();
        for entry in &sample.entries {
            if let Some(idx) = sample.rank_index(&entry.rank) {
                if idx < start_idx || idx > end_idx {
                    continue;
                }
                existing_by_key
                    .entry((idx, entry.taxid.clone()))
                    .or_insert_with(|| entry.clone());
                existing_by_rank
                    .entry(idx)
                    .or_default()
                    .push(entry.taxid.clone());
            }
        }

        let mut sums: HashMap<(usize, String), f64> = HashMap::new();
        let mut cache: HashMap<String, RankMap> = HashMap::new();

        for entry in &sample.entries {
            let Some(entry_rank_idx) = sample.rank_index(&entry.rank) else {
                continue;
            };
            if entry_rank_idx < start_idx || entry_rank_idx > base_idx {
                continue;
            }
            let fallback_entry = existing_by_key.get(&(entry_rank_idx, entry.taxid.clone()));
            let Some(rank_map) = rank_map_for(sample, taxonomy, &entry.taxid, &mut cache, || {
                fallback_entry.and_then(|e| fallback_rank_map(e, sample))
            }) else {
                continue;
            };

            for rank_idx in start_idx..=end_idx {
                if rank_idx == entry_rank_idx || rank_idx > entry_rank_idx {
                    continue;
                }
                if let Some(detail) = rank_map.get(&rank_idx) {
                    *sums.entry((rank_idx, detail.taxid.clone())).or_insert(0.0) +=
                        entry.percentage;
                }
            }
        }

        let mut new_entries: Vec<Entry> = Vec::new();
        for idx in 0..sample.ranks.len() {
            if idx < start_idx || idx > end_idx {
                continue;
            }
            let mut taxids: BTreeSet<String> = BTreeSet::new();
            if let Some(existing) = existing_by_rank.get(&idx) {
                for taxid in existing {
                    taxids.insert(taxid.clone());
                }
            }
            for ((s_idx, taxid), pct) in &sums {
                if *s_idx == idx && *pct > 0.0 {
                    taxids.insert(taxid.clone());
                }
            }

            for taxid in taxids {
                let existing_entry = existing_by_key.get(&(idx, taxid.clone()));
                let percentage = sums.get(&(idx, taxid.clone())).cloned().unwrap_or(0.0)
                    + existing_entry.map(|e| e.percentage).unwrap_or(0.0);
                if percentage <= 0.0 {
                    continue;
                }

                let Some(rank_map) = rank_map_for(sample, taxonomy, &taxid, &mut cache, || {
                    existing_entry.and_then(|e| fallback_rank_map(e, sample))
                }) else {
                    if let Some(entry) = existing_entry {
                        new_entries.push(entry.clone());
                    }
                    continue;
                };
                if !rank_map.contains_key(&idx) {
                    if let Some(entry) = existing_entry {
                        new_entries.push(entry.clone());
                    }
                    continue;
                }
                let detail = rank_map.get(&idx).unwrap();
                let (taxpath, taxpathsn) = build_paths(sample, &rank_map, idx);
                let (cami_genome_id, cami_otu, hosts) = existing_entry
                    .map(|e| {
                        (
                            e.cami_genome_id.clone(),
                            e.cami_otu.clone(),
                            e.hosts.clone(),
                        )
                    })
                    .unwrap_or((None, None, None));
                new_entries.push(Entry {
                    taxid: taxid.clone(),
                    rank: detail.rank.clone(),
                    taxpath,
                    taxpathsn,
                    percentage,
                    cami_genome_id,
                    cami_otu,
                    hosts,
                });
            }
        }
        sample.entries = new_entries;
    }
}

pub fn fill_up_default(samples: &mut [Sample], from_rank: Option<&str>, taxonomy: &Taxonomy) {
    for sample in samples.iter_mut() {
        if sample.entries.is_empty() {
            continue;
        }
        if let Some(target_rank) = sample.ranks.first().cloned() {
            fill_up_to(
                std::slice::from_mut(sample),
                from_rank,
                &target_rank,
                taxonomy,
            );
        }
    }
}

fn select_base_rank(sample: &Sample, from_rank: Option<&str>) -> Option<usize> {
    if let Some(rank) = from_rank {
        if let Some(idx) = sample.rank_index(rank) {
            if has_entries_at(sample, idx) {
                return Some(idx);
            }
        }
    }
    sample
        .ranks
        .iter()
        .enumerate()
        .rev()
        .find(|(idx, _)| has_entries_at(sample, *idx))
        .map(|(idx, _)| idx)
}

fn has_entries_at(sample: &Sample, idx: usize) -> bool {
    sample
        .entries
        .iter()
        .any(|e| sample.rank_index(&e.rank) == Some(idx))
}

fn rank_map_for<'a, F>(
    sample: &Sample,
    taxonomy: &Taxonomy,
    taxid: &str,
    cache: &'a mut HashMap<String, RankMap>,
    mut fallback: F,
) -> Option<&'a RankMap>
where
    F: FnMut() -> Option<RankMap>,
{
    if !cache.contains_key(taxid) {
        let mut map = build_rank_map(sample, taxonomy, taxid).unwrap_or_default();
        let needs_fallback = map.is_empty()
            || sample
                .ranks
                .iter()
                .enumerate()
                .any(|(idx, _)| map.get(&idx).is_none());
        if needs_fallback {
            if let Some(fallback_map) = fallback() {
                if map.is_empty() {
                    map = fallback_map;
                } else {
                    for (idx, value) in fallback_map {
                        map.entry(idx).or_insert(value);
                    }
                }
            }
        }
        cache.insert(taxid.to_string(), map);
    }
    let map = cache.get(taxid)?;
    if map.is_empty() { None } else { Some(map) }
}

fn build_rank_map(sample: &Sample, taxonomy: &Taxonomy, taxid: &str) -> Option<RankMap> {
    let tid = taxonomy.resolve_taxid_str(taxid)?;
    let lineage = taxonomy.lineage(tid);
    let mut map: RankMap = HashMap::new();
    for (tid_u32, rank, name) in lineage.iter() {
        let mut effective_rank = rank.clone();
        if rank.eq_ignore_ascii_case("no rank") {
            if *tid_u32 != tid {
                continue;
            }
            if sample.rank_index("strain").is_some() {
                effective_rank = "strain".to_string();
            }
        }
        if let Some(idx) = sample.rank_index(&effective_rank) {
            map.entry(idx).or_insert(RankDetail {
                rank: effective_rank,
                taxid: tid_u32.to_string(),
                name: name.clone(),
            });
        }
    }
    if let Some(other_idx) = sample.rank_index("other entries") {
        if map.get(&other_idx).is_none() {
            if let Some((tid_u32, _, name)) =
                lineage.iter().find(|(tid_u32, _, _)| *tid_u32 == 2_787_854)
            {
                map.insert(
                    other_idx,
                    RankDetail {
                        rank: "other entries".to_string(),
                        taxid: tid_u32.to_string(),
                        name: name.clone(),
                    },
                );
            }
        }
    }
    if map.is_empty() { None } else { Some(map) }
}

fn fallback_rank_map(entry: &Entry, sample: &Sample) -> Option<RankMap> {
    let ids: Vec<&str> = entry.taxpath.split('|').collect();
    let names: Vec<&str> = entry.taxpathsn.split('|').collect();
    if ids.is_empty() || names.is_empty() {
        return None;
    }
    let mut map: RankMap = HashMap::new();
    let upto = ids.len().min(names.len()).min(sample.ranks.len());
    for idx in 0..upto {
        let taxid = ids[idx].trim();
        let name = names[idx].trim();
        if !taxid.is_empty() {
            map.insert(
                idx,
                RankDetail {
                    rank: sample.ranks[idx].clone(),
                    taxid: taxid.to_string(),
                    name: name.to_string(),
                },
            );
        }
    }
    if map.is_empty() { None } else { Some(map) }
}

fn build_paths(sample: &Sample, rank_map: &RankMap, upto_idx: usize) -> (String, String) {
    let mut taxids = Vec::new();
    let mut names = Vec::new();
    let limit = upto_idx.min(sample.ranks.len().saturating_sub(1));
    for idx in 0..=limit {
        if let Some(detail) = rank_map.get(&idx) {
            taxids.push(detail.taxid.clone());
            names.push(detail.name.clone());
        } else {
            taxids.push(String::new());
            names.push(String::new());
        }
    }
    (taxids.join("|"), names.join("|"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cami::Sample;

    #[test]
    fn build_paths_includes_placeholders_for_missing_ranks() {
        let sample = Sample {
            id: String::new(),
            version: None,
            taxonomy_tag: None,
            ranks: vec![
                "cellular root".to_string(),
                "domain".to_string(),
                "kingdom".to_string(),
                "phylum".to_string(),
            ],
            rank_groups: Vec::new(),
            rank_aliases: HashMap::new(),
            entries: Vec::new(),
        };

        let mut rank_map: RankMap = HashMap::new();
        rank_map.insert(
            0,
            RankDetail {
                rank: "cellular root".to_string(),
                taxid: "131567".to_string(),
                name: "cellular organisms".to_string(),
            },
        );
        rank_map.insert(
            1,
            RankDetail {
                rank: "domain".to_string(),
                taxid: "2759".to_string(),
                name: "Eukaryota".to_string(),
            },
        );
        rank_map.insert(
            3,
            RankDetail {
                rank: "phylum".to_string(),
                taxid: "5794".to_string(),
                name: "Apicomplexa".to_string(),
            },
        );

        let (taxpath, taxpathsn) = build_paths(&sample, &rank_map, 3);

        assert_eq!(taxpath, "131567|2759||5794");
        assert_eq!(taxpathsn, "cellular organisms|Eukaryota||Apicomplexa");
    }
}
