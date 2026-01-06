use crate::cami::load_samples;
use anyhow::Result;
use std::collections::HashMap;
use std::io::{self, Write};
use std::path::PathBuf;

pub struct ListConfig<'a> {
    pub input: Option<&'a PathBuf>,
}

pub fn run(cfg: &ListConfig) -> Result<()> {
    let samples = load_samples(cfg.input)?;
    let mut out = io::stdout();

    for sample in &samples {
        writeln!(out, "Sample: {}", sample.id)?;
        let rank_tokens = sample.header_rank_tokens();
        writeln!(out, "  Ranks: {}", rank_tokens.join("|"))?;
        writeln!(out, "  Total taxa: {}", sample.entries.len())?;

        let mut stats: HashMap<usize, (usize, f64)> = HashMap::new();
        for entry in &sample.entries {
            if let Some(idx) = sample.rank_index(&entry.rank) {
                let stat = stats.entry(idx).or_insert((0, 0.0));
                if entry.percentage > 0.0 {
                    stat.0 += 1;
                    stat.1 += entry.percentage;
                }
            }
        }

        for (idx, token) in rank_tokens.iter().enumerate() {
            let (count, total) = stats.get(&idx).cloned().unwrap_or((0, 0.0));
            writeln!(out, "    {}: taxa={} total={:.3}", token, count, total)?;
        }
    }

    Ok(())
}
