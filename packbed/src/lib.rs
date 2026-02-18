use std::collections::HashMap;
use std::fmt::Debug;
use std::path::Path;

use genepred::{Bed12, GenePred, Reader};
use rayon::prelude::*;

#[cfg(feature = "dashmap")]
use dashmap::DashMap;

#[cfg(feature = "dashmap")]
pub type Map<K, V> = DashMap<K, V>;

#[cfg(not(feature = "dashmap"))]
pub type Map<K, V> = HashMap<K, V>;

#[cfg(feature = "dashmap")]
fn init_map<K, V>() -> DashMap<K, V> {
    DashMap::new()
}

#[cfg(not(feature = "dashmap"))]
fn init_map<K, V>() -> HashMap<K, V> {
    HashMap::new()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Role {
    Reference,
    Query,
}

impl From<&str> for Role {
    fn from(value: &str) -> Self {
        match value {
            "R" | "r" | "reference" => Role::Reference,
            "Q" | "q" | "query" => Role::Query,
            _ => panic!("ERROR: Cannot parse role!"),
        }
    }
}

impl std::str::FromStr for Role {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "R" | "r" | "reference" => Ok(Role::Reference),
            "Q" | "q" | "query" => Ok(Role::Query),
            _ => Err(format!(
                "Invalid role: '{}'. Expected R/r/reference or Q/q/query",
                s
            )),
        }
    }
}

pub fn pack<T: AsRef<Path> + Debug + Send + Sync>(
    files: Vec<T>,
    modes: Vec<Role>,
    overlap_type: OverlapType,
) -> Result<Map<String, Vec<Vec<GenePred>>>, Box<dyn std::error::Error>> {
    let mut accumulator: Map<String, Vec<GenePred>> = init_map();

    for (file, mode) in files.iter().zip(modes) {
        #[cfg(feature = "dashmap")]
        {
            Reader::<Bed12>::from_mmap(file)?
                .par_records()?
                .for_each(|record| {
                    let mut record = record.expect("Error reading record");

                    let chrom = std::str::from_utf8(&record.chrom)
                        .expect("ERROR: could not convert chrom to str");

                    let strand = match record.strand() {
                        Some(genepred::Strand::Forward) => '+',
                        Some(genepred::Strand::Reverse) => '-',
                        Some(genepred::Strand::Unknown) | None => '?',
                    };

                    match mode {
                        Role::Reference => record.add_extra("role", b"reference"),
                        Role::Query => record.add_extra("role", b"query"),
                    }

                    let key = format!("{chrom}:{strand}");
                    accumulator.entry(key).or_default().push(record);
                });
        }

        #[cfg(not(feature = "dashmap"))]
        {
            // 1) Build per-thread HashMaps
            // 2) Merge them into a single per-file HashMap via reduce
            let per_file: HashMap<String, Vec<GenePred>> = Reader::<Bed12>::from_mmap(file)?
                .par_records()?
                .fold(
                    HashMap::new,
                    |mut local: HashMap<String, Vec<GenePred>>, record| {
                        let mut record = record.expect("Error reading record");

                        let bind = record.chrom.clone();
                        let chrom = std::str::from_utf8(&bind)
                            .expect("ERROR: could not convert chrom to str");

                        let strand = match record.strand() {
                            Some(genepred::Strand::Forward) => '+',
                            Some(genepred::Strand::Reverse) => '-',
                            Some(genepred::Strand::Unknown) | None => '?',
                        };

                        match mode {
                            Role::Reference => record.add_extra("role", b"reference"),
                            Role::Query => record.add_extra("role", b"query"),
                        }

                        let key = format!("{chrom}:{strand}");
                        local.entry(key).or_default().push(record);
                        local
                    },
                )
                .reduce(HashMap::new, merge_maps);

            for (k, mut v) in per_file {
                accumulator.entry(k).or_default().append(&mut v);
            }
        }
    }

    let tracks = buckerize(accumulator, overlap_type);
    Ok(tracks)
}

#[cfg(not(feature = "dashmap"))]
fn merge_maps(
    mut a: HashMap<String, Vec<GenePred>>,
    b: HashMap<String, Vec<GenePred>>,
) -> HashMap<String, Vec<GenePred>> {
    for (k, mut v) in b {
        a.entry(k).or_default().append(&mut v);
    }
    a
}

#[derive(Debug, Clone)]
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    #[inline(always)]
    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    #[inline(always)]
    fn union(&mut self, x: usize, y: usize) {
        let root_x = self.find(x);
        let root_y = self.find(y);
        if root_x != root_y {
            self.parent[root_y] = root_x;
        }
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum OverlapType {
    CDS,      // CDS-overlap
    Exon,     // exon-overlap
    Boundary, // boundary-overlap
}

impl From<&str> for OverlapType {
    fn from(value: &str) -> Self {
        match value {
            "cds" => OverlapType::CDS,
            "exon" => OverlapType::Exon,
            "bounds" => OverlapType::Boundary,
            _ => panic!("ERROR: Cannot parse overlap type!"),
        }
    }
}

fn buckerize(
    tracks: Map<String, Vec<GenePred>>,
    overlap_type: OverlapType,
) -> Map<String, Vec<Vec<GenePred>>> {
    let mut cmap = init_map();

    tracks.into_iter().for_each(|(chr, transcripts)| {
        let mut exons = Vec::new();
        let mut id_map = HashMap::new();
        let mut uf = UnionFind::new(transcripts.len());

        // if base mode, tx boundaries will behave as exons ranges
        for (i, transcript) in transcripts.iter().enumerate() {
            id_map.insert(i, transcript);

            match overlap_type {
                OverlapType::Exon => {
                    for &(start, end) in &transcript.exons() {
                        exons.push((start, end, i));
                    }
                }
                OverlapType::CDS => {
                    for &(start, end) in &transcript.coding_exons() {
                        exons.push((start, end, i));
                    }
                }
                OverlapType::Boundary => {
                    exons.push((transcript.start, transcript.end, i));
                }
            }
        }

        exons.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));

        let mut prev_end = exons[0].1;
        let mut prev_idx = exons[0].2;
        for &(start, end, idx) in &exons[1..] {
            if start < prev_end {
                uf.union(prev_idx, idx);
                prev_end = prev_end.max(end);
            } else {
                // no overlap, update prev_end and prev_idx
                prev_end = end;
                prev_idx = idx;
            }
        }

        let mut groups = HashMap::new();
        for i in 0..transcripts.len() {
            let root = uf.find(i);
            groups
                .entry(root)
                .or_insert_with(Vec::new)
                .push(id_map[&i].clone());
        }

        let comps = groups.into_iter().map(|(_, v)| v).collect();
        cmap.insert(chr, comps);
    });

    cmap
}
