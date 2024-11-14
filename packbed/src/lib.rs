use std::borrow::Borrow;
use std::cmp::PartialOrd;
use std::collections::BTreeSet;
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use dashmap::DashMap;
use flate2::read::MultiGzDecoder;
use hashbrown::HashMap;
use memmap2::Mmap;
use num_traits::{Num, NumCast};
use rand::Rng;
use rayon::prelude::*;
use rmp_serde::{decode, encode};

pub mod record;
pub use record::{Bed12, GenePred};

pub type GenePredMap = HashMap<String, Vec<GenePred>>;

pub const RGB: [&str; 10] = [
    "255,0,0",    // red
    "0,255,0",    // green
    "0,0,255",    // blue
    "58,134,47",  // dark-green
    "255,0,255",  // magenta
    "0,255,255",  // cyan
    "255,128,0",  // orange
    "51,153,255", // sky-blue
    "118,115,15", // dark-yellow
    "172,126,0",  // brown
];

fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
    match file.as_ref().extension() {
        Some(ext) => match ext.to_str() {
            Some("gz") => with_gz(&File::open(file)?),
            _ => {
                let mut file = File::open(file)?;
                let mut contents = String::new();
                file.read_to_string(&mut contents)?;
                Ok(contents)
            }
        },
        None => Err("No extension found".into()),
    }
}

fn with_gz(file: &File) -> Result<String, Box<dyn std::error::Error>> {
    let mmap = unsafe { Mmap::map(file)? };
    let mut decoder = MultiGzDecoder::new(&mmap[..]);

    let mut contents = String::new();
    decoder.read_to_string(&mut contents)?;

    Ok(contents)
}

fn par_reader<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
) -> Result<String, anyhow::Error> {
    let contents: Vec<String> = files
        .par_iter()
        .map(|path| reader(path).unwrap_or_else(|e| panic!("Error reading file: {:?}", e)))
        .collect();

    Ok(contents.concat())
}

fn unpack<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
    cds_overlap: bool,
) -> Result<GenePredMap, anyhow::Error> {
    let contents = par_reader(files)?;
    let tracks = parse_tracks(&contents, cds_overlap)?;

    Ok(tracks)
}

fn parse_tracks<'a>(contents: &'a str, cds_overlap: bool) -> Result<GenePredMap, anyhow::Error> {
    let mut tracks = contents
        .par_lines()
        .filter(|x| !x.starts_with("#"))
        .filter_map(|x| Bed12::parse(x, cds_overlap).ok())
        .fold(
            || HashMap::new(),
            |mut acc: GenePredMap, record| {
                acc.entry(record.chrom.clone()).or_default().push(record);
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (k, v) in map {
                    let acc_v = acc.entry(k).or_insert(Vec::new());
                    acc_v.extend(v);
                }
                acc
            },
        );

    // sort by start/end in descending order
    tracks.par_iter_mut().for_each(|(_, v)| {
        v.par_sort_unstable_by(|a, b| a.start.cmp(&b.start).then(b.end.cmp(&a.end)));
    });

    let mut count = 0;
    count += tracks.values().map(|x| x.len()).sum::<usize>();
    dbg!(count);

    Ok(tracks)
}

#[allow(dead_code)]
#[inline(always)]
fn exonic_overlap<N, I>(exons_a: &BTreeSet<(N, N)>, exons_b: I) -> bool
where
    N: Num + NumCast + Copy + PartialOrd,
    I: IntoIterator,
    I::Item: Borrow<(N, N)>,
{
    let mut iter_a = exons_a.iter();
    let mut iter_b = exons_b.into_iter();

    let mut exon_a = iter_a.next();
    let mut exon_b = iter_b.next();

    loop {
        match (exon_a, exon_b.as_ref()) {
            (Some(&(start_a, end_a)), Some(exon_b_ref)) => {
                let (start_b, end_b) = exon_b_ref.borrow();

                if start_a < *end_b && *start_b < end_a {
                    return true;
                }

                if end_a < *end_b {
                    exon_a = iter_a.next();
                } else {
                    exon_b = iter_b.next();
                }
            }
            _ => break,
        }
    }

    false
}

fn buckerize(
    tracks: GenePredMap,
    overlap_cds: bool,
    overlap_exon: bool,
    colorize: bool,
) -> DashMap<String, Vec<Vec<GenePred>>> {
    let cmap = DashMap::new();

    tracks.into_par_iter().for_each(|(chr, transcripts)| {
        let mut exons = Vec::new();
        let mut id_map = HashMap::new();
        let mut uf = UnionFind::new(transcripts.len());

        // if base mode, tx boundaries will behave as exons ranges
        for (i, transcript) in transcripts.iter().enumerate() {
            id_map.insert(i, transcript);

            if !overlap_exon && !overlap_cds {
                exons.push((transcript.start, transcript.end, i));
            } else {
                for &(start, end) in &transcript.exons {
                    exons.push((start, end, i));
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

        let comps = groups
            .into_iter()
            .map(|(_, v)| {
                if colorize {
                    let color = choose_color();
                    v.into_iter().map(|gp| gp.colorline(color)).collect()
                } else {
                    v
                }
            })
            .collect();

        cmap.insert(chr, comps);
    });

    cmap
}

fn choose_color<'a>() -> &'a str {
    let mut rng = rand::thread_rng();
    let idx = rng.gen_range(0..RGB.len());
    RGB[idx]
}

pub fn packbed<T: AsRef<Path> + Debug + Send + Sync>(
    bed: Vec<T>,
    overlap_cds: bool,
    overlap_exon: bool,
    colorize: bool,
) -> Result<DashMap<String, Vec<Vec<GenePred>>>, anyhow::Error> {
    let tracks = unpack(bed, overlap_cds).unwrap();
    let buckets = buckerize(tracks, overlap_cds, overlap_exon, colorize);

    Ok(buckets)
}

pub fn binwriter<P: AsRef<Path> + Debug>(
    file: P,
    contents: DashMap<String, Vec<Vec<GenePred>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(file)?;
    let contents = contents.into_iter().collect::<HashMap<_, _>>();

    encode::write(&mut file, &contents)?;
    Ok(())
}

pub fn bedwriter<P: AsRef<Path> + Debug>(
    file: P,
    contents: DashMap<String, Vec<Vec<GenePred>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = BufWriter::new(File::create(file)?);

    for (_, components) in contents {
        for component in components {
            for tx in component {
                writeln!(file, "{}", tx.line())?;
            }
        }
    }

    Ok(())
}

pub fn get_component<T: AsRef<Path> + Debug + Send + Sync>(
    bed: Vec<T>,
    hint: Option<Vec<(String, Vec<usize>)>>,
    out: Option<T>,
    overlap_cds: Option<bool>,
    overlap_exon: Option<bool>,
    colorize: Option<bool>,
) {
    let buckets = packbed(
        bed,
        overlap_cds.unwrap_or(false),
        overlap_exon.unwrap_or(false),
        colorize.unwrap_or(false),
    )
    .expect("Error packing bed files");

    // [(chr, [1,2,3,4]), (chr, [5,6,7])] fmt to get components

    match hint {
        Some(hint) => {
            hint.into_par_iter().for_each(|(chr, comps)| {
                if let Some(bucket) = buckets.get(&chr) {
                    comps.into_par_iter().for_each(|comp| {
                        let filename = format!("{}_{}.bed", chr, comp);
                        let mut file = BufWriter::new(File::create(&filename).unwrap());

                        bucket[comp].iter().for_each(|x| {
                            writeln!(file, "{}", x.line()).unwrap();
                        });
                    });
                } else {
                    eprintln!("Chromosome {} not found in buckets", chr);
                }
            });
        }
        None => {
            let mut f_out = match out {
                Some(x) => BufWriter::new(File::create(x).unwrap()),
                None => BufWriter::new(File::create("comp.bed").unwrap()),
            };

            let bs = buckets.clone().into_read_only();
            let k = bs.keys().next().unwrap();

            buckets
                .get(k)
                .unwrap()
                .first()
                .unwrap()
                .iter()
                .for_each(|x| {
                    writeln!(f_out, "{}", x.line()).unwrap();
                });
        }
    }
}

pub fn binreader<P: AsRef<Path> + Debug>(
    file: P,
) -> Result<HashMap<String, Vec<Vec<GenePred>>>, Box<dyn std::error::Error>> {
    let file = File::open(file)?;
    let data: HashMap<String, Vec<Vec<GenePred>>> = decode::from_read(file)?;

    Ok(data)
}

pub fn compwriter<T: AsRef<Path> + Debug + Sync>(
    contents: DashMap<String, Vec<Vec<GenePred>>>,
    output: T,
    subdirs: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(&output)?;

    contents.iter().par_bridge().for_each(|comps| {
        let chr = comps.key();
        let buckets = comps.value();

        buckets
            .iter()
            .enumerate()
            .par_bridge()
            .for_each(|(i, bucket)| {
                let filename = if subdirs {
                    std::fs::create_dir_all(format!(
                        "{}/comp_{}_{}",
                        output.as_ref().display(),
                        chr,
                        i
                    ))
                    .expect("ERROR: Could not create directory");

                    format!(
                        "{}/comp_{}_{}/{}_{}.bed",
                        output.as_ref().display(),
                        chr,
                        i,
                        chr,
                        i
                    )
                } else {
                    format!("{}/{}_{}.bed", output.as_ref().display(), chr, i)
                };

                let mut file =
                    BufWriter::new(File::create(&filename).expect("ERROR: Could not create file"));

                bucket.iter().for_each(|x| {
                    writeln!(file, "{}", x.line()).unwrap();
                });
            });
    });

    Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_packbed_with_simple_bed() {
        let mut file = NamedTempFile::with_suffix(".bed").unwrap();
        let path = file.path().to_path_buf();

        let _ = write!(
            file,
            "s8\t100\t200\tread1\t0\t-\t110\t190\t0\t3\t20,20,20,\t0,30,60,\ns8\t100\t200\tread2\t0\t+\t110\t190\t0\t3\t20,20,20,\t0,30,60,"
        );

        let bed = vec![path];
        let overlap_cds = false;
        let overlap_exon = false;
        let colorize = false;

        let res = packbed(bed, overlap_cds, overlap_exon, colorize).unwrap();

        assert_eq!(res.len(), 1);
    }

    #[test]
    fn test_packbed_with_simple_bed_gz() {
        use std::process::Command;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::with_suffix(".bed").unwrap();
        let path = file.path().to_path_buf();

        write!(
            file,
            "s8\t100\t200\tread1\t0\t-\t110\t190\t0\t3\t20,20,20,\t0,30,60,\ns8\t100\t200\tread2\t0\t+\t110\t190\t0\t3\t20,20,20,\t0,30,60,"
        ).unwrap();

        let output = Command::new("gzip")
            .arg(&path)
            .output()
            .expect("Failed to execute gzip command");

        if !output.status.success() {
            panic!(
                "gzip command failed with status: {:?}, stderr: {}",
                output.status,
                String::from_utf8_lossy(&output.stderr)
            );
        }

        let bed = vec![path.with_extension("bed.gz")];
        let overlap_cds = false;
        let overlap_exon = false;
        let colorize = false;

        let res = packbed(bed, overlap_cds, overlap_exon, colorize).unwrap();

        assert_eq!(res.len(), 1);
    }

    #[test]
    fn test_exonic_overlap_true() {
        let r = &BTreeSet::from([(10, 20), (30, 40), (50, 60)]);
        let q = &Vec::from([(15, 25), (41, 49), (90, 110)]);

        let res = exonic_overlap(r, q);

        assert_eq!(res, true);
    }

    #[test]
    fn test_exonic_overlap_false() {
        let r = &BTreeSet::from([(10, 20), (30, 40), (50, 60)]);
        let q = &Vec::from([(21, 25), (41, 49), (90, 110)]);

        let res = exonic_overlap(r, q);

        assert_eq!(res, false);
    }
}
