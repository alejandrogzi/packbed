use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use hashbrown::HashMap;
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
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
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
        v.par_sort_unstable_by_key(|x| (x.start, x.end));
    });

    let mut count = 0;
    count += tracks.values().map(|x| x.len()).sum::<usize>();
    dbg!(count);

    Ok(tracks)
}

fn exonic_overlap(exons_a: &Vec<(u64, u64)>, exons_b: &Vec<(u64, u64)>) -> bool {
    let mut i = 0;
    let mut j = 0;

    while i < exons_a.len() && j < exons_b.len() {
        let (start_a, end_a) = exons_a[i];
        let (start_b, end_b) = exons_b[j];

        if start_a < end_b && start_b < end_a {
            return true;
        }

        if end_a < end_b {
            i += 1;
        } else {
            j += 1;
        }
    }

    false
}

fn buckerize(
    tracks: GenePredMap,
    overlap_cds: bool,
    colorize: bool,
) -> HashMap<String, Vec<Vec<Arc<GenePred>>>> {
    let cmap = Mutex::new(HashMap::new());

    tracks.into_par_iter().for_each(|(chr, records)| {
        let mut acc: Vec<(u64, u64, Vec<Arc<GenePred>>, &str)> = Vec::new();

        for tx in records {
            let group = acc.iter_mut().any(
                |(ref mut group_start, ref mut group_end, txs, group_color)| {
                    let (tx_start, tx_end) = if overlap_cds {
                        (tx.cds_start, tx.cds_end)
                    } else {
                        (tx.start, tx.end)
                    };

                    if tx_start >= *group_start && tx_start <= *group_end {
                        // loop over txs exons and see if they overlap
                        let exon_overlap = txs
                            .iter()
                            .any(|group| exonic_overlap(&group.exons, &tx.exons));

                        if exon_overlap {
                            *group_start = (*group_start).min(tx_start);
                            *group_end = (*group_end).max(tx_end);

                            let tx = Arc::new(tx.clone());

                            if colorize {
                                txs.push(tx.clone().colorline(*group_color));
                            } else {
                                txs.push(tx.clone());
                            }

                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        false
                    }
                },
            );

            if !group {
                let color = choose_color();
                let tx = Arc::new(tx);

                if colorize {
                    acc.push((tx.start, tx.end, vec![tx.clone().colorline(color)], color));
                } else {
                    acc.push((tx.start, tx.end, vec![tx.clone()], color));
                }
            }
        }

        let acc_map: Vec<Vec<Arc<GenePred>>> = acc.into_iter().map(|(_, _, txs, _)| txs).collect();
        cmap.lock().unwrap().insert(chr.to_string(), acc_map);
    });

    cmap.into_inner().unwrap()
}

fn choose_color<'a>() -> &'a str {
    let mut rng = rand::thread_rng();
    let idx = rng.gen_range(0..RGB.len());
    RGB[idx]
}

pub fn packbed<T: AsRef<Path> + Debug + Send + Sync>(
    bed: Vec<T>,
    overlap_cds: bool,
    colorize: bool,
) -> Result<HashMap<String, Vec<Vec<Arc<GenePred>>>>, anyhow::Error> {
    let tracks = unpack(bed, overlap_cds).unwrap();
    let buckets = buckerize(tracks, overlap_cds, colorize);

    Ok(buckets)
}

pub fn binwriter<P: AsRef<Path> + Debug>(
    file: P,
    contents: HashMap<String, Vec<Vec<Arc<GenePred>>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(file)?;

    encode::write(&mut file, &contents)?;
    Ok(())
}

pub fn bedwriter<P: AsRef<Path> + Debug>(
    file: P,
    contents: HashMap<String, Vec<Vec<Arc<GenePred>>>>,
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
    colorize: Option<bool>,
) {
    let buckets = packbed(bed, overlap_cds.unwrap_or(false), colorize.unwrap_or(false))
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

            buckets
                .get(buckets.keys().next().unwrap())
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
    contents: HashMap<String, Vec<Vec<Arc<GenePred>>>>,
    output: T,
    subdirs: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(&output)?;

    contents.iter().par_bridge().for_each(|(chr, buckets)| {
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
