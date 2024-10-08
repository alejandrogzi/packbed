use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::Mutex;

use hashbrown::HashMap;
use rayon::prelude::*;
use rmp_serde::{decode, encode};

pub mod record;
pub use record::{Bed12, GenePred};

pub type GenePredMap = HashMap<String, Vec<GenePred>>;

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
) -> Result<GenePredMap, anyhow::Error> {
    let contents = par_reader(files)?;
    let tracks = parse_tracks(&contents)?;

    Ok(tracks)
}

fn parse_tracks<'a>(contents: &'a str) -> Result<GenePredMap, anyhow::Error> {
    let mut tracks = contents
        .par_lines()
        .filter(|x| !x.starts_with("#"))
        .filter_map(|x| Bed12::parse(x).ok())
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

fn exonic_overlap(target_exons: &Vec<(u64, u64)>, query_exons: &Vec<(u64, u64)>) -> bool {
    for (tstart, tend) in target_exons {
        for (qstart, _) in query_exons {
            if tstart <= qstart && tend >= qstart {
                return true;
            }
        }
    }
    false
}

fn buckerize(tracks: GenePredMap, overlap_cds: bool) -> HashMap<String, Vec<Vec<Arc<GenePred>>>> {
    let cmap = Mutex::new(HashMap::new());

    tracks.into_par_iter().for_each(|(chr, records)| {
        let mut acc: Vec<(u64, u64, Vec<Arc<GenePred>>)> = Vec::new();

        for tx in records {
            let tx = Arc::new(tx);

            let group = acc
                .iter_mut()
                .any(|(ref mut group_start, ref mut group_end, txs)| {
                    let (tx_start, tx_end) = if overlap_cds {
                        (tx.cds_start, tx.cds_end)
                    } else {
                        (tx.start, tx.end)
                    };

                    if tx_start >= *group_start && tx_start <= *group_end {
                        // loop over txs exons and see if they overlap
                        let exon_overlap = txs.iter().any(|x| exonic_overlap(&x.exons, &tx.exons));

                        if exon_overlap {
                            *group_start = (*group_start).min(tx_start);
                            *group_end = (*group_end).max(tx_end);

                            txs.push(tx.clone());
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        false
                    }
                });

            if !group {
                acc.push((tx.start, tx.end, vec![tx.clone()]));
            }
        }

        let acc_map: Vec<Vec<Arc<GenePred>>> = acc.into_iter().map(|(_, _, txs)| txs).collect();
        cmap.lock().unwrap().insert(chr.to_string(), acc_map);
    });
    cmap.into_inner().unwrap()
}

pub fn packbed<T: AsRef<Path> + Debug + Send + Sync>(
    bed: Vec<T>,
    overlap_cds: bool,
) -> Result<HashMap<String, Vec<Vec<Arc<GenePred>>>>, anyhow::Error> {
    let tracks = unpack(bed).unwrap();
    let buckets = buckerize(tracks, overlap_cds);

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

pub fn get_component<T: AsRef<Path> + Debug + Send + Sync>(
    bed: Vec<T>,
    hint: Option<Vec<(String, Vec<usize>)>>,
    out: Option<T>,
    overlap_cds: Option<bool>,
) {
    let buckets = packbed(bed, overlap_cds.unwrap_or(false)).expect("Error packing bed files");

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
