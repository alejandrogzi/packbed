use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use flate2::read::MultiGzDecoder;
use hashbrown::HashMap;
use memmap2::Mmap;
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

#[inline(always)]
fn exonic_overlap(exons_a: &Vec<(u64, u64)>, exons_b: &Vec<(u64, u64)>) -> bool {
    let mut i = 0;
    let mut j = 0;

    while i < exons_a.len() && j < exons_b.len() {
        let (start_a, end_a) = exons_a[i];
        let (start_b, end_b) = exons_b[j];

        if start_a < end_b && start_b <= end_a || start_b < end_a && start_a <= end_b {
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
    overlap_exon: bool,
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

                        if overlap_cds || overlap_exon {
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
                            *group_start = (*group_start).min(tx_start);
                            *group_end = (*group_end).max(tx_end);

                            let tx = Arc::new(tx.clone());

                            if colorize {
                                txs.push(tx.clone().colorline(*group_color));
                            } else {
                                txs.push(tx.clone());
                            }

                            return true;
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
    overlap_exon: bool,
    colorize: bool,
) -> Result<HashMap<String, Vec<Vec<Arc<GenePred>>>>, anyhow::Error> {
    let tracks = unpack(bed, overlap_cds).unwrap();
    let buckets = buckerize(tracks, overlap_cds, overlap_exon, colorize);

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
            "s8	100	200	read1	0	-	110	190	0	3	20,20,20,	0,30,60,\ns8	100	200	read2	0	+	110	190	0	3	20,20,20,	0,30,60,"
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
        let r = &Vec::from([(10, 20), (30, 40), (50, 60)]);
        let q = &Vec::from([(15, 25), (41, 49), (90, 110)]);

        let res = exonic_overlap(r, q);

        assert_eq!(res, true);
    }

    #[test]
    fn test_exonic_overlap_false() {
        let r = &Vec::from([(10, 20), (30, 40), (50, 60)]);
        let q = &Vec::from([(21, 25), (41, 49), (90, 110)]);

        let res = exonic_overlap(r, q);

        assert_eq!(res, false);
    }
}
