use hashbrown::HashSet;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

const SCALE: u64 = 100000000000; // 100Gb

#[derive(Debug, PartialEq, Clone)]
pub struct Bed12;

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GenePred {
    pub name: String,
    pub chrom: String,
    pub strand: char,
    pub start: u64,
    pub end: u64,
    pub cds_start: u64,
    pub cds_end: u64,
    pub exons: Vec<(u64, u64)>,
    pub introns: Vec<(u64, u64)>,
    pub exon_count: usize,
    pub line: String,
}

impl GenePred {
    pub fn line(&self) -> &String {
        &self.line
    }

    pub fn name(&self) -> &String {
        &self.name
    }

    pub fn colorline(self: Arc<Self>, color: &str) -> Arc<Self> {
        let nline = self.line.clone();
        let mut fields = nline.split('\t').collect::<Vec<_>>();
        fields[8] = color;
        let new_line = fields.join("\t");

        Arc::new(GenePred {
            line: new_line,
            name: self.name.clone(),
            chrom: self.chrom.clone(),
            strand: self.strand,
            start: self.start,
            end: self.end,
            cds_start: self.cds_start,
            cds_end: self.cds_end,
            exons: self.exons.clone(),
            introns: self.introns.clone(),
            exon_count: self.exon_count,
        })
    }
}

impl Bed12 {
    pub fn parse(line: &str, cds_overlap: bool) -> Result<GenePred, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let mut fields = line.split('\t');
        let (
            chrom,
            tx_start,
            tx_end,
            name,
            _,
            strand,
            cds_start,
            cds_end,
            _,
            _,
            exon_sizes,
            exon_starts,
        ) = (
            fields.next().ok_or("Cannot parse chrom")?,
            fields.next().ok_or("Cannot parse tx_start")?,
            fields.next().ok_or("Cannot parse tx_end")?,
            fields.next().ok_or("Cannot parse name")?,
            fields.next().ok_or("Cannot parse score")?,
            fields
                .next()
                .ok_or("Cannot parse strand")?
                .chars()
                .next()
                .ok_or("Cannot parse strand as char")?,
            fields.next().ok_or("Cannot parse cds_start")?,
            fields.next().ok_or("Cannot parse cds_end")?,
            fields.next().ok_or("Cannot parse rgb")?,
            fields.next().ok_or("Cannot parse block_count")?,
            fields.next().ok_or("Cannot parse exon_sizes")?,
            fields.next().ok_or("Cannot parse exon_starts")?,
        );

        if strand != '+' && strand != '-' {
            return Err("Strand is not + or -");
        }

        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");
        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, get)?;

        let (exons, introns) = get_coords(exon_starts, exon_sizes, tx_start, tx_end, strand)?;

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by_key(|x| (x.0, x.1));
        if cds_overlap {
            // modify first and last exon to only preserve CDS
            exons[0].0 = cds_start;
            exons.last_mut().unwrap().1 = cds_end;
        }

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by_key(|x| (x.0, x.1));

        let exon_count = exons.len();

        Ok(GenePred {
            name: name.into(),
            chrom: chrom.into(),
            strand,
            start: tx_start,
            end: tx_end,
            cds_start,
            cds_end,
            exons,
            introns,
            exon_count,
            line: line.to_string(),
        })
    }
}

fn get_coords(
    start: &str,
    size: &str,
    plus: u64,
    minus: u64,
    strand: char,
) -> Result<(HashSet<(u64, u64)>, HashSet<(u64, u64)>), &'static str> {
    let group = |field: &str| -> Result<Vec<u64>, &'static str> {
        field
            .split(',')
            .filter_map(|num| {
                if !num.is_empty() {
                    Some(num.parse::<u64>().expect("Cannot parse number"))
                } else {
                    None
                }
            })
            .map(|num| Ok(num))
            .collect()
    };

    let ss = group(start)?;
    let sz = group(size)?;

    if ss.len() != sz.len() {
        return Err("Exon start and end vectors have different lengths");
    }

    let offset = match strand {
        '+' => plus,
        '-' => minus,
        _ => return Err("Strand is not + or -"),
    };

    let exons = ss
        .iter()
        .zip(&sz)
        .map(|(&s, &z)| match strand {
            '+' => Ok((s + offset, s + z + offset)),
            '-' => Ok((offset - s - z, offset - s)),
            _ => return Err("Strand is not + or -"),
        })
        .collect::<Result<HashSet<(u64, u64)>, &'static str>>()?;

    let introns = gapper(&exons);

    Ok((exons, introns))
}

fn abs_pos(
    tx_start: &str,
    tx_end: &str,
    cds_start: &str,
    cds_end: &str,
    strand: char,
    get: impl Fn(&str) -> Result<u64, &'static str>,
) -> Result<(u64, u64, u64, u64), &'static str> {
    match strand {
        '+' => {
            let tx_start = get(tx_start)?;
            let tx_end = get(tx_end)?;
            let cds_start = get(cds_start)?;
            let cds_end = get(cds_end)?;

            Ok((tx_start, tx_end, cds_start, cds_end))
        }
        '-' => {
            let tx_start = get(tx_start)?;
            let tx_end = get(tx_end)?;
            let cds_start = get(cds_start)?;
            let cds_end = get(cds_end)?;

            Ok((
                SCALE - tx_end,
                SCALE - tx_start,
                SCALE - cds_end,
                SCALE - cds_start,
            ))
        }
        _ => Err("Strand is not + or -"),
    }
}

fn gapper(intervals: &HashSet<(u64, u64)>) -> HashSet<(u64, u64)> {
    let mut vintervals: Vec<(u64, u64)> = intervals.iter().copied().collect();
    vintervals.sort_by(|a, b| a.0.cmp(&b.0));

    let mut gaps = HashSet::with_capacity(vintervals.len());
    for window in vintervals.windows(2) {
        if let [prev, next] = window {
            let gap_start = prev.1 + 1;
            let gap_end = next.0 - 1;

            if gap_start < gap_end {
                gaps.insert((gap_start, gap_end));
            }
        }
    }

    gaps
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed12_abs_pos_plus() {
        let tx_start = "10";
        let tx_end = "20";
        let cds_start = "10";
        let cds_end = "20";
        let strand = '+';

        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, |x| {
                Ok(x.parse().unwrap())
            })
            .unwrap();

        assert_eq!(tx_start, 10);
        assert_eq!(tx_end, 20);
        assert_eq!(cds_start, 10);
        assert_eq!(cds_end, 20);
    }

    #[test]
    fn test_bed12_abs_pos_minus() {
        let tx_start = "10";
        let tx_end = "20";
        let cds_start = "10";
        let cds_end = "20";
        let strand = '-';

        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, |x| {
                Ok(x.parse().unwrap())
            })
            .unwrap();

        assert_eq!(tx_start, SCALE - 20);
        assert_eq!(tx_end, SCALE - 10);
        assert_eq!(cds_start, SCALE - 20);
        assert_eq!(cds_end, SCALE - 10);
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_plus() {
        let start = "0,30";
        let size = "10,10";
        let plus = 10;
        let minus = 50;
        let strand = '+';

        let (exons, introns) = get_coords(start, size, plus, minus, strand).unwrap();

        assert_eq!(exons, [(10, 20), (40, 50)].iter().cloned().collect());
        assert_eq!(introns, [(21, 39)].iter().cloned().collect());
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_minus() {
        let start = "0,10,20";
        let size = "5,5,80";
        let plus = 800;
        let minus = 900;
        let strand = '-';

        let (exons, introns) = get_coords(start, size, plus, minus, strand).unwrap();

        assert_eq!(
            exons,
            [(800, 880), (885, 890), (895, 900)]
                .iter()
                .cloned()
                .collect()
        );
        assert_eq!(introns, [(881, 884), (891, 894)].iter().cloned().collect());
    }
}
