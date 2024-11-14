use hashbrown::HashSet;
use serde::{Deserialize, Serialize};

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

    pub fn colorline(self: Self, color: &str) -> Self {
        let nline = self.line.clone();
        let mut fields = nline.split('\t').collect::<Vec<_>>();
        fields[8] = color;
        let new_line = fields.join("\t");

        GenePred {
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
        }
    }
}

impl Bed12 {
    #[inline(always)]
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

        let (exons, introns) = get_coords(
            exon_starts,
            exon_sizes,
            tx_start,
            tx_end,
            cds_start,
            cds_end,
            strand,
            cds_overlap,
        )?;

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

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

#[inline(always)]
fn get_coords(
    starts: &str,
    sizes: &str,
    tx_start: u64,
    tx_end: u64,
    cds_start: u64,
    cds_end: u64,
    strand: char,
    cds_overlap: bool,
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

    let ss = group(starts)?;
    let sz = group(sizes)?;

    if ss.len() != sz.len() {
        return Err("Exon start and end vectors have different lengths");
    }

    let offset = match strand {
        '+' => tx_start,
        '-' => tx_end,
        _ => return Err("Strand is not + or -"),
    };

    let exons = ss
        .iter()
        .zip(&sz)
        .map(|(&s, &z)| match strand {
            '+' => {
                if cds_overlap {
                    if s + z + offset < cds_start || s + offset > cds_end {
                        return Err("UTRs are not allowed in CDS exons");
                    } else if s + offset < cds_start {
                        if s + z + offset > cds_end {
                            return Ok((cds_start, cds_end));
                        }
                        return Ok((cds_start, s + z + offset));
                    } else if s + z + offset > cds_end {
                        return Ok((s + offset, cds_end));
                    }
                }

                Ok((s + offset, s + z + offset))
            }
            '-' => {
                if cds_overlap {
                    if offset - s < cds_start || offset - s - z > cds_end {
                        return Err("UTRs are not allowed in CDS exons");
                    } else if offset - s - z < cds_start {
                        if offset - s > cds_end {
                            return Ok((cds_start, cds_end));
                        }
                        return Ok((cds_start, offset - s));
                    } else if offset - s > cds_end {
                        return Ok((offset - s - z, cds_end));
                    }
                };

                Ok((offset - s - z, offset - s))
            }
            _ => return Err("Strand is not + or -"),
        })
        .filter_map(Result::ok)
        .collect::<HashSet<_>>();

    let introns = gapper(&exons);

    Ok((exons, introns))
}

#[inline(always)]
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

#[inline(always)]
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
        let tx_start = 10;
        let tx_end = 50;
        let cds_start = 15;
        let cds_end = 45;
        let strand = '+';
        let cds_overlap = true;

        let (exons, introns) = get_coords(
            start,
            size,
            tx_start,
            tx_end,
            cds_start,
            cds_end,
            strand,
            cds_overlap,
        )
        .unwrap();

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

        assert_eq!(
            exons,
            [(15, 20), (40, 45)].iter().cloned().collect::<Vec<_>>()
        );
        assert_eq!(introns, [(21, 39)].iter().cloned().collect::<Vec<_>>());
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_minus() {
        let start = "0,20,40,60,80";
        let size = "10,10,10,10,10";
        let tx_start = "10";
        let tx_end = "100";
        let cds_start = "30";
        let cds_end = "80";
        let strand = '-';
        let cds_overlap = true;

        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");
        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, get).unwrap();

        let (exons, introns) = get_coords(
            start,
            size,
            tx_start,
            tx_end,
            cds_start,
            cds_end,
            strand,
            cds_overlap,
        )
        .unwrap();

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

        assert_eq!(
            exons,
            [
                (99999999920, 99999999930),
                (99999999940, 99999999950),
                (99999999960, 99999999970)
            ]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
        );
        assert_eq!(
            introns,
            [(99999999931, 99999999939), (99999999951, 99999999959)]
                .iter()
                .cloned()
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_plus_nested_utr() {
        let start = "0,20,40,60,80";
        let size = "10,10,10,10,10";
        let tx_start = 10;
        let tx_end = 100;
        let cds_start = 15;
        let cds_end = 95;
        let strand = '+';
        let cds_overlap = true;

        let (exons, introns) = get_coords(
            start,
            size,
            tx_start,
            tx_end,
            cds_start,
            cds_end,
            strand,
            cds_overlap,
        )
        .unwrap();

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

        assert_eq!(
            exons,
            [(15, 20), (30, 40), (50, 60), (70, 80), (90, 95)]
                .iter()
                .cloned()
                .collect::<Vec<_>>()
        );
        assert_eq!(
            introns,
            [(21, 29), (41, 49), (61, 69), (81, 89)]
                .iter()
                .cloned()
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_minus_nested_utr() {
        let start = "0,20,40,60";
        let size = "10,10,10,10";
        let tx_start = "10";
        let tx_end = "80";
        let cds_start = "15";
        let cds_end = "75";
        let strand = '-';
        let cds_overlap = true;

        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");
        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, get).unwrap();

        let (exons, introns) = get_coords(
            start,
            size,
            tx_start,
            tx_end,
            cds_start,
            cds_end,
            strand,
            cds_overlap,
        )
        .unwrap();

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

        assert_eq!(
            exons,
            [
                (99999999925, 99999999930),
                (99999999940, 99999999950),
                (99999999960, 99999999970),
                (99999999980, 99999999985)
            ]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
        );
        assert_eq!(
            introns,
            [
                (99999999931, 99999999939),
                (99999999951, 99999999959),
                (99999999971, 99999999979)
            ]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
        );
    }
}
