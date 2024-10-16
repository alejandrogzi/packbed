use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use hashbrown::HashMap;
use packbed::{get_component, packbed, GenePred};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rayon::prelude::*;
use rmp_serde::decode;
use serde::{Deserialize, Serialize};

#[pyfunction]
#[pyo3(signature = (bed, overlap_cds=true, colorize=true))]
fn pack(
    py: Python,
    bed: PyObject,
    overlap_cds: bool,
    colorize: bool,
) -> PyResult<Bound<'_, PyDict>> {
    let bed = bed
        .extract::<Vec<String>>(py)
        .expect("ERROR: failed to extract bed files");
    let buckets = packbed(bed, overlap_cds, colorize).expect("ERROR: failed to pack bed files");
    convert_map_to_pydict(py, buckets)
}

#[pyfunction]
fn binreader(py: Python, path: PyObject) -> PyResult<Bound<'_, PyDict>> {
    let f = File::open(
        path.extract::<PathBuf>(py)
            .expect("ERROR: failed to extract file path"),
    )?;

    let contents = decode::from_read(f).expect("ERROR: failed to read file");
    convert_map_to_pydict(py, contents)
}

#[pyfunction]
#[pyo3(signature = (bed, hint,overlap_cds=true, out=None, colorize=true))]
fn to_component(
    py: Python,
    bed: PyObject,
    hint: PyObject,
    overlap_cds: Option<bool>,
    out: Option<PyObject>,
    colorize: Option<bool>,
) -> PyResult<()> {
    let hint = hint.extract::<Vec<(String, Vec<usize>)>>(py).ok();
    let bed = bed.extract::<Vec<String>>(py)?;

    if let Some(out) = out {
        get_component(
            bed,
            hint,
            Some(out.extract::<String>(py)?),
            overlap_cds,
            colorize,
        );
    } else {
        get_component(bed, hint, None, overlap_cds, colorize);
    }

    Ok(())
}

#[pyfunction]
#[pyo3(signature = (contents, output="comps.bed", subdirs=false, out_type="bed"))]
fn write_components(
    py: Python,
    contents: PyObject,
    output: Option<&str>,
    subdirs: Option<bool>,
    out_type: &str,
) -> PyResult<()> {
    let mut map: HashMap<String, Vec<Vec<Arc<PyGenePred>>>> = HashMap::new();
    let py_dict = contents.downcast_bound::<PyDict>(py)?;
    let out_type = TypeChoice::from_str(out_type).expect("ERROR: invalid output type");

    for (chr, buckets) in py_dict.iter() {
        let chr = chr.extract::<String>()?;
        let buckets = buckets.extract::<Vec<PyObject>>()?;
        let mut new_buckets: Vec<Vec<Arc<PyGenePred>>> = Vec::with_capacity(buckets.len());

        for bucket in buckets.iter() {
            let bucket = bucket.downcast_bound::<PyList>(py)?;
            let mut new_bucket: Vec<Arc<PyGenePred>> = Vec::with_capacity(bucket.len());

            for py_gene_pred in bucket.iter() {
                let gene_pred: PyGenePred = py_gene_pred.extract()?;
                new_bucket.push(Arc::new(gene_pred));
            }
            new_buckets.push(new_bucket);
        }

        map.insert(chr, new_buckets);
    }

    match out_type {
        TypeChoice::Comp => {
            let output = Path::new(output.unwrap().trim_end_matches(".bed"));
            let _ = compwriter(map, output, subdirs.unwrap());
        }
        TypeChoice::Bed => {
            let output = Path::new(output.unwrap_or("comps.bed"));
            let _ = bedwriter(output, map);
        }
    }

    Ok(())
}

enum TypeChoice {
    Comp,
    Bed,
}

impl TypeChoice {
    fn from_str(s: &str) -> Option<Self> {
        match s {
            "comp" => Some(Self::Comp),
            "bed" => Some(Self::Bed),
            _ => None,
        }
    }
}

pub fn bedwriter<P: AsRef<Path> + Debug>(
    file: P,
    contents: HashMap<String, Vec<Vec<Arc<PyGenePred>>>>,
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

#[pymodule]
#[pyo3(name = "packbed")]
#[allow(unused_variables)]
fn py_chromsize(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pack, m)?)?;
    m.add_function(wrap_pyfunction!(binreader, m)?)?;
    m.add_function(wrap_pyfunction!(to_component, m)?)?;
    m.add_function(wrap_pyfunction!(write_components, m)?)?;
    Ok(())
}

#[pyclass]
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct PyGenePred {
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

#[pymethods]
impl PyGenePred {
    #[getter]
    pub fn line(&self) -> &String {
        &self.line
    }

    #[getter]
    pub fn name(&self) -> &String {
        &self.name
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "PyGenePred(name='{}', chrom='{}', strand='{}', start={}, end={}, cds_start={}, cds_end={}, exon_count={}, exons={:?}, introns={:?})",
            self.name,
            self.chrom,
            self.strand,
            self.start,
            self.end,
            self.cds_start,
            self.cds_end,
            self.exon_count,
            self.exons,
            self.introns
        ))
    }
}

impl From<GenePred> for PyGenePred {
    fn from(gp: GenePred) -> Self {
        PyGenePred {
            name: gp.name,
            chrom: gp.chrom,
            strand: gp.strand,
            start: gp.start,
            end: gp.end,
            cds_start: gp.cds_start,
            cds_end: gp.cds_end,
            exons: gp.exons,
            introns: gp.introns,
            exon_count: gp.exon_count,
            line: gp.line,
        }
    }
}

pub fn convert_map_to_pydict(
    py: Python,
    map: HashMap<String, Vec<Vec<Arc<GenePred>>>>,
) -> PyResult<Bound<'_, PyDict>> {
    let py_dict = PyDict::new_bound(py);

    for (key, gene_preds) in map.into_iter() {
        let py_list_of_lists = PyList::new_bound(
            py,
            gene_preds
                .into_iter()
                .map(|vec| {
                    PyList::new_bound(
                        py,
                        vec.into_iter()
                            .map(|arc_gp| PyGenePred::from((*arc_gp).clone()).into_py(py))
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>(),
        );
        py_dict.set_item(key, py_list_of_lists)?;
    }

    Ok(py_dict)
}

pub fn compwriter<T: AsRef<Path> + Debug + Sync>(
    contents: HashMap<String, Vec<Vec<Arc<PyGenePred>>>>,
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
