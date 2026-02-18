use std::fmt::Debug;

use genepred::GenePred;
use packbed::{Map, OverlapType, Role};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use serde::{Deserialize, Serialize};

#[pyfunction]
#[pyo3(signature = (bed, roles, overlap_type="exon", colorize=true))]
fn pack<'a>(
    py: Python<'a>,
    bed: PyObject,
    roles: PyObject,
    overlap_type: &'a str,
    colorize: bool,
) -> PyResult<Bound<'a, PyDict>> {
    let bed = bed
        .extract::<Vec<String>>(py)
        .expect("ERROR: failed to extract bed files");
    let roles = roles
        .extract::<Vec<String>>(py)
        .expect("ERROR: failed to extract roles");
    let roles: Vec<Role> = roles.into_iter().map(|x| Role::from(x.as_str())).collect();

    let overlap_type = OverlapType::from(overlap_type);
    let buckets = packbed::pack(bed, roles, overlap_type).expect("ERROR: failed to pack bed files");

    convert_map_to_pydict(py, buckets)
}

#[pymodule]
#[pyo3(name = "packbed")]
#[allow(unused_variables)]
fn py_chromsize(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pack, m)?)?;
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
    pub rgb: String,
}

#[pymethods]
impl PyGenePred {
    #[getter]
    pub fn name(&self) -> &String {
        &self.name
    }

    #[getter]
    pub fn rgb(&self) -> &String {
        &self.rgb
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "PyGenePred(name='{}', chrom='{}', strand='{}', start={}, end={}, cds_start={}, cds_end={}, exon_count={}, exons={:?}, introns={:?}, color={:?})",
            self.name,
            self.chrom,
            self.strand,
            self.start,
            self.end,
            self.cds_start,
            self.cds_end,
            self.exon_count,
            self.exons,
            self.introns,
            self.rgb
        ))
    }
}

impl From<GenePred> for PyGenePred {
    fn from(gp: GenePred) -> Self {
        PyGenePred {
            name: std::str::from_utf8(
                gp.name()
                    .unwrap_or_else(|| panic!("ERROR: failed to get name")),
            )
            .unwrap()
            .to_string(),
            chrom: std::str::from_utf8(gp.chrom()).unwrap().to_string(),
            strand: match gp
                .strand()
                .unwrap_or_else(|| panic!("ERROR: failed to get strand"))
            {
                genepred::Strand::Forward => '+',
                genepred::Strand::Reverse => '-',
                genepred::Strand::Unknown => '?',
            },
            start: gp.start(),
            end: gp.end(),
            cds_start: gp
                .thick_start()
                .unwrap_or_else(|| panic!("ERROR: failed to get thick_start")),
            cds_end: gp
                .thick_end()
                .unwrap_or_else(|| panic!("ERROR: failed to get thick_end")),
            exons: gp.exons(),
            introns: gp.introns(),
            exon_count: gp.exon_count(),
            rgb: "0,0,0".to_string(),
        }
    }
}

pub fn convert_map_to_pydict(
    py: Python,
    map: Map<String, Vec<Vec<GenePred>>>,
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
                            .map(|arc_gp| PyGenePred::from((arc_gp).clone()).into_py(py))
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>(),
        );
        py_dict.set_item(key, py_list_of_lists)?;
    }

    Ok(py_dict)
}
