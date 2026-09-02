#![allow(clippy::useless_conversion)]

use std::path::PathBuf;

use packbed_core::{GenePred, OverlapType, PackError, PackedComponents, Role};
use pyo3::exceptions::{PyOSError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

#[pyfunction]
#[pyo3(signature = (bed, roles, overlap_type="exon"))]
fn pack<'py>(
    py: Python<'py>,
    bed: Vec<PathBuf>,
    roles: Vec<String>,
    overlap_type: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let roles = roles
        .into_iter()
        .map(|role| {
            role.parse::<Role>()
                .map_err(|source| PyValueError::new_err(source.to_string()))
        })
        .collect::<PyResult<Vec<_>>>()?;
    let overlap_type = overlap_type
        .parse::<OverlapType>()
        .map_err(|source| PyValueError::new_err(source.to_string()))?;

    let buckets = packbed_core::pack(bed, roles, overlap_type).map_err(pack_error_to_pyerr)?;

    convert_map_to_pydict(py, buckets)
}

#[pymodule]
#[pyo3(name = "packbed")]
fn py_packbed(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGenePred>()?;
    m.add_function(wrap_pyfunction!(pack, m)?)?;
    Ok(())
}

#[pyclass(name = "GenePred", module = "packbed")]
#[derive(Debug, PartialEq, Clone)]
pub struct PyGenePred {
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub strand: char,
    #[pyo3(get)]
    pub start: u64,
    #[pyo3(get)]
    pub end: u64,
    #[pyo3(get)]
    pub cds_start: u64,
    #[pyo3(get)]
    pub cds_end: u64,
    #[pyo3(get)]
    pub exons: Vec<(u64, u64)>,
    #[pyo3(get)]
    pub introns: Vec<(u64, u64)>,
    #[pyo3(get)]
    pub exon_count: usize,
    #[pyo3(get)]
    pub rgb: String,
}

#[pymethods]
impl PyGenePred {
    fn __repr__(&self) -> String {
        format!(
            "GenePred(name={:?}, chrom={:?}, strand={:?}, start={}, end={}, cds_start={}, cds_end={}, exon_count={})",
            self.name,
            self.chrom,
            self.strand,
            self.start,
            self.end,
            self.cds_start,
            self.cds_end,
            self.exon_count,
        )
    }
}

impl From<GenePred> for PyGenePred {
    fn from(gp: GenePred) -> Self {
        PyGenePred {
            name: gp
                .name()
                .map(bytes_to_string)
                .unwrap_or_else(|| ".".to_string()),
            chrom: bytes_to_string(gp.chrom()),
            strand: gp
                .strand()
                .map(|strand| strand.to_string().chars().next().unwrap_or('?'))
                .unwrap_or('?'),
            start: gp.start(),
            end: gp.end(),
            cds_start: gp.thick_start().unwrap_or_else(|| gp.start()),
            cds_end: gp.thick_end().unwrap_or_else(|| gp.end()),
            exons: gp.exons(),
            introns: gp.introns(),
            exon_count: gp.exon_count(),
            rgb: "0,0,0".to_string(),
        }
    }
}

pub fn convert_map_to_pydict(py: Python<'_>, map: PackedComponents) -> PyResult<Bound<'_, PyDict>> {
    let py_dict = PyDict::new_bound(py);

    for (key, components) in map.into_iter() {
        let py_components = PyList::empty_bound(py);
        for component in components {
            let py_component = PyList::empty_bound(py);
            for gene in component {
                py_component.append(Py::new(py, PyGenePred::from(gene))?)?;
            }
            py_components.append(py_component)?;
        }
        py_dict.set_item(key, py_components)?;
    }

    Ok(py_dict)
}

fn bytes_to_string(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn pack_error_to_pyerr(error: PackError) -> PyErr {
    let message = error.to_string();
    match error {
        PackError::InputCountMismatch { .. } | PackError::UnsupportedFieldCount { .. } => {
            PyValueError::new_err(message)
        }
        PackError::Io { .. } => PyOSError::new_err(message),
        PackError::Read(_) => PyRuntimeError::new_err(message),
    }
}
