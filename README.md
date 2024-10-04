<p align="center">
  <h1 align="center">
    packbed
  </h1>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.0.1dev-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.0.1-green">
    </a>
    <a href="https://crates.io/crates/packbed" target="_blank">
      <img alt="Crates.io Version" src="https://img.shields.io/crates/v/packbed">
    </a>
    <a href="https://github.com/alejandrogzi/packbed" target="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/alejandrogzi/packbed?color=blue">
    </a>
    <a href="https://crates.io/crates/packbed" target="_blank">
      <img alt="Crates.io Total Downloads" src="https://img.shields.io/crates/d/packbed">
    </a>
  </p>


  <p align="center">
    pack a .bed into overlapping components
  </p>

</p>


<p align="center">
    <img width=700 align="center" src="https://imgur.com/WUrvVM8.gif">
</p>


## Features
- pack any number of .bed files into overlapping components through a binary, Rust library or Python module
- split components into separate .bed files through a binary, Rust library or Python module
- write serialized components to a binary file through a binary and Rust library
- read serialized components from a binary file through a Rust library or Python module
- write specific components each to a different .bed file through a Rust library or Python module

## Usage
### Binary
``` rust
Usage: packbed [OPTIONS] --bed <PATHS>... --output <PATH>

Arguments:
    -b, --bed <PATHS>...     Paths to BED12 files delimited by comma
    -o, --output <PATH>      Path to output BED12 file [not required if -c flag is set]

Options:
    -t, --threads <THREADS>  Number of threads [default: 8]
    -c, --comp               Flag to split components into separate BED files
    -h, --help               Print help
    --version: print version
```

### Library
``` rust
use packbed::packbed;

fn main() {
    let bed1 = PathBuf::new("/path/to/b1.bed");
    let bed2 = PathBuf::new("/path/to/b2.bed");
    let beds = vec![bed1, bed2];

    let comps: HashMap<String, Vec<Vec<Arc<GenePred>>>> = packbed(beds).unwrap();
}
```
### Python
build the port to install it as a pkg:
```bash
git clone https://github.com/alejandrogzi/packbed.git && cd packbed/py-packbed
hatch shell
maturin develop --release
```
use it:
``` python
from packbed import pack

beds = ["path/to/bed1.bed", "path/to/bed2.bed"]
comps = pack(beds)
```

#### crate: [https://crates.io/crates/packbed](https://crates.io/crates/packbed)
