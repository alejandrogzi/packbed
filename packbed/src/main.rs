use packbed::*;

use clap::{self, Parser};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[clap(
    name = "packbed",
    version = env!("CARGO_PKG_VERSION"),
    author = env!("CARGO_PKG_AUTHORS"),
    about = "pack a .bed into overlapping components"
)]
struct Args {
    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to BED12 files delimited by comma"
    )]
    pub bed: Vec<PathBuf>,

    #[arg(
        short = 'o',
        long = "output",
        required = true,
        value_name = "PATH",
        help = "Path to output BED12 file [not required if -c flag is set]"
    )]
    pub output: PathBuf,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'c',
        long = "comp",
        help = "Flag to split components into separate BED files",
        value_name = "COMPONENTS",
        default_value = "false"
    )]
    pub comp: bool,
}

impl Args {
    pub fn check(&self) -> anyhow::Result<()> {
        self.validate_args()
    }

    fn validate_args(&self) -> anyhow::Result<()> {
        self.check_dbs()?;

        Ok(())
    }

    fn check_dbs(&self) -> anyhow::Result<()> {
        if self.bed.is_empty() {
            let err = "No reference files provided".to_string();
            return Err(anyhow::anyhow!(err));
        }
        for db in &self.bed {
            validate(db)?;
        }
        Ok(())
    }
}

pub fn validate(arg: &PathBuf) -> anyhow::Result<()> {
    if !arg.exists() {
        return Err(anyhow::anyhow!("file {:?} does not exist", arg));
    }

    if !arg.is_file() {
        return Err(anyhow::anyhow!("file {:?} is not a file", arg));
    }

    match arg.extension() {
        Some(ext) if ext == "bed" => (),
        _ => {
            return Err(anyhow::anyhow!("file {:?} is not a BED file", arg));
        }
    }

    match std::fs::metadata(arg) {
        Ok(metadata) if metadata.len() == 0 => {
            return Err(anyhow::anyhow!("file {:?} is empty", arg));
        }
        Ok(_) => Ok(()),
        Err(e) => Err(e.into()),
    }
}

fn main() {
    let st = std::time::Instant::now();

    let args = Args::parse();
    args.check().unwrap_or_else(|e| {
        eprintln!("{}", e);
        std::process::exit(1);
    });

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let buckets = packbed(args.bed).expect("Error packing BED files");

    if args.comp {
        compwriter(buckets).expect("ERROR: Failed writing components to BED files")
    } else {
        binwriter(&args.output, buckets).expect("ERROR: Failed writing binary of components");
    }

    dbg!(st.elapsed());
}
