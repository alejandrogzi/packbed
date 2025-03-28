use packbed::{record::OverlapType, *};

use clap::{self, Parser, ValueEnum};
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
        help = "Path to output BED12 file [will interpret as dir if -t flag is set to comp]"
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
        long = "type",
        help = "Type of output",
        value_name = "TYPE",
        value_enum,
        default_value = "bed"
    )]
    pub out_type: TypeChoice,

    #[arg(
        long = "overlap_type",
        help = "Type of overlap",
        value_name = "TYPE",
        default_value = "exon",
        conflicts_with = "overlap_exon"
    )]
    pub overlap_type: OverlapType,

    #[arg(
        short = 's',
        long = "subdirs",
        help = "Flag to split components into separate BED files in subdirectories",
        value_name = "FLAG",
        default_value = "false"
    )]
    pub subdirs: bool,

    #[arg(
        long = "colorize",
        help = "Flag to colorize components in output BED file",
        value_name = "FLAG",
        default_value = "false"
    )]
    pub colorize: bool,
}

#[derive(ValueEnum, Debug, Clone)]
enum TypeChoice {
    Bin,
    Comp,
    Bed,
}

impl Args {
    pub fn check(&self) -> anyhow::Result<()> {
        self.validate_args()
    }

    fn validate_args(&self) -> anyhow::Result<()> {
        self.check_dbs()?;

        match self.out_type {
            TypeChoice::Bed => {
                if !self.colorize {
                    anyhow::bail!("ERROR: --colorize flag must be set for bed output");
                }
            }
            _ => (),
        }

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
        Some(ext) if ext == "bed" || ext == "gz" => (),
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
        .build()
        .unwrap();

    let buckets =
        packbed(args.bed, args.overlap_type, args.colorize).expect("Error packing BED files");

    match args.out_type {
        TypeChoice::Bin => {
            binwriter(&args.output, buckets).expect("ERROR: Failed writing binary of components");
        }
        TypeChoice::Comp => compwriter(buckets, &args.output, args.subdirs)
            .expect("ERROR: Failed writing components to BED files"),
        TypeChoice::Bed => {
            bedwriter(&args.output, buckets).expect("ERROR: Failed writing components to BED files")
        }
    }

    dbg!(st.elapsed());
}
