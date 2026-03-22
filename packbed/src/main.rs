use clap::{self, Parser};
use std::path::PathBuf;

use packbed::{OverlapType, Role};

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
        help = "Paths to BED files delimited by comma"
    )]
    pub bed: Vec<PathBuf>,

    #[arg(
        short = 'm',
        long = "mode",
        help = "Mode of the file",
        value_name = "MODE",
        value_delimiter = ',',
        num_args = 1..,
        required = true,
    )]
    pub modes: Vec<Role>,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        long = "overlap_type",
        help = "Type of overlap",
        value_name = "TYPE",
        default_value = "exon"
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
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    simple_logger::init_with_level(log::Level::Warn)?;
    let args = Args::parse();

    let components = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .map_err(Box::<dyn std::error::Error>::from)?
        .install(|| packbed::pack(args.bed, args.modes, args.overlap_type))?;

    drop(components);
    Ok(())
}
