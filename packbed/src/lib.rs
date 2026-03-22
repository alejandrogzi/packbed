use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use genepred::reader::ReaderError;
use genepred::{
    Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, BedFormat, GenePred, Reader, ReaderMode,
    ReaderOptions,
};
use log::warn;
use rayon::prelude::*;

#[cfg(feature = "dashmap")]
use dashmap::DashMap;
#[cfg(feature = "dashmap")]
use std::cmp::Eq;
#[cfg(feature = "dashmap")]
use std::hash::Hash;

/// Thread-safe or standard map type depending on compilation features.
///
/// When compiled with `dashmap` feature, uses `DashMap` for concurrent access.
/// Otherwise uses `HashMap`.
#[cfg(feature = "dashmap")]
pub type Map<K, V> = DashMap<K, V>;

/// Thread-safe or standard map type depending on compilation features.
///
/// When compiled with `dashmap` feature, uses `DashMap` for concurrent access.
/// Otherwise uses `HashMap`.
#[cfg(not(feature = "dashmap"))]
pub type Map<K, V> = HashMap<K, V>;

/// Creates a new thread-safe DashMap (requires `dashmap` feature).
#[cfg(feature = "dashmap")]
fn init_map<K: Hash + Eq, V>() -> DashMap<K, V> {
    DashMap::new()
}

/// Creates a new standard HashMap.
#[cfg(not(feature = "dashmap"))]
fn init_map<K, V>() -> HashMap<K, V> {
    HashMap::new()
}

/// Role of the input file in the packing process.
///
/// # Variants
/// - `Reference`: The BED file contains reference transcripts
/// - `Query`: The BED file contains query transcripts to compare
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    Reference,
    Query,
}

/// Parses a role from a string slice (panics on failure).
impl From<&str> for Role {
    fn from(value: &str) -> Self {
        value.parse().expect("ERROR: Cannot parse role!")
    }
}

/// Parses a role from a string.
///
/// Accepts "R", "r", "reference" for Reference, and "Q", "q", "query" for Query.
///
/// # Arguments
/// * `s` - The string to parse
///
/// # Returns
/// `Ok(Role)` on success, or `Err(String)` with an error message
impl std::str::FromStr for Role {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "R" | "r" | "reference" => Ok(Role::Reference),
            "Q" | "q" | "query" => Ok(Role::Query),
            _ => Err(format!(
                "Invalid role: '{}'. Expected R/r/reference or Q/q/query",
                s
            )),
        }
    }
}

/// Type of overlap detection for transcript grouping.
///
/// # Variants
/// - `CDS`: Overlap based on coding sequence (CDS) regions
/// - `Exon`: Overlap based on exon regions
/// - `Boundary`: Overlap based on transcript boundaries (start/end)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverlapType {
    CDS,
    Exon,
    Boundary,
}

/// Parses overlap type from a string slice (panics on failure).
impl From<&str> for OverlapType {
    fn from(value: &str) -> Self {
        value.parse().expect("ERROR: Cannot parse overlap type!")
    }
}

/// Parses overlap type from a string.
///
/// Accepts "cds", "exon", "bounds" or "boundary" (case-insensitive).
///
/// # Arguments
/// * `s` - The string to parse
///
/// # Returns
/// `Ok(OverlapType)` on success, or `Err(String)` with an error message
impl std::str::FromStr for OverlapType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "cds" => Ok(OverlapType::CDS),
            "exon" => Ok(OverlapType::Exon),
            "bounds" | "boundary" => Ok(OverlapType::Boundary),
            _ => Err(format!(
                "Invalid overlap type: '{}'. Expected cds, exon, or bounds",
                s
            )),
        }
    }
}

/// Error type for pack operations.
///
/// # Variants
/// - `InputCountMismatch`: Number of files and modes don't match
/// - `Io`: File read error with path and underlying error
/// - `UnsupportedFieldCount`: BED file has unsupported field count
/// - `Read`: Error reading BED records
#[derive(Debug)]
pub enum PackError {
    InputCountMismatch {
        files: usize,
        modes: usize,
    },
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    UnsupportedFieldCount {
        path: PathBuf,
        field_count: usize,
    },
    Read(ReaderError),
}

/// Formats the error for display.
impl Display for PackError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PackError::InputCountMismatch { files, modes } => write!(
                f,
                "expected the same number of BED paths and modes, got {files} paths and {modes} modes"
            ),
            PackError::Io { path, source } => {
                write!(f, "failed to read '{}': {source}", path.display())
            }
            PackError::UnsupportedFieldCount { path, field_count } => write!(
                f,
                "unsupported BED field count {field_count} in '{}'; supported BED widths are 3, 4, 5, 6, 8, 9, and 12+",
                path.display()
            ),
            PackError::Read(source) => Display::fmt(source, f),
        }
    }
}

/// Provides error source chain for `PackError`.
impl std::error::Error for PackError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            PackError::Io { source, .. } => Some(source),
            PackError::Read(source) => Some(source),
            PackError::InputCountMismatch { .. } | PackError::UnsupportedFieldCount { .. } => None,
        }
    }
}

/// Converts a `ReaderError` into a `PackError`.
impl From<ReaderError> for PackError {
    fn from(value: ReaderError) -> Self {
        PackError::Read(value)
    }
}

/// BED format variant based on field count.
///
/// # Variants
/// - `Bed3` to `Bed12`: Standard BED formats with corresponding field counts
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BedKind {
    Bed3,
    Bed4,
    Bed5,
    Bed6,
    Bed8,
    Bed9,
    Bed12,
}

/// Formats BedKind as string (e.g., "BED3").
impl Display for BedKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BedKind::Bed3 => f.write_str("BED3"),
            BedKind::Bed4 => f.write_str("BED4"),
            BedKind::Bed5 => f.write_str("BED5"),
            BedKind::Bed6 => f.write_str("BED6"),
            BedKind::Bed8 => f.write_str("BED8"),
            BedKind::Bed9 => f.write_str("BED9"),
            BedKind::Bed12 => f.write_str("BED12"),
        }
    }
}

impl BedKind {
    /// Determines BED format from field count.
    ///
    /// # Arguments
    /// * `path` - Path to the BED file (for error reporting)
    /// * `field_count` - Number of fields in the first data line
    ///
    /// # Returns
    /// `Ok(DetectedBedLayout)` or `Err(PackError::UnsupportedFieldCount)`
    fn from_field_count(path: &Path, field_count: usize) -> Result<DetectedBedLayout, PackError> {
        match field_count {
            3 => Ok(DetectedBedLayout::new(BedKind::Bed3, 0)),
            4 => Ok(DetectedBedLayout::new(BedKind::Bed4, 0)),
            5 => Ok(DetectedBedLayout::new(BedKind::Bed5, 0)),
            6 => Ok(DetectedBedLayout::new(BedKind::Bed6, 0)),
            8 => Ok(DetectedBedLayout::new(BedKind::Bed8, 0)),
            9 => Ok(DetectedBedLayout::new(BedKind::Bed9, 0)),
            12.. => Ok(DetectedBedLayout::new(BedKind::Bed12, field_count - 12)),
            _ => Err(PackError::UnsupportedFieldCount {
                path: path.to_path_buf(),
                field_count,
            }),
        }
    }

    /// Returns true if this BED format doesn't encode CDS boundaries.
    fn lacks_cds_bounds(self) -> bool {
        matches!(
            self,
            BedKind::Bed3 | BedKind::Bed4 | BedKind::Bed5 | BedKind::Bed6
        )
    }
}

/// Detected BED file layout with format type and additional fields.
///
/// # Fields
/// * `kind` - The BED format variant (Bed3-Bed12)
/// * `additional_fields` - Number of fields beyond the standard 12
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DetectedBedLayout {
    kind: BedKind,
    additional_fields: usize,
}

impl DetectedBedLayout {
    /// Creates a new DetectedBedLayout.
    ///
    /// # Arguments
    /// * `kind` - The BED format variant
    /// * `additional_fields` - Number of extra fields beyond standard 12
    const fn new(kind: BedKind, additional_fields: usize) -> Self {
        Self {
            kind,
            additional_fields,
        }
    }
}

/// Packs BED files into overlapping transcript groups.
///
/// Groups transcripts that overlap based on the specified overlap type.
/// Results are organized by chromosome and strand (e.g., "chr1:+", "chr1:-").
///
/// # Arguments
/// * `files` - Paths to BED files (supports BED3, BED4, BED5, BED6, BED8, BED9, BED12+)
/// * `modes` - Role for each file (Reference or Query)
/// * `overlap_type` - How to determine overlap (CDS, Exon, or Boundary)
///
/// # Returns
/// Map keyed by "chr:strand" containing vectors of overlapping transcript groups,
/// where each group contains vectors of GenePred records
///
/// # Example
/// ```rust, ignore
/// use packbed::{pack, Role, OverlapType};
///
/// let files = vec!["reference.bed", "query.bed"];
/// let modes = vec![Role::Reference, Role::Query];
/// let result = pack(files, modes, OverlapType::Exon).unwrap();
/// ```
pub fn pack<T: AsRef<Path> + Debug + Send + Sync>(
    files: Vec<T>,
    modes: Vec<Role>,
    overlap_type: OverlapType,
) -> Result<Map<String, Vec<Vec<GenePred>>>, PackError> {
    if files.len() != modes.len() {
        return Err(PackError::InputCountMismatch {
            files: files.len(),
            modes: modes.len(),
        });
    }

    let mut accumulator: Map<String, Vec<GenePred>> = init_map();

    for (file, mode) in files.iter().zip(modes.into_iter()) {
        let per_file = read_file(file.as_ref(), mode, overlap_type)?;
        for (key, mut values) in per_file {
            accumulator.entry(key).or_default().append(&mut values);
        }
    }

    Ok(buckerize(accumulator, overlap_type))
}

/// Reads and parses a BED file, returning track records grouped by key.
///
/// # Arguments
/// * `path` - Path to the BED file
/// * `role` - Role of the file (Reference or Query)
/// * `overlap_type` - Type of overlap detection
///
/// # Returns
/// HashMap of track records keyed by chromosome:strand
fn read_file(
    path: &Path,
    role: Role,
    overlap_type: OverlapType,
) -> Result<HashMap<String, Vec<GenePred>>, PackError> {
    let Some(layout) = sniff_bed_layout(path)? else {
        return Ok(HashMap::new());
    };

    if overlap_type == OverlapType::CDS && layout.kind.lacks_cds_bounds() {
        warn!(
            "CDS overlap requested for '{}' detected as {}. Falling back to transcript start/end because this BED format does not encode CDS bounds.",
            path.display(),
            layout.kind
        );
    }

    match layout.kind {
        BedKind::Bed3 => read_typed_file::<Bed3>(path, role, layout),
        BedKind::Bed4 => read_typed_file::<Bed4>(path, role, layout),
        BedKind::Bed5 => read_typed_file::<Bed5>(path, role, layout),
        BedKind::Bed6 => read_typed_file::<Bed6>(path, role, layout),
        BedKind::Bed8 => read_typed_file::<Bed8>(path, role, layout),
        BedKind::Bed9 => read_typed_file::<Bed9>(path, role, layout),
        BedKind::Bed12 => read_typed_file::<Bed12>(path, role, layout),
    }
}

/// Reads a BED file with a specific format type using parallel processing.
///
/// # Arguments
/// * `path` - Path to the BED file
/// * `role` - Role of the file (Reference or Query)
/// * `layout` - Detected BED format layout
///
/// # Returns
/// HashMap of track records keyed by chromosome:strand
fn read_typed_file<R>(
    path: &Path,
    role: Role,
    layout: DetectedBedLayout,
) -> Result<HashMap<String, Vec<GenePred>>, PackError>
where
    R: BedFormat + Into<GenePred> + Send,
{
    let options = ReaderOptions::new().additional_fields(layout.additional_fields);
    let reader = Reader::<R>::builder()
        .from_path(path)
        .mode(ReaderMode::Mmap)
        .options(options)
        .build()?;

    reader
        .par_records()?
        .fold(
            || Ok(HashMap::new()),
            |local_result, record| -> Result<HashMap<String, Vec<GenePred>>, PackError> {
                let mut local = local_result?;
                let mut gene = record.map_err(PackError::from)?;
                annotate_role(&mut gene, role);
                insert_gene(&mut local, gene);
                Ok(local)
            },
        )
        .reduce(
            || Ok(HashMap::new()),
            |left, right| -> Result<HashMap<String, Vec<GenePred>>, PackError> {
                let mut left = left?;
                merge_hash_maps(&mut left, right?);
                Ok(left)
            },
        )
}

/// Adds role annotation to a GenePred record.
///
/// # Arguments
/// * `record` - GenePred record to annotate
/// * `role` - Role to add ("reference" or "query")
fn annotate_role(record: &mut GenePred, role: Role) {
    match role {
        Role::Reference => record.add_extra("role", b"reference"),
        Role::Query => record.add_extra("role", b"query"),
    }
}

/// Inserts a GenePred into the tracks map.
///
/// For unstranded records, inserts into both strand keys.
/// For stranded records, inserts into the appropriate strand key.
///
/// # Arguments
/// * `tracks` - HashMap to insert into
/// * `gene` - GenePred to insert
fn insert_gene(tracks: &mut HashMap<String, Vec<GenePred>>, gene: GenePred) {
    let keys = track_keys(&gene);
    if keys.len() == 1 {
        tracks.entry(keys[0].clone()).or_default().push(gene);
        return;
    }

    for key in keys {
        tracks.entry(key).or_default().push(gene.clone());
    }
}

/// Generates track keys for a GenePred record.
///
/// # Arguments
/// * `record` - GenePred record
///
/// # Returns
/// Vector of keys: one for stranded ("chr:+"), two for unstranded ("chr:+", "chr:-")
fn track_keys(record: &GenePred) -> Vec<String> {
    let chrom = String::from_utf8_lossy(&record.chrom);
    match record.strand() {
        Some(genepred::Strand::Forward) => vec![format!("{chrom}:+")],
        Some(genepred::Strand::Reverse) => vec![format!("{chrom}:-")],
        Some(genepred::Strand::Unknown) | None => {
            vec![format!("{chrom}:+"), format!("{chrom}:-")]
        }
    }
}

/// Merges two HashMaps, appending values with matching keys.
///
/// # Arguments
/// * `left` - Target HashMap (mutated in place)
/// * `right` - Source HashMap (consumed)
fn merge_hash_maps<V>(left: &mut HashMap<String, Vec<V>>, right: HashMap<String, Vec<V>>) {
    for (key, mut values) in right {
        left.entry(key).or_default().append(&mut values);
    }
}

/// Detects BED format by reading the first non-comment, non-empty line.
///
/// Skips comment lines (starting with #) and empty lines.
///
/// # Arguments
/// * `path` - Path to the BED file
///
/// # Returns
/// `Ok(Some(DetectedBedLayout))` with format details, or `Ok(None)` if file is empty/comment-only
fn sniff_bed_layout(path: &Path) -> Result<Option<DetectedBedLayout>, PackError> {
    let file = File::open(path).map_err(|source| PackError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();

    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|source| PackError::Io {
                path: path.to_path_buf(),
                source,
            })?;

        if bytes == 0 {
            return Ok(None);
        }

        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let field_count = trimmed.split_ascii_whitespace().count();
        return BedKind::from_field_count(path, field_count).map(Some);
    }
}

/// Union-Find (Disjoint Set Union) data structure for grouping overlapping transcripts.
///
/// Uses path compression and union by rank for near-constant time operations.
#[derive(Debug, Clone)]
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    /// Creates a new UnionFind with n isolated elements.
    ///
    /// # Arguments
    /// * `n` - Number of elements
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    /// Finds the root of element x with path compression.
    ///
    /// # Arguments
    /// * `x` - Element index
    ///
    /// # Returns
    /// Root element index
    #[inline(always)]
    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    /// Unites elements x and y into the same set.
    ///
    /// # Arguments
    /// * `x` - First element index
    /// * `y` - Second element index
    #[inline(always)]
    fn union(&mut self, x: usize, y: usize) {
        let root_x = self.find(x);
        let root_y = self.find(y);
        if root_x != root_y {
            self.parent[root_y] = root_x;
        }
    }
}

/// Groups overlapping transcripts into buckets using Union-Find.
///
/// `tracks` must already be grouped by `chrom:strand`, for example `chr1:+` or `chr1:-`.
/// In `CDS` mode this function uses coding exons when present and falls back to the
/// transcript span when a record has no coding intervals.
pub fn buckerize<I>(tracks: I, overlap_type: OverlapType) -> Map<String, Vec<Vec<GenePred>>>
where
    I: IntoIterator<Item = (String, Vec<GenePred>)>,
{
    let mut cmap = init_map();

    tracks.into_iter().for_each(|(chr, transcripts)| {
        let mut spans = Vec::new();
        let mut uf = UnionFind::new(transcripts.len());

        for (idx, transcript) in transcripts.iter().enumerate() {
            for (start, end) in overlap_segments(transcript, overlap_type) {
                spans.push((start, end, idx));
            }
        }

        if !spans.is_empty() {
            spans.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
            let (_, first_end, first_idx) = spans[0];

            let mut prev_end = first_end;
            let mut prev_idx = first_idx;

            for &(start, end, idx) in &spans[1..] {
                if start < prev_end {
                    uf.union(prev_idx, idx);
                    prev_end = prev_end.max(end);
                } else {
                    prev_end = end;
                    prev_idx = idx;
                }
            }
        }

        let mut groups: HashMap<usize, Vec<GenePred>> = HashMap::new();
        for (idx, transcript) in transcripts.into_iter().enumerate() {
            let root = uf.find(idx);
            groups.entry(root).or_default().push(transcript);
        }

        cmap.insert(chr, groups.into_values().collect());
    });

    cmap
}

/// Extracts genomic segments for overlap detection from a `GenePred`.
///
/// # Arguments
/// * `record` - GenePred to extract segments from
/// * `overlap_type` - Type of segments to extract
///
/// # Returns
/// Vector of (start, end) tuples representing genomic intervals
fn overlap_segments(record: &GenePred, overlap_type: OverlapType) -> Vec<(u64, u64)> {
    match overlap_type {
        OverlapType::Exon => record.exons(),
        OverlapType::Boundary => vec![(record.start, record.end)],
        OverlapType::CDS => {
            let coding = record.coding_exons();
            if coding.is_empty() {
                vec![(record.start, record.end)]
            } else {
                coding
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use std::fs;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static TMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

    fn init_logger() {
        // simple_logger::init_with_level(log::Level::Info).unwrap_or_else(|source| {
        //     panic!("failed to initialize logger: {source}");
        // });
    }

    fn write_temp_bed(name: &str, contents: &str) -> PathBuf {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-packbed");
        fs::create_dir_all(&dir).unwrap();
        let id = TMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = dir.join(format!("{name}-{id}.bed"));
        fs::write(&path, contents).unwrap();
        path
    }

    fn canonicalize<I>(map: I) -> BTreeMap<String, Vec<Vec<String>>>
    where
        I: IntoIterator<Item = (String, Vec<Vec<GenePred>>)>,
    {
        let mut out = BTreeMap::new();

        for (key, components) in map {
            let mut canonical_components = components
                .into_iter()
                .map(|component| {
                    let mut genes = component
                        .into_iter()
                        .map(|gene| gene_signature(&gene))
                        .collect::<Vec<_>>();
                    genes.sort();
                    genes
                })
                .collect::<Vec<_>>();
            canonical_components.sort();
            out.insert(key, canonical_components);
        }

        out
    }

    fn gene_signature(gene: &GenePred) -> String {
        let name = gene
            .name()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .unwrap_or_else(|| ".".to_string());
        let strand = gene
            .strand()
            .map(|strand| strand.to_string())
            .unwrap_or_else(|| "?".to_string());
        let block_starts = gene
            .block_starts()
            .map(|starts| {
                starts
                    .iter()
                    .map(u64::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default();
        let block_ends = gene
            .block_ends()
            .map(|ends| {
                ends.iter()
                    .map(u64::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default();

        format!(
            "{}:{}-{}:{}:{}:{}:{}:{}:{}",
            String::from_utf8_lossy(gene.chrom()),
            gene.start(),
            gene.end(),
            name,
            strand,
            gene.thick_start()
                .map_or(".".to_string(), |value| value.to_string()),
            gene.thick_end()
                .map_or(".".to_string(), |value| value.to_string()),
            block_starts,
            block_ends
        )
    }

    fn legacy_bed12_pack(
        records: Vec<GenePred>,
        overlap_type: OverlapType,
    ) -> HashMap<String, Vec<Vec<GenePred>>> {
        let mut tracks = HashMap::new();

        for record in records {
            let chrom = String::from_utf8_lossy(record.chrom()).into_owned();
            let strand = match record.strand() {
                Some(genepred::Strand::Forward) => '+',
                Some(genepred::Strand::Reverse) => '-',
                Some(genepred::Strand::Unknown) | None => '?',
            };

            tracks
                .entry(format!("{chrom}:{strand}"))
                .or_insert_with(Vec::new)
                .push(record);
        }

        legacy_buckerize(tracks, overlap_type)
    }

    fn legacy_buckerize(
        tracks: HashMap<String, Vec<GenePred>>,
        overlap_type: OverlapType,
    ) -> HashMap<String, Vec<Vec<GenePred>>> {
        let mut cmap = HashMap::new();

        for (chr, transcripts) in tracks {
            let mut spans = Vec::new();
            let mut id_map = HashMap::new();
            let mut uf = UnionFind::new(transcripts.len());

            for (idx, transcript) in transcripts.iter().enumerate() {
                id_map.insert(idx, transcript);

                match overlap_type {
                    OverlapType::Exon => {
                        for &(start, end) in &transcript.exons() {
                            spans.push((start, end, idx));
                        }
                    }
                    OverlapType::CDS => {
                        for &(start, end) in &transcript.coding_exons() {
                            spans.push((start, end, idx));
                        }
                    }
                    OverlapType::Boundary => spans.push((transcript.start, transcript.end, idx)),
                }
            }

            spans.sort_unstable_by(|a, b| a.0.cmp(&b.0));

            let mut prev_end = spans[0].1;
            let mut prev_idx = spans[0].2;
            for &(start, end, idx) in &spans[1..] {
                if start < prev_end {
                    uf.union(prev_idx, idx);
                    prev_end = prev_end.max(end);
                } else {
                    prev_end = end;
                    prev_idx = idx;
                }
            }

            let mut groups = HashMap::new();
            for idx in 0..transcripts.len() {
                let root = uf.find(idx);
                groups
                    .entry(root)
                    .or_insert_with(Vec::new)
                    .push(id_map[&idx].clone());
            }

            cmap.insert(chr, groups.into_values().collect());
        }

        cmap
    }

    #[test]
    fn bed12_regression_matches_legacy_behavior() {
        init_logger();

        let bed12 = write_temp_bed(
            "bed12-regression",
            "chr1\t100\t260\ttx1\t0\t+\t120\t240\t0,0,0\t2\t50,40,\t0,120,\n\
             chr1\t130\t280\ttx2\t0\t+\t140\t250\t0,0,0\t2\t40,60,\t0,90,\n\
             chr1\t500\t600\ttx3\t0\t-\t520\t580\t0,0,0\t2\t30,20,\t0,80,\n\
             chr1\t540\t640\ttx4\t0\t-\t550\t620\t0,0,0\t2\t20,40,\t0,60,\n\
             chr2\t100\t150\ttx5\t0\t+\t110\t140\t0,0,0\t1\t50,\t0,\n",
        );

        let reference = {
            let reader = Reader::<Bed12>::builder()
                .from_path(&bed12)
                .mode(ReaderMode::Mmap)
                .build()
                .unwrap();
            reader
                .par_records()
                .unwrap()
                .map(|record| {
                    let mut gene = record.unwrap();
                    gene.add_extra("role", b"reference");
                    gene
                })
                .collect::<Vec<_>>()
        };

        for overlap_type in [OverlapType::Exon, OverlapType::CDS, OverlapType::Boundary] {
            let expected = canonicalize(legacy_bed12_pack(reference.clone(), overlap_type));
            let actual = canonicalize(
                pack(vec![bed12.clone()], vec![Role::Reference], overlap_type).unwrap(),
            );
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn automatically_detects_mixed_bed_formats() {
        init_logger();

        let bed3 = write_temp_bed("mixed-bed3", "chr1\t100\t220\n");
        let bed12_plus = write_temp_bed(
            "mixed-bed12-plus",
            "chr1\t150\t260\tplus_tx\t0\t+\t170\t240\t0,0,0\t2\t30,40,\t0,70,\n",
        );
        let bed12_minus = write_temp_bed(
            "mixed-bed12-minus",
            "chr1\t160\t250\tminus_tx\t0\t-\t180\t230\t0,0,0\t2\t20,30,\t0,60,\n",
        );

        let packed = canonicalize(
            pack(
                vec![bed3, bed12_plus, bed12_minus],
                vec![Role::Reference, Role::Query, Role::Query],
                OverlapType::Exon,
            )
            .unwrap(),
        );

        assert_eq!(packed["chr1:+"].len(), 1);
        assert_eq!(packed["chr1:-"].len(), 1);
        assert_eq!(packed["chr1:+"][0].len(), 2);
        assert_eq!(packed["chr1:-"][0].len(), 2);
        assert!(packed["chr1:+"]
            .iter()
            .flatten()
            .any(|gene| gene.contains("plus_tx")));
        assert!(packed["chr1:-"]
            .iter()
            .flatten()
            .any(|gene| gene.contains("minus_tx")));
    }

    #[test]
    fn cds_falls_back_to_transcript_bounds_for_bed6() {
        init_logger();

        let bed6 = write_temp_bed(
            "bed6-cds-fallback",
            "chr1\t100\t220\tgeneA\t0\t+\nchr1\t150\t260\tgeneB\t0\t+\n",
        );

        let packed =
            canonicalize(pack(vec![bed6], vec![Role::Reference], OverlapType::CDS).unwrap());

        assert_eq!(packed["chr1:+"].len(), 1);
        assert_eq!(packed["chr1:+"][0].len(), 2);
    }

    #[test]
    fn bed8_uses_thick_bounds_in_cds_mode() {
        init_logger();

        let bed8 = write_temp_bed(
            "bed8-cds",
            "chr1\t100\t220\tgeneA\t0\t+\t100\t140\nchr1\t120\t240\tgeneB\t0\t+\t180\t220\n",
        );

        let packed =
            canonicalize(pack(vec![bed8], vec![Role::Reference], OverlapType::CDS).unwrap());

        assert_eq!(packed["chr1:+"].len(), 2);
        assert_eq!(packed["chr1:+"][0].len(), 1);
        assert_eq!(packed["chr1:+"][1].len(), 1);
    }

    #[test]
    fn bed12_without_cds_no_longer_panics() {
        init_logger();

        let bed12 = write_temp_bed(
            "bed12-no-cds",
            "chr1\t100\t220\ttxA\t0\t+\t100\t100\t0,0,0\t2\t40,30,\t0,90,\n\
             chr1\t150\t260\ttxB\t0\t+\t150\t150\t0,0,0\t2\t20,40,\t0,70,\n",
        );

        let packed =
            canonicalize(pack(vec![bed12], vec![Role::Reference], OverlapType::CDS).unwrap());

        assert_eq!(packed["chr1:+"].len(), 1);
        assert_eq!(packed["chr1:+"][0].len(), 2);
    }

    #[test]
    fn buckerize_accepts_genepred_tracks_directly() {
        let mut gene_a = GenePred::from_coords(b"chr1".to_vec(), 100, 220, HashMap::new());
        gene_a.strand = Some(genepred::Strand::Forward);

        let mut gene_b = GenePred::from_coords(b"chr1".to_vec(), 150, 260, HashMap::new());
        gene_b.strand = Some(genepred::Strand::Forward);

        let tracks = HashMap::from([("chr1:+".to_string(), vec![gene_a, gene_b])]);
        let packed = canonicalize(buckerize(tracks, OverlapType::CDS));

        assert_eq!(packed["chr1:+"].len(), 1);
        assert_eq!(packed["chr1:+"][0].len(), 2);
    }

    #[test]
    fn mismatched_file_and_mode_counts_error() {
        init_logger();

        let bed3 = write_temp_bed("mismatch", "chr1\t0\t10\n");
        let err = pack(
            vec![bed3.clone(), bed3],
            vec![Role::Reference],
            OverlapType::Exon,
        )
        .unwrap_err();

        assert!(matches!(
            err,
            PackError::InputCountMismatch { files: 2, modes: 1 }
        ));
    }

    #[test]
    fn unsupported_field_count_errors() {
        init_logger();

        let invalid = write_temp_bed("invalid-bed7", "chr1\t10\t20\tgeneA\t0\t+\t15\n");
        let err = pack(vec![invalid], vec![Role::Reference], OverlapType::Exon).unwrap_err();

        assert!(matches!(
            err,
            PackError::UnsupportedFieldCount { field_count: 7, .. }
        ));
    }

    #[test]
    fn empty_files_are_ignored() {
        init_logger();

        let empty = write_temp_bed("empty", "# comment only\n\n");
        let packed = pack(vec![empty], vec![Role::Reference], OverlapType::Exon).unwrap();
        assert!(packed.is_empty());
    }
}
