# Changelog

All notable changes to **packbed** (Rust library + CLI) and **py-packbed**
(Python bindings) are documented in this file.

The changelog is maintained per release, derived from the repository's git
history (version-labeled commits). When cutting a new release, rename the
`[Unreleased]` section to the new version, add the release date, and leave an
empty `[Unreleased]` section on top.

## [Unreleased]

## [0.0.14] - 2026-09-02

### Changed

- `packbed` 0.0.14 published to crates.io (d47b609).
- `py-packbed` version aligned to 0.0.14 so the `v0.0.14` tag passes
  publish-workflow validation.
- README updates (3c3c982) and CI fix installing `clippy`/`rustfmt` components
  so the nightly matrix is green again (1f39c98).

## [0.0.13] - 2026-09-02

### Breaking

- BED-agnostic packing algorithm (099b2a4): automatic detection of BED3,
  BED4, BED5, BED6, BED8, BED9 and BED12+ inputs; overlap grouping by exon,
  CDS or transcript boundaries; reference/query role annotation; library-only
  release (CLI retains the packing flags).
- `py-packbed` API reworked: `pack(bed, roles, overlap_type="exon")` delegates
  entirely to the `packbed` crate, `GenePred` result objects, errors surface
  as Python exceptions (`ValueError`, `OSError`, `TypeError`) instead of
  panics; the `colorize` argument was removed (conversion-only package).

### Added

- Rust algorithm test suite: brute-force interval-graph oracle compared
  against the Union-Find sweep for all overlap modes on seeded random BED12
  inputs, plus edge fixtures — bridge chains, strict touch adjacency,
  intron-only overlap (boundary vs exon), CDS-vs-exon divergence, unstranded
  records on both strands, unsorted input, duplicate records,
  multi-chromosome separation, and missing-file errors.
- Python test suite (`py-packbed/tests/`): BED-width detection, unstranded
  records, role and overlap-type aliases, multi-file components, error
  contracts, and `GenePred` attribute/type checks.
- Type stubs (`__init__.pyi`) and `py.typed` marker for the Python module.
- Changelog, kept under `assets/changelog/` for every release.
- CI hardening: formatting checks, `dashmap`-backend test run, clippy with
  `--all-features` for `packbed`, clippy for `py-packbed`, wheel build +
  pytest in CI; new `python-release.yml` publishing `py-packbed` wheels
  (Linux, macOS, Windows, sdist) to PyPI on `v*` tags.

### Changed

- `genepred` dependency 0.0.11 → 0.0.16.
- `py-packbed` version aligned to 0.0.13 (was 0.0.7); dependency footprint
  reduced to `packbed` + `pyo3`; packaging metadata (`pyproject.toml`) and
  READMEs updated; version badge and image moved to `assets/`.
- `.gitignore` now covers `target/`, build artifacts (`dist/`, `build/`,
  `__pycache__/`, `*.pyc`, `*.egg-info/`) and `.venv/`.

### Fixed

- README image reference (`assets/img.png`) and stale version badge.
- Panicking `expect()` paths in the Python bindings replaced with typed
  exceptions.

## [0.0.12] - 2026-03-20

- Version bump (da47d92).

## [0.0.11] - 2026-03-05

### Fixed

- `dashmap` feature configuration so the crate builds under both the default
  `hashmap` and the `dashmap` backends (b0c85e4).

## [0.0.10] - 2025-12-22

### Breaking

- Record parsing migrated to the `genepred` crate (3cdf2e6).

### Added

- Stable `py-packbed` release (67d7f5c, 2026-02-18).

## [0.0.9] - 2025-03-28

### Breaking

- Overlap mode exposed as `OverlapType` across library, CLI
  (`--overlap_type`) and Python module (076a23c).

## [0.0.8] - 2024-12-02

### Added

- RGB (`rbg`) field exposed in the Python port (3000fc0).

## [0.0.7] - 2024-11-14

### Breaking

- Range-query packing through a Union-Find approach, replacing the previous
  gapping strategy (da59ab4, d46f3ed).

### Added

- Python port v0.0.4 compatibility (808e7ee).

## [0.0.6] - 2024-11-06

### Fixed

- `buckerize` recovers transcripts missed by the sorting step (2de0981).

## [0.0.5] - 2024-11-04

### Fixed

- Sorting boundaries in the overlap sweep and corrected non-overlapping CDS
  logic (18947e0).

### Added

- Python exposure of the RGB field (7f09dc2).

## [0.0.4] - 2024-10-30

### Added

- Gzip reader support (a174c42, 58ccca7).

### Fixed

- Strand-aware scaling of CDS coordinates (5ad6c48).
- Coordinates updated when `--overlap_cds` is used (71e4cf5).

## [0.0.3] - 2024-10-23

### Added

- CDS flag propagation, additional colors and faster exonic overlap
  (0b46a6e).

### Fixed

- CDS position updated when `--overlap_cds` is used (4de55ad).

## [0.0.2] - 2024-10-16

### Added

- `--colorize` and `--type` output selection with a BED writer (b078c89).

### Breaking

- `write_components` now takes an `out_type` argument (898353a).

## [0.0.1] - 2024-10-04

### Added

- Initial release: pack BED files into overlapping components (c964010).
- Exonic overlap and flags (f624ab3) and the first Python port (44da078).
