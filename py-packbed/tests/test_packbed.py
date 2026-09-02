from pathlib import Path

import pytest

import packbed


def write_bed(tmp_path: Path, name: str, contents: str) -> Path:
    path = tmp_path / name
    path.write_text(contents)
    return path


def iter_genes(components):
    for per_track in components.values():
        for component in per_track:
            yield from component


def test_pack_returns_converted_genepreds(tmp_path):
    reference = write_bed(
        tmp_path,
        "reference.bed",
        "chr1\t100\t260\tref1\t0\t+\t120\t240\t0,0,0\t2\t50,40,\t0,120,\n",
    )
    query = write_bed(
        tmp_path,
        "query.bed",
        "chr1\t130\t200\tquery1\t0\t+\t140\t180\t0,0,0\t1\t70,\t0,\n",
    )

    components = packbed.pack([reference, query], ["reference", "query"], "exon")

    assert set(components) == {"chr1:+"}
    assert len(components["chr1:+"]) == 1
    genes = sorted(iter_genes(components), key=lambda gene: gene.name)
    assert [gene.name for gene in genes] == ["query1", "ref1"]

    ref = genes[1]
    assert isinstance(ref, packbed.GenePred)
    assert ref.chrom == "chr1"
    assert ref.strand == "+"
    assert ref.start == 100
    assert ref.end == 260
    assert ref.cds_start == 120
    assert ref.cds_end == 240
    assert ref.exon_count == 2
    assert ref.exons == [(100, 150), (220, 260)]
    assert ref.introns == [(150, 220)]
    assert ref.rgb == "0,0,0"
    assert "GenePred(name=" in repr(ref)


def test_pack_rejects_invalid_python_inputs(tmp_path):
    bed = write_bed(tmp_path, "valid.bed", "chr1\t0\t10\n")

    with pytest.raises(ValueError, match="Invalid role"):
        packbed.pack([bed], ["target"], "exon")

    with pytest.raises(ValueError, match="Invalid overlap type"):
        packbed.pack([bed], ["reference"], "span")

    with pytest.raises(ValueError, match="expected the same number"):
        packbed.pack([bed], ["reference", "query"], "exon")

    with pytest.raises(TypeError):
        packbed.pack([bed], ["reference"], colorize=True)  # pyright: ignore[reportCallIssue]


def test_pack_maps_file_and_bed_errors_to_python_exceptions(tmp_path):
    invalid_bed = write_bed(tmp_path, "invalid.bed", "chr1\t10\t20\tgene\t0\t+\t15\n")

    with pytest.raises(ValueError, match="unsupported BED field count 7"):
        packbed.pack([invalid_bed], ["reference"], "exon")

    with pytest.raises(OSError, match="failed to read"):
        packbed.pack([tmp_path / "missing.bed"], ["reference"], "exon")


def test_all_bed_widths_are_detected(tmp_path):
    bed3 = write_bed(tmp_path, "w3.bed", "chr1\t0\t100\nchr1\t150\t250\n")
    bed6 = write_bed(tmp_path, "w6.bed", "chr1\t50\t120\tg6\t0\t+\n")
    bed8 = write_bed(tmp_path, "w8.bed", "chr1\t40\t200\tg8\t0\t+\t40\t90\n")
    bed9 = write_bed(tmp_path, "w9.bed", "chr1\t60\t220\tg9\t0\t+\t70\t210\t0,0,0\n")
    bed12_plus = write_bed(
        tmp_path,
        "w12.bed",
        "chr1\t10\t160\tg12\t0\t-\t20\t150\t0,0,0\t2\t50,40,\t0,110,\textra_col\n",
    )

    components = packbed.pack(
        [bed3, bed6, bed8, bed9, bed12_plus],
        ["reference"] * 5,
        "exon",
    )

    assert set(components) == {"chr1:+", "chr1:-"}
    # Unstranded BED3 records land on both strands; stranded ones on their own.
    assert len([g for g in iter_genes(components) if g.name == "g12"]) == 1
    plus = [g for c in components["chr1:+"] for g in c]
    minus = [g for c in components["chr1:-"] for g in c]
    assert len(plus) == 5
    assert len(minus) == 3
    assert {g.name for g in plus} == {".", "g6", "g8", "g9"}
    assert {g.name for g in minus} == {".", "g12"}


def test_unstranded_records_go_to_both_strands(tmp_path):
    bed = write_bed(
        tmp_path, "un.bed", "chr1\t0\t100\tu\t0\t.\t0\t100\t0,0,0\t1\t100,\t0,\n"
    )

    components = packbed.pack([bed], ["reference"], "exon")

    assert set(components) == {"chr1:+", "chr1:-"}
    for key in ("chr1:+", "chr1:-"):
        assert [g.name for g in iter_genes({key: components[key]})] == ["u"]


def test_empty_and_comment_only_files_yield_empty_result(tmp_path):
    bed = write_bed(tmp_path, "empty.bed", "# header line\n\n# only comments\n")

    assert packbed.pack([bed], ["reference"], "exon") == {}


def test_role_and_overlap_type_aliases(tmp_path):
    bed = write_bed(
        tmp_path, "aliases.bed", "chr1\t0\t100\tg\t0\t+\t0\t100\t0,0,0\t1\t100,\t0,\n"
    )

    for role in ("R", "r", "reference"):
        assert "chr1:+" in packbed.pack([bed], [role], "exon")

    for overlap in ("cds", "CDS", "Cds", "exon", "Exon", "bounds", "boundary"):
        assert "chr1:+" in packbed.pack([bed], ["reference"], overlap)


def test_multi_file_components_respect_roles_and_boundaries(tmp_path):
    reference = write_bed(
        tmp_path,
        "ref.bed",
        "chr1\t0\t100\tr1\t0\t+\t0\t100\t0,0,0\t1\t100,\t0,\n"
        "chr1\t200\t300\tr2\t0\t+\t200\t300\t0,0,0\t1\t100,\t0,\n",
    )
    query = write_bed(
        tmp_path, "query.bed", "chr1\t50\t150\tq1\t0\t+\t50\t150\t0,0,0\t1\t100,\t0,\n"
    )

    components = packbed.pack([reference, query], ["R", "q"], "boundary")

    assert set(components) == {"chr1:+"}
    assert len(components["chr1:+"]) == 2
    comp_sets = sorted(
        tuple(sorted(g.name for g in comp)) for comp in components["chr1:+"]
    )
    assert comp_sets == [("q1", "r1"), ("r2",)]


def test_gene_attributes_and_types(tmp_path):
    bed = write_bed(
        tmp_path,
        "attrs.bed",
        "chr1\t100\t260\tg\t0\t-\t120\t240\t0,0,0\t2\t50,40,\t0,120,\n",
    )

    gene = next(iter_genes(packbed.pack([bed], ["reference"], "exon")))

    assert isinstance(gene, packbed.GenePred)
    assert gene.name == "g"
    assert gene.chrom == "chr1"
    assert gene.strand == "-"
    assert gene.start == 100
    assert gene.end == 260
    assert gene.cds_start == 120
    assert gene.cds_end == 240
    assert gene.exon_count == 2
    assert gene.exons == [(100, 150), (220, 260)]
    assert gene.introns == [(150, 220)]
    assert gene.rgb == "0,0,0"
    assert all(isinstance(exon, tuple) and len(exon) == 2 for exon in gene.exons)
