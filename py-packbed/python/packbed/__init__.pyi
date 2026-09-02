import os

class GenePred:
    """A single transcript record as exposed by the compiled `packbed` module."""

    name: str
    chrom: str
    strand: str
    start: int
    end: int
    cds_start: int
    cds_end: int
    exon_count: int
    exons: list[tuple[int, int]]
    introns: list[tuple[int, int]]
    rgb: str

    def __repr__(self) -> str:
        """Return a developer-facing string representation."""


def pack(
    bed: list[os.PathLike],
    roles: list[str],
    overlap_type: str = "exon",
) -> dict[str, list[list[GenePred]]]:
    """Pack BED files into overlapping components keyed by ``chrom:strand``.

    :param bed: paths to BED files (BED3-BED12+)
    :param roles: one of ``R``/``r``/``reference`` or ``Q``/``q``/``query`` per file
    :param overlap_type: ``exon``, ``cds``, or ``bounds``/``boundary``
    """
