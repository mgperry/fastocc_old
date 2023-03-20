import pysam
from fastocc.range import Range
from contextlib import contextmanager
import os


def fetch_safe(bam: pysam.AlignmentFile, r: Range):
    """
    Wraps pysam.AlignmentFile.fetch so it doesn't fail on out-of-bounds ranges.
    """
    max_len = bam.header.get_reference_length(r.chr)
    return bam.fetch(r.chr, max(r.start, 1), min(r.end, max_len))


def fetch_chrom_sizes(h: pysam.AlignmentHeader) -> list[tuple[str, int]]:
    """
    Extract chromosome sizes from BAM header (via pysam).

    h["SQ"] => [{"SN": "chr1", "LN", 248956422}, ...]
    output => [("chr1", 248956422), ...]
    """
    return [(seq["SN"], seq["LN"]) for seq in h['SQ']]


def write_chrom_sizes(chrom_sizes: list[tuple[str, int]], f):
    """
    Write chromosome sizes to a file in UCSC format. Format is the same as
    the output of `fetch_chrom_sizes`, ie [("chr1", 248956422), ...]
    """
    for chr, l in chrom_sizes:
        f.write(f"{chr}\t{l}\n")


@contextmanager
def temp_chrom_info(chrs, fn=".chrom.info"):
    """
    Context manager which creates a temporary .chromInfo file for use by called CLI programs.
    """
    with open(fn, 'w') as tempfile:
        write_chrom_sizes(chrs, tempfile)

    yield fn

    os.unlink(fn)
