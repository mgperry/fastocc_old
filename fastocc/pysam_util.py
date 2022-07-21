import pysam
from fastocc.range import Range

def fetch_safe(bam: pysam.AlignmentFile, r: Range):
    """
    Wraps pysam.AlignmentFile.fetch so it doesn't fail on out-of-bounds ranges.
    """
    max_len = bam.header.get_reference_length(r.chr)
    return bam.fetch(r.chr, max(r.start, 1), min(r.end, max_len))


def fetch_chrom_sizes(h: pysam.AlignmentHeader) -> dict[str, tuple[int, int]]:
    """
    Extract chromosome sizes from BAM header (via pysam) into the format expected by
    pybedtools.

    header.to_dict() => {"SQ": [{"SN": "chr1", "LN", 248956422}, ...], ...}
    output => {"chr1": (1, 248956422), ...}
    """
    SQ = h.to_dict()['SQ']
    return {seq["SN"]: (0, seq["LN"]) for seq in SQ}
