import sys
from typing import Iterator

import pysam
import numpy as np

from .range import Range, BedFile
from .options import get_options, Options
from .fragments import FragmentLoader, Fragments


def fetch_chrom_sizes(h: pysam.AlignmentHeader) -> dict[str, int]:
    """
    Extract chromosome sizes from BAM header (via pysam).

    h["SQ"] => [{"SN": "chr1", "LN", 248956422}, ...]
    output => {"chr1": 248956422, ...}
    """
    return {seq["SN"]: seq["LN"] for seq in h['SQ']}


def process_bam_reads(reads, opts: Options) -> Fragments:
    """
    Reconstruct full-length ATAC-seq fragments from a bam file. Filters for the first
    read of proper pairs, corrects the Tn5 insertion distance and returns the start and
    length (size) of a fragment in a named tuple.
    """
    # q = r.extend(self.opts.max_insert//2) # this is to capture long fragments

    midpoint = []
    size = []

    for read in reads:
        if read.is_proper_pair and not read.is_reverse:
            #get left position
            left_insertion = read.pos + 4
            # insert size corrected by 8 base pairs to be insertion to insertion
            frag_size = read.template_length - 8

            if frag_size > opts.max_insert:
                continue

            midpoint.append(left_insertion + frag_size // 2)
            size.append(frag_size)

    return np.array(midpoint, dtype=int), np.array(size, dtype=int)


class Bam:
    """
    Wrap pysam.AlignmentFile to make a few things easier.
    """
    def __init__(self, bam: pysam.AlignmentFile, index=True):
        if index and not bam.has_index():
            print("No BAM index found, indexing BAM file...", file=sys.stderr)
            fn = bam.filename.decode()
            bam.close()

            pysam.index(fn)
            bam = pysam.AlignmentFile(fn) # reload index
        
        self.bam = bam
        self.chrom_info = fetch_chrom_sizes(bam.header)

    def fetch(self, r: Range):
        """
        Wrap pysam so that extended ranges don't fail due to out-of-bounds
        """
        return self.bam.fetch(
            r.chr,
            max(r.start, 1),
            min(r.end, self.chrom_info[r.chr])
        )
    
    @classmethod
    def open(cls, bam_file: str, index=True):
        return cls(pysam.AlignmentFile(bam_file), index)


class BamFragments(FragmentLoader):
    def __init__(self, bam: Bam, ranges: Iterator[Range], opts: Options=None) -> None:
        self.opts = get_options(opts)
        self.ranges = ranges
        self.bam = bam

    def __iter__(self):
        for r in self.ranges:
            r_ext = r.extend(2*self.opts.shift*self.opts.bin_size) # extend for windowing // smoothing
            reads = self.bam.fetch(r_ext.extend(self.opts.max_insert//2)) # capture long fragments
            midpoint, size = process_bam_reads(reads, self.opts)

            yield Fragments(
                r_ext,
                midpoint,
                size,
                self.opts
            )
    
    @classmethod
    def load(cls, bam: str, bed: str, opts: Options, index=True):
        bam_file = Bam.open(bam, index=index)

        ranges = BedFile.from_bed(bed, bam_file.chrom_info, opts)

        return cls(bam_file, ranges, opts)

