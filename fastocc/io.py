
import csv
import sys
from dataclasses import astuple

import numpy as np

import pysam
import pyBigWig

from .range import Range


def fetch_safe(bam: pysam.AlignmentFile, r: Range):
    """
    Wraps pysam.AlignmentFile.fetch so it doesn't fail on out-of-bounds ranges.
    """
    max_len = bam.header.get_reference_length(r.chr)
    return bam.fetch(r.chr, max(r.start, 1), min(r.end, max_len))


def fetch_chrom_sizes(h: pysam.AlignmentHeader) -> dict[str, int]:
    """
    Extract chromosome sizes from BAM header (via pysam).

    h["SQ"] => [{"SN": "chr1", "LN", 248956422}, ...]
    output => {"chr1": 248956422, ...}
    """
    return {seq["SN"]: seq["LN"] for seq in h['SQ']}


def fetch_bw(r: Range, bw):
    bins = bw.intervals(r.chr, r.start, r.end)
    return np.array([b[2] for b in bins])


def write_chrom_sizes(chrom_sizes: dict[str, int], f):
    """
    Write chromosome sizes to a file in UCSC format. Format is the same as
    the output of `fetch_chrom_sizes`, ie [("chr1", 248956422), ...]
    """
    for chr, size in chrom_sizes.items():
        f.write(f"{chr}\t{size}\n")


class BigWigWriter:
    def __init__(self, f: str, header: dict[str, int]):
         self.file = pyBigWig.open(f, 'w')
         self.file.addHeader(list(header.items()))

    def write(self, occ):
        self.file.addEntries(
            occ.range.chr,
            occ.range.start-1,
            values=occ.data["smooth"].values.tolist(),
            span=occ.options.bin_size,
            step=occ.options.bin_size,
        )


class BedGraphWriter:
    header = ["chr", "start", "end", "nfr", "nuc", "window", "smooth"]

    def __init__(self, f):
        if f == "-":
            self.fh = sys.stdout
        else:
            self.fh = open(f, 'w')

        self.fh.write("\t".join(self.header) + "\n")

    def write(self, occ):
        occ.to_data_frame().to_csv(
            self.fh,
            sep="\t",
            float_format="%.2f",
            header=False,
            index=False
        )


class BEDWriter:
    def __init__(self, f):
        if f == "-":
            fh = sys.stdout
        else:
            fh = open(f, 'w')

        self.writer = csv.writer(fh, delimiter="\t")

    def write(self, r: Range):
        self.writer.writerow(astuple(r))

    def writerows(self, ranges: list[Range]):
        for r in ranges:
            self.write(r)
