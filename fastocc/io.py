import csv
import sys

import numpy as np

import pyBigWig

from .range import Range
from .occupancy import Occupancy


class BigWigWriter:
    def __init__(self, f: str, chrom_info: dict[str, int], bin_size=10):
         self.file = pyBigWig.open(str(f), 'w')
         self.file.addHeader(list(chrom_info.items()))
         self.bin_size = bin_size

    def write(self, r: Range, data: np.ndarray | list[float]):
        self.file.addEntries(
            r.chr,
            r.start-1,
            values=data,
            span=self.bin_size,
            step=self.bin_size,
        )


class BedGraphWriter:
    header = ["chr", "start", "end", "nfr", "nuc", "occ"]

    def __init__(self, f):
        if f == "-":
            self.fh = sys.stdout
        else:
            self.fh = open(f, 'w')

        self.fh.write("\t".join(self.header) + "\n")

    def write(self, occ: Occupancy):
        occ.to_data_frame().to_csv(
            self.fh,
            sep="\t",
            float_format="%.2f",
            header=False,
            index=False
        )


class BEDWriter:
    """
    Class to write Ranges to a BED file.
    """
    def __init__(self, fh):
        self.writer = csv.writer(fh, delimiter="\t")

    def write(self, r: Range):
        self.writer.writerow(r.to_bed())

    def writerows(self, ranges: list[Range]):
        for r in ranges:
            self.write(r)

    @classmethod
    def open(cls, f: str):
        if f == "-":
            fh = sys.stdout
        else:
            fh = open(f, 'w')

        return cls(fh)
