import csv
import sys

import pyBigWig

from .range import Range


class BigWigWriter:
    def __init__(self, f: str, header: dict[str, int]):
         self.file = pyBigWig.open(f, 'w')
         self.file.addHeader(list(header.items()))

    def write(self, occ):
        self.file.addEntries(
            occ.range.chr,
            occ.range.start-1,
            values=occ.data["occ"].values.tolist(),
            span=occ.options.bin_size,
            step=occ.options.bin_size,
        )


class BedGraphWriter:
    header = ["chr", "start", "end", "nfr", "nuc", "raw_occ", "occ"]

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
