import csv
import sys
from dataclasses import astuple

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
