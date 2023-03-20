import pyBigWig
from fastocc.range import Range
from fastocc.options import Options
import csv
import sys
import numpy as np


class BigWigWriter:
    def __init__(self, f, header):
         self.file = pyBigWig.open(f, 'w')
         self.file.addHeader(header)

    def write(self, r: Range, data: dict[np.ndarray]):
        self.file.addEntries(r.chr, r.start-1, values=data["smooth"].tolist(), span=Options.bin_size, step=Options.bin_size)


class BedGraphWriter:
    def __init__(self, f):
        if f == "-":
            fh = sys.stdout
        else:
            fh = open(f, 'w')

        self.file = csv.writer(fh, delimiter="\t")
        self.file.writerow(["chr", "start", "end", "nfr", "nuc", "window", "smooth"])

    def write(self, r: Range, data: dict[np.ndarray]):
        for i, bin in enumerate(r.make_bins(Options.bin_size)):
            self.file.writerow([
                    bin.chr,
                    str(bin.start-1),
                    str(bin.end),
                    "%d" % data["nuc"][i],
                    "%d" % data["nfr"][i],
                    "%0.2f" % data["window"][i],
                    "%0.2f" % data["smooth"][i]
            ])
