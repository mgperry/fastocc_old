import csv
import pysam
from dataclasses import dataclass


@dataclass
class Range:
    chr: str
    start: int
    end: int
    name: str = ""
    score: float = 0
    strand: str = "."

    def extend(self, bp: int):
        return Range(self.chr, self.start - bp, self.end + bp, self.name)

    def window(self, bp: int):
        midpoint = self.start + ((self.end - self.start) // 2)
        return Range(self.chr, midpoint - (bp // 2), midpoint + (bp // 2), self.name)

    def __len__(self):
        return self.end - self.start + 1

    def __str__(self):
        return f"{self.chr}\t{self.start}\t{self.end}\t{self.name}\t{self.score}\t{self.strand}"


def bin_range(r: Range) -> list[Range]:
    start = list(range(round(r.start, 10), round(r.end, 10) - 1, 10))
    return [Range(r.chr, s, s + 9) for s in start]


def read_bed(f: str) -> list[Range]:
    with open(f, 'r') as fh:
        bed = []
        for row in csv.reader(fh, delimiter="\t"):
            bed.append(Range(row[0], int(row[1]), int(row[2]), row[3]))

    return bed
