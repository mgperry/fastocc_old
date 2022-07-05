import csv
from dataclasses import dataclass
from pydantic import validate_arguments

@validate_arguments
@dataclass
class Range:
    """
    NB this uses "granges" style inclusive ranges, not BED-style insanity.
    """
    chr: str
    start: int
    end: int
    name: str = ""
    score: float = 0
    strand: str = "."

    def extend(self, bp: int):
        return Range(self.chr, self.start - bp, self.end + bp, self.name, self.score, self.strand)

    def align(self):
        s = self.start // 10 * 10
        e = self.end // 10 * 10 + 9
        return Range(self.chr, s, e, self.name, self.score, self.strand)

    # should be idempotent if already aligned
    def bin(self):
        r = self.align()
        starts = range(r.start, r.end, 10)
        return [Range(self.chr, s, s+10-1, self.name, self.score, self.strand) for s in starts]

    def __len__(self):
        return self.end - self.start + 1

    def bed_record(self):
        return f"{self.chr}\t{self.start-1}\t{self.end}\t{self.name}\t{self.score}\t{self.strand}"

    def bed_graph(self):
        return f"{self.chr}\t{self.start-1}\t{self.end}\t{self.score}"


def read_bed(f: str) -> list[Range]:
    """
    Read a BED6 file into a list of Ranges. Expects correct 0-based index.
    """
    with open(f, 'r') as fh:
        bed = [Range(row[0], int(row[1]) + 1, int(row[2]), row[3], float(row[4]), row[5]) for row in csv.reader(fh, delimiter="\t")]

    return bed
