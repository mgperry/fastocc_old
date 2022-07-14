from dataclasses import dataclass
from pydantic import validate_arguments

@validate_arguments
@dataclass
class Range:
    """
    NB this uses "GRanges" style inclusive ranges, not BED-style insanity.
    """
    chr: str
    start: int
    end: int
    name: str = "."
    score: float = 0
    strand: str = "."

    def extend(self, bp: int):
        return Range(self.chr, self.start - bp, self.end + bp, self.name, self.score, self.strand)

    def align(self):
        s = (self.start - 1) // 10 * 10 + 1
        e = (self.end - 1) // 10 * 10 + 10
        return Range(self.chr, s, e, self.name, self.score, self.strand)

    def make_bins(self):
        r = self.align()
        starts = range(r.start, r.end, 10)
        return [Range(self.chr, s, s+9, self.name, self.score, self.strand) for s in starts]

    def __len__(self):
        return self.end - self.start + 1

