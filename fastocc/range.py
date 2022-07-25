from dataclasses import dataclass
from typing_extensions import Self
from pybedtools import Interval

@dataclass
class Range:
    """
    Class representing a genomic range or interval.

    NB this uses "GRanges" style inclusive ranges, not BED-style insanity.
    """

    chr: str
    start: int
    end: int
    name: str = "."
    score: float = 0
    strand: str = "."

    def extend(self, bp: int) -> Self:
        """
        Return a range extended by *bp* on either end. Does not check bounds.
        """

        start = self.start - bp
        end = self.end + bp

        return Range(self.chr, start, end, self.name, self.score, self.strand)

    def align(self) -> Self:
        """
        Return a window aligned to 10bp bins.

        e.g. chr2:54-89 => chr2:51-90
        """
        start = (self.start - 1) // 10 * 10 + 1
        end = (self.end - 1) // 10 * 10 + 10
        return Range(self.chr, start, end, self.name, self.score, self.strand)

    def make_bins(self) -> list[Self]:
        """
        Return a list of ranges
        """
        r = self.align()
        starts = range(r.start, r.end, 10)
        return [
            Range(self.chr, s, s + 9, self.name, self.score, self.strand)
            for s in starts
        ]

    def __len__(self) -> int:
        """
        Return the number of base pairs included in the range. Behaviour is similar to
        GRanges, so an interval "chr1:25-25" has length 1.
        """
        return self.end - self.start + 1


def interval_to_range(i: Interval) -> Range:
    """
    Convenience function to create Range from pybedtools.Interval. Assumes 0-based input
    (ie correct BED format).
    """
    score = i.score if i.score is not "." else 0
    return Range(i.chrom, i.start+1, i.end, i.name, score, i.strand)
