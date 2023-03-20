from dataclasses import dataclass, replace
from typing_extensions import Self
import pybedtools


@dataclass(frozen=True)
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

        s = self.start - bp
        e = self.end + bp

        return replace(self, start=s, end=e)

    def align(self, bin_size: int) -> Self:
        """
        Return a window aligned to 10bp bins.

        e.g. chr2:54-89 => chr2:51-90
        """

        s = (self.start - 1) // bin_size * bin_size + 1
        e = (self.end - 1) // bin_size * bin_size + bin_size

        return replace(self, start=s, end=e)

    def make_bins(self, bin_size: int) -> list[Self]:
        """
        Return a list of Ranges representing bins.

        e.g. make_bins(chr2:54-84) => [Range(chr2, 51, 60), ..., Range(chr2, 81, 90)]
        """

        r = self.align(bin_size)
        starts = range(r.start, r.end, bin_size)

        return [
            replace(self, start=s, end=s + bin_size - 1)
            for s in starts
        ]

    def __len__(self) -> int:
        """
        Return the number of base pairs included in the range. Behaviour is similar to
        GRanges, so an interval "chr1:25-25" has length 1.
        """
        return self.end - self.start + 1


def interval_to_range(i: pybedtools.Interval) -> Range:
    """
    Convenience function to create Range from pybedtools.Interval. Assumes 0-based input
    (ie correct BED format).
    """
    score = i.score if i.score != "." else 0
    return Range(i.chrom, i.start+1, i.end, i.name, score, i.strand)
