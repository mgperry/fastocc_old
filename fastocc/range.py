class Range:
    """
    Class representing a genomic range or interval.

    NB this uses "GRanges" style inclusive ranges, not BED-style insanity.
    """

    def __init__(
        self,
        chr: str,
        start: int,
        end: int,
        name: str = ".",
        score: float = 0,
        strand: str = ".",
    ):
        # explicit init with coercion for easier BED-file reading
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = float(score)
        self.strand = strand

    def extend(self, bp: int):
        """
        Return a range extended by *bp* on either end.
        """
        return Range(
            self.chr, self.start - bp, self.end + bp, self.name, self.score, self.strand
        )

    def align(self):
        """
        Return a window aligned to 10bp bins.

        e.g. chr2:54-89 => chr2:51-90
        """
        s = (self.start - 1) // 10 * 10 + 1
        e = (self.end - 1) // 10 * 10 + 10
        return Range(self.chr, s, e, self.name, self.score, self.strand)

    def make_bins(self):
        """
        Return a list of ranges
        """
        r = self.align()
        starts = range(r.start, r.end, 10)
        return [
            Range(self.chr, s, s + 9, self.name, self.score, self.strand)
            for s in starts
        ]

    def __len__(self):
        """
        Return the number of base pairs included in the range. Behaviour is similar to
        GRanges, so an interval "chr1:25-25" has length 1.
        """
        return self.end - self.start + 1
