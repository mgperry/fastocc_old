import abc
from dataclasses import dataclass

import pysam
import numpy as np

from .range import Range
from .options import Options
from .io import fetch_safe


@dataclass
class Fragments:
    range: Range
    midpoint: np.ndarray
    size: np.ndarray
    opts: Options = None


class FragmentSource:
    pass


class BamFragments(FragmentSource):
    def __init__(self, ranges: list[Range], bam: pysam.AlignmentFile, opts: Options=None) -> None:
        self.opts = Options.check(opts)
        self.ranges = ranges
        self.bam = bam

    def __iter__(self) -> list[Fragments]:
        for r in self.ranges:
            r_ext = r.extend(2*self.opts.shift*self.opts.bin_size)
            yield self.fetch(r_ext)

    def fetch(self, r: Range) -> Fragments:
        """
        Reconstruct full-length ATAC-seq fragments from a bam file. Filters for the first
        read of proper pairs, corrects the Tn5 insertion distance and returns the start and
        length (size) of a fragment in a named tuple.
        """
        q = r.extend(self.opts.max_insert//2) # this is to capture long fragments

        midpoint = []
        size = []

        for read in fetch_safe(self.bam, q):
            if read.is_proper_pair and not read.is_reverse:
                #get left position
                left_insertion = read.pos + 4
                # insert size corrected by 8 base pairs to be insertion to insertion
                frag_size = read.template_length - 8

                if frag_size > self.opts.max_insert:
                    continue

                midpoint.append(left_insertion + frag_size // 2)
                size.append(frag_size)

        return Fragments(
            r,
            np.array(midpoint, dtype=int),
            np.array(size, dtype=int),
            self.opts
        )


@dataclass
class BinnedFragments:
    range: Range
    nfr: np.ndarray
    nuc: np.ndarray
    opts: Options = None

    @classmethod
    def from_fragments(cls, frags: Fragments, opts: Options=None):
        """
        Read ATAC-seq fragments into a matrix of 10bp bins based on their midpoints, and
        classifies them based on insertion length (insertions less than 120 are 
        nucleosome-free. The insertion length is capped at 1000bp to improve efficiency (NB 
        these are rare in any case).
        """
        opts = Options.check(frags.opts, opts)

        if not frags.range.align(opts.bin_size) == frags.range:
            raise ValueError("Fragments.range must be aligned to provided bin_size")

        n_bins = len(frags.range) // opts.bin_size

        nfr = np.zeros(n_bins, dtype=np.intc)
        nuc = np.zeros(n_bins, dtype=np.intc)

        bin = (frags.midpoint - frags.range.start + 1) // opts.bin_size

        in_bounds = (bin >= 0) & (bin < n_bins)

        is_nfr = frags.size < opts.cutoff

        is_nuc = frags.size >= opts.cutoff

        np.add.at(nfr, bin[in_bounds & is_nfr], 1)
        np.add.at(nuc, bin[in_bounds & is_nuc], 1)

        return cls(
            frags.range,
            nfr,
            nuc,
            opts
        )

