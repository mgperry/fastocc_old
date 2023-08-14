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


class FragmentSource:
    pass


class BamFragments(FragmentSource):
    def __init__(self, bam: pysam.AlignmentFile) -> None:
        self.bam = bam

    def fetch(self, r: Range) -> Fragments:
        """
        Reconstruct full-length ATAC-seq fragments from a bam file. Filters for the first
        read of proper pairs, corrects the Tn5 insertion distance and returns the start and
        length (size) of a fragment in a named tuple.
        """
        q = r.extend(Options.max_insert//2) # this is to capture long fragments

        midpoint = []
        size = []

        for read in fetch_safe(self.bam, q):
            if read.is_proper_pair and not read.is_reverse:
                #get left position
                left_insertion = read.pos + 4
                # insert size corrected by 8 base pairs to be insertion to insertion
                frag_size = read.template_length - 8

                if frag_size > Options.max_insert:
                    continue

                midpoint.append(left_insertion + frag_size // 2)
                size.append(frag_size)

        return Fragments(
            r,
            np.array(midpoint, dtype=int),
            np.array(size, dtype=int),
        )
