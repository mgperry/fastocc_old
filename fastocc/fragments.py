from .range import Range
from .options import Options
from .io import fetch_safe
from collections import namedtuple
import pysam
from collections.abc import Iterable

Fragment = namedtuple("Fragment", ["start", "size"])

class Fragments:
    pass


class BamFragments:
    def __init__(self, bam: pysam.AlignmentFile) -> None:
        self.bam = bam

    def fetch(self, r: Range) -> Iterable[Fragment]:
        """
        Reconstruct full-length ATAC-seq fragments from a bam file. Filters for the first
        read of proper pairs, corrects the Tn5 insertion distance and returns the start and
        length (size) of a fragment in a named tuple.
        """
        q = r.extend(Options.max_insert//2) # this is to capture long fragments

        for read in fetch_safe(self.bam, q):
            if read.is_proper_pair and not read.is_reverse:
                #get left position
                left_insertion = read.pos + 4
                # insert size corrected by 8 base pairs to be insertion to insertion
                fragment_size = read.template_length - 8

                if fragment_size > Options.max_insert:
                    continue
                
                yield Fragment(left_insertion, fragment_size)
