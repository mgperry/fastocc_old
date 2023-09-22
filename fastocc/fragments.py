from dataclasses import dataclass
import numpy as np
import pandas as pd

from .range import Range
from .options import Options, get_options


@dataclass
class BinnedFragments:
    range: Range
    nfr: np.ndarray
    nuc: np.ndarray
    opts: Options = None

    def to_data_frame(self) -> pd.DataFrame:
        bins = self.range.make_bins(self.opts.bin_size)

        bed = pd.DataFrame(
            [(r.chr, r.start, r.end) for r in bins],
            columns=["chr", "start", "end"]
        )

        bed["nfr"] = self.nfr
        bed["nuc"] = self.nuc

        return bed


@dataclass
class Fragments:
    """
    Class to represent raw ATAC fragments from a single genomic interval.
    """
    range: Range
    midpoint: np.ndarray
    size: np.ndarray
    opts: Options = None

    def bin(self, opts=None):
        """
        Sort ATAC-seq fragments into fixed bins (default 10bp) based on their midpoints, and
        classify them based on insertion length (insertions less than 120 are 
        nucleosome-free. The insertion length is capped at 1000bp to improve efficiency (NB 
        these are rare in any case).
        """
        opts = get_options(self.opts, opts)

        if not self.range.align(opts.bin_size) == self.range:
            raise ValueError("Fragments.range must be aligned to provided bin_size")

        n_bins = len(self.range) // opts.bin_size

        nfr = np.zeros(n_bins, dtype=np.intc)
        nuc = np.zeros(n_bins, dtype=np.intc)

        bin = (self.midpoint - self.range.start + 1) // opts.bin_size

        in_bounds = (bin >= 0) & (bin < n_bins)

        is_nfr = self.size < opts.nfr_cutoff

        is_nuc = self.size >= opts.nfr_cutoff

        np.add.at(nfr, bin[in_bounds & is_nfr], 1)
        np.add.at(nuc, bin[in_bounds & is_nuc], 1)

        return BinnedFragments(
            self.range,
            nfr,
            nuc,
            opts
        )

    def __len__(self):
        return len(self.nfr)


class FragmentLoader:
    """
    Just needs iter -> Range method.
    """
    pass
