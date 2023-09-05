from dataclasses import dataclass

import numpy as np

from .range import Range
from .options import Options, get_options


@dataclass
class Fragments:
    range: Range
    midpoint: np.ndarray
    size: np.ndarray
    opts: Options = None


class FragmentLoader:
    """
    Just needs iter -> Range method.
    """
    pass


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
        opts = get_options(frags.opts, opts)

        if not frags.range.align(opts.bin_size) == frags.range:
            raise ValueError("Fragments.range must be aligned to provided bin_size")

        n_bins = len(frags.range) // opts.bin_size

        nfr = np.zeros(n_bins, dtype=np.intc)
        nuc = np.zeros(n_bins, dtype=np.intc)

        bin = (frags.midpoint - frags.range.start + 1) // opts.bin_size

        in_bounds = (bin >= 0) & (bin < n_bins)

        is_nfr = frags.size < opts.nfr_cutoff

        is_nuc = frags.size >= opts.nfr_cutoff

        np.add.at(nfr, bin[in_bounds & is_nfr], 1)
        np.add.at(nuc, bin[in_bounds & is_nuc], 1)

        return cls(
            frags.range,
            nfr,
            nuc,
            opts
        )

