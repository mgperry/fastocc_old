from dataclasses import dataclass, replace
from typing import Iterator

import numpy as np
import pandas as pd

import pyBigWig

from .options import Options, get_options
from .fragments import Fragments, BinnedFragments
from .range import Range

def nonzero_intervals(x: np.ndarray) -> Iterator[tuple]:
    """
    Calculate continuous regions in a numpy array, returning
    a generator of (start, stop) tuples.

    e.g. [0, 0, 1, 1, 1, 0, 1, 1, 0] => (2,5); (6,8)

    NB If the interval starts or ends with 1, this is handled, as a result the indexing
    is basically BED-like, ie (0, 1) represents an interval containing just the first
    item. This should be very rare for our exact use case, it means that the peak caller
    has missed a nucleosome free region.

    e.g. [1, 0, 0, 1, 1] => (0, 1); (3, 5)
    """
    d = np.diff(x)
    idx = d.nonzero()[0] + 1

    if x[0]:
        idx = np.r_[0, idx]

    if x[-1]:
        idx = np.r_[idx, x.size]

    for start, stop in np.reshape(idx, (-1, 2)):
        yield start, stop

# should be a rolling window function and a separate calculation
def rolling_window(arr: np.ndarray, window_size: int):
    """
    Sum counts over overlapping 110bp (or 11 bin) windows. The resulting array is 5 bins
    shorter on either side than the input.
    """
    # create matrix representing 11-bin windows into counts
    index = np.arange(window_size)[None, :] + np.arange(arr.shape[0] - window_size + 1)[:, None]

    # sum over windows
    return np.sum(arr[index], 1)


def occupancy(nfr: np.ndarray, nuc: np.ndarray, opts):
    """
    Calculates nucleosome occupancy from a vector of nucleosome-free and nucleosomal 
    reads, based on the formula (NFR*1.22) / (NFR + NUC).

    This function does not add a pseudocount, so if used outside the main script
    it can cointain NaNs in the output.
    """
    nfr_w = rolling_window(nfr, opts.window_size)
    nuc_w = rolling_window(nuc, opts.window_size)
    
    occ = (nfr_w * opts.nfr_ratio + opts.pseudocount_nuc) / (nfr_w + nuc_w + opts.pseudocount)

    # corrected values can't exceed 1
    return occ.clip(max=1)


def smooth(x: np.ndarray, gaussian, gaussian_norm):
    """
    Smoothes input 1D numpy array using hardcoded gaussian convolution. Output is
    normalised using a hardcoded constant derived from convolving the kernel with a
    vector of 1s. NaNs in x will result in errors. Returns an array with 5 fewer values
    on each side so that only 'valid' results are returned.
    """
    
    smoothed = np.convolve(gaussian, x, mode = 'valid')

    return smoothed / gaussian_norm


def fetch_bw(r: Range, bw):
    bins = bw.intervals(r.chr, r.start, r.end)
    return np.array([b[2] for b in bins])


@dataclass
class Occupancy:
    range: Range
    nfr: np.ndarray
    nuc: np.ndarray
    occ: np.ndarray
    opts: Options

    def to_data_frame(self) -> pd.DataFrame:
        bins = self.range.make_bins(self.opts.bin_size)

        bed = pd.DataFrame(
            [(r.chr, r.start, r.end) for r in bins],
            columns=["chr", "start", "end"]
        )

        bed["nfr"] = self.nfr
        bed["nuc"] = self.nuc
        bed["occ"] = self.occ

        return bed

    def threshold(self, th=None):
        if th is None:
            th = self.opts.nfr_threshold

        coords = nonzero_intervals(self.occ > th)

        r = self.range

        for i, (start, stop) in enumerate(coords):
            name = f"{r.chr}:{r.start}-{r.end}" if r.name == "." else r.name

            yield Range(
                r.chr,
                r.start + self.opts.bin_size*(start),
                r.start + self.opts.bin_size*stop - 1,
                f"{name}:NFR{i+1}",
                np.max(self.occ[start:stop]),
                "."
            )

    @classmethod
    def from_fragments(cls, fragments: Fragments, opts: Options=None):
        return cls.from_binned_fragments(fragments.bin(opts))

    @classmethod
    def from_binned_fragments(cls, fragments: BinnedFragments, opts: Options=None):
        opts = get_options(fragments.opts, opts)

        occ = occupancy(fragments.nfr, fragments.nuc, opts)

        if opts.smooth:
            occ = smooth(occ, opts.gaussian, opts.gaussian_norm)
        else:
            occ = occ[opts.shift:-opts.shift]

        s2 = opts.shift * 2

        r = replace(
            fragments.range.extend(-s2 * opts.bin_size),
            score = np.max(occ)
        )

        return cls(
            r,
            fragments.nfr[s2:-s2],
            fragments.nuc[s2:-s2],
            occ,
            opts,
        )
    
    @classmethod
    def from_bigwig(cls, bw, r: Range, opts: Options=None):
        opts = get_options(opts)
    
        # TODO should check that intervals are the right size
        bins = bw.intervals(r.chr, r.start, r.end)
        occ = np.array([b[2] for b in bins])
        
        return Occupancy(
            replace(r, score=np.max(occ)),
            None,
            None,            
            occ,
            opts
        )

