from dataclasses import dataclass, astuple

import numpy as np
import pandas as pd

from .options import Options
from .fragments import Fragments
from .range import Range

def bin_and_split_fragments(frags: Fragments, bin_size, cutoff) -> tuple[np.ndarray, np.ndarray]:
    """
    Read ATAC-seq fragments into a matrix of 10bp bins based on their midpoints, and
    classifies them based on insertion length (insertions less than 120 are 
    nucleosome-free. The insertion length is capped at 1000bp to improve efficiency (NB 
    these are rare in any case).
    """

    if not frags.range.align(bin_size) == frags.range:
        raise ValueError("Fragments.range must be aligned to provided bin_size")

    n_bins = len(frags.range) // bin_size

    nfr = np.zeros(n_bins, dtype=np.intc)
    nuc = np.zeros(n_bins, dtype=np.intc)

    bin = (frags.midpoint - frags.range.start + 1) // bin_size

    in_bounds = (bin >= 0) & (bin < n_bins)

    is_nfr = frags.size < cutoff

    is_nuc = frags.size >= cutoff

    np.add.at(nfr, bin[in_bounds & is_nfr], 1)
    np.add.at(nuc, bin[in_bounds & is_nuc], 1)

    return nfr, nuc

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


def calculate_occupancy(nfr: np.ndarray, nuc: np.ndarray, nfr_ratio: int):
    """
    Calculates nucleosome occupancy from a vector of nucleosome-free and nucleosomal 
    reads, based on the formula (NFR*1.22) / (NFR + NUC).

    This function does not add a pseudocount, so if used outside the main script
    it can cointain NaNs in the output.
    """
    occ = (nfr * nfr_ratio) / (nfr + nuc)

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


@dataclass
class Occupancy:
    range: Range
    data: pd.DataFrame
    options: Options

    def to_data_frame(self) -> pd.DataFrame:
        bins = self.range.make_bins(self.options.bin_size)

        bed = pd.DataFrame(
            [(r.chr, r.start, r.end) for r in bins],
            columns=["chr", "start", "end"]
        )

        return pd.concat([bed, self.data], axis=1)


def occ_from_fragments(fragments: Fragments, options: Options = None) -> Occupancy:
    if options is None:
        options = Options()

    nfr, nuc = bin_and_split_fragments(fragments, options.bin_size, options.nfr_cutoff)

    nfr_windows = rolling_window(nfr, options.window_size)
    nuc_windows = rolling_window(nuc + options.pseudocount, options.window_size)

    occ = calculate_occupancy(nfr_windows, nuc_windows, options.nfr_ratio)
    smoothed = smooth(occ, options.gaussian, options.gaussian_norm)

    s = options.shift

    data = pd.DataFrame({
        "nfr": nfr[2*s:-2*s],
        "nuc": nuc[2*s:-2*s],
        "occ": occ[s:-s],
        "smooth": smoothed,
    })

    shift_bp = 2 * s * options.bin_size

    # if we are clipping data, we should clip range
    return Occupancy(
        fragments.range.extend(-shift_bp),
        data,
        options,
    )