from dataclasses import dataclass

import numpy as np
from copy import copy

def gaussian(M, std):
    """
    Adapted from scipy.signal.gaussian
    """
    n = np.arange(0, M) - (M - 1.0) / 2.0
    sig2 = 2 * std * std
    w = np.exp(-n ** 2 / sig2)

    return w

@dataclass
class Options:
    """
    These aren't meant to be changed after release (or by users), but instead to avoid
    magic numbers and allow us to test alternatives if needed. Be aware that as it stands
    changing these numbers is untested, there may be subtle bugs (particularly around
    windowing calculations).
    """
    window_size: int = 5 # number of bins for window estimation
    bin_size: int = 10 # width of bins
    extend: int = 500 # base pairs to extend input ranges
    merge: int = 200 # merge input ranges closer than this (after extension)
    max_insert: int = 300
    smooth: bool = True
    gaussian_sd: float = 1.5
    nfr_cutoff: int = 119 # threshold above which reads are considered to come from nucleosomes
    nfr_ratio: float = 2 # ratio of overall NFR reads to reads < nfr_cutoff
    method: str = "linear"
    p_short: tuple[float] = (0.2, 0.95)
    pseudocount: int = 10 # can be float probably
    pseudocount_nuc: int = 1
    nfr_threshold: float = 0.6
    shift_atac: bool = True
    filter_duplicates: bool = True

    def __post_init__(self):
        self.shift = self.window_size // 2
        self.gaussian = gaussian(self.window_size, self.gaussian_sd)
        self.gaussian_norm = np.convolve(self.gaussian, np.ones(self.window_size), mode="valid")


def get_options(*args):
    if not args:
        return Options()
    
    opts = [arg for arg in args if arg is not None]

    if len(opts) == 0:
        return Options()

    if len(opts) == 1:
        return opts[0]
    
    print("multiple options supplied, take care.")

    return copy(opts[-1])
