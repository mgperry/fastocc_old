from dataclasses import dataclass
import json

import numpy as np
from copy import copy

from .np_utils import gaussian

@dataclass
class Options:
    """
    These aren't meant to be changed after release (or by users), but instead to avoid
    magic numbers and allow us to test alternatives if needed. Be aware that as it stands
    changing these numbers is untested, there may be subtle bugs (particularly around
    windowing calculations).
    """
    # bins + windows
    window_size: int = 5 # number of bins for window estimation
    bin_size: int = 10 # width of bins
    extend: int = 500 # base pairs to extend input ranges
    merge: int = 200 # merge input ranges closer than this (after extension)
    # ATAC fragment options
    max_insert: int = 300
    shift_atac: bool = True
    filter_duplicates: bool = True
    # occupancy params
    nfr_cutoff: int = 119 # threshold above which reads are considered to come from nucleosomes
    nfr_ratio: float = 2.0 # ratio of overall NFR reads to reads < nfr_cutoff
    pseudocount: int = 20 
    pseudocount_nuc: int = 1
    # smoothing params
    smooth: bool = True
    gaussian_sd: float = 1.5
    # nfr params
    nfr_threshold: float = 0.6
    # nfr_min: float = 0.8
    # nfr_extend: int = 20

    def __post_init__(self):
        self.shift = self.window_size // 2
        self.gaussian = gaussian(self.window_size, self.gaussian_sd)
        self.gaussian_norm = np.convolve(self.gaussian, np.ones(self.window_size), mode="valid")

    @classmethod
    def load(cls, f: str):
        jsn = json.load(open(f))

        return cls(**jsn)



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
