import numpy as np
from copy import copy

class Options:
    """
    These aren't meant to be changed after release (or by users), but instead to avoid
    magic numbers and allow us to test alternatives if needed. Be aware that as it stands
    changing these numbers is untested, there may be subtle bugs (particularly around
    windowing calculations).
    """
    window_size = 11 # number of bins for window estimation
    bin_size = 10 # width of bins
    extend = 500 # base pairs to extend input ranges
    merge = 200 # merge input ranges closer than this (after extension)
    max_insert = 1000
    # scipy.signal.gaussian(bins, 2) 
    gaussian = np.array([0.04393693, 0.13533528, 0.32465247, 0.60653066, 0.8824969, 1., 0.8824969 , 0.60653066, 0.32465247, 0.13533528, 0.04393693])
    # gaussian_norm = 4.98590449 # np.convolve(gaussian, np.ones(bins), mode='valid')
    nfr_cutoff = 119 # threshold above which reads are considered to come from nucleosomes
    nfr_ratio = 1.22 # ratio of overall NFR reads to reads < nfr_cutoff
    p_short = [0.2, 0.95]
    pseudocount = 1
    nfr_threshold = 0.6

    def __init__(self, **kwargs):
        for attr, value in kwargs.items():
            setattr(self, attr, value)

        # could also initialise guassian here

        self.shift = self.bin_size // 2
        self.gaussian_norm = np.convolve(self.gaussian, np.ones(self.bins), mode="valid")

    def replace(self, **kwargs):
        opts = copy(self)

        for attr, value in kwargs.items():
            setattr(opts, attr, value)

        return opts

    @classmethod
    def check(cls, *args):
        if not args:
            return cls()
        
        opts =  [arg for arg in args if arg is not None]

        if len(opts) == 0:
            return cls()

        if len(opts) == 1:
            return opts[0]
        
        print("multiple options supplied, take care.")

        return copy(opts[-1])
