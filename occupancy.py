import numpy as np
import pysam
from bed import Range
import scipy.signal as signal

def round(x):
    return x - (x % 10)


def get_binned_reads(b: pysam.AlignmentFile, r: Range, max_len=1000):
    """
    Read ATAC-seq fragments into a matrix of 10bp bins.
    """

    r_bin = Range(r.chr, round(r.start), round(r.end + 10) - 1)

    bins = np.zeros((2, len(r_bin)//10 + 1), dtype=int)

    r_ext = r_bin.extend(max_len)

    for read in b.fetch(r_ext.chr, r_ext.start, r_ext.end):
        if read.is_proper_pair and not read.is_reverse:
            #get left position
            l_pos = read.pos + 4
            #get insert size
            #correct by 8 base pairs to be inserion to insertion
            ilen = abs(read.template_length) - 8
            
            midpoint = (l_pos + ilen)
            pos_in_read = midpoint - r.start + 1
            
            if pos_in_read > 0 and pos_in_read < len(r):
                if ilen <= 120:
                    bins[0, pos_in_read // 10] += 1
                else:
                    bins[1, pos_in_read // 10] += 1

    return r_bin, bins


def calc_occ_bins(bins, n):
    # fancy indexing into 1D arrays
    index = np.arange(n)[None, :] + np.arange(bins.shape[1] - n + 1)[:, None]
    nfr = bins[0, :][index].sum(1)
    nuc = bins[1, :][index].sum(1)
    return ((nfr / (nuc + nfr)) * 1.22).clip(max=1)


def smooth(x, window_len=12, sd=2):
    """
    smoothes input signal using gaussian window. NaNs in x will result in errors.
    """
    w = signal.gaussian(window_len, sd) # TODO avoid full scipy dep, hardcode matrix?
    smoothed = np.convolve(w, x, mode = 'valid')
    norm = np.convolve(w, np.ones(len(x)), mode = 'valid')

    return smoothed / norm