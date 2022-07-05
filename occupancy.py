import numpy as np
import pysam
from bed import Range
import scipy.signal as signal

def occupancy(b: pysam.AlignmentFile, r: Range, max_len=1000, window=120):
    """
    Helper function to get reads, calculated windowed + binned occupancy and smooth
    the result. It also wraps the windowing necessary to do this accurately so nothing
    else needs to know about this.
    """
    # bins in window
    wb = (window - 1) // 10 + 1
    # extend r to acconut for window + smoothing
    nfr, nuc = bin_and_split_reads(b, r.extend((wb-1)*10), max_len)
    binned_occ = calc_occ(nfr, nuc, wb)
    smooth_occ = smooth(binned_occ)

    return smooth_occ


def bin_and_split_reads(bam: pysam.AlignmentFile, r: Range, max_ilen=1000):
    """
    Read ATAC-seq fragments into a matrix of 10bp bins.
    """

    q = r.extend(max_ilen) # this is to capture long fragments

    nfr = np.zeros(len(r)//10, dtype=int)
    nuc = np.zeros(len(r)//10, dtype=int)

    for read in bam.fetch(q.chr, q.start, q.end):
        if read.is_proper_pair and not read.is_reverse:
            #get left position
            l_pos = read.pos + 4
            # insert size corrected by 8 base pairs to be inserion to insertion
            ilen = abs(read.template_length) - 8

            if ilen > max_ilen:
                continue
            
            midpoint = (l_pos + ilen)
            pos_in_read = midpoint - r.start + 1
            
            if pos_in_read > 0 and pos_in_read < len(r):
                if ilen <= 120:
                    nfr[pos_in_read // 10] += 1
                else:
                    nuc[pos_in_read // 10] += 1

    return nfr, nuc


def calc_occ(nfr, nuc, n):
    # fancy indexing into 1D arrays
    index = np.arange(n)[None, :] + np.arange(nfr.shape[0] - n + 1)[:, None]
    sum_nfr = np.sum(nfr[index])
    sum_nuc = np.sum(nuc[index])
    pc = 1 # pseudocount avoids division by zero
    return ((nfr / (nuc + nfr + pc)) * 1.22).clip(max=1)


def smooth(x, window_len=12, sd=2):
    """
    smoothes input signal using gaussian window. NaNs in x will result in errors.
    """
    w = signal.gaussian(window_len, sd) # TODO avoid full scipy dep, hardcode matrix?
    smoothed = np.convolve(w, x, mode = 'valid')
    norm = np.convolve(w, np.ones(len(x)), mode = 'valid')

    return smoothed / norm
