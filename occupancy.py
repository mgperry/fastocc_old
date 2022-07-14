import numpy as np
import pysam
from bed import Range
import scipy.signal as signal

def run_occupancy(b: pysam.AlignmentFile, r: Range, max_ilen=1000, window=120, pc=1):
    """
    Helper function to get reads, calculated windowed + binned occupancy and smooth
    the result. It also wraps the windowing necessary to do this accurately so nothing
    else needs to know about this.
    """

    nbin = 11

    r_ext = r.extend((nbin-1)*10)
    nfr, nuc = bin_and_split_reads(b, r_ext, max_ilen)

    nuc = nuc + pc

    binned_occ = occupancy(nfr, nuc, nbin, pc)

    return smooth(binned_occ)


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
            
            midpoint = (l_pos + (ilen-1)//2)
            pos_in_read = midpoint - r.start + 1
            
            if pos_in_read > 0 and pos_in_read < len(r):
                if ilen <= 120:
                    nfr[pos_in_read // 10] += 1
                else:
                    nuc[pos_in_read // 10] += 1

    return nfr, nuc


def occupancy(nfr, nuc):
    # fancy indexing into 1D arrays
    index = np.arange(11)[None, :] + np.arange(nfr.shape[0] - 10)[:, None]
    sum_nfr = np.sum(nfr[index], 1)
    sum_nuc = np.sum(nuc[index], 1)
    occ = (sum_nfr * 1.22) / (sum_nuc + sum_nfr)
    return occ.clip(max=1)

GAUSS = np.array([0.04393693, 0.13533528, 0.32465247, 0.60653066, 0.8824969, 1., 0.8824969 , 0.60653066, 0.32465247, 0.13533528, 0.04393693]) # signal.gaussian(11, 2)

def smooth(x):
    """
    smoothes input signal using gaussian window, hardcoded as 
    `scipy.signal.gaussian(11, 2)`. NaNs in x will result in errors.
    """
    
    smoothed = np.convolve(GAUSS, x, mode = 'valid')
    norm = np.convolve(GAUSS, np.ones(len(x)), mode = 'valid')

    return smoothed / norm
