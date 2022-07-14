import numpy as np
import pysam
from range import Range


def bin_and_split_reads(bam: pysam.AlignmentFile, r: Range):
    """
    Read ATAC-seq fragments into a matrix of 10bp bins.
    """

    max_ilen = 1000

    q = r.extend(max_ilen//2) # this is to capture long fragments

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
    # create matrix representing 11-bin windows into counts
    index = np.arange(11)[None, :] + np.arange(nfr.shape[0] - 10)[:, None]

    nfr_window_sums = np.sum(nfr[index], 1)
    nuc_window_sums = np.sum(nuc[index], 1)

    # super secret occupancy formula
    occ = (nfr_window_sums * 1.22) / (nfr_window_sums + nuc_window_sums)

    return occ.clip(max=1) # corrected values can't exceed 1


# scipy.signal.gaussian(11, 2)
GAUSS = np.array([0.04393693, 0.13533528, 0.32465247, 0.60653066, 0.8824969, 1., 0.8824969 , 0.60653066, 0.32465247, 0.13533528, 0.04393693])

# np.convolve(GAUSS, np.ones(11), mode='valid')
GAUSS_NORM = 4.98590449

def smooth(x):
    """
    smoothes input signal using hardcoded gaussian function (see source). NaNs in x will result in errors.

    Output is normalised using a hardcoded constant derived from convolving the kernel with a vector of 1s.
    """
    
    smoothed = np.convolve(GAUSS, x, mode = 'valid')

    return smoothed / GAUSS_NORM
