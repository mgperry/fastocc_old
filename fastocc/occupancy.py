import numpy as np
import pysam
from fastocc.range import Range
from fastocc.options import Options

def bin_and_split_reads(bam: pysam.AlignmentFile, r: Range):
    """
    Read ATAC-seq fragments into a matrix of 10bp bins based on their midpoints, and
    classifies them based on insertion length (insertions less than 120 are 
    nucleosome-free. The insertion length is capped at 1000bp to improve efficiency (NB 
    these are rare in any case).
    """

    q = r.extend(Options.max_insert//2) # this is to capture long fragments

    n_bins = len(r) // Options.bin_width

    nfr = np.zeros(n_bins)
    nuc = np.zeros(n_bins)

    for read in bam.fetch(q.chr, q.start, q.end):
        if read.is_proper_pair and not read.is_reverse:
            #get left position
            left_insertion = read.pos + 4
            # insert size corrected by 8 base pairs to be inserion to insertion
            fragment_size = read.template_length - 8

            if fragment_size > Options.max_insert:
                continue
            
            midpoint = left_insertion + (fragment_size-1) // 2
            
            bin = (midpoint - r.start + 1) // Options.bin_width
            
            if bin >= 0 and bin < n_bins:
                if fragment_size <= 120:
                    nfr[bin] += 1
                else:
                    nuc[bin] += 1

    return nfr, nuc


def occupancy(nfr: np.array, nuc: np.array):
    """
    Calculates nucleosome occupancy from a vector of nucleosome-free and nucleosomal 
    reads, based on the formula (NFR*1.22) / (NFR + NUC). This is calculated for
    overlapping 110bp (or 11 bin) windows. Output is clipped so occupancy values
    cannot exceed 1, and the resulting array is 5 bins shorter on either side than the
    input due to windowing.

    NB This function does not add a pseudocount, so if used outside the main script
    it can cointain NaNs in the output.
    """
    # create matrix representing 11-bin windows into counts
    index = np.arange(Options.bins)[None, :] + np.arange(nfr.shape[0] - Options.bins + 1)[:, None]

    # sum over windows
    nfr_window_sums = np.sum(nfr[index], 1)
    nuc_window_sums = np.sum(nuc[index], 1)

    # super secret occupancy formula
    occ = (nfr_window_sums * 1.22) / (nfr_window_sums + nuc_window_sums)

    # corrected values can't exceed 1
    return occ.clip(max=1)


def smooth(x: np.array):
    """
    Smoothes input 1D numpy array using hardcoded gaussian convolution. Output is
    normalised using a hardcoded constant derived from convolving the kernel with a
    vector of 1s. NaNs in x will result in errors. Returns an array with 5 fewer values
    on each side so that only 'valid' results are returned.
    """
    
    smoothed = np.convolve(Options.gaussian, x, mode = 'valid')

    return smoothed / Options.gaussian_norm
