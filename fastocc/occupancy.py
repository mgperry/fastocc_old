import numpy as np
import pysam
from collections import namedtuple
from collections.abc import Iterable

from .range import Range
from .options import Options
from .io import fetch_safe


Fragment = namedtuple("Fragment", ["start", "size"])


def fetch_fragments(bam: pysam.AlignmentFile, r: Range) -> Iterable[Fragment]:
    """
    Reconstruct full-length ATAC-seq fragments from a bam file. Filters for the first
    read of proper pairs, corrects the Tn5 insertion distance and returns the start and
    length (size) of a fragment in a named tuple.
    """
    q = r.extend(Options.max_insert//2) # this is to capture long fragments

    for read in fetch_safe(bam, q):
        if read.is_proper_pair and not read.is_reverse:
            #get left position
            left_insertion = read.pos + 4
            # insert size corrected by 8 base pairs to be insertion to insertion
            fragment_size = read.template_length - 8

            if fragment_size > Options.max_insert:
                continue
            
            yield Fragment(left_insertion, fragment_size)


def bin_and_split_reads(bam: pysam.AlignmentFile, r: Range) -> tuple[np.ndarray, np.ndarray]:
    """
    Read ATAC-seq fragments into a matrix of 10bp bins based on their midpoints, and
    classifies them based on insertion length (insertions less than 120 are 
    nucleosome-free. The insertion length is capped at 1000bp to improve efficiency (NB 
    these are rare in any case).
    """

    n_bins = len(r) // Options.bin_size

    nfr = np.zeros(n_bins)
    nuc = np.zeros(n_bins)

    for f in fetch_fragments(bam, r):           

        midpoint = f.start + (f.size - 1) // 2
        
        bin = (midpoint - r.start + 1) // Options.bin_size
        
        if bin >= 0 and bin < n_bins:
            if f.size <= Options.nfr_cutoff:
                nfr[bin] += 1
            else:
                nuc[bin] += 1

    return nfr, nuc


def raw_occ(nfr: np.ndarray, nuc: np.ndarray):
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
    occ = (nfr_window_sums * Options.nfr_ratio) / (nfr_window_sums + nuc_window_sums)

    # corrected values can't exceed 1
    return occ.clip(max=1)


def smooth(x: np.ndarray):
    """
    Smoothes input 1D numpy array using hardcoded gaussian convolution. Output is
    normalised using a hardcoded constant derived from convolving the kernel with a
    vector of 1s. NaNs in x will result in errors. Returns an array with 5 fewer values
    on each side so that only 'valid' results are returned.
    """
    
    smoothed = np.convolve(Options.gaussian, x, mode = 'valid')

    return smoothed / Options.gaussian_norm


def occupancy(r: Range, bam: pysam.AlignmentFile, pseudocount: int) -> dict[np.ndarray]:
    shift = Options.bins // 2

    # double shift for windowing and smoothing
    nfr, nuc = bin_and_split_reads(bam, r.extend(2*shift*Options.bin_size))
    occ = raw_occ(nfr, nuc + pseudocount)

    return {
        "nfr": nfr[2*shift:-2*shift],
        "nuc": nuc[2*shift:-2*shift],
        "window": occ[shift:-shift],
        "smooth": smooth(occ)
    }

