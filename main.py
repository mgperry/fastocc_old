import argparse
import pysam
from pybedtools import BedTool
import numpy as np

from bed import Range
from occupancy import bin_and_split_reads, occupancy, smooth

parser = argparse.ArgumentParser()
parser.add_argument("--bed", help="BED file containing ranges")
parser.add_argument("--bam", help="BAM file containing ATAC-seq reads")
parser.add_argument("--pseudocount", help="Pseudocount added to NUC reads (number is for each 10bp bin)", type=float, default=1)


if __name__=="__main__":
    args = parser.parse_args()

    # generate sorted, extended and merged bed w/ hacked 1-based coords
    temp_bed = BedTool(args.bed).sort().slop(l=499, r=500, genome="hg38").merge(d=200)

    # read in the resulting intervals and align to 10bp bins
    ranges = [Range(*i).align() for i in temp_bed]

    bam = pysam.AlignmentFile(args.bam)

    for r in ranges:

        nbin = 11

        r_ext = r.extend((nbin-1)*10)

        nfr, nuc = bin_and_split_reads(bam, r_ext)

        nuc = nuc + args.pseudocount 

        window_occ = occupancy(nfr, nuc)

        smooth_occ = smooth(window_occ)

        scores = np.vstack((
            nfr,
            nuc,
            np.pad(window_occ, 5),
            np.pad(smooth_occ, 10)
        ))
        
        for i, bin in enumerate(r_ext.make_bins()):
            cols = "%d\t%d\t%0.2f\t%0.2f" % tuple(scores[:,i])
            print(f"{bin.chr}\t{bin.start-1}\t{bin.end}\t" + cols)

