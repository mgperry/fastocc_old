import argparse
import pysam
from pybedtools import BedTool
import numpy as np

from range import Range
from occupancy import bin_and_split_reads, occupancy, smooth

parser = argparse.ArgumentParser()
parser.add_argument("--bed", help="BED file containing ranges")
parser.add_argument("--bam", help="BAM file containing ATAC-seq reads")
parser.add_argument("--pseudocount", help="Pseudocount added to NUC reads (number is for each 10bp bin)", type=float, default=1)


if __name__=="__main__":
    args = parser.parse_args()

    if not args.pseudocount > 0:
        args.pseudocount = 0.001 # avoid potential div by zero errors

    # generate sorted, extended and merged bed w/ hacked 1-based coords
    temp_bed = BedTool(args.bed).sort().slop(l=499, r=500, genome="hg38").merge(d=200)

    # read in the resulting intervals and align to 10bp bins
    ranges = [Range(*i).align() for i in temp_bed]

    bam = pysam.AlignmentFile(args.bam)

    nbin = 11

    # print column headings
    print("\t".join(["chr", "start", "end", "nfr", "nuc", "window", "smooth"]))

    for r in ranges:

        r_ext = r.extend((nbin-1)*10)

        nfr, nuc = bin_and_split_reads(bam, r_ext)

        window_occ = occupancy(nfr, nuc + args.pseudocount)

        smooth_occ = smooth(window_occ)
        
        window_occ = np.pad(window_occ, 5)
        smooth_occ = np.pad(smooth_occ, 10)
        
        for i, bin in enumerate(r_ext.make_bins()):
            row = (bin.chr, bin.start-1, bin.end, nfr[i], nuc[i], window_occ[i], smooth_occ[i])
            print("%s\t%d\t%d\t%d\t%d\t%0.2f\t%0.2f" % row)

