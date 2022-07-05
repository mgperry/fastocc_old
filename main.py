import argparse
import pysam
import sys
from pybedtools import BedTool

from bed import Range
from occupancy import occupancy

parser = argparse.ArgumentParser()
parser.add_argument("--bed", help="BED file containing ranges")
parser.add_argument("--bam", help="BAM file containing ATAC-seq reads")

if __name__=="__main__":
    args = parser.parse_args()

    # generate sorted, extended and merged bed w/ hacked 1-based coords
    temp_bed = BedTool(args.bed).sort().slop(l=499, r=500, genome="hg38").merge(d=100)
    # read in the resulting intervals and align to 10bp bins
    ranges = [Range(*i).align() for i in temp_bed]

    bam = pysam.AlignmentFile(args.bam)

    for r in ranges:
        occ = occupancy(bam, r)
        
        for i, b in enumerate(r.bin()):
            b.score = occ[i]
            print(b.bed_graph())

