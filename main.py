from bed import read_bed, bin_range
import argparse
import pysam
from occupancy import calc_occ_bins, get_binned_reads, smooth

parser = argparse.ArgumentParser()
parser.add_argument("--bed", help="BED file containing ranges")
parser.add_argument("--bam", help="BAM file containing ATAC-seq reads")

if __name__=="__main__":
    args = parser.parse_args()

    ranges = read_bed(args.bed)

    bam_file = pysam.AlignmentFile(args.bam)

    for r in ranges:
        # extend region and align to grid
        r1, binned_reads = get_binned_reads(bam_file, r)
        binned_occ = calc_occ_bins(binned_reads, 12)
        smooth_occ = smooth(binned_occ)
        
        for i, r_bin in enumerate(bin_range(r1)):
            try:
                r_bin.score = smooth_occ[i]
                print(r_bin)
            except:
                break # this is just because i can't count


