import click

import pysam
from pybedtools import BedTool

from fastocc.range import interval_to_range
from fastocc.occupancy import bin_and_split_reads, occupancy, smooth
from fastocc.options import Options
from fastocc.pysam_util import fetch_chrom_sizes


@click.command()
@click.option("--bam", help="BAM file containing ATAC-seq reads", required=True, type=str)
@click.option("--bed", help="BED file containing ranges. If not supplied, calulate occupancy genome-wide.", type=str)
@click.option("--pseudocount", help="Pseudocount added to NUC reads (number is for each 10bp bin)", type=float, default=1)
def cli(bam: str, bed: str=None, pseudocount: float=1):
    if not pseudocount > 0:
        pseudocount = 0.001 # avoid potential div by zero errors

    b = pysam.AlignmentFile(bam)

    # extract chrom sizes from bam header into pybedtools format
    chrom_sizes = fetch_chrom_sizes(b.header)

    if bed:
        intervals = BedTool(bed).sort().slop(l=Options.extend, r=Options.extend, g=chrom_sizes).merge(d=Options.merge)
    else:
        intervals = BedTool().window_maker(g=chrom_sizes, w=1000)

    # print column headings
    print("\t".join(["chr", "start", "end", "nfr", "nuc", "window", "smooth"]))

    shift = Options.bins // 2

    for i in intervals:
        r = interval_to_range(i).align()

        nfr, nuc = bin_and_split_reads(b, r.extend((Options.bins-1)*Options.bin_width))
        window_occ = occupancy(nfr, nuc + pseudocount)
        smooth_occ = smooth(window_occ)

        for i, bin in enumerate(r.make_bins()):
            row = (bin.chr, bin.start-1, bin.end, nfr[i+shift*2], nuc[i+shift*2], window_occ[i+shift], smooth_occ[i])
            print("%s\t%d\t%d\t%d\t%d\t%0.2f\t%0.2f" % row)

