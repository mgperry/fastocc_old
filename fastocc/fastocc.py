import click

import pysam
from pybedtools import BedTool

from fastocc.range import Range
from fastocc.occupancy import bin_and_split_reads, occupancy, smooth
from fastocc.options import Options

@click.command()
@click.option("--bam", help="BAM file containing ATAC-seq reads", required=True, type=str)
@click.option("--bed", help="BED file containing ranges", type=str)
@click.option("--pseudocount", help="Pseudocount added to NUC reads (number is for each 10bp bin)", type=float, default=1)
def cli(bam: str, bed: str=None, pseudocount: float=1, raw: bool=False):
    if not pseudocount > 0:
        pseudocount = 0.001 # avoid potential div by zero errors

    # generate sorted, extended and merged bed w/ hacked 1-based coords
    if bed:
        intervals = BedTool(bed).sort().slop(l=Options.extend-1, r=Options.extend, genome="hg38").merge(d=Options.merge)
    else:
        intervals = BedTool().window_maker(genome="hg38", w=1000).slop(l=-1, r=0, genome="hg38")

    b = pysam.AlignmentFile(bam)

    # print column headings
    print("\t".join(["chr", "start", "end", "nfr", "nuc", "window", "smooth"]))

    shift = Options.bins // 2

    for i in intervals:
        r = Range(*i).align()

        nfr, nuc = bin_and_split_reads(b, r.extend((Options.bins-1)*Options.bin_width))
        window_occ = occupancy(nfr, nuc + pseudocount)
        smooth_occ = smooth(window_occ)

        for i, bin in enumerate(r.make_bins()):
            row = (bin.chr, bin.start-1, bin.end, nfr[i+shift*2], nuc[i+shift*2], window_occ[i+shift], smooth_occ[i])
            print("%s\t%d\t%d\t%d\t%d\t%0.2f\t%0.2f" % row)

