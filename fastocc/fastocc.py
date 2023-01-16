import click

import pysam
from pybedtools import BedTool
import pyBigWig

from fastocc.range import interval_to_range
from fastocc.occupancy import bin_and_split_reads, occupancy, smooth
from fastocc.options import Options
from fastocc.pysam_util import fetch_chrom_sizes


@click.command()
@click.option("--bam", help="BAM file containing ATAC-seq reads", required=True, type=str)
@click.option("--bed", help="Sorted BED file containing ranges. If not supplied, calulate occupancy genome-wide.", type=str)
@click.option("--bigwig", help="Output smooth occupancy to bigwig file.", type=str)
@click.option("--pseudocount", help="Pseudocount added to NUC reads (number is for each 10bp bin)", type=float, default=1)
def cli(bam: str, bed: str=None, bigwig: str="", pseudocount: float=1):
    if not pseudocount > 0:
        pseudocount = 0.001 # avoid potential div by zero errors

    b = pysam.AlignmentFile(bam)

    chrom_sizes = fetch_chrom_sizes(b.header)

    if bed:
        intervals = BedTool(bed).slop(l=Options.extend, r=Options.extend, g=chrom_sizes).sort(g=chrom_sizes).merge(d=Options.merge)
    else:
        # TODO test to see if this creates missign edges, may have to extend by half bin width
        intervals = BedTool().window_maker(g=chrom_sizes, w=1000)

    # column headings
    if not bigwig:
        print("\t".join(["chr", "start", "end", "nfr", "nuc", "window", "smooth"]))

    if bigwig:
        bw = pyBigWig.open(bigwig, 'w')
        bw.addHeader([(chr, end) for chr, (_, end) in chrom_sizes.items()]) # different format req'd

    shift = Options.bins // 2

    for i in intervals:
        r = interval_to_range(i).align()

        nfr, nuc = bin_and_split_reads(b, r.extend(2*shift*Options.bin_size))
        window_occ = occupancy(nfr, nuc + pseudocount)
        smooth_occ = smooth(window_occ)

        if bigwig:
            bw.addEntries(r.chr, r.start-1, values=smooth_occ.tolist(), span=Options.bin_size, step=Options.bin_size)
        else:
            for j, bin in enumerate(r.make_bins()):
                row = (bin.chr, bin.start-1, bin.end, nfr[j+shift*2], nuc[j+shift*2], window_occ[j+shift], smooth_occ[j])
                print("%s\t%d\t%d\t%d\t%d\t%0.2f\t%0.2f" % row)

