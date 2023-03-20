import sys

import click

import pysam
from pybedtools import BedTool

from fastocc.range import Range
from fastocc.occupancy import occupancy
from fastocc.options import Options
from fastocc.pysam_util import fetch_chrom_sizes, temp_chrom_info
import fastocc.io as io


@click.command()
@click.option("--bam", help="BAM file containing ATAC-seq reads", required=True, type=str)
@click.option("--bed", help="BED file containing ranges. If not supplied, calulate occupancy genome-wide.", type=str)
@click.option("-o", "--output", help="Output file location (defaults to stdout).", type=str, default="-")
@click.option("--format", help="Output file format.", type=click.Choice(["bedgraph", "bigwig"]), default="bedgraph")
@click.option("--pseudocount", help="Pseudocount added to NUC reads (number is for each 10bp bin)", type=float, default=1)
def cli(bam: str, bed: str, output: str, format: str, pseudocount: float):
    if not pseudocount > 0:
        pseudocount = 0.001 # avoid potential div by zero errors

    if format == "bigwig" and output == "-":
        raise click.UsageError("Output must be a file if output format is bigwig")

    b = pysam.AlignmentFile(bam)

    if not b.has_index():
        print("No BAM index found, indexing BAM file...", file=sys.stderr)
        pysam.index(b.filename)
        b = pysam.AlignmentFile(bam) # reload index

    chrom_sizes = fetch_chrom_sizes(b.header)

    with temp_chrom_info(chrom_sizes) as genome:

        if bed:
            intervals = BedTool(bed).slop(l=Options.extend, r=Options.extend, g=genome).sort(g=genome).merge(d=Options.merge)
        else:
            intervals = BedTool().window_maker(g=genome, w=10000).slop(l=Options.extend, r=Options.extend, g=genome)

    ranges = (Range.from_interval(i).align(Options.bin_size) for i in intervals)

    if format == "bigwig":
        writer = io.BigWigWriter(output, header=chrom_sizes)
    else:
        writer = io.BedGraphWriter(output)
    
    for r in ranges:
        occ = occupancy(r, b, pseudocount)
        writer.write(r, occ)


