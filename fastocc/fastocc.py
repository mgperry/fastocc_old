import sys
import click

import pysam
import pyBigWig
from pybedtools import BedTool

from .range import Range
from .occupancy import calculate_occupancy
from .options import Options
from .split import threshold_occupancy
from .fragments import BamFragments
from . import io


@click.group()
def fastocc():
    pass


@click.command()
@click.option("--bam", help="BAM file containing ATAC-seq reads", required=True, type=str)
@click.option("--bed", help="BED file containing ranges. If not supplied, calulate occupancy genome-wide.", type=str)
@click.option("-o", "--output", help="Output file location (defaults to stdout).", type=str, default="-")
@click.option("--format", help="Output file format.", type=click.Choice(["bedgraph", "bigwig"]), default="bedgraph")
@click.option("--pseudocount", help="Pseudocount added to NUC reads (number is for each 10bp bin)", type=float, default=1)
def occ(bam: str, bed: str, output: str, format: str, pseudocount: float):
    if not pseudocount > 0:
        pseudocount = 0.001 # avoid potential div by zero errors

    if format == "bigwig" and output == "-":
        raise click.UsageError("Output must be a file if output format is bigwig")

    b = pysam.AlignmentFile(bam)

    fragments = BamFragments(b)

    if not b.has_index():
        print("No BAM index found, indexing BAM file...", file=sys.stderr)
        pysam.index(b.filename)
        b = pysam.AlignmentFile(bam) # reload index

    chrom_sizes = io.fetch_chrom_sizes(b.header)

    if bed:
        intervals = io.fix_bed(BedTool(bed), chrom_sizes, Options.extend, Options.merge)
    else:
        intervals = io.make_windows(chrom_sizes)

    ranges = (Range.from_interval(i).align(Options.bin_size) for i in intervals)

    if format == "bigwig":
        writer = io.BigWigWriter(output, header=chrom_sizes)
    else:
        writer = io.BedGraphWriter(output)
    
    for r in ranges:
        occ = calculate_occupancy(r, fragments, pseudocount)
        writer.write(r, occ)


@click.command
@click.option("--bigwig", help="Bigwig file containing `fastocc occ` output.", required=True, type=str)
@click.option("--bed", help="BED file containing ranges.", required=True, type=str)
@click.option("--output", help="Output file location (defaults to stdout).", type=str, default="-")
@click.option("--threshold", help="Threshold for calling nucleosome free regions.", type=float, default=0.6)
@click.option("--merge/--no-merge", help="Disable extending and merging input BED: user must guarantee non-overlapping ranges", default=True)
def nfr(bigwig: str, bed: str, output: str, threshold: float = 0.6, merge: bool=True):
    bw = pyBigWig.open(bigwig)

    chrom_sizes = bw.chroms()

    if merge:
        intervals = io.fix_bed(BedTool(bed), chrom_sizes, Options.extend, Options.merge)
    else:
        with io.temp_chrom_info(chrom_sizes) as genome:
            intervals = BedTool(bed).sort(g=genome)

    ranges = (Range.from_interval(i).align(Options.bin_size) for i in intervals)

    writer = io.BEDWriter(output)

    for r in ranges:
        reduced_peaks = threshold_occupancy(r, bw, threshold)
        writer.writerows(reduced_peaks)


fastocc.add_command(occ)
fastocc.add_command(nfr)
