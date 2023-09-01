import sys
import click

import pysam
import pyBigWig

from .range import BedFile
from .occupancy import Occupancy
from .options import Options
from .fragments import BamFragments
from .io import BigWigWriter, BedGraphWriter, BEDWriter, fetch_chrom_sizes


@click.group()
def fastocc():
    pass


@click.command()
@click.option("--bam", help="BAM file containing ATAC-seq reads", required=True, type=str)
@click.option("--bed", help="BED file containing ranges.", type=str)
@click.option("-o", "--output", help="Output file location (defaults to stdout).", type=str, default="-")
@click.option("--format", help="Output file format.", type=click.Choice(["bedgraph", "bigwig"]), default="bedgraph")
@click.option("--pseudocount", help="Pseudocount added to NUC reads (number is for each 10bp bin)", type=float, default=1)
def occ(bam: str, bed: str, output: str, format: str, pseudocount: float):
    if not pseudocount > 0:
        pseudocount = 0.001 # avoid potential div by zero errors

    opts = Options(pseudocount=pseudocount)

    if format == "bigwig" and output == "-":
        raise click.UsageError("Output must be a file if output format is bigwig")

    bam_file = pysam.AlignmentFile(bam)

    if not bam_file.has_index():
        print("No BAM index found, indexing BAM file...", file=sys.stderr)
        pysam.index(bam_file.filename)
        bam_file = pysam.AlignmentFile(bam) # reload index

    chrom_sizes = fetch_chrom_sizes(bam_file.header)

    ranges = BedFile(bed, chrom_sizes, opts)

    fragments = BamFragments(bam_file, ranges, opts)

    if format == "bigwig":
        writer = BigWigWriter(output, header=chrom_sizes)
    else:
        writer = BedGraphWriter(output)
    
    for frags in fragments:
        occupancy = Occupancy.from_fragments(frags)
        writer.write(occupancy)


@click.command
@click.option("--bigwig", help="Bigwig file containing `fastocc occ` output.", required=True, type=str)
@click.option("--bed", help="BED file containing ranges.", required=True, type=str)
@click.option("--output", help="Output file location (defaults to stdout).", type=str, default="-")
@click.option("--threshold", help="Threshold for calling nucleosome free regions.", type=float, default=0.6)
@click.option("--merge/--no-merge", help="Disable extending and merging input BED: user must guarantee non-overlapping ranges", default=True)
def nfr(bigwig: str, bed: str, output: str, threshold: float = 0.6, merge: bool=True):
    bw = pyBigWig.open(bigwig)

    chrom_sizes = bw.chroms()

    opts = Options(
        nfr_threshold = threshold
    )

    # I think this is required to avoid getting NAs from the bigwig file
    # TODO should just default to zero
    if not merge:
        opts.extend = 0
        opts.merge = 0

    ranges = BedFile(bed, chrom_sizes, opts)

    writer = BEDWriter(output)

    for r in ranges:
        occupancy = Occupancy.from_bigwig(bw, r, opts)
        nfrs = occupancy.threshold()
        writer.writerows(nfrs)


fastocc.add_command(occ)
fastocc.add_command(nfr)
