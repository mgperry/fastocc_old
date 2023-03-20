from fastocc.fastocc import cli
run = cli.callback

from fastocc.range import Range
from fastocc.occupancy import occupancy
from fastocc.pysam_util import fetch_chrom_sizes, write_chrom_sizes
