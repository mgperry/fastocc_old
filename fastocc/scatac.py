import numpy as np
import h5py
import pandas as pd
from pybedtools import BedTool
from collections.abc import Iterator

from .range import Range
from .fragments import Fragment
from . import io
from .options import Options
from .occupancy import occupancy


def load_fragments(chr: h5py.Group):
    return {k: v[:] for k, v in chr.items()}


def fragments_from_chr(chr, cell):
    
    is_correct_id = np.argwhere(chr["RGValues"][0,:] == cell.encode("ascii"))

    if len(is_correct_id) == 1:
        index = is_correct_id.flatten()[0]
    elif len(is_correct_id) == 0:
        # not all chromosomes contain fragments from all cells
        return pd.DataFrame({
            "start": np.array([], dtype=int),
            "length": np.array([], dtype=int),
            "cell": [],
        })
    else:
        raise Exception("Cell ID should be unique.")

    start = chr["RGLengths"][0,:index].sum()
    end = start + chr["RGLengths"][0,index]

    return pd.DataFrame({
        "start": chr["Ranges"][0, start:end],
        "length": chr["Ranges"][1, start:end],
        "cell": cell,
    })


class ArrowFragments:
    """
    Read all fragments from a subset of cells into memory, reading the file once. 
    """
    def __init__(self, arrow, cells):
        frags = {}

        # load all fragments for chosen cells into memory
        for chr in arrow["Fragments"]:
            data = load_fragments(arrow["Fragments"][chr])

            frags_by_cell = [fragments_from_chr(data, cell) for cell in cells]

            frags[chr] = self.sort_by_midpoint(pd.concat(frags_by_cell))

        self.fragments = frags
        self.cells = cells

    @staticmethod
    def sort_by_midpoint(df):
        df["midpoint"] = df["start"] + df["length"] // 2
        return df.sort_values("midpoint")
    

    def iter_ranges(self, ranges: list[Range]) -> Iterator[list[Fragment]]:
        for r in ranges:
            yield self.fetch(r)


    def fetch(self, r: Range) -> list[Fragment]:
        idx = np.searchsorted(self.fragments[r.chr]["midpoint"], [r.start, r.end])
        df = self.fragments[r.chr].iloc[idx[0]:idx[1]]
        return [Fragment(f.start, f.length) for f in df.itertuples()]


def sc_occ(arrow: str, bed: str, output: str, chrom_sizes: dict[str, int], pseudocount: float):
    if not pseudocount > 0:
        pseudocount = 0.001 # avoid potential div by zero errors

    fragments = ArrowFragments(arrow)

    intervals = io.fix_bed(BedTool(bed), chrom_sizes, Options.extend, Options.merge)

    ranges = (Range.from_interval(i).align(Options.bin_size) for i in intervals)

    writer = io.BedGraphWriter(output)
    
    for r in ranges:
        occ = occupancy(r, fragments, pseudocount)
        writer.write(r, occ)
