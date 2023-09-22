import numpy as np
import h5py
import pandas as pd

from .range import Range
from .fragments import Fragments
from .intersect import intersect_sorted


def load_cells(grp, ids=None):
    """
    NB cells are not order wrt `ids`, and can be missing.
    """

    # we cache this because indexing doesn't work with default archR compression
    cells = grp["RGValues"][0, :]
    n_frags = grp["RGLengths"][0, :]

    ranges = grp["Ranges"][:]
    
    start = 0

    for cell, N in zip(cells, n_frags):
        end = start + N

        id = cell.decode('UTF-8')

        if ids is None or id in ids:
            yield {
                "id": id,
                "start": ranges[0, start:end],
                "size": ranges[1, start:end],
            }
    
        start = end


def fragments(cell, ranges: list[Range]):
    df = pd.DataFrame({
        "start": cell["start"],
        "size": cell["size"]
    })

    df["midpoint"] = df["start"] + df["size"] // 2
    df = df.sort_values("midpoint")

    for r, (a, b) in intersect_sorted(df["midpoint"].values, ranges):
        fs = Fragments(
            r, 
            df["midpoint"].values[a:b],
            df["size"].values[a:b]
        )

        fs.cell = cell.get("id", None)

        yield fs


def merge_cells(cells: list[str], id=None):
    dfs = []

    for cell in cells:
        dfs.append(pd.DataFrame({
            "start": cell["start"],
            "size": cell["size"]
        }))

    cluster_df = pd.concat(dfs)

    return {
        "id": id,
        "start": cluster_df["start"].values,
        "size": cluster_df["size"].values,
    }
