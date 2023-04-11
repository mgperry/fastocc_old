import numpy as np

from .io import fetch_bw
from .range import Range
from .options import Options
from typing import Iterator


def contiguous_regions(x: np.ndarray) -> Iterator[tuple]:
    """
    Calculate continuous regions in a numpy array, returning
    a generator of (start, stop) tuples.

    e.g. [0, 0, 1, 1, 1, 0, 1, 1, 0] => (2,5); (6,8)

    NB If the interval starts or ends with 1, this is handled, as a result the indexing
    is basically BED-like, ie (0, 1) represents an interval containing just the first
    item. This should be very rare for our exact use case, it means that the peak caller
    has missed a nucleosome free region.

    e.g. [1, 0, 0, 1, 1] => (0, 1); (3, 5)
    """
    d = np.diff(x)
    idx = d.nonzero()[0] + 1

    if x[0]:
        idx = np.r_[0, idx]

    if x[-1]:
        idx = np.r_[idx, x.size]

    for start, stop in np.reshape(idx, (-1, 2)):
        yield start, stop


def threshold_occupancy(r: Range, bw, th: float, bin_size = Options.bin_size) -> Iterator[Range]:
    occ = fetch_bw(r, bw)

    coords = contiguous_regions(occ > th)

    for i, (start, stop) in enumerate(coords):
        name = f"{r.chr}:{r.start}-{r.end}" if r.name == "." else r.name

        yield Range(
            r.chr,
            r.start + bin_size*(start),
            r.start + bin_size*stop-1,
            f"{name}:NFR{i+1}",
            np.max(occ[start:stop]),
            "."
        )
