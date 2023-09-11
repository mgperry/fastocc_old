import h5py
import numpy as np

from itertools import groupby
from dataclasses import asdict

from .range import Range

def count_bins(r: Range, size):
    return (r.end - r.start + 1) // size


def write_datasets(grp: h5py.Group, datasets: dict[str, np.ndarray], **kwargs):
    for k, d in datasets.items():
        grp.create_dataset(k, data=d, **kwargs)


class HDF5Writer:
    compression = "gzip"
    chunks = 10000
    scaleoffset = 2 # decimal places for floats

    def __init__(self, file: h5py.File, schema: dict, opts=None, genome_info=None):
        self.file = file
        self.schema = schema

        self.file.create_group("chromosomes")

        meta = self.file.create_group("metadata")
        
        meta.create_group("schema")
        meta["schema"].attrs.update(schema)

        if opts:
            meta.create_group("opts")
            meta["opts"].attrs.update(asdict(opts))

        if genome_info:
            meta.create_group("genome_info")
            names, sizes = zip(*genome_info.items())
            write_datasets(meta["genome_info"], {"names": names, "sizes": sizes})

    def write(self, range_data):
        data_by_chr = groupby(range_data, key=lambda fs: fs.range.chr)

        for chr, windows in data_by_chr:
            grp = self.file["chromosomes"].create_group(chr)

            windows = list(windows) # don't exhaust iterator, mem in tiny anyway

            ranges, data = self.encode(windows)

            # fix for empty // small chrs:
            total = sum(ranges["bins"])

            chunks = min(total, self.chunks)

            write_datasets(grp, ranges)

            # write float arrays with scaleoffset
            for col, type in zip(self.schema["cols"], self.schema["types"]):
                if type[:5] == "float" and self.scaleoffset:
                    grp.create_dataset(col, data=data[col], chunks=chunks, scaleoffset=self.scaleoffset)
                    data.pop(col)

            write_datasets(grp, data, chunks=chunks, compression="gzip")

            print(f"{chr} written.")

    def encode(self, windows) -> (dict, dict):
        """
        Encode a list of binned range data (e.g. BinnedFragments, Occupancy) into
        numpy arrays representing the ranges and the data.
        """
        rs = (w.range for w in windows)
        start, end, bins = zip(
            *[(r.start, r.end, count_bins(r, self.schema["bin_size"])) for r in rs]
        )

        total = sum(bins)

        ranges = {
            "start": start,
            "end": end,
            "bins": bins
        }

        data = {}
        
        for col, typ in zip(self.schema["cols"], self.schema["types"]):
            data[col] = np.zeros(total, dtype=np.dtype(typ))

        start = 0

        for i, w in enumerate(windows):
            end = start + bins[i]
            for col in self.schema["cols"]:
                data[col][start:end] = getattr(w, col)
            start = end

        return ranges, data

    @classmethod
    def open(cls, f: str, *args, **kwargs):
        return cls(h5py.File(f, 'w'), *args, **kwargs)


class HDF5Loader:
    data_group = "chromosomes"

    def __init__(self, f: h5py.File, cls=None):
        self.file = f
        self.cls = cls

        meta = self.file["metadata"]

        self.schema = dict(meta["schema"].attrs)

        if "genome_info" in meta:
            self.chrom_info = {
                str(k, "utf-8"): v for k, v in
                zip(
                    meta["genome_info"]["names"],
                    meta["genome_info"]["sizes"]
                )
            }

        if "opts" in meta:
            self.opts = dict(meta["opts"].attrs)
        else:
            self.opts = None

    @classmethod
    def open(cls, f: str, *args, **kwargs):
        return cls(h5py.File(f, 'r'), *args, **kwargs)

    def load(self):
        for chr in self.file[self.data_group]:
            yield from self.decode(chr)

    def decode(self, chr):
        grp = self.file[self.data_group][chr]

        ranges = [
            Range(chr, s, e) for s, e in zip(grp["start"], grp["end"])
        ]

        # load info to avoid chunk nonsense
        datasets = {col: grp[col][:] for col in self.schema["cols"]}

        start = 0

        for i, n in enumerate(grp["bins"][:]):
            end = start + n

            data = {}

            for col in self.schema["cols"]:
                data[col] = datasets[col][start:end]

            if self.cls:
                args = {
                    "range": ranges[i],
                    **data,
                    "opts": self.opts
                }

                yield self.cls(**args)
            else:
                yield {
                    "range": ranges[i],
                    "data": data
                }

            start = end
