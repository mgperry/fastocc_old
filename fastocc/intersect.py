from .range import Range
import numpy as np

def find_next(A: np.ndarray, y: int):
    step = 1000

    for i in range(0, len(A), step):
        if i + step > len(A):
            step = len(A) - i

        if A[i + step - 1] <= y:
            continue

        return i + np.searchsorted(A[i:i+step], y)
    

def find_range(A: np.ndarray, r: Range):
    a = find_next(A, r.start)

    # r is beyond all values in A
    if a is None:
        return None

    b_rel = find_next(A[a:], r.end)

    
    if b_rel is None:
        # r includes final value in A
        b = len(A)
    else:
        b = a + b_rel

    return (a, b)


def intersect_sorted(A: np.ndarray, rs: list[Range]):
    idx = 0

    for r in rs:
        ab = find_range(A[idx:], r)

        if ab is None:
            break

        a, b = ab

        yield r, (idx+a, idx+b)

        idx = idx + b
