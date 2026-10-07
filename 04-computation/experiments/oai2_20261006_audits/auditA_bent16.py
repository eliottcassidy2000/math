"""Audit A: bent switchings of order 16 (tower vs Sylvester)."""
import numpy as np, itertools
def tower(k):
    H = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(k - 1): H = np.block([[H, H], [-H.T, H.T]])
    return H
def syl(k):
    H = np.array([[1]], dtype=np.int64)
    for _ in range(k): H = np.block([[H, H], [H, -H]])
    return H
D = np.array([[1] + list(t) for t in itertools.product((1, -1), repeat=15)], dtype=np.int64)
for name, H in (("tower", tower(4)), ("sylvester", syl(4))):
    P = D @ H.T
    print(f"{name}: #d (d0=1) with |Hd| = 4 entrywise: {int(np.all(np.abs(P) == 4, axis=1).sum())}; max_d min_i |(Hd)_i| = {int(np.abs(P).min(axis=1).max())}")
