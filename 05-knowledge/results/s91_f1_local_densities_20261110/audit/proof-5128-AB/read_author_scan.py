#!/usr/bin/env python3
"""Audit proof-5128-AB: read the author's scan_003_101.npz (data only, no pickle) and print the per-m maxima for
the primes up to 101, to compare with the band stated in THM-5128 (B)."""
import sys

import numpy as np

path = sys.argv[1]
d = np.load(path, allow_pickle=False)
print("arrays:", {k: d[k].shape for k in d.files})
m = d["m"]
rate = d["rate"]
est = d["est"]
primes = [p for p in range(3, 102) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
for p in primes:
    sel = m == p
    r = rate[sel]
    e = est[sel]
    rmax = np.nanmax(r) if np.any(~np.isnan(r)) else float("nan")
    print(f"m = {p:3d}: dense max rate {rmax:.6f}   (max estimate {np.max(e):.6f}, rows {int(sel.sum())})")
