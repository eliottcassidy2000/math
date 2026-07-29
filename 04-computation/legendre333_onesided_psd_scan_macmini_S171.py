#!/usr/bin/env python3
"""
psd_scan.py — one-sided obstruction scan for Legendre pairs of length 333:
if A is H-invariant (|H| >= 7, <= 25 orbits) and (A, B) is a Legendre pair
with B UNSTRUCTURED, then PSD_A(k) <= 2L+2 = 668 for all k != 0 is NECESSARY
(PSD_B = 668 - PSD_A >= 0).  Count PSD-admissible invariant candidates.
Vectorized: chunked batch FFT over all row-sum-±1 invariant sign patterns.
"""
import numpy as np
import sys
sys.argv = ["x", "333"]
import legendre333_mim as LM

L = 333
subs = LM.all_subgroups()
big = [H for H in subs if len(H) >= 7]
total_admissible = 0
rows = []
for H in big:
    orbs = LM.orbits_of(H)
    k = len(orbs)
    if k > 25:
        continue
    sizes = np.array([len(o) for o in orbs], dtype=np.int64)
    M = np.zeros((k, L))
    for i, o in enumerate(orbs):
        M[i, o] = 1.0
    nbits = k
    total = 1 << nbits
    admissible = 0
    classes = 0
    minmax = np.inf
    CH = 1 << 14
    for start in range(0, total, CH):
        cnt = min(CH, total - start)
        bits = (np.arange(start, start + cnt)[:, None] >> np.arange(nbits)[None, :]) & 1
        signs = 2 * bits - 1                       # cnt x k
        rowsum = signs @ sizes
        keep = np.abs(rowsum) == 1
        if not keep.any():
            continue
        S = signs[keep] @ M                        # sequences, cnt' x 333
        classes += S.shape[0]
        psd = np.abs(np.fft.rfft(S, axis=1)) ** 2  # cnt' x 167
        mx = psd[:, 1:].max(axis=1)
        minmax = min(minmax, float(mx.min()) if mx.size else np.inf)
        admissible += int((mx <= 668 + 1e-6).sum())
    rows.append((len(H), k, classes, admissible, minmax))
    print(f"|H|={len(H):3d} orbits={k:2d} rowsum-ok={classes:7d} "
          f"PSD-admissible={admissible:6d} min-maxPSD={minmax if minmax != np.inf else '—'}",
          flush=True)
    total_admissible += admissible
print(f"TOTAL one-sided PSD-admissible invariant candidates: {total_admissible}", flush=True)
