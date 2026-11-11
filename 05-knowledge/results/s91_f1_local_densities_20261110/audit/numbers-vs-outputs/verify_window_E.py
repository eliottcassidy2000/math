#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: independent recomputation (own enumeration kernel of verify_misc.py) of the measured
excess E of the distinguished class on the pre-registered window of odd N in [2^24, 2^27), for the primes quoted in
README section 1 and THM-5128 (G), against the committed window_split_prereg_27_24.out."""
import sys

import numpy as np
from verify_misc import build_q, window

ref = {}
for line in open(sys.argv[1]):
    f = line.split()
    if len(f) == 10 and f[0].isdigit() and f[1] in ("s=0", "s>=1"):
        ref[(int(f[0]), 0 if f[1] == "s=0" else 1)] = (float(f[2]), float(f[9]))
LO, K = 24, 27
A, B = 1 << LO, 1 << K
q = build_q(B // 2 + 2)
nodd = (B - A) // 2
r0 = np.zeros(nodd, np.int32)
r1 = np.zeros(nodd, np.int32)
t0, t1 = window(q, A, B, K, r0, r1)
del q
Nodd = np.arange(A + 1, B, 2, dtype=np.int64)
worst = 0.0
nprimes = 0
for p in sorted({k[0] for k in ref}):
    res = Nodd % p
    nclass = np.bincount(res, minlength=p)
    row = []
    for sheet, r in ((0, r0), (1, r1)):
        cs = (p - 1) if sheet == 0 else 0
        counts = np.bincount(res, weights=r, minlength=p)
        means = counts / nclass
        mo = np.delete(means, cs).mean()
        E = (means[cs] - mo) / mo
        var = float((r.astype(np.float64) ** 2).mean() - r.mean(dtype=np.float64) ** 2)
        se = (var / nclass[cs]) ** 0.5 / mo
        d = abs(E - ref[(p, sheet)][0])
        worst = max(worst, d, abs(se - ref[(p, sheet)][1]))
        row.append(f"E = {E:+.5f} (committed {ref[(p, sheet)][0]:+.5f}), se {se:.5f} (committed {ref[(p, sheet)][1]:.5f})")
    nprimes += 1
    if p in (3, 5, 7, 11, 13, 17, 19, 29, 31, 43, 257, 2731, 43691, 65537, 65521, 65539, 8191):
        print(f"p = {p:6d}: sheet 0 {row[0]}; sheet 1 {row[1]}")
print(f"{nprimes} primes, both sheets: largest |difference| from the committed E and se columns = {worst:.2e}")
