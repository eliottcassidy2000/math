#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (A). R_flat with more digits and in units of the sampling scale
# sqrt(2^(j-1) m / #T), for 12, 20, 36, 72 at X = 2^20 .. 2^28 (from the saved histograms).
import numpy as np, json, math
print("# rflat_detail.py (auditor hypotheses): R_flat = max over e mod m, odd a mod 2^j of |2^(j-1) #{N=a (2^j), N=e (m)}/#{N=e (m)} - 1|")
for (j, m) in ((2, 3), (2, 5), (2, 9), (3, 9), (3, 3), (3, 5)):
    Mod = (1 << j) * m
    row = []
    for E in (20, 22, 24, 26, 28):
        hist = np.array(json.load(open("even_moduli_hist_%d.json" % E))["hist"], dtype=np.int64).sum(axis=1)
        cells = []
        for sh in (0, 1):
            T = hist[sh].sum()
            cnt = np.zeros(Mod)
            cm = np.zeros(m)
            for c in range(2520):
                cnt[c % Mod] += hist[sh, c]
                cm[c % m] += hist[sh, c]
            r = max(abs(2 ** (j - 1) * cnt[c] / cm[c % m] - 1) for c in range(1, Mod, 2))
            cells.append("%.5f (%.1f)" % (r, r / math.sqrt(2 ** (j - 1) * m / T)))
        row.append("2^%d: %s / %s" % (E, cells[0], cells[1]))
    print("modulus %3d:  " % Mod + ";  ".join(row))
print("# in brackets: R_flat divided by sqrt(2^(j-1) m / #T)")
