#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (A). Signed flatness defects at 36 = 4 * 9, cell by cell:
#   delta(e) = 2 #{N = 1 mod 4, N = e mod 9} / #{N = e mod 9} - 1,   in units of 1e-4, and in units of
#   sd(e) = sqrt(kappa_disp / #{N = e mod 9}) with kappa_disp = 3 (index of dispersion of r, THM-5128 Remark 3).
# Is the defect systematic (same sign and size across X) or sampling noise?
import numpy as np, json, math
print("# rflat_cells_36.py (auditor hypotheses)")
for m, j in ((9, 2), (3, 2), (5, 2)):
    for sh in (0, 1):
        print("modulus %d, sheet %d: delta(e) * 1e4  [z-score with dispersion 3]" % (4 * m, sh))
        for E in (22, 24, 26, 28):
            hist = np.array(json.load(open("even_moduli_hist_%d.json" % E))["hist"], dtype=np.int64).sum(axis=1)
            N = np.arange(2520)
            cells = []
            for e in range(m):
                tot = hist[sh][N % m == e].sum()
                one = hist[sh][(N % m == e) & (N % 4 == 1)].sum()
                d = 2 * one / tot - 1
                z = d / math.sqrt(3.0 / tot)
                cells.append("%+6.2f [%+4.1f]" % (d * 1e4, z))
            print("   2^%d: " % E + "  ".join(cells))
