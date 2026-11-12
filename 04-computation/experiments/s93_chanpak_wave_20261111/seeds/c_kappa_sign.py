#!/usr/bin/env python3
"""S93 lane seeds: for all pairs {e, f} of all classes X of order <= NMAX: joint sign of kappa_X({e, f}) and of
2 H(X) - W_ef; the ratio kappa / (Pi_2/H * (2H - W))."""
import sys
from collections import Counter
import numpy as np
from seedlib import dp_tables, out_from_colcode
from c_neutral import gadget_tables

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
for n in range(2, nmax + 1):
    reps = np.atleast_1d(np.loadtxt(f"cls_{n}.txt", dtype=np.int64))
    cnt = Counter()
    rat = []
    for code in reps.tolist():
        out = np.zeros(n, np.int64)
        out_from_colcode(n, code, out)
        fe = np.zeros((1 << n, n), np.int64); gs = np.zeros((1 << n, n), np.int64); hS = np.zeros(1 << n, np.int64)
        UR = np.zeros((1 << n, n), np.int64); K = np.zeros((n, n), np.int64); Wm = np.zeros((n, n), np.int64)
        dp_tables(n, out, fe, gs, hS)
        pi2 = int(gadget_tables(n, fe, gs, hS, UR, K, Wm))
        H = int(hS[(1 << n) - 1])
        for e in range(n):
            for f in range(e + 1, n):
                kap = pi2 - int(K[e, f])
                d = 2 * H - int(Wm[e, f])
                cnt[((kap > 0) - (kap < 0), (d > 0) - (d < 0))] += 1
                if d != 0:
                    rat.append(kap * H / (pi2 * d))
    print(f"order {n}: (sign kappa, sign(2H - W)) counts: {dict(sorted(cnt.items()))}; kappa H / (Pi_2 (2H - W)) in [{min(rat):.4f}, {max(rat):.4f}]" if rat else f"order {n}: {dict(cnt)}")
