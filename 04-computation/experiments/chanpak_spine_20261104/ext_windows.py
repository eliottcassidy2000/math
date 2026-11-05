#!/usr/bin/env python3
"""ext_windows.py -- shape of the true one-vertex extension spectrum E(X): the smallest positive increment
(H(X+v) - h)/2 over cuts, and the longest window of consecutive odd values inside E(X), per tournament X of order n.
Also records every (X, v) with H(X + v) = h + 2 (increment 2)."""
import sys
import numpy as np
from ext_spectrum_true import tournaments, ext_spectrum, strong

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for n in range(3, NMAX + 1):
        inc2 = {}; stats = []; by_h = {}
        for adj in tournaments(n):
            H = ext_spectrum(n, adj); h = int(H[0]); E = sorted(set(int(x) for x in H)); M = E[-1]
            pos = [x for x in E if x > h]; minpos = (pos[0] - h) // 2 if pos else None
            if minpos == 1: inc2[h] = inc2.get(h, 0) + 1
            # longest window of consecutive odd values in E
            best = (0, None, None); start = E[0]; prev = E[0]
            for x in E[1:] + [None]:
                if x is not None and x == prev + 2: prev = x; continue
                length = (prev - start) // 2 + 1
                if length > best[0]: best = (length, start, prev)
                if x is not None: start = prev = x
            L, p, q = best
            st = strong(n, adj)
            stats.append((h, M, minpos, L, p, q, st))
            by_h.setdefault(h, []).append((L, p, q, M, st))
        hs = sorted(by_h)
        # for each h: the realiser with the longest window; report window position relative to h
        rep = []
        for h in hs:
            L, p, q, M, st = max(by_h[h], key=lambda r: r[0])
            rep.append((h, L, p, q, M))
        print("n=%d: tournaments %d; increment-2 cases by h: %s; min positive increment overall (h>=5): %s" % (
            n, len(stats), inc2, min((s[2] for s in stats if s[0] >= 5 and s[2] is not None), default=None)), flush=True)
        print("   per h, best window (h, L, p, q, M) for the 12 largest h and the first 8 h >= 23:")
        print("     ", [r for r in rep if r[0] >= 23][:8])
        print("     ", rep[-6:])
        # how many X have a window covering [2h, 3h]?
        cov = sum(1 for (h, M, mp, L, p, q, st) in stats if p is not None and p <= 2 * h + 1 and q >= 3 * h)
        covs = sum(1 for (h, M, mp, L, p, q, st) in stats if st and p is not None and p <= 2 * h + 1 and q >= 3 * h)
        print("   X whose longest window covers [2h+1, 3h]: %d of %d (strong: %d of %d); min over X of (q - p)/h: %.3f" % (
            cov, len(stats), covs, sum(1 for s in stats if s[6]), min(((q - p) / h) for (h, M, mp, L, p, q, st) in stats if p is not None)), flush=True)
