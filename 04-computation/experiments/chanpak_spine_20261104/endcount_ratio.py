#!/usr/bin/env python3
"""endcount_ratio.py -- HYP-9469 (4): max over tournaments X of order n and non-sink vertices a of e_a(X)/H(X);
also the increment-10/12 levels' h-sets (HYP-9469 (2)) at orders of X <= 8."""
import sys
from collections import defaultdict
sys.path.insert(0, '.')
from ext_spectrum_true import tournaments, end_start_tables, ext_spectrum
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for n in range(3, NMAX + 1):
    worst = (0, None); lev = defaultdict(set)
    for adj in tournaments(n):
        e, s = end_start_tables(n, adj); full = (1 << n) - 1
        h = int(e[full].sum())
        if h >= 5:
            for a in range(n):
                if adj[a] != 0:   # not a sink
                    r = int(e[full][a]) / h
                    if r > worst[0]: worst = (r, (h, int(e[full][a])))
        Hs = ext_spectrum(n, adj)
        for A in range(1 << n):
            k2 = int(Hs[A]) - h
            if 0 < k2 <= 12: lev[k2].add(h)
    print("n=%d: max e_a/h over non-sink a (h >= 5): %.4f at (h, e_a) = %s; h-sets by increment: %s" % (
        n, worst[0], worst[1], {k2: sorted(lev[k2]) for k2 in sorted(lev)}), flush=True)
