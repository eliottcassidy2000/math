#!/usr/bin/env python3
"""Vertex-wise balance in the strong transitive-plus-one family: does every vertex (z and every u_i) have a partner x with
1/3 <= P(v before x) <= 2/3?  Exact, all strong words up to m = MMAX (uses trans_plus_one.stats)."""
import sys
from fractions import Fraction as F
from itertools import product
src = open('trans_plus_one.py').read().split("M = int(sys.argv[1])")[0]
exec(src)
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 13
for m in range(2, MMAX + 1):
    fails = 0; worst = None
    for mid in product((0, 1), repeat=m - 2):
        w = (1,) + mid + (0,)
        H, inA, cnt = stats(w)
        third = lambda c: 3 * c >= H and 3 * c <= 2 * H
        okz = any(third(inA[i]) for i in range(m))
        oku = []
        for i in range(m):
            ok = third(inA[i])                                  # partner z
            for j in range(m):
                if j == i: continue
                a, b = min(i, j), max(i, j)
                c = cnt[a][b]                                  # #paths with u_b before u_a
                if third(c): ok = True
            oku.append(ok)
        if not (okz and all(oku)):
            fails += 1
            if worst is None: worst = (''.join(map(str, w)), okz, oku)
    print(f"m={m} (n={m+1}): {2**(m-2)} words; words where some vertex has no balanced partner: {fails}" + (f"  e.g. {worst}" if worst else ""))
