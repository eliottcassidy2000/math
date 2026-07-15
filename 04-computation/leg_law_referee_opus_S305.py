#!/usr/bin/env python3
"""THM-790 leg-law referee (opus-S305): on ALL tilings n=4..7 verify
(1) dx = 8(e1 - en); (2) the flip score action d -> -d with leg defect
(-2 at v=1, +2 at v=n); (3) the layer GF 2^{m-2n+5}(1+z)^{n-3}(1+1/z)^{n-3}(z+1/z);
(4) the blue GF 2^{(m+f)/2-(n-2)}(z+1/z)^{n-2}. See the session .out for the run."""
from collections import defaultdict
from math import comb
for n in (4, 5, 6, 7):
    V = list(range(1, n+1))
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]
    m = len(tiles); ti = {t: i for i, t in enumerate(tiles)}
    refl = [ti[(n-y+1, n-x+1)] for (x, y) in tiles]
    f = sum(1 for i in range(m) if refl[i] == i)
    ok = [True]*4
    hist = defaultdict(int); bhist = defaultdict(int)
    for b in range(1 << m):
        e = {v: 0 for v in V}; eb = {v: 0 for v in V}
        for i, (x, y) in enumerate(tiles):
            if b >> i & 1: e[x] += 1; eb[y] += 1
            else: e[y] += 1; eb[x] += 1
        d = {v: 2*((1 if v >= 2 else 0) + e[v]) - (n-1) for v in V}
        db = {v: 2*((1 if v >= 2 else 0) + eb[v]) - (n-1) for v in V}
        dx = sum(db[v]**2 for v in V) - sum(d[v]**2 for v in V)
        e1 = sum(1 for i,(x,y) in enumerate(tiles) if y == 1 and not (b >> i & 1))
        en = sum(1 for i,(x,y) in enumerate(tiles) if x == n and (b >> i & 1))
        if dx != 8*(e1-en): ok[0] = False
        for v in V:
            if db[v] != -d[v] + (-2 if v == 1 else (2 if v == n else 0)): ok[1] = False
        hist[e1-en] += 1
        if all((b >> i & 1) == (b >> refl[i] & 1) for i in range(m)): bhist[e1-en] += 1
    pred = defaultdict(int)
    for A in range(n-2):
        for B in range(n-2):
            for c in (0, 1):
                pred[A - B + (1-2*c)] += comb(n-3, A)*comb(n-3, B) * (1 << (m-2*n+5))
    ok[2] = dict(hist) == {k: v for k, v in pred.items() if v}
    bpred = {2*k-(n-2): comb(n-2, k) * (1 << ((m+f)//2 - (n-2))) for k in range(n-1)}
    ok[3] = dict(bhist) == {k: v for k, v in bpred.items() if v}
    print(f"n={n}: leg law {ok[0]}, score action {ok[1]}, layer GF {ok[2]}, blue GF {ok[3]}")
