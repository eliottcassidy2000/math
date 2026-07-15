#!/usr/bin/env python3
"""locker_tournament_and_corridor_kps_S128c20.py -- kind-pasteur S128 cont.20.
(1) THE DEEP-WELL CORRIDOR LAW referee: m({1..12}; lam) = (2 H_12 / 13)(1 - 13 lam) on [1/14, 1/13];
    value at the covering-min 14/183 = 2 H_12 / 2379. (2) THE LOCKER (DIVISIBILITY) TOURNAMENT D_n:
    tile (x,y) bit = [y | x] (binary divisibility conditions as edges); first invariants n=5..9."""
import sys
from fractions import Fraction as F
from math import comb
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
# (1) corridor law
def good_measure(speeds, lam):
    pieces = []
    for w in speeds:
        r = lam / w
        for k in range(w):
            c = F(k, w); lo, hi = c - r, c + r
            if lo < 0: pieces.append((F(0), hi)); pieces.append((lo + 1, F(1)))
            elif hi > 1: pieces.append((lo, F(1))); pieces.append((F(0), hi - 1))
            else: pieces.append((lo, hi))
    pieces.sort(); tot = F(0); cur = F(0)
    for lo, hi in pieces:
        if lo > cur: tot += lo - cur
        cur = max(cur, hi)
    if cur < 1: tot += 1 - cur
    return tot
H12 = sum(F(1, j) for j in range(1, 13))
ok = True
for lam in [F(1,14), F(14,183), F(15,196), F(1,13)]:
    pred = (2*H12/13)*(1-13*lam)
    got = good_measure(list(range(1,13)), lam)
    ok &= (got == pred)
    print("corridor lam=%s: m=%s pred=%s %s" % (lam, got, pred, "OK" if got==pred else "FAIL"))
print("covering-min margin: m({1..12}; 14/183) = %s = 2H_12/2379: %s" % (good_measure(list(range(1,13)), F(14,183)), 2*H12/2379))
print("CORRIDOR LAW:", "EXACT" if ok else "FAIL")
# (2) locker tournament
def ham(B, n):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for S in range(1<<n):
        for v in range(n):
            c = dp[S][v]
            if not c or not (S>>v)&1: continue
            for u in range(n):
                if not (S>>u)&1 and B[v][u]:
                    dp[S|1<<u][u] += c
    return sum(dp[(1<<n)-1][v] for v in range(n))
def iso(B1, B2, n):
    d1=[sum(B1[u]) for u in range(n)]; d2=[sum(B2[u]) for u in range(n)]
    if sorted(d1)!=sorted(d2): return False
    from collections import defaultdict as dd
    buckets=dd(list)
    for v in range(n): buckets[d2[v]].append(v)
    order=sorted(range(n), key=lambda u:(d1[u],u)); assign=[-1]*n; used=[False]*n
    def bt(i):
        if i==n: return True
        u=order[i]
        for v in buckets[d1[u]]:
            if used[v]: continue
            g=True
            for j in range(i):
                w=order[j]
                if B1[u][w]!=B2[v][assign[w]] or B1[w][u]!=B2[assign[w]][v]: g=False; break
            if g:
                assign[u]=v; used[v]=True
                if bt(i+1): return True
                used[v]=False; assign[u]=-1
        return False
    return bt(0)
print()
print("THE LOCKER TOURNAMENT D_n: tile (x,y) bit=1 iff y | x  (base path fixed)")
for n in range(5, 10):
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    tidx = {t: i for i, t in enumerate(tiles)}
    m = len(tiles)
    gmap = [tidx[(n - y + 1, n - x + 1)] for (x, y) in tiles]
    tv = [1 if x % y == 0 else 0 for (x, y) in tiles]
    B = [[False]*n for _ in range(n)]
    for k in range(2, n+1): B[k-1][k-2] = True
    for (x, y), i in tidx.items():
        B[(x-1) if tv[i] else (y-1)][(y-1) if tv[i] else (x-1)] = True
    d = [sum(r) for r in B]
    c3 = comb(n,3) - sum(comb(x,2) for x in d)
    H = ham(B, n)
    BT = [[B[u][v] for u in range(n)] for v in range(n)]
    sc = iso(B, BT, n)
    kv = [1-b for b in tv]
    Bk = [[False]*n for _ in range(n)]
    for k in range(2, n+1): Bk[k-1][k-2] = True
    for (x, y), i in tidx.items():
        Bk[(x-1) if kv[i] else (y-1)][(y-1) if kv[i] else (x-1)] = True
    qf = iso(B, Bk, n)
    gs = all(tv[i]==tv[gmap[i]] for i in range(m))
    print("n=%d: scores=%s c3=%2d H=%4d  SC:%s  quasi-fixed:%s  gridsym:%s  (H mod 4 = %d)"
          % (n, sorted(d), c3, H, sc, qf, gs, H % 4))
print("DONE")
