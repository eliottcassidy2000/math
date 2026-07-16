#!/usr/bin/env python3
"""oddcycle_functional_n67_kps_S128c31.py -- kind-pasteur S128 cont.31.
Verify C(n,2) E[dH|T] = a_n - b_n H - w5 c5 - w7 c7 at n=6,7; pin whether c7 enters."""
import sys
from math import comb
from fractions import Fraction as F
from itertools import permutations, combinations
sys.stdout.reconfigure(line_buffering=True)
def run(n):
    pairs = [(u, v) for u in range(n) for v in range(u+1, n)]
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]
    tidx = {t: i for i, t in enumerate(tiles)}
    m = len(tiles)
    m2 = comb(n, 2)
    def ham(B):
        dp = [[0]*n for _ in range(1<<n)]
        for v in range(n): dp[1<<v][v] = 1
        for S in range(1<<n):
            row = dp[S]
            for v in range(n):
                c = row[v]
                if not c or not (S>>v)&1: continue
                Bv = B[v]
                for u in range(n):
                    if not (S>>u)&1 and Bv[u]: dp[S|1<<u][u] += c
        full = (1<<n)-1
        return sum(dp[full][v] for v in range(n))
    def codd(B, L):
        tot = 0
        for S in combinations(range(n), L):
            u = S[0]
            for perm in permutations(S[1:]):
                prev = u; ok = True
                for w in perm:
                    if not B[prev][w]: ok = False; break
                    prev = w
                if ok and B[prev][u]: tot += 1
        return tot
    tuples = {}
    for t in range(1 << m):
        tv = [(t>>i)&1 for i in range(m)]
        B = [[False]*n for _ in range(n)]
        for k in range(2, n+1): B[k-1][k-2] = True
        for (x, y), i in tidx.items():
            if tv[i]: B[x-1][y-1] = True
            else: B[y-1][x-1] = True
        H0 = ham(B)
        c5 = codd(B, 5)
        c7 = codd(B, 7) if n >= 7 else 0
        s = 0
        for u in range(n):
            for v in range(n):
                if u != v and B[u][v]:
                    B[u][v], B[v][u] = False, True
                    s += ham(B) - H0
                    B[u][v], B[v][u] = True, False
        key = (H0, c5, c7)
        val = F(s, m2)
        if key in tuples:
            if tuples[key] != val:
                tuples[key] = 'SPLIT'   # (H,c5,c7) does NOT determine drift
        else:
            tuples[key] = val
        if n == 7 and t % 4096 == 0:
            print("  ...%d/32768, tuples %d" % (t, len(tuples)), flush=True)
    splits = [k for k, v in tuples.items() if v == 'SPLIT']
    print("n=%d: distinct (H,c5,c7) tuples %d ; fibers where (H,c5,c7) fails to determine drift: %d %s"
          % (n, len(tuples), len(splits), splits[:4] if splits else ""))
    clean = {k: v for k, v in tuples.items() if v != 'SPLIT'}
    # exact affine fit: m2*E = a - b H - w5 c5 - w7 c7 : solve from 4 tuples, verify all
    import itertools as it
    keys = list(clean)
    from fractions import Fraction
    def solve(sel):
        # unknowns a, b, w5, w7
        rows = []
        for (H, c5, c7) in sel:
            rows.append([1, -H, -c5, -c7, m2*clean[(H,c5,c7)]])
        # gaussian elim exact
        M = [r[:] for r in rows]
        cols = 4
        piv = []
        r0 = 0
        for c in range(cols):
            p = next((r for r in range(r0, len(M)) if M[r][c] != 0), None)
            if p is None: continue
            M[r0], M[p] = M[p], M[r0]
            M[r0] = [x / M[r0][c] for x in M[r0]]
            for r in range(len(M)):
                if r != r0 and M[r][c] != 0:
                    M[r] = [x - M[r][c]*y for x, y in zip(M[r], M[r0])]
            piv.append(c); r0 += 1
        if r0 < 4: return None
        sol = [None]*4
        for i, c in enumerate(piv): sol[c] = M[i][4]
        return sol
    found = None
    for sel in it.combinations(keys, 4):
        sol = solve(sel)
        if sol and all(x is not None for x in sol):
            ok = all(m2*clean[k] == sol[0] - sol[1]*k[0] - sol[2]*k[1] - sol[3]*k[2] for k in keys)
            if ok:
                found = sol
                break
    if found:
        print("n=%d LAW EXACT: C(n,2) E[dH] = %s - %s H - %s c5 - %s c7   (c7 %s)"
              % (n, found[0], found[1], found[2], found[3], "ENTERS" if found[3] != 0 else "ABSENT"))
    else:
        print("n=%d: NO exact affine law in (1, H, c5, c7) -- higher terms needed" % n)
run(6)
run(7)
print("DONE")
