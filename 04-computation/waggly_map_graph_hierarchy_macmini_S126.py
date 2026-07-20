#!/usr/bin/env python3
"""
THM-1390: the waggly filtration as a map-graph hierarchy (mac-mini-2026-07-20-S126)
===================================================================================
Map graphs generalize planar duals by letting faces meet at a POINT, not only an
edge -- and gain unbounded cliques. The waggly filtration does the same:
  d=1 (wiggly, flip one tile)  ~ the planar dual (edge contacts)
  d>=2 waggly layers           ~ point contacts
  G^(<=k)                      ~ a k-map graph truncation
RESULT: G^(1) is sparse and sparsifying (density .833/.455/.188 at n=4,5,6) while
G^(<=d) saturates to the COMPLETE graph -- the clique explosion.
NEW INVARIANT: d_sat(n) = least d with G^(<=d) complete = the S_n-quotient's
tile-flip diameter.  Computed exactly: 2, 3, 4, 7 for n = 4,5,6,7.
CAUTION: d_sat = n-2 fits n<=6 and is FALSE at n=7 (7, not 5) -- a small-n mirage.
"""
import numpy as np, collections
from itertools import combinations, permutations

def tiles(n): return [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]

def tour(n, mask, T):
    A = [[0]*n for _ in range(n)]
    fl = {T[i] for i in range(len(T)) if mask >> i & 1}
    for x in range(1, n+1):
        for y in range(1, x):
            if (x-y >= 2) and ((x, y) in fl): A[y-1][x-1] = 1
            else: A[x-1][y-1] = 1
    return A

def canon(n, A):
    """canonical form: refine by score, brute force only inside score-equal cells"""
    sc = [sum(A[i]) for i in range(n)]
    order = sorted(range(n), key=lambda i: sc[i]); cells = []; i = 0
    while i < len(order):
        j = i
        while j < len(order) and sc[order[j]] == sc[order[i]]: j += 1
        cells.append(order[i:j]); i = j
    def pg(c):
        if not c: yield []; return
        for p in permutations(c[0]):
            for r in pg(c[1:]): yield list(p)+r
    best = None
    for p in pg(cells):
        b = 0
        for a in range(n):
            for bb in range(n):
                if a != bb and A[p[a]][p[bb]]: b |= 1 << (a*n+bb)
        if best is None or b < best: best = b
    return best

def classes(n):
    T = tiles(n); m = len(T); N = 1 << m
    cid = {}; cls = np.empty(N, dtype=np.int32)
    for mask in range(N):
        c = canon(n, tour(n, mask, T))
        if c not in cid: cid[c] = len(cid)
        cls[mask] = cid[c]
    return m, N, cls, len(cid)

print("=== layer densities (small n, full layer enumeration) ===")
for n in (4, 5, 6):
    m, N, cls, V = classes(n)
    layer = collections.defaultdict(set)
    for mask in range(N):
        a = int(cls[mask])
        for d in range(1, m+1):
            for S in combinations(range(m), d):
                mm = mask
                for i in S: mm ^= 1 << i
                b = int(cls[mm])
                if a != b: layer[d].add((min(a, b), max(a, b)))
    tot = V*(V-1)//2; cum = set(); out = []
    for d in range(1, m+1):
        cum |= layer[d]; out.append(f"d<={d}:{len(cum)/tot:.3f}")
        if len(cum) == tot: break
    print(f"  n={n}: classes={V}, |E(d=1)|={len(layer[1])}, "
          f"density(d=1)={len(layer[1])/tot:.3f}  ->  " + "  ".join(out))

print()
print("=== saturation depth = S_n-quotient tile-flip diameter (exact) ===")
for n in (4, 5, 6, 7):
    m, N, cls, V = classes(n)
    masks = np.arange(N, dtype=np.int32)
    pc = np.zeros(N, dtype=np.int8)
    for i in range(m): pc += ((masks >> i) & 1).astype(np.int8)
    for D in range(1, m+1):
        pats = np.array([p for p in range(N) if pc[p] <= D], dtype=np.int32)
        seen = np.zeros((V, V), dtype=bool)
        step = max(1, 2**18 // max(1, len(pats)))
        for s in range(0, N, step):
            blk = masks[s:s+step]
            x = (blk[:, None] ^ pats[None, :]).ravel()
            seen[np.repeat(cls[blk], len(pats)), cls[x]] = True
        seen |= seen.T
        if not (~seen).any():
            print(f"  n={n}: m={m}, classes={V}, d_sat={D}   (n-2 = {n-2}, "
                  f"{'MATCH' if D == n-2 else 'MISMATCH -- the linear pattern is a mirage'})")
            break
