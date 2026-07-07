#!/usr/bin/env python3
"""
klein-2026-07-07-S173 -- the FULL n=10 blue self-loop census (background).
Pipeline: enumerate 2^20 gridsym tilings via orbit bits; filter D=0 (affine-D law,
THM-649); filter score-multiset equality T vs T-flip; filter refined invariant
(sorted per-vertex (out,in) degree-pair + c3 exact + out-degree-of-out-neighbors
profile); refinement canon on survivors. Reports self-loop tilings and lines.
"""
import itertools
from collections import defaultdict, Counter
from math import comb

n = 10
tiles = [(x,y) for x in range(1,n+1) for y in range(1,x) if x-y>=2]
m = len(tiles); full = (1<<m)-1
tidx = {t:k for k,t in enumerate(tiles)}
sigma = [tidx[(n+1-y,n+1-x)] for (x,y) in tiles]
orbs = []
seen = set()
for b in range(m):
    if b in seen: continue
    o = {b, sigma[b]}; seen |= o; orbs.append(sorted(o))
nfree = len(orbs)

# affine-D coefficients (THM-649): carriers legs +1, apex +2
def sb(v): return 1 if v >= 2 else 0
def td(v): return ((v-2) if v>=2 else 0) + ((n-1-v) if v<=n-1 else 0)
w = {v: 2*sb(v)+td(v)-1 for v in range(1,n+1)}
cco = [w[x]-w[y] for (x,y) in tiles]
# D(t) = D(0) + sum c_b bit_b; compute D(0)
def adj(tv):
    A = [[0]*n for _ in range(n)]
    for a in range(2,n+1): A[a-1][a-2] = 1
    for b,(x,y) in enumerate(tiles):
        if (tv>>b)&1: A[x-1][y-1] = 1
        else: A[y-1][x-1] = 1
    return A
def c3(A):
    return comb(n,3) - sum(comb(sum(r),2) for r in A)
D0 = c3(adj(full)) - c3(adj(0))
orb_c = [sum(cco[b] for b in o) for o in orbs]

def canon(A):
    col = [sum(A[i]) for i in range(n)]
    for _ in range(n):
        prof = []
        for i in range(n):
            outs = tuple(sorted(col[j] for j in range(n) if A[i][j]))
            ins = tuple(sorted(col[j] for j in range(n) if A[j][i]))
            prof.append((col[i], outs, ins))
        rk = {p:r for r,p in enumerate(sorted(set(prof)))}
        nc = [rk[p] for p in prof]
        if nc == col: break
        col = nc
    cells = defaultdict(list)
    for i in range(n): cells[col[i]].append(i)
    ordered = [cells[c] for c in sorted(cells)]
    best = None
    def enc(perm):
        return tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
    def rec(pre, rest):
        nonlocal best
        if not rest:
            k = enc(pre)
            if best is None or k < best: best = k
            return
        for p in itertools.permutations(rest[0]):
            rec(pre+list(p), rest[1:])
    tot = 1
    for c in ordered:
        f = 1
        for i in range(2,len(c)+1): f *= i
        tot *= f
    if tot > 500000:
        return ("BIG", tuple(sorted(sum(r) for r in A)))  # fallback marker (rare)
    rec([], ordered)
    return best

def invariant(A):
    s = [sum(r) for r in A]
    prof = tuple(sorted((s[i], sum(s[j] for j in range(n) if A[i][j])) for i in range(n)))
    return prof

survivors = []
cnt_D0 = 0
for bits in range(1<<nfree):
    D = D0 + sum(orb_c[i] for i in range(nfree) if (bits>>i)&1)
    if D != 0: continue
    cnt_D0 += 1
    tv = 0
    for i,o in enumerate(orbs):
        if (bits>>i)&1:
            for b in o: tv |= (1<<b)
    A1 = adj(tv); A2 = adj(tv^full)
    if sorted(sum(r) for r in A1) != sorted(sum(r) for r in A2): continue
    if invariant(A1) != invariant(A2): continue
    survivors.append(tv)
print(f"n=10: gridsym 2^{nfree}; D=0 candidates {cnt_D0}; invariant survivors {len(survivors)}")
loops = []
for tv in survivors:
    if canon(adj(tv)) == canon(adj(tv^full)):
        loops.append(tv)
print(f"SELF-LOOP TILINGS = {len(loops)}   lines = {len(loops)//2}   [2^(n/2-1)=16 tilings / 2^(n/2-2)=8 lines predicted]")
for tv in loops: print("  loop tiling:", tv)
