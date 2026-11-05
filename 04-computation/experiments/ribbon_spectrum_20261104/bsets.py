#!/usr/bin/env python3
"""bsets.py -- exact root-response sets of ribbon words rooted at their last vertex.
For a tournament T with last vertex t: h = H(T), e = e_t(T), R = R_t(T) = sum_X h(T[X]) s_t(T - X).
In the reduced span-2 coordinates these are read off the post-sink state A_0 y: H, R - e = (A_0 y)_6, 2R = (A_0 y)_4.
(1) Ribbon B-sets: for every odd h <= HMAX, the set B(h) = {R : (h, e, R) occurs, R <= CMAX*h} over all span-2 words
    of length <= LMAX (BFS on reduced states pruned by H <= HMAX), with the density of same-parity-of-e lattice points
    hit in [h, CMAX*h]; saturation check between LMAX-4 and LMAX.
(2) All-tournament B-sets from gentourng up to order 8 (brute-force subset DP) for comparison.
Usage: bsets.py HMAX LMAX CMAX"""
import sys, subprocess
from collections import defaultdict
import numpy as np
HMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 99
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20
CMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 4
def load_intrep(path):
    with open(path) as f:
        r, nl = map(int, f.readline().split())
        y0 = list(map(int, f.readline().split())); hc = list(map(int, f.readline().split()))
        A = {}
        for _ in range(nl):
            mask = int(f.readline().strip()); A[mask] = [list(map(int, f.readline().split())) for _ in range(r)]
    return r, y0, hc, A
r, y0, hc, A = load_intrep('span2_intrep.txt')
Am = {m: np.array(M, dtype=np.int64) for m, M in A.items()}
masks = sorted(Am); hidx = hc.index(1); A0 = Am[0]
def triple(y):
    z = A0 @ np.array(y, dtype=np.int64)
    H = int(z[hidx]); d = int(z[5]); q = int(z[3]); R = q // 2; return H, R - d, R
B = defaultdict(set)      # h -> set of (e, R)
Bby = {}
cur = {tuple(y0)}
for n in range(1, LMAX + 1):
    nxt = set()
    for y in cur:
        v = np.array(y, dtype=np.int64)
        for m in masks:
            if n >= 2 and (m >> (n - 1)): continue
            if n == 1 and m: continue
            z = tuple((Am[m] @ v).tolist())
            if z[hidx] > HMAX: continue
            nxt.add(z)
    cur = nxt
    for y in cur:
        H, e, R = triple(y)
        if R <= CMAX * H: B[H].add((e, R))
    if n in (LMAX - 4, LMAX):
        Bby[n] = {h: set(s) for h, s in B.items()}
    print(f"length {n}: {len(cur)} states with H <= {HMAX}", flush=True)
print(f"(1) ribbon B-sets, words of length <= {LMAX}, R <= {CMAX} h; density = |{{R}}| / #{{R in [h, {CMAX}h] with R = e mod 2 for some observed e}}")
sat = 0
for h in sorted(B):
    Rs = sorted({R for e, R in B[h]})
    es = sorted({e for e, R in B[h]})
    # lattice points: R in [h, CMAX h] with parity of some e; if both parities of e occur, all integers count
    pars = {e % 2 for e in es}
    cand = [R for R in range(h, CMAX * h + 1) if R % 2 in pars]
    dens = len(Rs) / max(1, len(cand))
    same = Bby[LMAX - 4].get(h, set()) == Bby[LMAX].get(h, set())
    sat += same
    if h <= 45 or h in (63, 65, 67, 75, 81, 99):
        gaps = [R for R in cand if R not in set(Rs)]
        print(f"    h={h:3d}: {len(B[h]):4d} (e,R) pairs, R-values {len(Rs):4d} of {len(cand):4d} candidates (density {dens:.2f}), e-values {len(es):3d}; unchanged from length {LMAX-4}: {same}; first R-gaps {gaps[:6]}")
print(f"    B-sets unchanged between length {LMAX-4} and {LMAX}: {sat} of {len(B)} values of h")
# (2) all tournaments of order <= 8, every root
def gentourng(n):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    res = []
    for bits in out:
        adj = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits[k] == '1': adj[i] |= 1 << j
                else: adj[j] |= 1 << i
                k += 1
        res.append(adj)
    return res
def tables(n, adj):
    # h[mask], s[mask][v] (paths of T[mask] starting at v), e[mask][v]
    full = 1 << n
    s = [[0] * n for _ in range(full)]; e = [[0] * n for _ in range(full)]
    # dp over (mask, last) counting paths covering mask ending at last
    dp = [[0] * n for _ in range(full)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, full):
        for v in range(n):
            c = dp[mask][v]
            if not c: continue
            out = adj[v] & ~mask
            while out:
                u = (out & -out).bit_length() - 1; out &= out - 1
                dp[mask | (1 << u)][u] += c
    h = [sum(dp[mask]) for mask in range(full)]
    for mask in range(full):
        for v in range(n): e[mask][v] = dp[mask][v]
    # starts: reverse graph
    radj = [0] * n
    for p in range(n):
        o = adj[p]
        while o:
            q = (o & -o).bit_length() - 1; o &= o - 1; radj[q] |= 1 << p
    dp2 = [[0] * n for _ in range(full)]
    for v in range(n): dp2[1 << v][v] = 1
    for mask in range(1, full):
        for v in range(n):
            c = dp2[mask][v]
            if not c: continue
            out = radj[v] & ~mask
            while out:
                u = (out & -out).bit_length() - 1; out &= out - 1
                dp2[mask | (1 << u)][u] += c
    for mask in range(full):
        for v in range(n): s[mask][v] = dp2[mask][v]
    return h, s, e
Ball = defaultdict(set)
for n in range(1, 9):
    for adj in gentourng(n):
        h, s, e = tables(n, adj)
        full = (1 << n) - 1
        for t in range(n):
            rest = full & ~(1 << t)
            R = 0
            sub = rest
            while True:
                X = sub
                comp = full & ~X
                R += (h[X] if X else 1) * s[comp][t]
                if sub == 0: break
                sub = (sub - 1) & rest
            Ball[h[full]].add((e[full][t], R))
print("(2) all tournaments of orders <= 8, every root: B-sets for small h (R <= %d h):" % CMAX)
for h in [3, 5, 9, 11, 13, 15, 17, 19, 21, 25, 27, 29, 31, 33, 35, 37, 41, 45]:
    S = {(e, R) for e, R in Ball.get(h, set()) if R <= CMAX * h}
    Rs = sorted({R for e, R in S}); es = sorted({e for e, R in S})
    pars = {e % 2 for e in es}
    cand = [R for R in range(h, CMAX * h + 1) if R % 2 in pars]
    rib = sorted({R for e, R in B.get(h, set())})
    print(f"    h={h:3d}: all-tournament R-values {len(Rs):3d} of {len(cand):3d} (density {len(Rs)/max(1,len(cand)):.2f}); ribbon R-values {len(rib):3d}; R-values of all tournaments not in ribbons: {sorted(set(Rs) - set(rib))[:8]}")
