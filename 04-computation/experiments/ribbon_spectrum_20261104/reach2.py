#!/usr/bin/env python3
"""reach2.py -- (1) the exact lattice coset of the reachable states: L_n = Z-span of differences of Y_n
(full set, flint HNF), its stabilisation in n, and the invariance of L under the letters (so Y_n subset
y_0 + L for all n once it holds for the generators: a proof); (2) cumulative (n <= NMAX) multiplicities of
the next-level values N = H(u a) on the solid interval, per letter; (3) the rooted-triple set
T = {(H, R - e, 2R)} cumulative: per-H slice sizes and the fraction of the (e, R)-box that is hit.
Usage: reach2.py NMAX"""
import sys, math
from collections import Counter, defaultdict
import numpy as np, flint
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
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
masks = sorted(Am); hidx = hc.index(1)
levels = []; cur = {tuple(y0)}
for n in range(1, NMAX + 1):
    nxt = set()
    for y in cur:
        v = np.array(y, dtype=np.int64)
        for m in masks:
            if n >= 2 and (m >> (n - 1)): continue
            if n == 1 and m: continue
            nxt.add(tuple((Am[m] @ v).tolist()))
    cur = nxt; levels.append(cur)
print("levels:", [len(l) for l in levels], flush=True)
def hnf_rows(vectors):
    H = flint.fmpz_mat([[int(x) for x in v] for v in vectors]).hnf()
    return [[int(H[i, j]) for j in range(H.ncols())] for i in range(H.nrows()) if any(int(H[i, j]) for j in range(H.ncols()))]
# (1) lattice of differences, per level, from the full level
base = None
for n in range(3, NMAX + 1):
    Yn = sorted(levels[n - 1])
    y_ref = Yn[0]
    diffs = [[a - b for a, b in zip(y, y_ref)] for y in Yn[1:]]
    rows = hnf_rows(diffs)
    piv = [next(x for x in rw if x) for rw in rows]
    print(f"(1) n={n}: lattice of differences rank {len(rows)}, pivots {piv}, index {np.prod(piv)}", flush=True)
    base = (rows, y_ref)
rows, y_ref = base
# invariance of L = Z-span(rows) under the nonempty letters: A_m L subset L  (so y_0 + L is invariant iff A_m y_0 - y_0 in L)
def in_lattice(v, rows):
    v = list(v)
    for rw in rows:
        p = next(c for c, x in enumerate(rw) if x)
        q, rem = divmod(v[p], rw[p])
        if rem: return False
        v = [a - q * b for a, b in zip(v, rw)]
    return not any(v)
for m in masks:
    inv = all(in_lattice((Am[m] @ np.array(rw, dtype=np.int64)).tolist(), rows) for rw in rows)
    shift = in_lattice((Am[m] @ np.array(y_ref, dtype=np.int64) - np.array(y_ref)).tolist(), rows)
    print(f"(1) letter mask {m}: A L subset L: {inv}; A y_ref - y_ref in L: {shift}")
# residues: describe the coset by congruences on coordinates (mod 2, 4, 8) over all levels >= 4
allY = set().union(*levels[3:])
mods = {}
for i in range(r):
    for M in (8, 4, 2):
        res = {y[i] % M for y in allY}
        if len(res) == 1: mods[i] = (M, res.pop()); break
print("(1) fixed residues per coordinate (mod M -> residue), coordinates 1-based:", {i + 1: v for i, v in mods.items()})
# (2) cumulative multiplicities of next-level values by letter
cum = {m: Counter() for m in (1, 2, 3)}
for n, lev in enumerate(levels, start=1):
    for y in lev:
        v = np.array(y, dtype=np.int64)
        for m in (1, 2, 3):
            cum[m][int((Am[m] @ v)[hidx])] += 1
tot = Counter()
for m in (1, 2, 3): tot.update(cum[m])
vals = [v for v in range(49, 3034, 2)]
print("(2) cumulative multiplicity (words of length <= %d+1) on [49, 3033]: per-letter min / total min, median:" % NMAX)
for m in (1, 2, 3):
    mins = [cum[m][v] for v in vals]
    print(f"    letter {m}: min {min(mins)}, zeros {sum(1 for x in mins if x == 0)}, median {sorted(mins)[len(mins)//2]}")
mins = [tot[v] for v in vals]
print(f"    total: min {min(mins)} (at {vals[mins.index(min(mins))]}), median {sorted(mins)[len(mins)//2]}, max {max(mins)}")
# growth of the minimum multiplicity on a fixed window [49, 639] with the length bound
for L in range(8, NMAX + 1):
    c = Counter()
    for lev in levels[:L]:
        for y in lev:
            v = np.array(y, dtype=np.int64)
            for m in (1, 2, 3): c[int((Am[m] @ v)[hidx])] += 1
    w = [c[v] for v in range(49, 640, 2)]
    print(f"    length <= {L + 1}: min multiplicity on [49, 639] = {min(w)}, median {sorted(w)[len(w)//2]}")
# (3) rooted triples (cumulative)
A0 = Am[0]
T = set()
for lev in levels:
    for y in lev:
        z = A0 @ np.array(y, dtype=np.int64)
        H = int(z[hidx]); d = int(z[5]); q = int(z[3])   # d = R - e, q = 2R
        R = q // 2; e = R - d
        T.add((H, e, R))
byH = defaultdict(set)
for H, e, R in T: byH[H].add((e, R))
print(f"(3) cumulative rooted triples (h, e, R) from ribbons of length <= {NMAX}: {len(T)}; distinct h: {len(byH)}")
for H in [9, 11, 13, 15, 17, 21, 25, 31, 45, 63, 99, 129, 257]:
    if H not in byH: continue
    S = byH[H]
    es = sorted({e for e, R in S}); Rs = sorted({R for e, R in S})
    Rset = sorted({R for e, R in S})
    gaps = sum(1 for x in range(Rset[0], Rset[-1] + 1, 2) if x not in set(Rset))
    print(f"    h={H}: {len(S)} (e,R) pairs; e in [{es[0]},{es[-1]}] ({len(es)} values); R in [{Rs[0]},{Rs[-1]}] ({len(Rs)} values, {gaps} same-parity gaps)")
