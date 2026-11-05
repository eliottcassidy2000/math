#!/usr/bin/env python3
"""reach.py -- structure of the reachable set Y_n of reduced states of span-2 ribbon words.
Reads span2_intrep.txt (rank-14 integer representation; y_1 = y_2 = H on nonempty words) and
span2_condensation.txt; builds Y_n for n <= NMAX (all four letters) and reports:
 (1) the Z-lattice spanned by differences of states (is it the full hyperplane lattice?), parity classes hit;
 (2) the cone: extreme ratios y_i / H;
 (3) 2D projection densities (H, y_i): hit fraction of the odd-H lattice points inside the convex hull;
 (4) the post-sink triples (H, y_3, q) = (h, R - e, 2R): per-H slice sizes and the density of (y_3, q);
 (5) the H-functionals of the three nonempty letters: H' = H + y_3 | H + y_3 + y_4 | H + y_7, and the
     coverage of [c, d] by the three sumsets, with multiplicities.
Usage: reach.py NMAX  (run inside 04-computation/experiments/ribbon_spectrum_20261104)"""
import sys, itertools, math
from collections import Counter, defaultdict
import numpy as np
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 11

def load_intrep(path):
    with open(path) as f:
        r, nl = map(int, f.readline().split())
        y0 = list(map(int, f.readline().split()))
        hc = list(map(int, f.readline().split()))
        A = {}
        for _ in range(nl):
            mask = int(f.readline().strip())
            A[mask] = [list(map(int, f.readline().split())) for _ in range(r)]
    return r, y0, hc, A
r, y0, hc, A = load_intrep('span2_intrep.txt')
Am = {m: np.array(M, dtype=np.int64) for m, M in A.items()}
masks = sorted(Am)
hidx = hc.index(1)
levels = []
cur = {tuple(y0)}
for n in range(1, NMAX + 1):
    nxt = set()
    for y in cur:
        v = np.array(y, dtype=np.int64)
        for m in masks:
            if n >= 2 and (m >> (n - 1)): continue
            if n == 1 and m: continue
            nxt.add(tuple((Am[m] @ v).tolist()))
    cur = nxt
    levels.append(cur)
    print(f"n={n}: {len(cur)} states", flush=True)
Y = np.array(sorted(levels[-1]), dtype=np.int64)
H = Y[:, hidx]
print("coordinate 1 == coordinate 2 on all states:", bool(np.all(Y[:, 0] == Y[:, 1])))
# (1) lattice of differences
D = Y[1:] - Y[0]
import flint
sub = D[np.random.RandomState(0).choice(len(D), size=min(len(D), 400), replace=False)]
Hn = flint.fmpz_mat([[int(x) for x in row] for row in sub]).hnf()
rows = [[int(Hn[i, j]) for j in range(Hn.ncols())] for i in range(Hn.nrows()) if any(int(Hn[i, j]) for j in range(Hn.ncols()))]
print(f"(1) difference lattice (from 400 random differences): rank {len(rows)}; pivot entries {[next(x for x in rw if x) for rw in rows]}; index in the saturated lattice = product of pivots = {np.prod([next(x for x in rw if x) for rw in rows])}")
par = Counter(tuple(int(x) % 2 for x in y) for y in Y)
print(f"(1) parity classes hit at n={NMAX}: {len(par)} of {2**13} (on the 13-dim hyperplane); H odd always: {bool(np.all(H % 2 == 1))}")
# (2) cone
print("(2) extreme ratios y_i/H over states with H >= 9:")
mask9 = H >= 9
for i in range(r):
    rat = Y[mask9, i] / H[mask9]
    print(f"    y_{i+1}/H in [{rat.min():+.3f}, {rat.max():+.3f}], mean {rat.mean():+.3f}")
# (3) 2D projection densities (H, y_i)
print("(3) 2D projection (H, y_i): distinct points vs lattice points of the bounding region (same parity as observed):")
for i in range(r):
    if i in (0, 1): continue
    pts = set(zip(H.tolist(), Y[:, i].tolist()))
    byH = defaultdict(list)
    for h, v in pts: byH[h].append(v)
    # hull proxy: for each H, the range [min, max] of y_i, step by observed gcd parity
    tot = 0; hit = len(pts)
    for h, vs in byH.items():
        lo, hi = min(vs), max(vs)
        step = 0
        for v in vs: step = math.gcd(step, v - lo)
        step = step or 1
        tot += (hi - lo) // step + 1
    print(f"    y_{i+1}: {hit} points, {tot} lattice points in per-H ranges, density {hit / tot:.3f}")
# (4) post-sink triples
A0 = Am[0]
trip = set()
for y in Y:
    z = A0 @ y
    trip.add((int(z[hidx]), int(z[5]), int(z[3])))  # (H, y_6 = y_3-image, y_4-image = q)
byH = defaultdict(set)
for h, e, q in trip: byH[h].add((e, q))
sizes = sorted(len(s) for s in byH.values())
print(f"(4) post-sink triples: {len(trip)} distinct (H, R-e, 2R); H-slices: {len(byH)}; slice sizes min/median/max {sizes[0]}/{sizes[len(sizes)//2]}/{sizes[-1]}")
# (5) sumset coverage by the three letters
cov = Counter()
for y in Y:
    h = int(y[hidx])
    for m in (1, 2, 3):
        cov[int((Am[m] @ y)[hidx])] += 1
vals = sorted(cov)
odd_lo, odd_hi = 49, max(vals)
missing = [v for v in range(odd_lo, odd_hi + 1, 2) if v not in cov]
solid_hi = odd_lo - 2
for v in range(odd_lo, odd_hi + 1, 2):
    if v in cov: solid_hi = v
    else: break
mult = [cov[v] for v in range(odd_lo, solid_hi + 1, 2)]
print(f"(5) next-level values from level {NMAX}: solid [49, {solid_hi}], multiplicity on it min/median/max {min(mult)}/{sorted(mult)[len(mult)//2]}/{max(mult)}; first misses above: {missing[:8]}")
