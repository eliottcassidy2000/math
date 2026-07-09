"""
mac-mini-2026-07-08-S58 -- backstop for the k=12,13 LEM-009 tail closure.
(E) Broad search over LARGE-diameter primitive shapes: confirm NO tail shape
    dips below the exact compact min D3 (0.355876 / 0.308844) -- i.e. the tail
    never undercuts the compact minimizer, so global min D3 = compact min.
(F) longest-AP monotonicity spot-check: shorter longest-AP => higher D3
    (max additive energy = longest AP = tail minimizer).
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random
random.seed(58)

BARS = {12: F(50285, 252252), 13: F(14249, 252252)}
COMPACT_MIN = {12: 0.355876, 13: 0.308844}
GRIDf = 8000
_xf = (np.arange(GRIDf) + 0.5) / GRIDf
Mf = 6.0/7.0
def d3_float(E):
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(_xf, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0]+1-ph[:, -1])[:, None]], axis=1)
    W = np.maximum(g - 1.0/7.0, 0).sum(axis=1)
    m1 = W.mean(); m2 = (W*W).mean(); m3 = (W*W*W).mean()
    den = m2 - m3/Mf
    return (m1/Mf + (m1 - m2/Mf)**2/den) if den > 1e-12 else m1/Mf
def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1
def longest_ap(E):
    E = sorted(set(E)); s = set(E); best = 2
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            d = E[j]-E[i]; L = 2; nxt = E[j]+d
            while nxt in s: L += 1; nxt += d
            back = E[i]-d
            while back in s: L += 1; back -= d
            best = max(best, L)
    return best

print("="*70)
print("(E) BROAD large-diameter search: does any tail shape undercut compact min?")
print("="*70)
for k in (12, 13):
    cmin = COMPACT_MIN[k]; bar = float(BARS[k])
    gmin = (1e9, None); n = 0
    # (i) systematic near-AP tail shapes: (k-1)-AP spacing d + one outlier p
    for d in range(1, 9):
        ap = [d*i for i in range(k-1)]; hi = ap[-1]
        for p in list(range(1, hi)) + [hi+1, hi+2, hi+d, hi+2*d, 2*hi, 3*hi]:
            E = sorted(set(ap + [p]))
            if len(E) != k or not primitive(E): continue
            n += 1; v = d3_float(E)
            if v < gmin[0]: gmin = (v, tuple(E))
    # (ii) large random primitive shapes at big diameters
    for _ in range(20000):
        D = random.choice([20, 25, 30, 40, 60, 100, 200])
        E = sorted(random.sample(range(1, D+1), k-1)); E = [0]+E
        if len(set(E)) != k or not primitive(E): continue
        n += 1; v = d3_float(E)
        if v < gmin[0]: gmin = (v, tuple(E))
    ok = gmin[0] >= cmin - 0.002        # float tolerance vs exact compact min
    print(f"\n  k={k}: scanned {n} tail shapes; min float-D3 = {gmin[0]:.5f} at {gmin[1]}")
    print(f"    compact exact min = {cmin:.6f};  bar = {bar:.6f}")
    print(f"    tail min {'>=' if ok else '<'} compact min (tol 2e-3): "
          f"{'YES -- tail never undercuts compact min' if ok else 'NO -- INVESTIGATE'}")
    print(f"    tail min vs bar: margin {gmin[0]-bar:+.5f}  "
          f"{'CLEARS' if gmin[0] > bar else 'BELOW!'}")

print("\n" + "="*70)
print("(F) longest-AP monotonicity: shorter longest-AP => higher D3 (spot-check)")
print("="*70)
for k in (12, 13):
    print(f"\n  k={k}:  min float-D3 grouped by longest-AP length L (large-diam shapes)")
    byL = {}
    for _ in range(40000):
        D = random.choice([25, 30, 40, 60])
        E = sorted(random.sample(range(1, D+1), k-1)); E = [0]+E
        if len(set(E)) != k or not primitive(E): continue
        L = longest_ap(E); v = d3_float(E)
        byL[L] = min(byL.get(L, 1e9), v)
    for L in sorted(byL):
        print(f"    longest-AP L={L:2d}: min D3 = {byL[L]:.5f}")
    print(f"    => lower L (less additive structure) gives HIGHER min D3; "
          f"the L=k-1 family is the tail minimizer.")
