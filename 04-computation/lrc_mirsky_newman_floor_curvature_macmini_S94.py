#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S94 -- HYP-3850 parts (a),(d),(e).

(a) CONTINUOUS MIRSKY-NEWMAN FLOOR (candidate THM-594, the tower-floor input):
    For a finite arc system D_w(r) with DISTINCT speeds, pick any divisor-minimal
    w0 in F (no other element of F divides it).  Then the coverage Fourier
    coefficient at m=w0 is EXACTLY sin(2*pi*r)/pi (single term), so
        int (C - A)^2 >= 2 sin^2(2 pi r)/pi^2,
    and with u = |{C=0}|, A = 2r|F|, Cmax = max coverage:
        u >= [ 2 sin^2(2 pi r)/pi^2 - (Cmax-1)(A-1) + (A-1)^2 ] / Cmax.        (*)
    At the critical mass A=1 (j = 1/(2r) elements -- EXACTLY where the union bound
    dies, opus-S32) this is u >= 2 sin^2(2 pi r)/(pi^2 Cmax) > 0.
    COROLLARY (rigidity): no finite distinct-speed arc system tiles the circle
    exactly at any r in (0,1/2) -- tiling would need C == 1, contradicting the
    nonzero divisor-minimal coefficient.  Exact-tiling limits (D7(k/7)=0, opus)
    must be INFINITE chain limits = the fixed locus of the scale action.
    VERIFY: u exactly (Fraction) for many 7-clusters at r=1/14 (A=1) incl.
    adversarial divisor-chains; check (*); find the observed min u.

(b=d) SECOND FAREY LAYER CURVATURE: AP {1..13} vs GW {1..11,13,24} tie at first
    order (slope 1666/6435 on [1/15,1/14]).  Does the tie break on the SECOND
    cell (below 1/15, the q+q'=15 layer)?  Exact profiles via lonely_profile.

(c=e) TIGHT => NO MULTIPLE OF n: exhaustive small search for tight sets
    CONTAINING a multiple of n (deep-witness branch), n=5..8, elements <= cap.
"""
from fractions import Fraction as F
from math import gcd, sin, pi
import itertools, sys
sys.path.insert(0, '04-computation')
from lonely_profile import profile

def dist(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def coverage_stats(S, r):
    """Exact uncovered measure u, overlap profile, Cmax for arc system at radius r."""
    events = []  # (position, +1/-1)
    for v in S:
        half = F(r) / v
        for k in range(v):
            c = F(k, v)
            a, b = c - half, c + half
            if a < 0:
                events += [(a + 1, 1), (F(1), -1), (F(0), 1), (b, -1)]
            elif b > 1:
                events += [(a, 1), (F(1), -1), (F(0), 1), (b - 1, -1)]
            else:
                events += [(a, 1), (b, -1)]
    events.sort()
    u = F(0); depth = 0; last = F(0); cmax = 0; l2 = F(0)  # l2 = int C^2
    intCm1sq = F(0)
    for pos, d in events:
        if pos > last:
            seg = pos - last
            if depth == 0:
                u += seg
            l2 += seg * depth * depth
            intCm1sq += seg * (depth - 1) ** 2
            cmax = max(cmax, depth)
        depth += d
        last = pos
    if last < 1:
        seg = 1 - last
        if depth == 0: u += seg
        l2 += seg * depth * depth
        intCm1sq += seg * (depth - 1) ** 2
    return u, cmax, l2, intCm1sq

def mn_bound(S, r):
    """The (*) lower bound; also returns the exact Parseval mass int(C-A)^2."""
    k = len(S)
    A = 2 * F(r) * k
    # divisor-minimal element exists always; coefficient is sin(2 pi r)/pi
    coef2 = 2 * (sin(2 * pi * float(r)) / pi) ** 2
    u, cmax, l2, intCm1sq = coverage_stats(S, r)
    parseval = float(l2) - float(A) ** 2   # int C^2 - A^2 = int (C-A)^2
    bound = (coef2 - (cmax - 1) * (float(A) - 1) + (float(A) - 1) ** 2) / cmax
    return u, cmax, A, parseval, coef2, bound

print("=" * 78)
print("(a) CONTINUOUS MIRSKY-NEWMAN FLOOR at the union-bound death j=7, r=1/14 (A=1)")
print("=" * 78)
r = F(1, 14)
clusters = {
    "consecutive {1..7}": [1, 2, 3, 4, 5, 6, 7],
    "shifted {2..8}": [2, 3, 4, 5, 6, 7, 8],
    "shifted {8..14}": [8, 9, 10, 11, 12, 13, 14],
    "odds {1,3,5,7,9,11,13}": [1, 3, 5, 7, 9, 11, 13],
    "divisor chain {1,2,4,8,16,32,64}": [1, 2, 4, 8, 16, 32, 64],
    "divisor web {2,3,4,6,8,12,24}": [2, 3, 4, 6, 8, 12, 24],
    "coprimes {1,2,3,5,7,11,13}": [1, 2, 3, 5, 7, 11, 13],
    "geometric-ish {3,6,12,24,48,96,192}": [3, 6, 12, 24, 48, 96, 192],
    "near-tiling try {7,14,21,28,35,42,49}": [7, 14, 21, 28, 35, 42, 49],
    "big spread {1,20,41,83,167,335,671}": [1, 20, 41, 83, 167, 335, 671],
}
min_u = None
for name, S in clusters.items():
    u, cmax, A, pv, coef2, bound = mn_bound(S, r)
    ok = float(u) >= bound - 1e-12
    print(f"  {name:38s} u={float(u):.6f}  Cmax={cmax}  A={A}  "
          f"Parseval={pv:.5f}>={coef2:.5f}? {pv >= coef2 - 1e-9}  bound(*)={bound:.6f}  u>=(*)? {ok}")
    if min_u is None or u < min_u[0]:
        min_u = (u, name)
print(f"  min observed u: {float(min_u[0]):.6f} ({min_u[1]}); "
      f"universal floor 2sin^2(2pi r)/(pi^2 Cmax) with Cmax=7: {2*(sin(2*pi/14)/pi)**2/7:.6f}")

# random adversarial search for small u at A=1
import random
random.seed(93)
worst = (F(1), None)
for trial in range(400):
    S = sorted(random.sample(range(1, 120), 7))
    u, cmax, A, pv, coef2, bound = mn_bound(S, r)
    if float(u) < bound - 1e-12:
        print(f"  !! VIOLATION {S}: u={float(u):.6f} < bound {bound:.6f}")
    if u < worst[0]:
        worst = (u, S)
print(f"  400 random 7-clusters: min u = {float(worst[0]):.6f} at {worst[1]}  (no violations printed above)")

print()
print("=" * 78)
print("(d) SECOND FAREY LAYER: AP vs GW curvature below 1/15")
print("=" * 78)
AP = list(range(1, 14))
GW = list(range(1, 12)) + [13, 24]
for name, S in [("AP", AP), ("GW-24", GW)]:
    p = profile(S, F(1, 14))
    cells = [(a, b, sl) for a, b, sl, _ in p._cells if b > F(1, 17)]
    print(f"  {name}: cells above 1/17:")
    for a, b, sl in cells:
        print(f"    [{a}, {b}]  slope {sl} = {float(sl):.6f}   m({b})={float(p.measure(b)):.8f}")
pAP = profile(AP, F(1, 14)); pGW = profile(GW, F(1, 14))
for rr in [F(1, 15), F(1, 16), F(2, 31), F(1, 17)]:
    dAP, dGW = pAP.measure(rr), pGW.measure(rr)
    cmp = "AP<GW" if dAP < dGW else ("AP>GW" if dAP > dGW else "TIE")
    print(f"  m at r={rr}: AP={dAP} = {float(dAP):.8f}   GW={dGW} = {float(dGW):.8f}   {cmp}")

print()
print("=" * 78)
print("(e) TIGHT => NO MULTIPLE OF n: search the deep-witness branch")
print("=" * 78)
def M_exact(S):
    Sl = sorted(set(S)); dens = set()
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            dens.add(v + w)
            if w > v: dens.add(w - v)
    best = F(0)
    for den in dens:
        for m in range(1, den):
            t = F(m, den)
            mn = min(dist(v * t) for v in Sl)
            if mn > best: best = mn
    return best

for n, cap in [(5, 40), (6, 36), (7, 30)]:
    hits = []
    mults = [m for m in range(n, cap + 1, n)]
    others = [v for v in range(1, cap + 1) if v % n != 0]
    count = 0
    for w in mults:
        for rest in itertools.combinations(others, n - 2):
            S = sorted(rest + (w,))
            count += 1
            if M_exact(S) == F(1, n):
                hits.append(S)
    print(f"  n={n} (|S|={n-1}, one multiple of n, elements<={cap}): "
          f"{count} sets tested, tight-with-multiple found: {len(hits)} {hits[:3]}")
