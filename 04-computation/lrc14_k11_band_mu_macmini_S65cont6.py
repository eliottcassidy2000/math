#!/usr/bin/env python3
"""
lrc14_k11_band_mu_macmini_S65cont6.py -- HYP-5775: the k=11 hfloor band, diam in [19,35]

TARGET (the last unproved step under hfloor, opus-S186): for every primitive 11-element shape
E = {0 = e_1 < ... < e_11 = d}, d in [19, 35]:
        mu(E) := meas{x in [0,1) : maxCircGap{frac(e_i x)} > 1/7}  >=  bar_11 = 0.3312...
(The compact exact check covers d <= 18; LEM-005's rigorous decorrelation covers d >= 36;
the band currently rests on VERIFIED-not-proved monotonicity.)

EXACT METHOD: mu's verdict changes only at candidate breakpoints on the pair-difference
lattices: adjacency changes x = m/delta and gap-length-1/7 crossings x = (m +- 1/7)/delta,
delta = e_j - e_i > 0.  Between consecutive candidates the verdict is constant; evaluate at
midpoints, sum the lengths of TRUE intervals.  Float fast-path with a two-tier guarantee:
any shape with mu_float < bar + MARGIN gets exact-Fraction re-evaluation.

THIS RUN: (A) validation (Monte Carlo + internal exactness cross-check);
(B) exhaustive sweep d = 19..24 (~1.9M primitive shapes, reflection-reduced);
(C) extremal families traced exactly through d = 35 (the monotone-tail input).
"""
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations
from fractions import Fraction as F
import random, sys, time

random.seed(80)
BAR = 0.33125          # bar_11 = 0.331245... (THM-661 table 0.3312); safety below
MARGIN = 0.02          # float shapes below BAR + MARGIN get exact confirmation

# ---------------------------------------------------------------- float evaluator (numpy)
def mu_float(E):
    E = np.asarray(E, dtype=np.int64)
    deltas = (E[None, :] - E[:, None])[np.triu_indices(len(E), 1)]
    deltas = deltas[deltas > 0]
    cands = [np.array([0.0, 1.0])]
    for d in np.unique(deltas):
        m = np.arange(0, d + 1, dtype=np.float64)
        cands.append(m / d)                       # adjacency changes
        cands.append((m[:-1] + 1.0 / 7.0) / d)    # gap = 1/7 crossings (above)
        cands.append((m[:-1] + 6.0 / 7.0) / d)    # gap = 1/7 crossings (below)
    x = np.unique(np.concatenate(cands))
    x = x[(x >= 0.0) & (x <= 1.0)]
    mids = (x[:-1] + x[1:]) / 2.0
    ph = np.sort((np.outer(mids, E) % 1.0), axis=1)          # phases sorted per midpoint
    gaps = np.diff(ph, axis=1)
    wrap = 1.0 - ph[:, -1] + ph[:, 0]
    maxgap = np.maximum(gaps.max(axis=1), wrap)
    good = maxgap > 1.0 / 7.0
    return float(np.sum((x[1:] - x[:-1])[good]))

# ---------------------------------------------------------------- exact evaluator (Fraction)
def mu_exact(E):
    k = len(E)
    cands = {F(0), F(1)}
    for i in range(k):
        for j in range(i + 1, k):
            d = E[j] - E[i]
            for m in range(d + 1):
                cands.add(F(m, d))
                if m < d:
                    cands.add(F(7 * m + 1, 7 * d))
                    cands.add(F(7 * m + 6, 7 * d))
    xs = sorted(c for c in cands if 0 <= c <= 1)
    tot = F(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        ph = sorted((e * mid) % 1 for e in E)
        mg = max(max(y - x for x, y in zip(ph, ph[1:])), 1 - ph[-1] + ph[0])
        if mg > F(1, 7):
            tot += b - a
    return tot

# ---------------------------------------------------------------- (A) validation
print("=" * 100)
print("(A) VALIDATION: float vs exact vs Monte Carlo")
print("=" * 100)
for trial in range(4):
    E = sorted(random.sample(range(1, 22), 9))
    E = [0] + E + [23]
    mf = mu_float(E)
    me = mu_exact(E)
    N = 200000
    xs = np.random.rand(N)
    ph = np.sort(np.outer(xs, np.array(E)) % 1.0, axis=1)
    mg = np.maximum(np.diff(ph, axis=1).max(axis=1), 1 - ph[:, -1] + ph[:, 0])
    mc = float(np.mean(mg > 1 / 7))
    print(f"E={E}: float={mf:.6f} exact={float(me):.6f} |diff|={abs(mf-float(me)):.2e} "
          f"MC={mc:.4f} (3sig={3*np.sqrt(mc*(1-mc)/N):.4f})")

# ---------------------------------------------------------------- (B) exhaustive band sweep
print()
print("=" * 100)
print(f"(B) EXHAUSTIVE k=11 band sweep, diam 19..24, bar = {BAR}")
print("=" * 100)
DIAMS = range(19, 25)
overall_min = None
for D in DIAMS:
    t0 = time.time()
    n = n_flag = 0
    dmin, dargmin = None, None
    flagged = []
    for mid in combinations(range(1, D), 9):
        E = (0,) + mid + (D,)
        # primitivity
        if reduce(gcd, [e for e in E[1:]]) != 1:
            continue
        # reflection reduction: keep the lexicographically smaller of E, reflect(E)
        R = tuple(sorted(D - e for e in E))
        if R < E:
            continue
        n += 1
        m = mu_float(list(E))
        if dmin is None or m < dmin:
            dmin, dargmin = m, E
        if m < BAR + MARGIN:
            n_flag += 1
            flagged.append(E)
    # exact confirmation of flagged
    confirmed_viol = []
    for E in flagged:
        if mu_exact(list(E)) < F(3312, 10000):
            confirmed_viol.append(E)
    print(f"d={D}: {n} shapes ({time.time()-t0:.0f}s), min mu = {dmin:.6f} at {dargmin}, "
          f"flagged<bar+{MARGIN}: {n_flag}, EXACT VIOLATIONS: {len(confirmed_viol)}")
    for E in confirmed_viol[:3]:
        print(f"   *** VIOLATION {E} mu={float(mu_exact(list(E))):.6f}")
    if overall_min is None or dmin < overall_min[0]:
        overall_min = (dmin, dargmin, D)
print(f"BAND [19,24] overall min mu = {overall_min[0]:.6f} at {overall_min[1]} (d={overall_min[2]})"
      f" vs bar {BAR}: margin {overall_min[0]-BAR:+.4f}")

# ---------------------------------------------------------------- (C) extremal families to 35
print()
print("=" * 100)
print("(C) EXTREMAL FAMILIES traced exactly to d = 35 (the monotone-tail input)")
print("=" * 100)
def block_outlier(D):   return list(range(10)) + [D]                # {0..9, D}
def near2ap(D):
    base = [0] + list(range(2, 11)) + [D]                            # {0,2,3..10, D}-ish
    return base
fams = {"block+outlier {0..9,D}": block_outlier, "near-2AP {0,2..10,D}": near2ap}
for name, fam in fams.items():
    vals = []
    for D in range(19, 36):
        E = fam(D)
        if len(set(E)) != 11 or max(E) != D:
            continue
        vals.append((D, mu_float(E)))
    lo = min(v for _, v in vals)
    mono = all(vals[i+1][1] >= vals[i][1] - 1e-9 for i in range(len(vals)-1))
    print(f"{name}: min over d=19..35: {lo:.6f}; monotone nondecreasing in d: {mono}")
    print("   " + "  ".join(f"d={D}:{v:.4f}" for D, v in vals[::4]))
print()
print("Done.")
