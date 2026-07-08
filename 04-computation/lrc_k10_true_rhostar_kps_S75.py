#!/usr/bin/env python3
r"""
lrc_k10_true_rhostar_kps_S75.py  (kind-pasteur-2026-07-07-S75, HYP-5247)

Is the k=10 (A') leg TRUE?  Directly minimize the ACTUAL rho*(P,E) = meas(G_P cap
{cluster maxgap > 1/7}) over 10-element primitive clusters E and 3-element slow parts P,
independent of any tent bound.  If min true rho* >> m_P = 14249/252252, k=10 holds and the
only task is a proof (the composed conditional tent); if it dips near/below m_P, find where.

rho* computed by exact interval integration is expensive; here use a fine deterministic
grid over G_P (adaptive per interval) -- accurate to ~1e-4, enough to locate the true
minimum and confirm the margin.  Cross-check the worst finds at higher resolution.
"""
from fractions import Fraction as F
from math import gcd
import numpy as np
import random

M_P = 14249 / 252252  # ~0.056487

def GP_intervals(P):
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            bad.append((F(j, p) - w, F(j, p) + w))
    bad = [(max(l, F(0)), min(h, F(1))) for l, h in bad if h > 0 and l < 1]
    bad.sort()
    merged = []
    for l, h in bad:
        if merged and l <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], h))
        else:
            merged.append((l, h))
    good = []; prev = F(0)
    for l, h in merged:
        if l > prev: good.append((float(prev), float(l)))
        prev = max(prev, h)
    if prev < 1: good.append((float(prev), 1.0))
    return good

def rhostar(E, P, npts=60000):
    """meas(G_P cap {maxgap(E) > 1/7}), numeric."""
    iv = GP_intervals(P)
    meas = sum(h - l for l, h in iv)
    Ev = np.array(E, dtype=np.float64)
    total = 0.0
    for l, h in iv:
        m = max(50, int(npts * (h - l)))
        xs = l + (h - l) * (np.arange(m) + 0.5) / m
        ph = np.mod(np.outer(xs, Ev), 1.0)
        ph.sort(axis=1)
        gaps = np.diff(ph, axis=1)
        wrap = (ph[:, 0] + 1.0 - ph[:, -1])[:, None]
        mg = np.maximum(gaps.max(axis=1), wrap[:, 0])
        total += (h - l) * np.mean(mg > 1.0 / 7.0)
    return total, meas

# the binding |P|=3 slow parts (min-meas + hard-core + mixed)
PSET = [(1, 12, 13), (11, 12, 13), (1, 2, 13), (10, 11, 12), (1, 2, 3), (2, 12, 13), (1, 11, 13)]
print("=" * 88)
print("k=10 TRUE rho* minimization over primitive 10-clusters x binding slow parts")
print(f"m_P = {M_P:.6f}")
print("=" * 88)

def rand_primitive(k, dmin, dmax, rng):
    for _ in range(200):
        pts = sorted({0, rng.randint(dmin, dmax)} | set(rng.sample(range(1, dmax), k - 2)))
        if len(pts) == k and pts[-1] - pts[0] >= dmin:
            g = 0
            for a in pts[1:]:
                g = gcd(g, a)
            if g == 1:
                return pts
    return None

rng = random.Random(75)
overall = (1.0, None, None)
# structured seeds known to be dangerous: blocks, 2-APs, dilated-AP+bump
seeds = [
    [0,1,2,3,4,5,6,7,8,9],            # AP_10 (compact)
    [0,2,4,6,8,10,12,14,16,17],       # near-2AP primitive
    [0,1,2,3,4,15,16,17,18,19],       # two blocks
    [0,1,2,3,4,5,6,7,8,20],           # block + far
    [0,2,4,6,8,10,12,14,16,25],       # 2AP + far primitive
    [0,3,6,9,12,15,18,21,24,25],      # 3AP + bump
    [0,1,2,3,4,5,15,16,17,18],
    [0,5,10,15,20,25,30,35,40,41],
]
for k in (10,):
    for P in PSET:
        best = (1.0, None)
        for s in seeds:
            if len(s) == k:
                r, meas = rhostar(s, P, npts=40000)
                if r < best[0]:
                    best = (r, list(s))
        # random search across diam ranges
        for dmax in (10, 14, 18, 25, 40, 70):
            for _ in range(90):
                E = rand_primitive(k, max(9, dmax - 8), dmax, rng)
                if E is None: continue
                r, meas = rhostar(E, P, npts=25000)
                if r < best[0]:
                    best = (r, E)
        # local hill-descent on the best
        cur_r, cur_E = best
        for _ in range(250):
            E = list(cur_E)
            idx = rng.randrange(1, k)
            E[idx] += rng.choice([-2, -1, 1, 2])
            E = sorted(set(E))
            if len(E) != k or E[0] < 0: continue
            g = 0
            for a in E[1:]: g = gcd(g, a - E[0])
            if g != 1: continue
            r, _ = rhostar([e - E[0] for e in E], P, npts=30000)
            if r < cur_r:
                cur_r, cur_E = r, [e - E[0] for e in E]
        r_hi, meas = rhostar(cur_E, P, npts=120000)  # high-res recheck
        flag = "OK" if r_hi >= M_P else "*** BELOW m_P ***"
        print(f"  P={str(P):>12} meas(G_P)={meas:.4f}: min rho* ~ {r_hi:.5f} at E={cur_E} (diam {cur_E[-1]-cur_E[0]})  [{flag}]")
        if r_hi < overall[0]:
            overall = (r_hi, cur_E, P)
print()
print(f"OVERALL min true rho* ~ {overall[0]:.5f} at E={overall[1]}, P={overall[2]}")
print(f"  vs m_P = {M_P:.5f}  =>  margin {overall[0] - M_P:+.5f}  "
      f"({overall[0]/M_P:.2f}x)")
print("  (numeric; the true minimizer + margin guide the composed-bound proof)")
