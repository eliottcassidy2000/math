#!/usr/bin/env python3
"""
opus-2026-07-05-S86 -- HYP-4146: cluster desert localization + resonance immunity.

A c-cluster C = {v_1 < ... < v_c} (distinct integers, ratio v_c/v_1 < 26/11, scale W)
at band rho = 1/13.  DESERT = maximal interval where every t has some v_i within 1/13
of an integer (the union of the c combs covers it).  The window/descent machinery dies
iff the base window lands inside a desert.

THE RESONANCE PICTURE (to verify): within a short window at t0, phases split as a fast
common rotation theta = v_1 t plus slow offsets {d_j t}, d_j = v_j - v_1 <= (15/11)W.
Covering near t0 <=> the offset set {d_j t0 mod 1} is a 2/13-net of the circle.  At
t0 = p/q the offsets sit on the q-grid: a 2/13-net needs q >= 7 (1/q <= 2/13).  So
deserts should LOCALIZE near rationals p/q, q >= 7, satisfying the net condition, and
their total measure should be small.

MEASUREMENTS:
 1. exact desert census for random and adversarial clusters (c = 7..12, W = 150..3000):
    max desert length x W (the K(c) constant), desert count, total desert measure;
 2. desert LOCATIONS vs rationals with q in [7, 13]: localization check;
 3. immunity check: distance from t = 1/6, 1/4, 1/3, 5/6 (q <= 6 points) to the
    nearest desert, in units of 1/W (should be >> 1);
 4. the worst structure: consecutive blocks vs random vs two-scale clusters.

Exact interval arithmetic (Fractions): the union of teeth on [0,1] is computed exactly;
deserts = gaps of the complement.
"""
from fractions import Fraction as F
import random, sys

RHO = F(1, 13)

def deserts(cluster, lo=F(0), hi=F(1)):
    """exact: merge all teeth intervals of the cluster in [lo,hi]; return the maximal
    covered intervals (deserts of the GOOD set = components of the union)."""
    ivs = []
    for w in cluster:
        m0 = int(w * lo)
        m1 = int(w * hi) + 1
        for m in range(m0, m1 + 1):
            a, b = (F(m) - RHO) / w, (F(m) + RHO) / w
            if b > lo and a < hi:
                ivs.append((max(a, lo), min(b, hi)))
    ivs.sort()
    merged = []
    for a, b in ivs:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged

def analyze(cluster, label):
    W = min(cluster)
    cov = deserts(cluster)
    lens = [(b - a) for a, b in cov]
    mx = max(lens)
    mi = max(range(len(lens)), key=lambda i: lens[i])
    tot = sum(lens)
    # localization: nearest rational p/q (q <= 13) to the longest desert's midpoint
    mid = (cov[mi][0] + cov[mi][1]) / 2
    best = None
    for q in range(1, 14):
        p = round(mid * q)
        d = abs(mid - F(p, q))
        if best is None or d < best[2]:
            best = (p, q, d)
    # immunity: distance from q<=6 grid points to nearest desert
    imm = []
    for t0 in (F(1,6), F(1,4), F(1,3), F(5,6), F(1,2)):
        dmin = min(min(abs(t0 - a), abs(t0 - b)) for a, b in cov)
        inside = any(a <= t0 <= b for a, b in cov)
        imm.append((t0, "INSIDE" if inside else float(dmin * W)))
    print(f"{label}: c={len(cluster)} W={W} maxdesert*W={float(mx*W):.2f} "
          f"count={len(lens)} total={float(tot):.4f} "
          f"longest@{float(mid):.4f}~{best[0]}/{best[1]} (dist {float(best[2]):.2e})")
    print(f"    immunity (dist to desert in 1/W units): "
          + ", ".join(f"{t0}={v if isinstance(v,str) else round(v,1)}" for t0, v in imm))
    return mx, tot

random.seed(86)
print("=" * 100)
print("1. CONSECUTIVE clusters (the S85-extremal shape)")
for c in (7, 9, 12):
    for W in (150, 400, 1000):
        analyze(list(range(W, W + c)), f"consec c={c}")
print()
print("2. RANDOM clusters (ratio < 26/11)")
for c in (7, 9):
    for W in (150, 400, 1000):
        hi = W * 26 // 11 - 1
        cl = sorted(random.sample(range(W, hi), c))
        analyze(cl, f"random")
print()
print("3. ADVERSARIAL: arithmetic clusters d_j multiples of q (resonance-seeking, q=7)")
for W in (150, 400, 1000):
    # differences all multiples chosen so offsets at p/7 form the full 7-grid
    base = W
    cl = [base + j * (W // 8) + j for j in range(7)]  # spread, differences ~W/8, distinct
    cl = sorted(set(cl))[:7]
    analyze(cl, "adversarial-ish")
print()
print("4. TWO-SCALE: near-coincident pairs inside a 7-cluster")
for W in (400, 1000):
    cl = [W, W + 1, W + W // 3, W + W // 3 + 1, W + 2 * W // 3, W + 2 * W // 3 + 1, W + W - W // 10]
    cl = sorted(set(cl))[:7]
    analyze(cl, "paired")
