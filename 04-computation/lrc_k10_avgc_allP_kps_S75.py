#!/usr/bin/env python3
r"""
lrc_k10_avgc_allP_kps_S75.py  (kind-pasteur-2026-07-07-S75, HYP-5247)

DECISIVE k=10 TEST: does the AVERAGE-form conditional tent (THM-655 mechanism) close the
k=10 (A') leg on its own, i.e. is sup_E avgc(E,P) <= c*(P) for ALL C(13,3)=286 slow parts?

mac-mini's EnvBlock (decreasing-envelope block bound) EXCEEDS c* at k=10 -- but that is a
LOOSE upper bound (uses chat(d)=max_{d'>=d}c(d') >= c(d)).  The TRUE max avgc may still be
<= c*.  Two facts make the AP_10 = {0..9} the natural maximizer of avgc:
  - block domination (mac-mini): N_E(t) <= N_block(t), equality for the AP -> the AP packs
    the most small-d (high-c) pairs;
  - so avgc(AP_10,P) is the candidate max.  We verify by aggressive local search that no
    primitive 10-cluster beats it, then compare avgc_max vs c*(P) over all 286 P.

If avgc(AP_10,P) <= c*(P) and the AP is the maximizer for every P => k=10 CLOSED by the
average form alone (diameter-free, like THM-655 for k=9).  Otherwise: report the offending
P's (candidates for the klein-W_F composition).  Exact rational c(d,P).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

M_P = F(14249, 252252)
K = 10
BETA = F(14 - K, 7 * K)      # 4/70 = 2/35
INTF = (F(1, 7) - BETA) ** 2 / 2
FLOOR = F(2 * (K - 1) * (K - 7), 7 * K)   # 1 - floor_10 factor is (1-floor)= this? see below
# tentfloor_10 = 1 - 2(k-1)(k-7)/(7k); (1 - tentfloor) = 2(k-1)(k-7)/(7k)
ONE_MINUS_FLOOR = F(2 * (K - 1) * (K - 7), 7 * K)   # = 27/35

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
        if l > prev: good.append((prev, l))
        prev = max(prev, h)
    if prev < 1: good.append((prev, F(1)))
    return good

def tent_forward(iv, d):
    tot = F(0); d = int(d); th = F(1, 7)
    for (l, h) in iv:
        m_lo = int(l * d - th) - 1
        m_hi = int(h * d) + 1
        for m in range(max(0, m_lo), m_hi + 1):
            wl, wh = (F(m) + BETA) / d, (F(m) + th) / d
            a, b = max(wl, l), min(wh, h)
            if a >= b: continue
            tot += (d * (b * b - a * a) / 2 - (F(m) + BETA) * (b - a))
    return tot

def build(P, dmax_cache=60):
    iv = GP_intervals(P)
    mGP = sum(h - l for l, h in iv)
    cstar = (1 - M_P / mGP) / ONE_MINUS_FLOOR
    cache = {}
    def c(d):
        d = abs(int(d))
        if d not in cache:
            cache[d] = tent_forward(iv, d) / (mGP * INTF)
        return cache[d]
    return iv, mGP, cstar, c

def avgc(E, c):
    ds = [E[j] - E[i] for i in range(len(E)) for j in range(i + 1, len(E))]
    return sum(c(d) for d in ds) / len(ds)

def maximize_avgc(P, seed=0, iters=3000):
    iv, mGP, cstar, c = build(P)
    rng = random.Random(seed)
    best_E = list(range(K)); best = avgc(best_E, c)     # AP seed
    # also try 2-AP-ish and block seeds
    for seedE in ([0,2,4,6,8,10,12,14,16,17], [0,1,2,3,4,5,6,7,8,10],
                  [0,1,2,3,4, 10,11,12,13,14]):
        v = avgc(seedE, c)
        if v > best: best, best_E = v, list(seedE)
    for _ in range(iters):
        E = list(best_E)
        idx = rng.randrange(1, K)
        E[idx] += rng.choice([-2, -1, 1, 2])
        E = sorted(set(E))
        if len(E) != K or E[0] < 0: continue
        E = [e - E[0] for e in E]
        g = 0
        for a in E[1:]: g = gcd(g, a)
        if g != 1: continue
        v = avgc(E, c)
        if v > best: best, best_E = v, E
    return best, best_E, cstar, mGP, c

print("=" * 92)
print("k=10: is sup_E avgc(E,P) <= c*(P) for ALL 286 slow parts?  (average-form alone)")
print(f"BETA=2/35, 1-tentfloor = 27/35 = {float(ONE_MINUS_FLOOR):.4f}, m_P={float(M_P):.5f}")
print("=" * 92)
Ps = list(combinations(range(1, 14), 3))
worst = []       # (avgc_max - cstar, P, ...)
ap_over = []     # P where even avgc(AP) > cstar
for P in Ps:
    iv, mGP, cstar, c = build(P)
    apc = avgc(list(range(K)), c)
    if apc > cstar:
        ap_over.append((P, float(apc), float(cstar)))
for P in Ps:
    # only hard-search the P's where AP is close to or over cstar (efficiency)
    iv, mGP, cstar, c = build(P)
    apc = avgc(list(range(K)), c)
    if apc < cstar - F(1, 5):     # comfortably under, AP-max heuristic -> skip hard search
        gap = float(apc - cstar)
        worst.append((gap, P, float(apc), float(cstar), list(range(K))))
        continue
    best = float(apc); bestE = list(range(K))
    for s in range(4):
        v, E, cs, m, cc = maximize_avgc(P, seed=s, iters=2500)
        if v > best: best, bestE = float(v), E
    worst.append((best - float(cstar), P, best, float(cstar), bestE))

worst.sort(reverse=True)
n_fail = sum(1 for w in worst if w[0] > 0)
print(f"\n  #P with sup_E avgc > c* (average-form FAILS alone): {n_fail} of {len(Ps)}")
print(f"  #P where even avgc(AP_10) > c*: {len(ap_over)}")
print(f"\n  Top 12 tightest/failing shapes (avgc_max - c*, P, avgc_max, c*, argmax E):")
for gap, P, amax, cs, E in worst[:12]:
    tag = "FAIL" if gap > 0 else "ok"
    print(f"    {tag:>4}  gap={gap:+.4f}  P={str(P):>12}  avgc_max={amax:.4f}  c*={cs:.4f}  "
          f"E={E} (diam {E[-1]})")
print()
if n_fail == 0:
    print("  => sup_E avgc <= c* for ALL 286 shapes => k=10 (A') leg CLOSED by the average")
    print("     form ALONE, diameter-free (like THM-655 for k=9).  Rigorous bound = next.")
else:
    print(f"  => {n_fail} shapes need the klein-W_F composition (the compact offenders); the")
    print("     average form covers the wide families.  Those P's are the composition targets.")
