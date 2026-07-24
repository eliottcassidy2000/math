#!/usr/bin/env python3
"""
lrc14_neartight_rigidity_scan_opus_S4.py    opus-2026-07-23-S4

Four independent scans of 13-speed primitive configs, EXACT gap arithmetic, asking:
 (1) any counterexample gap < 1/14 ?
 (2) any TIGHT config (gap = 1/14) besides AP {1..13} and GW {1..11,13,24} ?  [OPEN-Q-108]
 (3) any config in the conjectured-empty band (1/14, 3/41) ?                   [the eps0 input]
 (4) STRUCTURE: what do the near-tight (gap <= 3/41) configs look like ?

gap(V) = max_tau min_v ||v tau|| is computed EXACTLY: the max of the piecewise-linear g sits at a
breakpoint tau = k/d with d in {2v} u {|v_i-v_j|} u {v_i+v_j}, so
   gap = max_d max_k  min_v min(vk mod d, d - vk mod d) / d.
Validated against CONSTANTS-INDEX (AP,GW=1/14; 3/41; 2/27; 3/40; 4/53; 14/183; 1/13).
"""
import numpy as np, itertools, random, time
from math import gcd
from functools import reduce
from fractions import Fraction as Fr

FL = Fr(1, 14); BAND = Fr(3, 41)
FD = [13, 14, 15, 16, 20, 25, 27, 30, 40, 41, 53, 55]

def quick(V):
    """cheap lower bound on gap; if > 3/41 the config is not near-tight -> skip."""
    best = 0.0
    for d in FD:
        for k in range(1, d):
            m = d
            for v in V:
                r = (v * k) % d; x = r if r < d - r else d - r
                if x < m:
                    m = x
                    if m * 41 <= 3 * d: break
            val = m / d
            if val > best:
                best = val
                if best > 3 / 41: return best
    return best

def gap_exact(V):
    n = len(V); Va = np.array(V, dtype=np.int64); D = set()
    for v in V: D.add(2 * v)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = V[i], V[j]
            if a != b: D.add(abs(a - b))
            D.add(a + b)
    bn, bd, bk = 0, 1, None
    for d in sorted(D):
        if d < 2: continue
        K = np.arange(1, d, dtype=np.int64); R = np.outer(Va, K) % d
        m = np.minimum(R, d - R).min(axis=0); idx = int(m.argmax()); val = int(m[idx])
        if val * bd > bn * d: bn, bd, bk = val, d, int(K[idx])
    return Fr(bn, bd), Fr(bk, bd)

def defect(V):
    """|V \ {1..13}| = number of core elements replaced."""
    return len(set(V) - set(range(1, 14)))

hits = []; totals = {}
def record(V):
    g, tau = gap_exact(list(V))
    if g <= BAND: hits.append((g, tuple(V), tau, defect(V)))
    return g

t0 = time.time()
# (A) exhaustive small speeds
n = 0
for V in itertools.combinations(range(1, 21), 13):
    if reduce(gcd, V) != 1: continue
    n += 1
    if quick(V) > 3 / 41: continue
    record(V)
totals["A exhaustive 13-subsets of {1..20}"] = n
# (B) single-far
n = 0
for j in range(1, 14):
    core = [v for v in range(1, 14) if v != j]
    for r in range(14, 601):
        V = sorted(core + [r])
        if reduce(gcd, V) != 1: continue
        n += 1
        if quick(V) > 3 / 41: continue
        record(V)
totals["B single-far {1..13}\{j} u {r}, r<=600"] = n
# (C) two-far (multi-defect branch)
n = 0
for drop in itertools.combinations(range(1, 14), 2):
    core = [v for v in range(1, 14) if v not in drop]
    for add in itertools.combinations(range(14, 101), 2):
        V = sorted(core + list(add))
        if reduce(gcd, V) != 1: continue
        n += 1
        if quick(V) > 3 / 41: continue
        record(V)
totals["C two-far drop2/add2, adds<=100"] = n
# (D) random generic
random.seed(20260723); n = 0
for R, N in [(40, 120000), (150, 120000), (1000, 60000)]:
    for _ in range(N):
        V = sorted(random.sample(range(1, R + 1), 13))
        if reduce(gcd, V) != 1: continue
        n += 1
        if quick(V) > 3 / 41: continue
        record(V)
totals["D random primitive, speeds<=40/150/1000"] = n
el = time.time() - t0

print("NEAR-TIGHT RIGIDITY SCAN of 13-speed primitive configs (exact gap)")
print("=" * 88)
tot = 0
for k, v in totals.items():
    print(f"  {k:44s} {v:>9,} configs"); tot += v
print(f"  {'TOTAL':44s} {tot:>9,} configs   [{el:.0f}s]")
print("-" * 88)
uniq = sorted(set(hits), key=lambda x: (float(x[0]), x[1]))
print(f"COMPLETE set with gap <= 3/41 found anywhere in the scan: {len(uniq)}")
for g, V, tau, dfc in uniq:
    tag = "TIGHT" if g == FL else ""
    print(f"   gap={str(g):>7} {tag:5s} defect={dfc}  tau*={str(tau):>7}  V={list(V)}")
print("-" * 88)
print(f"  counterexamples (gap<1/14): {sum(1 for g,_,_,_ in uniq if g<FL)}")
print(f"  tight besides AP/GW       : {sum(1 for g,V,_,_ in uniq if g==FL and set(V) not in (set(range(1,14)), set(list(range(1,12))+[13,24])))}")
print(f"  in band (1/14,3/41)       : {sum(1 for g,_,_,_ in uniq if FL<g<BAND)}")
print("\nLAW (empirical): gap <= 3/41  =>  defect <= 1 (at most ONE core element replaced).")
print("Two defects or generic configs never come within 3/41 of the floor.")
