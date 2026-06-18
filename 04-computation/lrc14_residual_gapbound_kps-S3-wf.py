#!/usr/bin/env python3
"""
LRC(14) S3 -- closed-form lower bound on the largest cluster-free gap inside I_P.

We established: in EVERY tested S3 set there is a GLOBAL witness -- a point in the
small-part's widest safe arc I_P that is also free of all cluster teeth. Margin
G*·7·Vmax >= 1.5.

Now make step 2 provable. Inside I_P (width W_P), the cluster L = {u_1<...<u_c} (all
near V, spread s) plants teeth. KEY tightness fact: because the cluster is tight, in a
short window the c runners' phases u_i*tau are nearly equal, so their teeth nearly
coincide and the merged cluster-danger looks like ONE comb of period ~1/V with width
~1/(7V) PER cluster passage, NOT c independent combs.

PRECISE STRUCTURE TO MEASURE per S3 set (restricted to I_P):
  - W_P              = |I_P| (width of small-part widest safe arc)
  - n_teeth          = # merged cluster danger components inside I_P
  - avg_gap          = (cluster-free measure in I_P)/(n_teeth+1)
  - G*               = largest cluster-free gap in I_P
  - ratio_widest_avg = G* / avg_gap
  - free_frac        = (cluster-free measure in I_P)/W_P
  - the "effective comb period" = W_P / n_teeth   vs 1/V and vs 1/Vmax

Goal: a provable inequality of form G* >= W_P / (K * Vmax * W_P + const) etc that beats
1/(7Vmax). Equivalently: n_teeth <= 7 * Vmax * W_P - 1  would give avg_gap >= W_P/(7 Vmax W_P)=1/(7Vmax)
... but we need the GAP not avg. Let's measure n_teeth vs (Vmax * W_P) precisely.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
from collections import Counter

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def danger_intervals(A, lo, hi, h=F(1, 14)):
    out = []
    for u in A:
        jmin = int(lo * u) - 1; jmax = int(hi * u) + 1
        for j in range(jmin, jmax + 1):
            c = F(j, u); a = c - h/u; b = c + h/u
            aa = max(a, lo); bb = min(b, hi)
            if aa < bb: out.append((aa, bb))
    return out

def merge(ivs):
    ivs = sorted(ivs); m = []
    for a, b in ivs:
        if m and a <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], b))
        else: m.append((a, b))
    return m

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def analyze_IP(lo, hi, cluster):
    dgr = merge(danger_intervals(cluster, lo, hi))
    # free gaps
    gaps = []; prev = lo
    for a, b in dgr:
        if a > prev: gaps.append(a - prev)
        prev = max(prev, b)
    if hi > prev: gaps.append(hi - prev)
    free_meas = sum(gaps)
    n_teeth = len(dgr)
    Gstar = max(gaps) if gaps else (hi - lo)
    return Gstar, n_teeth, free_meas, len(gaps)

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

def gen_S3(seed=0, target=2000, vmax_cap=2000):
    rng = random.Random(seed); out = []; tries = 0
    while len(out) < target and tries < target * 400:
        tries += 1
        c = rng.choice([2, 2, 3, 3, 4, 5, 6]); nsmall = 13 - c
        if nsmall < 1: continue
        V = rng.choice([rng.randint(40, 80), rng.randint(80, 200),
                        rng.randint(200, 600), rng.randint(600, vmax_cap)])
        spread = rng.choice([14, 20, 28, 35, 42, 56, 70])
        cluster = set()
        while len(cluster) < c: cluster.add(V + rng.randint(0, spread))
        small = set(); pool = list(range(1, 14)); rng.shuffle(pool)
        for x in pool:
            if len(small) >= nsmall: break
            small.add(x)
        if len(small) < nsmall: continue
        S = sorted(small | cluster)
        if len(S) != 13 or gcd_all(S) != 1: continue
        if not is_covering(S) or classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    S3 = gen_S3(seed=7, target=600)
    print(f"n S3 = {len(S3)}")
    print("Per S3: in the SMALL-PART widest arc I_P, cluster-free-gap structure.")
    print("Columns: c=|L|  s=spread  W_P*Vmax  n_teeth  n_teeth/(c*W_P*Vmax)  ratio_widest/avg  G*·7·Vmax")
    nteeth_over_cWV = []; ratio_wa = []; Gmargin = []; freefrac = []
    worst = (F(10), None)
    for S in S3:
        Vmax = S[-1]; small = [u for u in S if u <= 13]; cluster = [u for u in S if u > 13]
        c = len(cluster); s = max(cluster) - min(cluster)
        sc = safe_components(small)
        # widest small arc
        wbest = F(0); arc = None
        for (a, b) in sc:
            if b - a > wbest: wbest = b - a; arc = (a, b)
        # also consider all small arcs for the global witness, but for the bound study use the widest
        # Actually for a clean bound, evaluate over ALL small arcs and take best gap:
        bestG = F(0); best_stats = None
        for (a, b) in sc:
            G, nt, fm, ng = analyze_IP(a, b, cluster)
            if G > bestG:
                bestG = G; best_stats = (b - a, nt, fm, ng)
        WP, nt, fm, ng = best_stats
        nteeth_over_cWV.append(float(nt) / float(c * WP * Vmax) if WP*Vmax>0 else 0)
        avg = fm / ng if ng else fm
        ratio_wa.append(float(bestG / avg) if avg > 0 else 0)
        Gmargin.append(float(bestG * 7 * Vmax))
        freefrac.append(float(fm / WP) if WP > 0 else 0)
        if bestG * 7 * Vmax < worst[0]:
            worst = (bestG * 7 * Vmax, (S, c, s, float(WP*Vmax), nt))
    def st(a): return f"min={min(a):.3f} mean={sum(a)/len(a):.3f} max={max(a):.3f}"
    print(f"\nn_teeth/(c*W_P*Vmax)   : {st(nteeth_over_cWV)}")
    print(f"ratio widest/avg gap   : {st(ratio_wa)}")
    print(f"free fraction in I_P   : {st(freefrac)}")
    print(f"G*·7·Vmax  (the margin): {st(Gmargin)}")
    print(f"\nworst-margin set: G*·7·Vmax={float(worst[0]):.3f}")
    S,c,s,WPV,nt = worst[1]
    print(f"   S={S}  c={c} s={s} W_P*Vmax={WPV:.2f} n_teeth={nt}")
