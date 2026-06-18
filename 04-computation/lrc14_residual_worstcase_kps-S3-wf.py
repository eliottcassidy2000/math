#!/usr/bin/env python3
"""
LRC(14) S3 -- pin the WORST global-witness margin exactly, and test a clean lemma.

We minimize G*·7·Vmax (the global-witness margin) over a large adversarial S3 sample
including extreme spreads / cluster sizes / large Vmax, to find the true infimum and
whether it stays > 1. Also we test the cleanest candidate provable lemma:

CANDIDATE LEMMA (single-cluster-runner suffices for a global witness):
  For S3, consider the two-element reduction P ∪ {Vmin}. The small part P safe set has
  arcs; the SINGLE cluster runner Vmin plants teeth of width 1/(7 Vmin) spaced 1/Vmin.
  Inside a small safe arc I_P of width W_P, the largest Vmin-free gap is
     >= (a clean function of W_P, Vmin).
  THEN add the remaining cluster runners: they lie within spread s of Vmin, so on a
  window of width w their phases differ from Vmin's by <= s*w/?... We test whether the
  global witness can be CONSTRUCTED from P ∪ {Vmin} plus a spread correction.

Mainly: get the exact infimum margin, and the exact worst S3 sets, to either (a) prove
a clean bound or (b) honestly bound the residual to a finite check (Vmax-bounded).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
from collections import Counter

def merge(ivs):
    ivs = sorted(ivs); m = []
    for a, b in ivs:
        if m and a <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], b))
        else: m.append((a, b))
    return m

def teeth_in(u, lo, hi, h=F(1, 14)):
    out = []
    jmin = int(lo * u) - 1; jmax = int(hi * u) + 1
    for j in range(jmin, jmax + 1):
        c = F(j, u); a = c - h/u; b = c + h/u
        aa = max(a, lo); bb = min(b, hi)
        if aa < bb: out.append((aa, bb))
    return out

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

def largest_free_gap(lo, hi, cluster):
    dgr = []
    for u in cluster: dgr.extend(teeth_in(u, lo, hi))
    dgr = merge(dgr)
    gaps = []; prev = lo
    for a, b in dgr:
        if a > prev: gaps.append((a - prev, prev, a))
        prev = max(prev, b)
    if hi > prev: gaps.append((hi - prev, prev, hi))
    if not gaps: return (F(0), None)
    return max(gaps)[0], max(gaps)

def global_witness_margin(S):
    small = [u for u in S if u <= 13]; cluster = [u for u in S if u > 13]
    Vmax = max(cluster)
    sc = safe_components(small)
    best = F(0); arc = None
    for (a, b) in sc:
        g, gg = largest_free_gap(a, b, cluster)
        if g > best: best = g; arc = (a, b, gg)
    return best * 7 * Vmax, best, arc

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

def gen_S3_adv(seed=0, target=2000):
    """Adversarial: bias toward small spreads (tight cluster) + large Vmax + c=2,3
       which empirically minimize the margin; also include drop-1 small parts
       (the worst small part is {1..13}\\{j})."""
    rng = random.Random(seed); out = []; tries = 0
    while len(out) < target and tries < target * 800:
        tries += 1
        c = rng.choice([2, 2, 2, 3, 3, 4]); nsmall = 13 - c
        if nsmall < 1: continue
        V = rng.choice([rng.randint(50, 100), rng.randint(100, 400),
                        rng.randint(400, 1500), rng.randint(1500, 6000)])
        spread = rng.choice([2, 3, 4, 5, 7, 10, 14, 20, 28, 40])
        cluster = set()
        while len(cluster) < c: cluster.add(V + rng.randint(0, spread))
        # bias small part toward drop-1 from {1..13}
        pool = list(range(1, 14)); rng.shuffle(pool)
        small = set(pool[:nsmall])
        S = sorted(small | cluster)
        if len(S) != 13 or gcd_all(S) != 1: continue
        if not is_covering(S) or classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    S3 = gen_S3_adv(seed=21, target=2500)
    print(f"n adversarial S3 = {len(S3)}")
    margins = []
    worst = []
    for S in S3:
        m, g, arc = global_witness_margin(S)
        margins.append((m, S))
    margins.sort(key=lambda x: x[0])
    print(f"min global-witness margin G*·7·Vmax = {float(margins[0][0]):.4f}")
    print(f"  (all > 1 ?  {'YES' if margins[0][0] > 1 else 'NO -- FAILURE'})")
    print(f"mean = {float(sum(m for m,_ in margins)/len(margins)):.3f}")
    print("\n10 worst (smallest-margin) S3 sets:")
    for m, S in margins[:10]:
        cl = [u for u in S if u > 13]; sm = [u for u in S if u <= 13]
        print(f"   margin={float(m):.4f}  small={sm} cluster={cl} (c={len(cl)},s={max(cl)-min(cl)})")
