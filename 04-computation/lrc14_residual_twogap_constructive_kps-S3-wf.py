#!/usr/bin/env python3
r"""
LRC(14) S3 -- the actual "two-gap-lemma2plus": a CONSTRUCTIVE r-gap band-fit lemma,
and how much rigorous coverage it adds over the single-gap Lemma A.

r-GAP LEMMA (constructive multi-band-fit). Fix integers K_1 < K_2 < ... and a split of
the cluster L into contiguous sub-bands by speed: L_1=[Vmin..b_1], L_2=[b_1+1..b_2], ...
We seek a common t such that:
   (i) for each sub-band L_i = [lo_i, hi_i] (speed range) and chosen index K_i:
          K_i + 1/14 <= lo_i * t   and   hi_i * t <= K_i + 13/14,
       i.e. t in J_i := ( (K_i+1/14)/lo_i , (K_i+13/14)/hi_i ),  AND
   (ii) the small part P is safe at t.
If the intersection of all J_i over i, intersected with a small-safe arc, is nonempty,
then t is a global witness => M(S) >= 1/14.  (each J_i is a per-subband cluster window.)

This is exactly Lemma A iterated per sub-band. The freedom is the SPLIT and the indices
K_i. We search over splits (by sorting L and cutting into r contiguous pieces) and over
index assignments forced by t (K_i = floor(lo_i * t)). To make it constructive WITHOUT
knowing t: for r=2, sweep the cut point and a candidate base index, solve the two-window
intersection in closed form, check small-safe.

We measure: rigorous coverage of (single-gap) vs (single + two-gap) vs (single+two+three-gap),
to quantify how much the multi-band generalization buys in PROVABLE terms.
"""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
from itertools import combinations
import random

def merge(ivs):
    ivs = sorted(ivs); m = []
    for a, b in ivs:
        if m and a <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], b))
        else: m.append((a, b))
    return m
def small_safe_arcs(P, h=F(1, 14)):
    iv = []
    for u in P:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    m = merge(iv); safe = []; prev = F(0)
    for a, b in m:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def subband_window(lo_speed, hi_speed, K):
    """J = ((K+1/14)/lo_speed, (K+13/14)/hi_speed). Safe for all speeds in [lo,hi]
    at gap index K. (For a sub-band, lo_speed=min, hi_speed=max of the sub-band.)"""
    lo = F(14*K+1, 14)/lo_speed; hi = F(14*K+13, 14)/hi_speed
    return (lo, hi) if lo < hi else None

def intersect(a, b):
    if a is None or b is None: return None
    lo = max(a[0], b[0]); hi = min(a[1], b[1])
    return (lo, hi) if lo < hi else None

def small_safe_in(J, safe):
    if J is None: return False
    for (a, b) in safe:
        if max(a, J[0]) < min(b, J[1]): return True
    return False

def rgap_closable(S, rmax=1):
    """Try r-gap band fit for r=1..rmax. Returns True if any closes (rigorous proof)."""
    P = [u for u in S if u <= 13]; L = sorted(u for u in S if u > 13)
    safe = small_safe_arcs(P); c = len(L); Vmax = L[-1]
    # candidate base index range: K up to about 13*Vmin/Vmax windows; cap
    Kcap = 14  # gap indices we'll try as the LOWEST sub-band index; t<1 keeps it small
    # r=1
    Vmin = L[0]
    for K in range(0, Kcap+1):
        J = subband_window(Vmin, Vmax, K)
        if small_safe_in(J, safe): return True
    if rmax >= 2:
        # r=2: cut L into [L[0..i]] and [L[i+1..]] for each i; pick base index K0 for
        # band1; band2 index can be K0, K0+1, ..., K0+span. We let t determine K via
        # the constraint; simplest: try band1 index k1, band2 index k2 with k2>=k1, and
        # intersect windows.
        for i in range(0, c-1):
            b1 = L[:i+1]; b2 = L[i+1:]
            lo1, hi1 = b1[0], b1[-1]; lo2, hi2 = b2[0], b2[-1]
            for k1 in range(0, Kcap+1):
                for k2 in range(k1, k1+ Kcap+1):
                    J1 = subband_window(lo1, hi1, k1)
                    J2 = subband_window(lo2, hi2, k2)
                    J = intersect(J1, J2)
                    if J is None: continue
                    if J[0] >= 1: continue
                    if small_safe_in(J, safe): return True
                    if J2 is not None and J2[0] >= 1: break
    if rmax >= 3:
        for i in range(0, c-2):
            for j in range(i+1, c-1):
                bands = [L[:i+1], L[i+1:j+1], L[j+1:]]
                rng_b = [(bb[0], bb[-1]) for bb in bands]
                for k1 in range(0, 10):
                    for k2 in range(k1, k1+10):
                        for k3 in range(k2, k2+10):
                            Js = [subband_window(rng_b[t][0], rng_b[t][1], kk)
                                  for t, kk in enumerate([k1,k2,k3])]
                            J = Js[0]
                            for x in Js[1:]: J = intersect(J, x)
                            if J is None or J[0] >= 1: continue
                            if small_safe_in(J, safe): return True
    return False

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'
def missing_q(P): return [q for q in range(2, 15) if not any(v % q == 0 for v in P)]
def gen_constructive(seed=0, target=2000, Vrange=(50, 2000)):
    rng = random.Random(seed); out = []; tries = 0
    smalls = []
    base = list(range(1, 14))
    for sz in (11, 10, 9, 8, 7):
        for P in combinations(base, sz): smalls.append(list(P))
    while len(out) < target and tries < target * 300:
        tries += 1
        P = rng.choice(smalls); c = 13 - len(P)
        if c < 2: continue
        miss = missing_q(P); V = rng.randint(*Vrange)
        spread = rng.choice([2,3,5,7,10,14,20,28,40,56,70,90])
        window = list(range(V, V + spread + 1))
        if len(window) < c: continue
        cluster = rng.sample(window, c)
        if not all(any(v % q == 0 for v in cluster) for q in miss): continue
        S = sorted(set(P) | set(cluster))
        if len(S) != 13 or gcd_all(S) != 1: continue
        if not is_covering(S) or classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    allS = []
    for (lo, hi) in [(50, 200), (200, 800), (800, 3000)]:
        allS += gen_constructive(seed=lo+5, target=1000, Vrange=(lo, hi))
    n = len(allS)
    print(f"n S3 = {n}")
    c1 = sum(1 for S in allS if rgap_closable(S, rmax=1))
    print(f"r=1 (single-gap Lemma A):        {c1}/{n} = {100*c1/n:.1f}%")
    c2 = sum(1 for S in allS if rgap_closable(S, rmax=2))
    print(f"r<=2 (single + TWO-gap):         {c2}/{n} = {100*c2/n:.1f}%")
    c3 = sum(1 for S in allS if rgap_closable(S, rmax=3))
    print(f"r<=3 (single + two + three-gap): {c3}/{n} = {100*c3/n:.1f}%")
    print(f"\nGain from two-gap over single: +{100*(c2-c1)/n:.1f}%")
    print(f"Gain from three-gap over two:  +{100*(c3-c2)/n:.1f}%")
