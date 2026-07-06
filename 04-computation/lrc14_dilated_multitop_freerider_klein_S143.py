#!/usr/bin/env python3
"""
klein-2026-07-05-S143 (HYP-4192) - THE DILATED MULTI-TOP FREE-RIDER BY COUNTING.

Base c*{1..k} (c >= 2), tops w_1..w_l (l = 12-k), target margin 2/25 at t = a/Q,
Q = (k+1)c. Base-safety: a must be a unit mod (k+1) (then base dist >= 1/(k+1) > 2/25
for k <= 11 -- wait: base runner cj at a/Q: dist(cj a/((k+1)c)) = dist(ja/(k+1)) >=
1/(k+1); need >= 2/25: 1/(k+1) >= 2/25 <=> k <= 11.5 OK all k <= 11).
Top-safety: R_j = (w_j * a) mod Q in [ceil(2Q/25), floor(23Q/25)].

H2' (the counting hypothesis): the unsafe a's per top are ~(4/25) of the steering set
{a in [1,Q): gcd(a, k+1) = 1}; union over l tops < 1 for l <= 6 => a good a EXISTS.
Test: over dilated-AP bases (c = 2..12, k = 6..11) with random valid tops (nonzero,
not making the family non-primitive... just primitive families), count:
  - families where NO good a exists (counterexamples to the free-rider),
  - the min good-a fraction (vs the predicted 1 - 4l/25).
"""
from math import gcd
from fractions import Fraction as F
import random

random.seed(23)
worst_frac = 1.0
fails = []
tested = 0
for trial in range(4000):
    k = random.randint(6, 11)
    l = 12 - k
    c = random.randint(2, 12)
    Q = (k+1)*c
    base = [c*j for j in range(1, k+1)]
    # tops: mid-scale, avoid (k+1)-divisibility only in the c=1 sense? test raw first
    tops = []
    while len(tops) < l:
        w = random.randint(int(11.5*c)+1, 40*Q)
        tops.append(w)
    V = base + tops
    if gcd(*V) != 1: continue
    tested += 1
    lo = -(-2*Q // 25)      # ceil(2Q/25)
    hi = (23*Q) // 25       # floor(23Q/25)
    good = 0; total = 0
    for a in range(1, Q):
        if gcd(a, k+1) != 1: continue
        total += 1
        ok = True
        for w in tops:
            R = (w*a) % Q
            if not (lo <= R <= hi):
                ok = False; break
        if ok: good += 1
    if total == 0: continue
    frac = good/total
    if frac < worst_frac: worst_frac = frac
    if good == 0:
        fails.append((k, c, tops))
print(f"tested {tested} dilated-AP-base families (c=2..12, l=1..6 tops)")
print(f"families with NO good steering a: {len(fails)}")
for k, c, tops in fails[:8]:
    print(f"  FAIL k={k} c={c} Q={(k+1)*c} tops={tops}")
print(f"worst good-a fraction: {worst_frac:.3f} (union-bound prediction >= {1-4*6/25:.3f} for l=6)")
