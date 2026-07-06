#!/usr/bin/env python3
"""
klein-2026-07-05-S141b (HYP-4161 cont.) - THE SECOND-VALUE INDUCTION: verify the
architecture [gap+uniq](n-1) + AP-base(n) + drag(n) => [gap+uniq](n), and probe the
MIDDLE STRIP that neither the midrange witness nor the one-window drag closes.

Key numerology (to verify):
  rho_n = 2/(2n+1) (target = level-n second value), beta_n = 2/(2n-1) (induction floor
  = level-(n-1) second value). One-window drag threshold = rho/(beta-rho) = (2n-1)/2.
  Midrange closes ratio v_max/v_min <= (2n-1)/2. THE SAME CONSTANT both sides.
  => remaining strip at level n: v_max/v_min > (2n-1)/2 AND v_max/secondmax <= (2n-1)/2
     (spread overall but compressed at top), all-loose.
Probe: enumerate structured 12-family middle-strip members, exact M, min over the strip.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

def dist(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def exact_M(V):
    qs = set()
    Vl = sorted(V)
    for i, u in enumerate(Vl):
        qs.add(2*u)
        for v in Vl[i+1:]:
            qs.add(u+v); qs.add(v-u)
    qs.discard(0)
    best = F(0)
    for q in qs:
        for a in range(1, q):
            if gcd(a, q) > 1: continue
            m = min(dist(x*F(a,q)) for x in Vl)
            if m > best: best = m
    return best

n = 12
rho = F(2, 2*n+1)    # 2/25
beta = F(2, 2*n-1)   # 2/23
thr = rho/(beta-rho)
print(f"n={n}: rho={rho}, beta={beta}, drag threshold rho/(beta-rho) = {thr} = {float(thr)}")
print(f"midrange boundary: m/(m+Mx) >= {rho} <=> Mx <= {(1-rho)/rho} m -> {float((1-rho)/rho)}")
print(f"MEETING POINT check: (2n-1)/2 = {F(2*n-1,2)} -- both = 11.5: {thr == F(2*n-1,2) - 0 and (1-rho)/rho == F(2*n-1,2)}")

# middle strip probe: v_max/v_min > 11.5, v_max/secondmax <= 11.5
print("\nmiddle strip probe (structured + random, 12-families, primitive):")
random.seed(9)
worst = None; count = 0
cands = []
# structured: small base + two-or-more large near-equal tops
for lo_top in range(2, 8):
    base = list(range(1, 13 - lo_top + 1))[:12 - lo_top]
    for T in range(int(11.5*max(base))+2, 40*max(base), max(base)):
        tops = [T + j for j in range(lo_top)]
        V = base + tops
        if len(V) != 12: continue
        cands.append(V)
for _ in range(400):
    k = random.randint(2, 5)
    base = sorted(random.sample(range(1, 30), 12 - k))
    T = random.randint(int(11.5*min(base))+1, 60*min(base))
    tops = sorted(random.sample(range(T, T + int(10.5*T)), k))
    V = sorted(set(base + tops))
    if len(V) == 12: cands.append(V)
for V in cands:
    V = sorted(V)
    if gcd(*V) != 1: continue
    ratio_full = F(V[-1], V[0])
    ratio_top = F(V[-1], V[-2])
    if not (ratio_full > F(23,2) and ratio_top <= F(23,2)): continue
    count += 1
    M = exact_M(V)
    if worst is None or M < worst[0]: worst = (M, V)
print(f"strip members tested: {count}")
if worst:
    print(f"strip minimum M = {worst[0]} = {float(worst[0]):.6f} at {worst[1]}")
    print(f"vs 2/25 = {float(F(2,25)):.6f}: {'ALL CLEAR >= 2/25' if worst[0] >= F(2,25) else 'VIOLATION'}")
