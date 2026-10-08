#!/usr/bin/env python3
"""The contracting translation-only rank-3 map on Z_5 with multipliers (1, 1, 6, 11, 16) at residues 0..4:
T(x) = x/5, (x + 4)/5, (6x + 3)/5, (11x + 2)/5, (16x + 1)/5   [r_i = -m_i i mod 5].
Lambda = (1/5) ln(1056) - ln 5 = -0.217 (contracting).  THM-4609 predicts a transient debt walk: Haar y and y+e
merge with probability < 1.  (a) chain estimate of q(T) from (1, 1); (b) actual integers: n uniform 2000-digit (base 5)
and n+1, merge within T steps; (c) orbits of small integers all reach cycles (contraction), yet consecutive integers
rarely merge at equal time."""
import random, math, sys
from fractions import Fraction as Fr
m = [1, 1, 6, 11, 16]; d = 5
r = [(-m[i] * i) % d for i in range(d)]
for i in range(d): assert (m[i] * i + r[i]) % d == 0
def T(x):
    i = x % d; return (m[i] * x + r[i]) // d
print("map:", [f"({m[i]}x+{r[i]})/5" for i in range(d)], " Lambda =", round(sum(math.log(x / d) for x in m) / d, 4))
rnd = random.Random(6)
# (b) actual integers
for T_ in (100, 1000, 4000):
    merged = 0; N = 400
    for _ in range(N):
        n = rnd.randrange(5 ** 2000, 2 * 5 ** 2000); a, b = n, n + 1
        for t in range(T_):
            a, b = T(a), T(b)
            if a == b: merged += 1; break
    print(f"integers (2000 base-5 digits): fraction of n merging with n+1 within {T_} steps: {merged/N:.3f}", flush=True)
# (c) small integers: do orbits reach cycles?  find the cycles reached from n <= 20000
cycles = {}
for n in range(1, 20001):
    x = n; seen = {}
    for t in range(5000):
        if x in seen: break
        seen[x] = t; x = T(x)
    if x in seen:
        cyc = []; y = x
        while True:
            cyc.append(y); y = T(y)
            if y == x: break
        key = min(cyc); cycles[key] = cycles.get(key, 0) + 1
print("cycles reached from n <= 20000 (min element: count):", dict(sorted(cycles.items())[:12]), "...", len(cycles), "cycles")
same = 0; tot = 0
for n in range(1, 20000):
    a, b = n, n + 1; ok = False
    for t in range(3000):
        a, b = T(a), T(b)
        if a == b: ok = True; break
    same += ok; tot += 1
print(f"n <= 20000: fraction with T^t(n) = T^t(n+1) for some t <= 3000: {same/tot:.3f}")
