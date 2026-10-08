#!/usr/bin/env python3
"""The +1-barrier child as an integer process.  c = 3^m (v + 1) with v the child's Terras orbit and m the remaining
3-adic denominator exponent: c even -> c/2, m -> m-1; c odd -> (c + 3^m)/2, m unchanged.  Start c = 2, m = D
(v_0 = 2/3^D - 1); landing at m = 0 gives N_D = c - 1.  (Check against the direct rational orbit.)
L(c, m) := landing value c at m = 0 from state (c, m).
(1) exactness check; (2) distribution of L(c, m) for all c in [1, 2*3^m] (r = c/3^m in (0, 2]), m <= 13;
(3) i.i.d.-parity model: r -> (r+1)/2 or 3r/2 (fair coin), stop after m '3r/2' steps; tail of the landing r."""
import random, math
from collections import Counter
def land_direct(D):
    a, m = 2 - 3**D, D
    while m > 0:
        if a & 1: a = (a + 3**(m-1)) // 2; m -= 1
        else: a //= 2
    return a
def L(c, m):
    # returns None if the child hits v = 0 (c = 3^m, the fixed point 0 of T) before landing
    while m > 0:
        p = 3**m
        if c == p: return None
        if c & 1: c = (c + p) // 2
        else: c //= 2; m -= 1
    return c
assert all(L(2, D) - 1 == land_direct(D) for D in range(3, 400))
print("(1) c-process reproduces N_D = L(2, D) - 1 for D < 400: OK")
print("(2) exhaustive landing values L(c, m), c in [1, 2*3^m]:")
for m in range(2, 12):
    allv = [L(c, m) for c in range(1, 2 * 3**m + 1)]
    vals = [v for v in allv if v is not None]; zero = len(allv) - len(vals)
    n = len(vals); s = sorted(vals)
    tail = {x: sum(1 for v in vals if v > x) / n for x in (4, 16, 64, 256)}
    print(f"   m={m:2d}: max {s[-1]:7d}  median {s[n//2]:3d}  stuck at 0: {zero}  P(>4,16,64,256) = " + ", ".join(f"{tail[x]:.4f}" for x in (4, 16, 64, 256)), flush=True)
rnd = random.Random(5)
def iid(m):
    r = 2 / 3**m * 3**m  # start r = c/3^m = 2/3^m * ... use generic start r=1
    r = 1.0; k = m
    while k > 0:
        if rnd.random() < 0.5: r = (r + 1) / 2
        else: r = 1.5 * r; k -= 1
    return r
print("(3) i.i.d.-parity model, landing r after m expansions (20000 runs each):")
for m in (6, 10, 13, 20, 40):
    vals = [iid(m) for _ in range(20000)]
    tail = {x: sum(1 for v in vals if v > x) / len(vals) for x in (4, 16, 64, 256)}
    print(f"   m={m:2d}: max {max(vals):9.1f}  P(>4,16,64,256) = " + ", ".join(f"{tail[x]:.4f}" for x in (4, 16, 64, 256)))
