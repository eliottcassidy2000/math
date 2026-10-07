#!/usr/bin/env python3
"""Audit A, item 6(e): HYP-9214's sampler event meet(n, (n-1)/2) (equal U-index, stopped at 1) vs absorption of the
pair chain from (0, 1) on (u, v) = (n, n - 1) driven by v's Terras parities. Also: the first-reset-2 condition."""
import random
from a3_drift import step, v2
def U(x):
    y = 3 * x + 1
    return y >> v2(y)
def T(x): return x >> 1 if x % 2 == 0 else (3 * x + 1) >> 1
def meet(n, m):   # verbatim logic of sevens_20261006_debt_resolution_trend.py
    j = 0
    while n != 1 and m != 1:
        n, m, j = U(n), U(m), j + 1
        if n == m:
            return j if n != 1 else None
    return None
def chain(n):
    k, N = 0, 1; u, v = n, n - 1; t = 0
    while True:
        if (k, N) == (0, 0): return t, u
        if u <= 2 or v <= 2: return None, None
        k, N = step(k, N, v & 1); u, v = T(u), T(v); t += 1

# (1) first-reset-2 condition: n's reset exponent == 2  <=>  v2(3^r t - 1) >= 2  (n = 2^(r+1) t - 1, r >= 1, t odd)
bad = 0
for n in range(3, 200001, 4):
    r = v2(n + 1) - 1; t = (n + 1) >> (r + 1)
    x = n
    while v2(3 * x + 1) == 1: x = U(x)
    reset = v2(3 * x + 1)
    if (reset == 2) != (v2(3 ** r * t - 1) >= 2): bad += 1
print(f"[reset-2 condition] odd n = 3 mod 4 below 2e5: mismatches {bad}")
# (2) sampler event vs chain absorption
rng = random.Random(9214)
for B in (32, 64, 256, 1024):
    tot = hit = ab = both = ab_only = hit_only = 0; Kmax_ok = 0
    while tot < (3000 if B <= 256 else 800):
        n = rng.getrandbits(B) | (1 << (B - 1)) | 1
        r = v2(n + 1) - 1
        if r < 1: continue
        t = (n + 1) >> (r + 1)
        if 1 + v2(3 ** r * t - 1) < 3: continue
        tot += 1
        j = meet(n, (n - 1) // 2)
        ta, z = chain(n)
        hit += j is not None; ab += ta is not None
        if j is not None and ta is not None: both += 1
        elif ta is not None: ab_only += 1
        elif j is not None: hit_only += 1
    print(f"[B={B:5d}] sources {tot}: sampler merge {hit}, chain absorbed (values > 2) {ab}, both {both}, absorbed-but-no-sampler-merge {ab_only}, sampler-merge-but-not-absorbed {hit_only}")
