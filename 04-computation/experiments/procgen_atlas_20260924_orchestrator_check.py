#!/usr/bin/env python3
"""Orchestrator's independent audit of Theorem M' (THM-4469), written without reading the lane's code.

Checks, each raising on failure:
  O1  R_B = 4726, R_B' = 4727 for B = 0111101110, B' = 1101100111; p = 3^7 > 2q = 2^11.
  O2  cylinder fact: {x mod 2^10 : first 10 T-parities of x = B} = {x : 2^10 | 3^7 x + R_B}, stable on 49 lifts.
  O3  census of adjacent pairs (R_B' = R_B + 1, equal weight a, 3^a > 2^(L+1)) for L <= 20:
      exactly (10,7):4, (16,11):5, (19,13):1, (20,14):8.
  O4  T-side counts on [1, 2^22]: starts whose first d blocks lie in {B,B'}: 8192, 16, 0 for d = 1, 2, 3.
  O5  xi-side, exact rationals: for 56 sampled integer parts x0 (40 random one-block starts plus all 16
      two-block starts), "some t0 in I with (x0+t0) alpha^j in Z+I for j <= 2" holds iff x0 has two blocks;
      and 200 random non-cylinder x0 fail already at j = 1.
Run: python3 04-computation/experiments/procgen_atlas_20260924_orchestrator_check.py
"""
import math, random
from fractions import Fraction as F
import numpy as np

def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
def R_of(w):
    R = 0
    for i, ch in enumerate(w):
        e = int(ch); R = (3 ** e) * R + e * 2 ** i
    return R
def word(x, L):
    s = ''
    for _ in range(L): s += str(x % 2); x = T(x)
    return s

B, Bp = '0111101110', '1101100111'
p, q = 3 ** 7, 2 ** 10
# O1
assert (R_of(B), R_of(Bp), p - q) == (4726, 4727, 1163) and p > 2 * q
print("O1  R_B = 4726, R_B' = 4727, p - q = 1163, p = 2187 > 2q = 2048: ok")
# O2
for blk in (B, Bp):
    cyl = [x for x in range(q) if word(x, 10) == blk]
    div = [x for x in range(q) if (p * x + R_of(blk)) % q == 0]
    assert cyl == div and len(cyl) == 1
    assert all(word(x + q * t, 10) == blk for x in cyl for t in range(1, 50))
    print(f"O2  cylinder of {blk} = class {cyl[0]} mod 1024 = divisibility class: ok")
# O3
R = np.array([0], dtype=np.int64); W = np.array([0], dtype=np.int64); found = {}
for L in range(1, 21):
    R = np.concatenate([R, 3 * R + 2 ** (L - 1)]); W = np.concatenate([W, W + 1])
    for a in range(1, L):
        if 3 ** a > 2 ** (L + 1):
            v = np.sort(R[W == a]); k = int((np.diff(v) == 1).sum())
            if k: found[(L, a)] = k
assert found == {(10, 7): 4, (16, 11): 5, (19, 13): 1, (20, 14): 8}, found
print("O3  adjacent-pair census L <= 20:", found, ": ok")
# O4
Bset = {B, Bp}; N = 2 ** 22
cres = {x for x in range(q) if word(x, 10) in Bset}
d1 = [x for x in range(1, N + 1) if x % q in cres]
def blocks_ok(x, d):
    for _ in range(d):
        if word(x, 10) not in Bset: return False
        for _ in range(10): x = T(x)
    return True
d2 = [x for x in d1 if blocks_ok(x, 2)]; d3 = [x for x in d2 if blocks_ok(x, 3)]
assert (len(d1), len(d2), len(d3)) == (8192, 16, 0)
print("O4  T-side counts d = 1, 2, 3 on [1, 2^22]: 8192, 16, 0: ok; two-block starts:", d2)
# O5
alpha = F(p, q); I = (F(4726, 1163), F(4727, 1163))
def xi_ok(x0, d):
    ivs = [I]
    for j in range(1, d + 1):
        aj = alpha ** j; new = []
        for lo, hi in ivs:
            ylo, yhi = (x0 + lo) * aj, (x0 + hi) * aj
            for n in range(math.floor(ylo - I[1]) - 1, math.ceil(yhi - I[0]) + 2):
                a_, b_ = max(ylo, n + I[0]), min(yhi, n + I[1])
                if a_ <= b_: new.append((a_ / aj - x0, b_ / aj - x0))
        ivs = new
        if not ivs: return False
    return True
random.seed(1)
samp = random.sample(d1, 40) + d2
mism = [x for x in samp if xi_ok(x, 2) != (x in d2)]
assert not mism, mism
d1s = set(d1); nonc = [x for x in random.sample(range(1, N), 200) if x not in d1s]
assert all(not xi_ok(x, 1) for x in nonc)
print(f"O5  xi-side = T-side at depth 2 on {len(samp)} samples; {len(nonc)} non-cylinder starts fail at j = 1: ok")
print("ALL CHECKS PASSED")
