#!/usr/bin/env python3
"""kps-S4: for |L|=3,4 covering S3 sets, does the via-max criterion W(S\\max)*7*Vmax>1 hold
above a Vmax threshold (extending the k=2 Vmax>=63 result)? Find the failure region (finite core)
and whether large-Vmax always passes. Also record M directly (never <1/14). Exact."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random, sys
random.seed(2)
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe
def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1: ws.append(sc[0][1] + (1 - sc[-1][0]))
    return max(ws)
def is_cov(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def klarge(S): return sum(1 for v in S if v > 13)

# build |L|=3 and |L|=4 covering S3 sets: drop 3 or 4 from {1..13}, add large that complete covering
smalls = list(range(1, 14))
for kk in [3, 4]:
    print("=== |L| = %d ===" % kk); sys.stdout.flush()
    fails = []  # criterion-via-max fails
    maxfail_vmax = 0; passcount = 0; total = 0
    Dlist = list(combinations(smalls, kk)); random.shuffle(Dlist)
    for D in Dlist[:40]:
        P = [s for s in smalls if s not in D]
        for _ in range(700):
            CAP = random.choice([60, 100, 160, 260, 400])
            L = random.sample(range(14, CAP+1), kk)
            S = sorted(P + L)
            if len(S) != 13 or reduce(gcd, S) != 1 or not is_cov(S): continue
            if klarge(S) != kk or max(S) < 13*min(S): continue
            total += 1
            V = max(S); ratio = Wwidth([u for u in S if u != V]) * 7 * V
            if ratio > 1: passcount += 1
            else:
                fails.append((V, ratio, S))
                if V > maxfail_vmax: maxfail_vmax = V
    print("  tested %d sets; via-max criterion passed %d, FAILED %d" % (total, passcount, len(fails)))
    print("  largest Vmax where via-max criterion FAILED: %d" % maxfail_vmax)
    if fails:
        fails.sort(key=lambda x: -x[0])
        print("  highest-Vmax criterion-failures (still need finite check / other v):")
        for V, r, S in fails[:4]:
            print("     Vmax=%d ratio=%.3f S=%s" % (V, float(r), S))
    sys.stdout.flush()
