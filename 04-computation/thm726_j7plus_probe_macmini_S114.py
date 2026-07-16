#!/usr/bin/env python3
"""thm726_j7plus_probe_macmini_S114.py -- j >= 7 (|P| <= 6) adversarial probe:
can 7+ outliers >= 13 push a covering 13-set below M = 1/13? (The fragmentation
inequality is vacuous for j >= 7; this regime routes to the loose/density tile.)"""
import sys, math, random
from fractions import Fraction as Fr
sys.stdout.reconfigure(line_buffering=True)
LAM = 1.0 / 13
LAMx = Fr(1, 13)

def good_components_f(speeds):
    ex = []
    for u in speeds:
        r = LAM / u
        for a in range(u + 1):
            lo, hi = a / u - r, a / u + r
            ex.append((max(lo, 0.0), min(hi, 1.0)))
            if lo < 0: ex.append((lo + 1, 1.0))
            if hi > 1: ex.append((0.0, hi - 1))
    ex.sort()
    comps, cur = [], 0.0
    for lo, hi in ex:
        if lo > cur + 1e-15: comps.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1 - 1e-15: comps.append((cur, 1.0))
    return comps

def exact_covered(S):
    ex = []
    for u in S:
        r = LAMx / u
        for a in range(u + 1):
            c = Fr(a, u); lo, hi = c - r, c + r
            ex.append((max(lo, Fr(0)), min(hi, Fr(1))))
            if lo < 0: ex.append((lo + 1, Fr(1)))
            if hi > 1: ex.append((Fr(0), hi - 1))
    ex.sort()
    cur = Fr(0)
    for lo, hi in ex:
        if lo > cur: return False
        cur = max(cur, hi)
    return cur >= 1

rng = random.Random(20260716)
found = []; tested = 0
for trial in range(8000):
    szP = rng.choice([4, 5, 6])
    P = sorted(rng.sample([3, 4, 5, 6, 7, 8, 9, 10, 11, 12], szP))
    j = 13 - szP
    comps = good_components_f(P)
    if not comps: continue
    W = []
    need = [q for q in range(2, 15) if not any(u % q == 0 for u in P)]
    for q in need:
        if len(W) >= j: break
        W.append(q * rng.randint(1, 15))
    ok = True
    while len(W) < j:
        cc = good_components_f(P + W)
        if not cc: break
        x, y = max(cc, key=lambda c: c[1] - c[0])
        mid = (x + y) / 2
        best = min(range(13, 500), key=lambda ww: abs(ww * mid - round(ww * mid)) / ww * 13)
        W.append(best if best not in W else rng.randint(13, 500))
    S = sorted(set(P + W))
    if len(S) != 13: continue
    tested += 1
    if not good_components_f(S) and exact_covered(S):
        found.append(S)
print(f"j>=7 probe: trials {tested}; M < 1/13 configs found: {found[:4] if found else 'NONE'}"
      + (f" ({len(found)} total)" if found else ""))
