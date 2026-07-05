#!/usr/bin/env python3
"""
mac-mini-2026-07-04-S49 -- HYP-4095: the LOOSE-BASE SWEEP through THM-619's
alignment-band pipeline + the TWO-PIN data for the generic lemma.

Base families swept (structured, per MISTAKE-102's lesson -- no sampling):
  (F1) {1..11, x},        x = 12..60
  (F2) {1..10, x, y},     12 <= x < y <= 40
  (F3) drop-families: {1..13}\{a,b} u {x},  a<b<=13, x = 14..40  (12 elts)
  (F4) pinless ten-cover extensions: {5..14} u {x, y}, 15 <= x < y <= 45
  (F5) dilation-mixed: c*{1..11} u {x} for c = 2,3, x chosen primitive
For each PRIMITIVE base with M(B) > 1/13 (loose; tight bases go to the CRT
case): run bands + pins + window; report candidates; exact-check survivors.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import sys, itertools
sys.path.insert(0, '04-computation')
from lonely_profile import profile

R13 = F(1, 13)

def dist(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def components_at(S, r):
    iv = []
    for v in S:
        half = F(r) / v
        for k in range(v):
            c = F(k, v); a, b = c - half, c + half
            if a < 0: iv.append((a + 1, F(1))); iv.append((F(0), b))
            elif b > 1: iv.append((a, F(1))); iv.append((F(0), b - 1))
            else: iv.append((a, b))
    iv.sort()
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            if b > merged[-1][1]: merged[-1] = (merged[-1][0], b)
        else: merged.append((a, b))
    comps = []
    for i in range(len(merged)):
        b_cur = merged[i][1]
        a_next = merged[(i + 1) % len(merged)][0] + (1 if i + 1 == len(merged) else 0)
        if a_next > b_cur:
            comps.append((b_cur, a_next))
    return comps

def M_exact_prof(S):
    for cap in (8, 5, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m
    return None

def pipeline(B):
    """Returns (status, ncand, nbad). status: 'tight'|'swept'."""
    B = sorted(B)
    if reduce(gcd, B) != 1:
        return ('nonprim', 0, 0)
    MB = M_exact_prof(B)
    if MB is None or MB == R13:
        return ('tight', 0, 0)
    if MB < R13:
        return ('below', 0, 0)   # impossible by LRC(13); flag if seen
    comps = components_at(B, R13)
    b = max(B)
    wmax = 13 * b
    pins = [q for q in range(2, 15) if not any(v % q == 0 for v in B)]
    # pin lcm shortcut
    L = 1
    for q in pins: L = L * q // gcd(L, q)
    cands = []
    w0 = ((b // L) + 1) * L if L > 1 else b + 1
    step = L if L > 1 else 1
    w = w0
    while w <= wmax:
        ok = True
        for (a, bb) in comps:
            cJ = (a + bb) / 2; h = (bb - a) / 2
            beta = R13 - h * w
            if beta < 0 or dist(w * cJ) > beta:
                ok = False; break
        if ok and reduce(gcd, B + [w]) == 1 and all(any(v % q == 0 for v in B + [w]) for q in range(2, 15)):
            cands.append(w)
        w += step
    bad = []
    for w in cands:
        MV = M_exact_prof(B + [w])
        if MV is not None and MV < R13:
            bad.append((w, MV))
    return ('swept', len(cands), len(bad)), cands, bad

results = {'tight': 0, 'swept_empty': 0, 'swept_cand': 0, 'bad': 0, 'nonprim': 0, 'below': 0}
survivors = []
def run(B, fam):
    global results, survivors
    r = pipeline(B)
    if isinstance(r[0], tuple):
        (st, nc, nb), cands, bad = r
        if nc == 0: results['swept_empty'] += 1
        else:
            results['swept_cand'] += 1
            survivors.append((fam, sorted(B), cands, bad))
        results['bad'] += nb
    else:
        results[r[0]] += 1

import sys as _s
# F1
for x in range(12, 61):
    run(list(range(1, 12)) + [x], 'F1')
# F2 (subsample the grid for time: y-x pattern full)
for x in range(12, 31):
    for y in range(x + 1, 31):
        run(list(range(1, 11)) + [x, y], 'F2')
# F3
for a, b2 in itertools.combinations(range(1, 14), 2):
    base12 = [v for v in range(1, 14) if v not in (a, b2)]
    for x in range(14, 28):
        run(base12 + [x], 'F3')
# F4 pinless ten-cover extensions
core = list(range(5, 15))
for x in range(15, 36):
    for y in range(x + 1, 36):
        run(core + [x, y], 'F4')
# F5
for c in (2, 3):
    for x in range(1, 40):
        Bc = [c * k for k in range(1, 12)] + [x]
        if len(set(Bc)) == 12:
            run(Bc, 'F5')

print("SWEEP RESULTS (loose-base space, structured families F1-F5):")
for k, v in results.items():
    print(f"  {k}: {v}")
print(f"  survivors (bases with nonempty band intersection): {len(survivors)}")
for fam, B, cands, bad in survivors[:12]:
    print(f"   [{fam}] {B}: candidates {cands[:6]}{'...' if len(cands)>6 else ''} "
          f"bad={[(w, str(m)) for w, m in bad]}")
print(f"\nTOTAL M(V) < 1/13 violations: {results['bad']}")
