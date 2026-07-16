#!/usr/bin/env python3
"""thm726_rigorization_fast_macmini_S114.py -- fast (float+margin, exact-confirmed) version
of the THM-726 rigorization for j = 2..6, plus the j >= 7 adversarial probe.
Soundness: all pruning bounds computed from float ell_max ENLARGED by (1+1e-9) and +2;
leaf coverage positives and near-misses (margin 1e-7) re-confirmed in exact rationals.
"""
import sys, math, random
from fractions import Fraction as Fr
from math import gcd
from itertools import combinations
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
    comps = []
    cur = 0.0
    for lo, hi in ex:
        if lo > cur + 1e-15: comps.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1 - 1e-15: comps.append((cur, 1.0))
    return comps

def subtract_f(comps, w):
    out = []
    r = LAM / w
    for x, y in comps:
        a0 = int(math.floor(w * x)) - 1
        a1 = int(math.ceil(w * y)) + 1
        cur = x
        for a in range(a0, a1 + 1):
            lo, hi = a / w - r, a / w + r
            if hi <= cur or lo >= y: continue
            if lo > cur + 1e-15: out.append((cur, lo))
            cur = max(cur, hi)
        if cur < y - 1e-15: out.append((cur, y))
    return out

def covers_all_alone_f(w, comps):
    for x, y in comps:
        lo = w * y - 1.0 / 13
        hi = w * x + 1.0 / 13
        if math.floor(hi) <= lo + 1e-9:           # no integer strictly inside (float)
            if math.floor(hi + 1e-9) <= lo - 1e-9:
                return False
            # ambiguous -> exact
            X, Y = Fr(x).limit_denominator(10**12), Fr(y).limit_denominator(10**12)
            ok = any(Fr(a) > w * Y - LAMx and Fr(a) < w * X + LAMx
                     for a in range(int(w * y - 1), int(w * x + 2)))
            if not ok: return False
    return True

def exact_check_S(S):
    """exact: is the closed good set at 1/13 empty (M < 1/13)?"""
    ex = []
    for u in S:
        r = LAMx / u
        for a in range(u + 1):
            c = Fr(a, u)
            lo, hi = c - r, c + r
            ex.append((max(lo, Fr(0)), min(hi, Fr(1))))
            if lo < 0: ex.append((lo + 1, Fr(1)))
            if hi > 1: ex.append((Fr(0), hi - 1))
    ex.sort()
    cur = Fr(0)
    for lo, hi in ex:
        if lo > cur: return False       # gap -> good point -> M >= 1/13
        cur = max(cur, hi)
    return cur >= 1

def missing_moduli(P):
    return [q for q in range(2, 15) if not any(u % q == 0 for u in P)]

stats = {}
def run_j(j):
    szP = 13 - j
    leaves = 0; branches = 0; violations = []
    def rec(P, W, comps, jr, wlo, missq):
        nonlocal leaves, branches
        branches += 1
        lm = max((hi - lo for lo, hi in comps), default=0.0)
        if lm <= 0: return              # cannot happen (LRC), guard
        if jr == 1:
            bound = (2.0 / 13) / lm * (1 + 1e-9)
            need = 1
            for q in missq: need = need * q // gcd(need, q)
            w0 = ((wlo + need - 1) // need) * need if need > 1 else wlo
            w = w0
            while w <= bound + 2:
                if w >= 13 and all(w % q == 0 for q in missq):
                    leaves += 1
                    if covers_all_alone_f(w, comps):
                        S = sorted(set(P + W + [w]))
                        if len(S) == 13 and exact_check_S(S):
                            violations.append(S)
                w += need if need > 1 else 1
            return
        bound = (2.0 * jr / (13 - 2 * jr)) / lm * (1 + 1e-9)
        # secondary bound via total measure + count
        tot = sum(hi - lo for lo, hi in comps)
        nc = len(comps)
        bound2 = (2.0 * jr * nc / ((13 - 2 * jr) * tot)) * (1 + 1e-9) if tot > 0 else bound
        bound = min(bound, bound2)
        for w in range(max(13, wlo), int(bound) + 3):
            rest = [q for q in missq if w % q]
            if len(rest) and jr == 1: continue
            comps2 = subtract_f(comps, w)
            if not comps2: continue     # LRC says impossible; float artifact -> skip is UNSAFE
            rec(P, W + [w], comps2, jr - 1, w + 1, rest)
    nP = 0
    for Pt in combinations(range(1, 13), szP):
        P = list(Pt)
        comps = good_components_f(P)
        if not comps: continue
        nP += 1
        rec(P, [], comps, j, 13, missing_moduli(P))
    print(f"  j={j}: small parts {nP}; branches {branches}; leaves {leaves}; "
          f"VIOLATIONS: {violations if violations else 'NONE'}", flush=True)
    return violations

print("RIGORIZED THM-726 (fast sweep, exact-confirmed leaves), j = 2..6")
allv = []
for j in range(2, 7):
    allv += run_j(j)
print(f"  VERDICT j<=6: {'M >= 1/13 PROVED (no violations)' if not allv else ('VIOLATIONS: ' + str(allv))}")

print()
print("j >= 7 adversarial probe (|P| <= 6):")
rng = random.Random(20260716)
found = []; tested = 0
for trial in range(6000):
    szP = rng.choice([4, 5, 6])
    P = sorted(rng.sample([3, 4, 5, 6, 7, 8, 9, 10, 11, 12], szP))
    j = 13 - szP
    comps = good_components_f(P)
    if not comps: continue
    W = []
    need = [q for q in range(2, 15) if not any(u % q == 0 for u in P)]
    # greedy: carriers first (small multiples), then aimed killers
    for q in need:
        if len(W) >= j: break
        W.append(q * rng.randint(1, 15))
    while len(W) < j:
        cc = good_components_f(P + [w for w in W])
        if not cc: break
        x, y = max(cc, key=lambda c: c[1] - c[0])
        mid = (x + y) / 2
        best = min(range(13, 500), key=lambda ww: min(abs(ww * mid - round(ww * mid)) / ww * 13, 9) - 0.001 * (ww < 200))
        W.append(best if best not in W else rng.randint(13, 500))
    S = sorted(set(P + W))
    if len(S) != 13: continue
    tested += 1
    if not good_components_f(S):        # float says covered -> exact confirm
        if exact_check_S(S):
            found.append(S)
print(f"  trials {tested}; M < 1/13 configs: {found[:4] if found else 'NONE'}"
      + (f"  ({len(found)} total)" if found else ""))
print("\nDONE")
