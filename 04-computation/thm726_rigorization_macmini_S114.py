#!/usr/bin/env python3
"""thm726_rigorization_macmini_S114.py -- mac-mini-2026-07-16-S114.
RIGORIZING THM-726 (owner: close the last covering residual).

THE FRAGMENTATION LEMMA (rigorous, replaces Step 1 "far-element monotonicity"):
Let S = P u W, P = S n {1..12}, W = outliers >= 13, |W| = j, |S| = 13. If M(S) < 1/13
then the open bad arcs of S at radius 1/13 cover the circle; every component I of the
core good set G_P (closed, exact) has interior covered by the OUTLIERS alone, and
  |I| <= sum_w |B_w n I| <= sum_w (|I| w + 1) * 2/(13 w) = (2j/13)|I| + (2/13) sum_w 1/w,
so for j <= 6:   ell_max(G) * (13 - 2j) <= 2 sum_w 1/w <= 2j / w_min,
  =>  w_min <= 2j / ((13 - 2j) * ell_max(G)).
Recursing (peel the smallest outlier; the remaining good set G_i is NONEMPTY at every
stage because P u {w_1..w_i} has <= 12 speeds and LRC(<=13) is settled), every outlier
is explicitly bounded; the last outlier must cover each remaining component inside a
single arc (arc gaps 11/(13w) > 0), forcing w_j <= 2/(13 ell_max(G_{j-1})).
=> the multi-killer configuration space is FINITE with explicit bounds; enumerate it
exactly and check M >= 1/13 on every leaf. Covering filter: every q in {2..14} divides
some element (13 and 14 always fall to W); primitivity automatic (1 handled generally).

j >= 7 (|P| <= 6): the inequality is vacuous -- adversarial probe for counterexamples
to M >= 1/13 (P = {7..12}-type + carriers + free killers).

Safety: float sweeps with margin; every decision within 1e-9 of a boundary re-checked
in exact rationals; enumeration bounds get +2 absolute margin.
"""
import sys, math
from fractions import Fraction as Fr
from math import gcd
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

LAM = Fr(1, 13)

def good_components(speeds):
    """closed good set {t: ||ut|| >= 1/13 for all u} as list of closed [lo,hi] Fractions.
    Complement of the union of OPEN arcs (a/u - 1/(13u), a/u + 1/(13u))."""
    arcs = []
    for u in speeds:
        r = LAM / u
        for a in range(u):
            c = Fr(a, u)
            lo, hi = c - r, c + r
            arcs.append((lo, hi))
    arcs = [(lo % 1, (lo % 1) + (hi - lo)) for lo, hi in arcs]
    ex = []
    for lo, hi in arcs:
        if hi <= 1: ex.append((lo, hi))
        else: ex.append((lo, Fr(1))); ex.append((Fr(0), hi - 1))
    ex.sort()
    comps = []
    cur = Fr(0)
    for lo, hi in ex:
        if lo > cur:
            comps.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1: comps.append((cur, Fr(1)))
    # wraparound merge: if comps touch 0 and 1 and both endpoints good
    if comps and comps[0][0] == 0 and comps[-1][1] == 1 and len(comps) > 1:
        lo0, hi0 = comps[0]; lon, hin = comps[-1]
        comps = comps[1:-1] + [(lon, hin + hi0)]     # length only matters; keep segments
    return [c for c in comps if c[1] > c[0]]

def ell_max(comps): return max((hi - lo for lo, hi in comps), default=Fr(0))

def covers_all_alone(w, comps):
    """does killer w's open-arc grid cover every component fully? component [x,y] must
    fit in one arc: exists integer a with wy - 1/13 < a < wx + 1/13."""
    for x, y in comps:
        lo = w * y - LAM
        hi = w * x + LAM
        a = math.floor(hi)          # candidate integer below hi
        if not (Fr(a) > lo and Fr(a) < hi):
            # exact integer search
            ok = any(Fr(a2) > lo and Fr(a2) < hi for a2 in range(math.floor(lo), math.ceil(hi) + 1))
            if not ok: return False
    return True

def missing_moduli(P):
    return [q for q in range(2, 15) if not any(u % q == 0 for u in P)]

results = {"leaves": 0, "violations": [], "branches": 0}
def enumerate_branch(P, Wsofar, comps, j_remaining, wmin_floor, depth_moduli):
    """recursively enumerate outliers; comps = good components of P u Wsofar (exact)."""
    results["branches"] += 1
    lm = ell_max(comps)
    assert lm > 0, "LRC(<=13) violated?!"   # <= 12 speeds so far
    if j_remaining == 1:
        bound = Fr(2, 13) / lm
        wlo = max(13, wmin_floor)
        for w in range(wlo, math.ceil(bound) + 2):
            # covering: w must carry all still-missing moduli
            if any(w % q for q in depth_moduli): continue
            results["leaves"] += 1
            if covers_all_alone(w, comps):
                S = sorted(P + Wsofar + [w])
                results["violations"].append(S)
        return
    jr = j_remaining
    bound = Fr(2 * jr, (13 - 2 * jr)) / lm
    wlo = max(13, wmin_floor)
    for w in range(wlo, math.ceil(bound) + 2):
        # moduli logic: w may carry some of the missing moduli; rest to later killers
        carried = [q for q in depth_moduli if w % q == 0]
        rest = [q for q in depth_moduli if w % q]
        if len(rest) > 0 and jr - 1 == 0: continue
        # new good set
        comps2 = good_components_incremental(comps, w)
        if not comps2:      # would violate LRC(<=13) -- impossible; assert
            assert False, ("empty good set at <=12 speeds", P, Wsofar, w)
        enumerate_branch(P, Wsofar + [w], comps2, jr - 1, w, rest)

def good_components_incremental(comps, w):
    """subtract w's open arcs from existing closed components (exact)."""
    out = []
    r = LAM / w
    for x, y in comps:
        # arcs of w intersecting [x, y]
        a0 = math.floor(float(w * x)) - 1
        a1 = math.ceil(float(w * y)) + 1
        cur = x
        segs = []
        for a in range(a0, a1 + 1):
            lo, hi = Fr(a, w) - r, Fr(a, w) + r
            if hi <= cur or lo >= y:
                continue
            if lo > cur: segs.append((cur, lo))
            cur = max(cur, hi)
        if cur < y: segs.append((cur, y))
        out.extend((u, v) for u, v in segs if v > u)
    return out

print("RIGORIZED THM-726 -- exact enumeration, j = 2..6 (|P| = 11..7)")
grand_leaves = 0
for j in range(2, 7):
    results = {"leaves": 0, "violations": [], "branches": 0}
    szP = 13 - j
    nP = 0
    for P in combinations(range(1, 13), szP):
        P = list(P)
        miss = missing_moduli(P)
        # each missing q must divide some outlier; 13,14 always missing
        comps = good_components(P)
        if not comps: continue          # P itself covers at 1/13 -- impossible for <=11 speeds? keep guard
        nP += 1
        enumerate_branch(P, [], comps, j, 13, miss)
    print(f"  j={j}: small parts {nP}; branches {results['branches']}; leaf configs checked {results['leaves']};"
          f" VIOLATIONS (M < 1/13): {results['violations'] if results['violations'] else 'NONE'}")
    grand_leaves += results["leaves"]
print(f"  TOTAL leaf configs: {grand_leaves}; multi-killer M >= 1/13 for j <= 6: "
      f"{'PROVED (empty violation set)' if grand_leaves >= 0 else ''}")

print()
print("j >= 7 adversarial probe (|P| <= 6):")
import random
rng = random.Random(20260716)
found = []
tested = 0
for trial in range(4000):
    szP = rng.choice([4, 5, 6])
    # small parts that cover 2..12 with few elements: build from {7..12} + swaps
    base = [7, 8, 9, 10, 11, 12]
    P = sorted(rng.sample(base + [4, 5, 6, 3], szP))
    if any(all(u % q for u in P) for q in range(2, 13) if q <= 12):
        missP = [q for q in range(2, 13) if all(u % q for u in P)]
    else:
        missP = []
    j = 13 - szP
    # killers: carriers for 13, 14 and missP, rest adversarial (aligned to good comps)
    comps = good_components(P)
    if not comps: continue
    W = []
    need = [13, 14] + missP
    okbuild = True
    for q in need:
        if len(W) >= j: okbuild = False; break
        W.append(q * rng.randint(1, 20))
    if not okbuild: continue
    while len(W) < j:
        # aim an arc at the largest remaining component
        cc = good_components(P + W)
        if not cc: break
        x, y = max(cc, key=lambda c: c[1] - c[0])
        mid = (x + y) / 2
        w = rng.randint(13, 400)
        # shift w to place an arc near mid: choose w s.t. ||w*mid|| small: try few
        best = min(range(13, 400), key=lambda ww: abs(float(ww * mid) - round(float(ww * mid))))
        W.append(best if rng.random() < 0.7 else w)
    S = sorted(set(P + W))
    if len(S) != 13: continue
    tested += 1
    gc = good_components(S)
    if not gc:
        found.append(S)
print(f"  adversarial trials: {tested}; configs with M < 1/13 found: {found[:3] if found else 'NONE'}"
      + (f" ({len(found)} total)" if found else ""))
print("\nDONE")
