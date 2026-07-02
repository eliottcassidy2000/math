#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S94 -- HYP-3850 parts (b2),(f),(g).

(1) PAIR-OVERLAP LAW (THM-594 Part B) exact verification, both branches, r=1/14
    and one off-critical radius (r=1/10) for the general-r algebraic form.
(2) OVERTAKING-DEFECT SIEVE (THM-592(v) consumer): exact concave-kink defect
    K(r0, 1/14) for covering-type 13-sets.  If K is uniformly small relative to
    the anchor slope, the tangent ladder closes covering floors from sieve-
    certified anchors.  Census over structured + random covering sets.
(3) ARC x RADIUS TWO-VARIABLE LP (demo): the layer-cake matrix
    M[c,i] = |{t in I_c : f_S(t) >= r_i}| refines BOTH HYP-3824's Q-localized LP
    (a row-sum) and the global profile (a column-sum).  Demo: AP vs GW tie in
    the last-cell column sums (slope) but differ in the matrix; and the matrix
    at Q=14 pins WHERE the lonely mass sits (the unit fractions).
"""
from fractions import Fraction as F
from math import gcd
import sys, random, itertools
sys.path.insert(0, '04-computation')
from lonely_profile import profile, _measure

def overlap_exact(p, q, r):
    iv = []
    for v in (p, q):
        half = F(r) / v
        for k in range(v):
            c = F(k, v); a, b = c - half, c + half
            if a < 0: iv.append((a + 1, F(1), v)); iv.append((F(0), b, v))
            elif b > 1: iv.append((a, F(1), v)); iv.append((F(0), b - 1, v))
            else: iv.append((a, b, v))
    P = [x for x in iv if x[2] == p]; Q = [x for x in iv if x[2] == q]
    tot = F(0)
    for a1, b1, _ in P:
        for a2, b2, _ in Q:
            lo, hi = max(a1, a2), min(b1, b2)
            if hi > lo: tot += hi - lo
    return tot

print("=" * 78)
print("(1) THM-594 Part B: two-branch pair law, exhaustive coprime pairs p<q<=13")
print("=" * 78)
r = F(1, 14)
viol = 0
for p in range(1, 14):
    for q in range(p + 1, 14):
        if gcd(p, q) != 1: continue
        ov = overlap_exact(p, q, r)
        if p + q <= 14:
            pred = 2 * r / q
        else:
            pred = F(q + 2 * p - 14, 7 * p * q)
        if ov != pred:
            viol += 1
            print(f"  !! ({p},{q}): overlap={ov} pred={pred}")
print(f"  all coprime pairs p<q<=13 at r=1/14: violations = {viol} (law EXACT)")
# general-r algebraic form: |D∩D| = 4r^2 - 2(1-2rp)(1-rq)/(pq) on branch 2
r2 = F(1, 10)
viol2 = 0
for p in range(1, 10):
    for q in range(p + 1, 12):
        if gcd(p, q) != 1: continue
        ov = overlap_exact(p, q, r2)
        if F(r2) * (p + q) <= 1:
            pred = 2 * r2 / q
        else:
            pred = 4 * r2 * r2 - 2 * (1 - 2 * r2 * p) * (1 - r2 * q) / (p * q)
        if ov != pred:
            viol2 += 1
            print(f"  !! r=1/10 ({p},{q}): overlap={ov}={float(ov):.6f} pred={pred}={float(pred):.6f}")
print(f"  general-r form at r=1/10, all coprime p<q (p<10,q<12): violations = {viol2}")

print()
print("=" * 78)
print("(2) OVERTAKING-DEFECT SIEVE: K(r0, 1/14) for covering-type 13-sets")
print("=" * 78)
r0 = F(1, 16)
def covering_sets_sample():
    sets = {
        "deep well {1..12,182}": list(range(1, 13)) + [182],
        "single-drop {1..11,13,24?} cover14 {1,2,..,11,13,28}": list(range(1, 12)) + [13, 28],
        "cover {1..10,12,14,22}": list(range(1, 11)) + [12, 14, 22],
        "cover {2..13,14}": list(range(2, 14)) + [14],
        "cover {1,2,3,4,5,6,7,8,9,10,22,26,14}": [1,2,3,4,5,6,7,8,9,10,22,26,14],
    }
    random.seed(94)
    tries = 0
    while len(sets) < 25 and tries < 4000:
        tries += 1
        S = sorted(random.sample(range(1, 60), 13))
        # covering predicate: multiple of every q in 2..14
        if all(any(v % q == 0 for v in S) for q in range(2, 15)):
            sets[f"rand{len(sets)}: {S}"] = S
    return sets

maxK = (0.0, None)
print(f"  anchor r0={r0}; window ({r0}, 1/14); K = sum of concave slope-jumps in window")
for name, S in covering_sets_sample().items():
    p = profile(S, F(1, 14))
    K = p.defect(r0, F(1, 14))
    m0 = p.measure(r0); sl = p.slope(r0 + F(1, 10**9))
    ladder = m0 + (sl - K) * (F(1, 14) - r0)   # THM-592(v) lower bound at 1/14
    ok = "LADDER>0 ✓" if ladder > 0 else "ladder fails (needs later anchor)"
    if float(K) > maxK[0]: maxK = (float(K), name)
    print(f"  K={float(K):7.4f}  m(r0)={float(m0):.4f}  slope={float(sl):8.3f}  "
          f"ladder@1/14={float(ladder):+.4f}  {ok}   {name[:46]}")
print(f"  max K over sample: {maxK[0]:.4f} ({maxK[1][:40]})")

print()
print("=" * 78)
print("(3) ARC x RADIUS LAYER-CAKE MATRIX (two-variable localization demo)")
print("=" * 78)
def layer_matrix(S, Q, radii):
    """M[c][i] = |{t in [c/Q,(c+1)/Q) : min_v||vt|| >= r_i}| exact."""
    rows = []
    for c in range(Q):
        lo, hi = F(c, Q), F(c + 1, Q)
        row = []
        for rr in radii:
            # measure of lonely set within [lo,hi): complement of danger arcs
            iv = []
            for v in S:
                half = F(rr) / v
                for k in range(v):
                    ctr = F(k, v); a, b = ctr - half, ctr + half
                    for (aa, bb) in ([(a + 1, F(1)), (F(0), b)] if a < 0 else
                                     [(a, F(1)), (F(0), b - 1)] if b > 1 else [(a, b)]):
                        na, nb = max(aa, lo), min(bb, hi)
                        if nb > na: iv.append((na, nb))
            iv.sort(); tot = F(0); ca = cb = None
            for a, b in iv:
                if cb is None: ca, cb = a, b
                elif a <= cb: cb = max(cb, b)
                else: tot += cb - ca; ca, cb = a, b
            if cb is not None: tot += cb - ca
            row.append((hi - lo) - tot)
        rows.append(row)
    return rows

radii = [F(1, 16), F(1, 15), F(27, 392)]  # 27/392 is inside (1/15, 1/14)
AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
Q = 14
MAP = layer_matrix(AP, Q, radii); MGW = layer_matrix(GW, Q, radii)
print(f"  Q={Q} arcs x radii {[str(x) for x in radii]}; showing nonzero rows (arc c: values)")
for name, M in [("AP", MAP), ("GW", MGW)]:
    tot = [sum(M[c][i] for c in range(Q)) for i in range(len(radii))]
    print(f"  {name}: column sums (global profile) = {[f'{float(x):.5f}' for x in tot]}")
    for c in range(Q):
        if any(x > 0 for x in M[c]):
            print(f"    arc {c:2d} [{c}/14,{c+1}/14): " +
                  "  ".join(f"{float(x):.5f}" for x in M[c]))
print("  => same last-column structure test: does the matrix separate AP/GW where")
print("     the final-cell slope ties?  Compare per-arc values at r=27/392 (last cell):")
diff = [(c, MAP[c][2], MGW[c][2]) for c in range(Q) if MAP[c][2] != MGW[c][2]]
if diff:
    for c, a, g in diff:
        print(f"    arc {c:2d}: AP={float(a):.6f}  GW={float(g):.6f}")
else:
    print("    identical per-arc at the last cell (tie extends to Q=14 localization);")
    print("    separation appears at the second layer (r=1/16,1/15 columns) instead:")
    d2 = [(c, MAP[c][0], MGW[c][0]) for c in range(Q) if MAP[c][0] != MGW[c][0]]
    for c, a, g in d2[:6]:
        print(f"    arc {c:2d} @1/16: AP={float(a):.6f}  GW={float(g):.6f}")
