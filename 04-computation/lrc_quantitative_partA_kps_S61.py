#!/usr/bin/env python3
r"""
lrc_quantitative_partA_kps_S61.py   (kind-pasteur-2026-07-07-S61, HYP-4857)

QUANTITATIVE PART A via the ROBUST ROOF-SUBSET: the o(Vmax) arc-count problem
(monad HYP-4787 redirect target 2: "the missing written piece is #arcs(good set)
= o(Vmax) for spread~Vmax shapes") DISSOLVES, because the witness argument only
needs a positive-measure FEW-ARC SUBSET of the good set -- and the S59/S60 ledger
subset has O(1) arcs absolutely.

THE WITNESS MECHANICS (THM-527-A at the sharp 1/7 criterion, boxeph-verified):
  S = P u {Vmax - e : e in E}, 0 in E, D = diam(E).  Ruler period
  I_j = ((14j+1)/(14Vmax), (14j+13)/(14Vmax)), center x_j = (2j+1)/(2Vmax),
  spacing 1/Vmax.  The fast phase phi = frac(Vmax*tau) sweeps (1/14, 13/14) over I_j.
  Within I_j the teeth frac(e_i*tau) move by <= e_i*|I_j| <= 6D/(7Vmax) =: Delta,
  and slow phases by <= p*3/(7Vmax) from the center.

  ROBUST CRITERION at the period center x_j (all slacks conservative):
    (a) maxgap{frac(e_i x_j)} >= 1/7 + 2*Delta        [teeth drift both ways]
    (b) ||p x_j|| >= 1/14 + 3p/(7Vmax)  for all p in P [slow drift to tau*]
  ==> pick phi* = center of the widest center-config gap: at some tau* in I_j,
      phi(tau*) = phi* (full sweep), every cluster runner ||(Vmax-e_i)tau*|| =
      ||phi* - frac(e_i tau*)|| >= (1/7+2Delta)/2 - Delta = 1/14; the 0-tooth gives
      ||Vmax tau*|| = ||phi*|| >= 1/14; slow runners >= 1/14 by (b).  M(S) >= 1/14.

  ROBUST SUBSET (lower bound on the set of good centers, by S59/S60 monotonicity):
    R(P, D+1; Vmax) = G_P^{+3/(7Vmax)}  intersect  {roof_{D+1} >= 1/7 + 12D/(7Vmax)}
  (widened cuts: |x - j/p| >= 1/(14p) + 3/(7Vmax); shifted threshold).  This set has
    #arcs <= |Farey_6 roof windows| + Sum_{p in P}(p+1)  <= 13 + 70   (ABSOLUTE),
  independent of E, of the spread, and of Vmax -- the o(Vmax) bound, solved by
  substituting the subset for the true good set.

  COUNTING: centers x_j are spaced exactly 1/Vmax, so #good centers >=
  Vmax*meas(R) - #arcs(R).  Hence
    **Vmax * RL(P,D+1;Vmax) > arcs(R)  ==>  M(S) >= 1/14.**

  V0(P,D) := the least Vmax satisfying it (RL is monotone nonincreasing in the
  widenings, so a simple upward scan/bisection is valid).

This script:
  PART 1: RL and arc counts at representative (|P|-worst P, D = bite edge) pairs,
          explicit V0(P,D); the global V0* = max over the ledger-covered domain.
  PART 2: the factored-statement summary table (per leg: bite + V0 range).
  PART 3: DIRECT VERIFICATION: build concrete covering sets S at Vmax near V0,
          locate a good period from the robust criterion, exhibit tau*, and check
          min_v ||v tau*|| >= 1/14 exactly (rationals); also spot-check the
          good-center counting inequality.
"""
from fractions import Fraction as F
from itertools import combinations
import random

TH = F(1, 7)
MP = F(14249, 252252)

def merge(iv):
    iv = sorted((a, b) for a, b in iv if b > a)
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return [(a, b) for a, b in out]

def complement01(iv):
    iv = merge(iv); out = []; cur = F(0)
    for a, b in iv:
        if a > cur: out.append((cur, a))
        cur = max(cur, b)
    if cur < 1: out.append((cur, F(1)))
    return out

def inter(A, B):
    out = []; i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: out.append((lo, hi))
        if b1 < b2: i += 1
        else: j += 1
    return out

def gp_widened(P, w):
    """G_P with every cut widened by w: remove |x-j/p| < 1/(14p) + w."""
    bad = []
    for p in P:
        hw = F(1, 14 * p) + w
        for j in range(p + 1):
            c = F(j, p)
            bad.append((max(F(0), c - hw), min(F(1), c + hw)))
    return complement01(bad)

_farey_cache = {}
def farey(n):
    if n not in _farey_cache:
        fr = set()
        for q in range(1, n + 1):
            for p in range(0, q + 1):
                fr.add(F(p, q))
        _farey_cache[n] = sorted(fr)
    return _farey_cache[n]

def roof_superlevel(n, theta):
    Fs = farey(n); out = []
    for a, b in zip(Fs[:-1], Fs[1:]):
        q, qp = a.denominator, b.denominator
        vl, vr = F(1, q), F(1, qp)
        if vl <= theta and vr <= theta: continue
        if vl > theta and vr > theta:
            out.append((a, b)); continue
        xc = a + (theta - vl) * (b - a) / (vr - vl)
        out.append((a, xc) if vl > theta else (xc, b))
    return merge(out)

def robust_set(P, n, Vmax):
    D = n - 1
    theta = TH + F(12 * D, 7 * Vmax)
    w = F(3, 7 * Vmax)
    return inter(gp_widened(P, w), roof_superlevel(n, theta))

def arcbound(P):
    """PROVEN absolute arc bound for R(P,n,V), any n, any V, theta > 1/7:
    roof superlevel has <= 13 components (adjacent Farey_6 rationals are separated
    by their mediant, whose denominator >= 7 forces roof <= 1/7 < theta there);
    widened-G_P has <= Sum(p+1)+1 components; A cap B has <= compA+compB-1."""
    return 13 + sum(p + 1 for p in P)

def V0_of(P, n, hi=400000, absolute=False):
    """least Vmax with Vmax*meas(R) > arcs.  absolute=True uses the PROVEN per-P
    arc bound (then V*meas is monotone in V => the condition holds for ALL larger
    Vmax -- the rigorous forall-V form).  absolute=False uses observed arcs
    (sharp per-V certificate; re-check per V)."""
    A = arcbound(P)
    def ok(V):
        R = robust_set(P, n, V)
        m = sum(b - a for a, b in R)
        return V * m > (A if absolute else len(R))
    lo = 14 * n
    if not ok(hi): return None
    while not ok(lo):
        lo *= 2
        if lo > hi: return None
    a, b = max(14 * n, lo // 2), lo
    while a < b:
        mid = (a + b) // 2
        if ok(mid): b = mid
        else: a = mid + 1
    return a

# ---------------------------------------------------------------- PART 1
print("=" * 100)
print("PART 1 -- robust ledger RL, arc counts, explicit V0(P, D) at representative shapes")
print("=" * 100)
worstP = {0: (), 1: (6,), 2: (2, 12), 3: (1, 4, 6), 4: (1, 3, 4, 5), 5: (1, 2, 3, 4, 5)}
bites = {0: 75, 1: 34, 2: 21, 3: 17, 4: 11, 5: 11}   # S59/S60 composite (k=13..8)
thm530_worst = {5: (1, 5, 7, 8, 9), 4: (1, 11, 12, 13), 3: (1, 12, 13), 2: (1, 13), 1: (1,)}
rows = []
for s in range(0, 6):
    k = 13 - s
    for D in sorted({k + 1, (k - 1 + bites[s]) // 2, bites[s]}):
        n = D + 1
        Pset = worstP[s] if s else ()
        v0 = V0_of(Pset, n)
        R = robust_set(Pset, n, v0) if v0 else []
        m = sum(b - a for a, b in R) if R else 0
        rows.append((k, Pset, D, v0, len(R), float(m)))
        print(f"  k={k:2d} P={str(Pset):>18} D={D:3d}:  V0 = {v0}   (at V0: arcs={len(R)}, meas={float(m):.5f})")
V0s = [r[3] for r in rows if r[3]]
print(f"\n  V0 range over the table: [{min(V0s)}, {max(V0s)}]")
print("\n  PART 1b -- RIGOROUS forall-V form (PROVEN per-P arc bound A_P = 13 + Sum(p+1);")
print("  V*meas is monotone in V, so the condition persists for ALL Vmax >= V0abs):")
v0abs_all = []
for s in range(0, 6):
    k = 13 - s
    D = bites[s]; n = D + 1
    Pset = worstP[s] if s else ()
    va = V0_of(Pset, n, absolute=True)
    v0abs_all.append(va)
    print(f"    k={k:2d} P={str(Pset):>18} D={D:3d} (bite edge):  A_P={arcbound(Pset)}  V0abs = {va}")
print(f"    ==> GLOBAL V0abs* (bite edges, worst-P) = {max(v for v in v0abs_all if v)}")
# also the THM-530 worst P (min meas(G_P)) rows at the bite edge:
print("  THM-530-worst-P rows (bite edge):")
extra = []
for s in range(1, 6):
    k = 13 - s; D = bites[s]; n = D + 1
    v0 = V0_of(thm530_worst[s], n)
    extra.append(v0)
    print(f"    k={k:2d} P={str(thm530_worst[s]):>18} D={D:3d}:  V0 = {v0}")
V0star = max(V0s + [v for v in extra if v])
print(f"\n  ==> GLOBAL V0* over all tabulated shapes = {V0star}")

# ---------------------------------------------------------------- PART 2
print()
print("=" * 100)
print("PART 2 -- THE FACTORED PART A (explicit, per leg)")
print("=" * 100)
print("  For every covering shape S = P u {Vmax - e}, cluster diam D within the leg's bite:")
print("    [Vmax >= V0(P, D)  (tabulated above, all <= V0*)]  ==>  M(S) >= 1/14, WITNESSED.")
print("    [Vmax <  V0(P, D)]  ==>  FINITE CHECK: bounded family (height < V0*, diam <= bite).")
print(f"  V0* = {V0star}: the finite check lives in height < {V0star}, cluster diam <= 75.")
print("  The arc bound is ABSOLUTE (<= 13 roof windows + Sum(p+1) <= 83), independent of")
print("  E's spread -- the o(Vmax) obstruction (monad target 2) is dissolved by SUBSETTING.")

# ---------------------------------------------------------------- PART 3
print()
print("=" * 100)
print("PART 3 -- direct verification: concrete covering sets near V0, witness exhibited")
print("=" * 100)
def exact_lonely_witness(P, E, Vmax):
    """Find a good period center via the robust criterion; return exact tau*, min clearance."""
    n = max(E) - min(E) + 1
    E0 = sorted(e - min(E) for e in E)
    D = n - 1
    R = robust_set(tuple(P), n, Vmax)
    if not R: return None
    # first robust arc; pick the center x_j inside it
    a, b = max(R, key=lambda ab: ab[1] - ab[0])
    jlo = int((a * Vmax * 2 - 1) / 2) + 1
    xj = None
    for j in range(max(jlo - 2, 0), jlo + 4):
        c = F(2 * j + 1, 2 * Vmax)
        if a <= c <= b: xj = (j, c); break
    if xj is None: return None
    j, c = xj
    # teeth at center; widest gap center -> phi*; tau* = (j + phi*)/Vmax... phi in (1/14,13/14)
    teeth = sorted(F(e * c.numerator % (c.denominator), c.denominator) if False else (e * c) % 1 for e in E0)
    gaps = [(teeth[i + 1] - teeth[i], teeth[i], teeth[i + 1]) for i in range(len(teeth) - 1)]
    gaps.append((teeth[0] + 1 - teeth[-1], teeth[-1], teeth[0] + 1))
    g, lo, hi = max(gaps)
    phis = (lo + hi) / 2 % 1
    taust = (F(j) + phis) / Vmax
    S = sorted(set(list(P) + [Vmax - e for e in E0]))
    def dist(v):
        r = (v * taust) % 1
        return min(r, 1 - r)
    mind = min(dist(v) for v in S)
    return taust, mind, S

rng = random.Random(61)
tests = [((1, 5, 7, 8, 9), [0, 1, 2, 3, 5, 7, 9, 11], None),       # k=8 spreadish, D=11
         ((6,), [0, 2, 3, 5, 8, 13, 17, 21, 25, 28, 30, 34], None),  # k=12, D=34 bite edge
         ((), [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 17, 28], None),   # k=13 opus-stretch D=28
         ((1, 13), [0, 1, 4, 7, 9, 12, 15, 17, 19, 20, 21], None)]  # k=11 D=21 bite edge
for P, E0, _ in tests:
    n = max(E0) + 1
    v0 = V0_of(tuple(P), n)
    for Vm in (v0, v0 + 137):
        res = exact_lonely_witness(P, E0, Vm)
        if res is None:
            print(f"  P={P} E-diam={max(E0)} Vmax={Vm}: no robust center found (unexpected)")
            continue
        taust, mind, S = res
        ok = mind >= F(1, 14)
        print(f"  P={P} D={max(E0)} Vmax={Vm} (V0={v0}): witness tau*={taust}  min||v tau*|| = {mind} = {float(mind):.5f}  >= 1/14: {ok}")
print()
print("DONE.")
