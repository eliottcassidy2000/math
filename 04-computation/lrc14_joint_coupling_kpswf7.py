#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 2 -- the FINITE LAYER-3 aggregate wall: the JOINT coupling.

Goal: prove consec maximizes measS7 on the full-residue stratum via a coupling
across ALL 6 cells simultaneously (mac-mini ROUTE-3 proved per-cell maximality
FALSE; consec trades single cells to win the SUM).

measS7(E) = sum_{a=1..6} W_a(E),  W_a(E)= meas{ y in cell a : E*y mod 1 hits all 7 sectors }.
cell a = [a/7 - 1/14, a/7 + 1/14].

KEY REFRAME (this script): the 6 cells a=1..6 plus the y=0 cell (a=0, which is the
short / always-covered piece) TILE the circle y in [0,1) into 7 arcs of length 1/7
each centered at a/7. So

    measS7(E) = meas{ y in [0,1) \ cell_0 : all 7 sectors hit }   (the a=0 cell is
    the trivial covered piece around y=0; excluded from S7 by the LRC structure).

So measS7 is ONE coverage measure on the whole circle, NOT 6 separate problems.
The JOINT coupling should live on this single circle.

ATTACKS TESTED:
 (1) WHOLE-CIRCLE coupling: is there a measure-preserving map phi_E : [0,1) -> [0,1)
     carrying consec's NON-covered set INTO E's non-covered set (so covered(E) subset
     phi(covered(consec)) => meas covered(E) <= meas covered(consec))?  Test the
     candidate: phi = identity refined by the conductance/binding structure. We test
     whether the SORTED gap-length vector of E's uncovered components is MAJORIZED by
     consec's (Schur route, Huffer-Shepp style on the WHOLE circle).
 (2) CONDUCTANCE-VECTOR majorization: c(E)=(c_1,...,c_6), c_r=sum_{e=r mod7}1/|e|.
     Does c(consec) MAJORIZE c(E) for every full-residue rival, AND is measS7 a
     Schur-convex function of c?  (HYP-2760: measS7 = function of c exactly.)
 (3) WIN harmonic-sum vector: does consec's per-sector binding-speed vector majorize?
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ---------- exact measS7 machinery (ported from route3 v2) ----------
def W_a(E, a):
    E = sorted(set(E))
    lo = F(a, 7) - F(1, 14); hi = F(a, 7) + F(1, 14)
    bps = {lo, hi}
    for e in E:
        if e == 0: continue
        d = 7 * abs(e)
        j0 = math.floor(lo * d); j1 = math.ceil(hi * d)
        for j in range(j0 - 1, j1 + 2):
            x = F(j, d)
            if lo <= x <= hi: bps.add(x)
    bps = sorted(bps); tot = F(0)
    for l, h in zip(bps, bps[1:]):
        if h <= l: continue
        xm = (l + h) / 2; hit = set()
        for e in E:
            v = e * xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if len(hit) == 7: tot += h - l
    return tot

def Wvec(E): return [W_a(E, a) for a in range(1, 7)]
def measS7(E): return sum(W_a(E, a) for a in range(1, 7))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def conductance_vec(E):
    """c_r = sum_{e != 0, e = r mod 7} 1/|e|, for r=1..6 (r=0 is the short)."""
    c = [F(0)]*7
    for e in E:
        if e == 0: continue
        c[e % 7] += F(1, abs(e))
    return c  # index 0..6; c[0] is residue-0 (the 'short')

def full_residue(E):
    return set(e % 7 for e in E) == set(range(7))

def majorizes(u, v):
    """Does vector u majorize v? (both same sum check + partial-sum dominance)
       Returns (is_major, equal_sum)."""
    su = sorted(u, reverse=True); sv = sorted(v, reverse=True)
    n = max(len(su), len(sv))
    su += [F(0)]*(n-len(su)); sv += [F(0)]*(n-len(sv))
    eq = sum(su) == sum(sv)
    pu = 0; pv = 0; ok = True
    for i in range(n):
        pu += su[i]; pv += sv[i]
        if pu < pv: ok = False
    return ok, eq

# ---------- build full-residue stratum ----------
def stratum(k, box):
    """All primitive, full-residue, min=0, sorted k-subsets within [0,box]."""
    out = []
    for combo in itertools.combinations(range(1, box+1), k-1):
        E = (0,) + combo
        if not full_residue(E): continue
        if not primitive(E): continue
        out.append(E)
    return out

if __name__ == "__main__":
    print("="*88)
    print("THREAD 2: JOINT COUPLING for consec-maximizes-measS7 on full-residue stratum")
    print("="*88)

    for k, box in [(8, 14), (9, 14), (10, 13)]:
        print(f"\n{'#'*84}\n# k={k}, box=[0,{box}]\n{'#'*84}")
        C = consec(k)
        mC = measS7(C)
        cC = conductance_vec(C)
        print(f"consec={C}")
        print(f"  measS7(consec) = {mC} = {float(mC):.6f}")
        print(f"  conductance c(consec) (r=0..6) = {[str(x) for x in cC]}")
        print(f"  c-vec(r=1..6) sorted desc = {[str(x) for x in sorted(cC[1:],reverse=True)]}")

        S = stratum(k, box)
        print(f"  full-residue primitive stratum size (box {box}) = {len(S)}")

        # ============ ATTACK 2: conductance-vector majorization ============
        # Hypothesis: c(consec) restricted to r=1..6 majorizes c(E) restricted to r=1..6,
        # and measS7 is Schur-convex in this 6-vector.
        print("\n[ATTACK 2] Conductance-vector majorization (r=1..6 active sectors):")
        cmaj_fail = 0; cmaj_ok = 0; eqsum_fail = 0
        worst_examples = []
        for E in S:
            if E == tuple(C): continue
            cE = conductance_vec(E)
            uC = cC[1:]   # r=1..6
            uE = cE[1:]
            ismaj, eq = majorizes(uC, uE)
            if not eq: eqsum_fail += 1
            if ismaj: cmaj_ok += 1
            else:
                cmaj_fail += 1
                if len(worst_examples) < 5:
                    worst_examples.append((E, [str(x) for x in sorted(uE,reverse=True)]))
        print(f"  consec c-vec majorizes rival c-vec (r=1..6): OK={cmaj_ok}  FAIL={cmaj_fail}")
        print(f"  (equal-sum holds for? FAIL count = {eqsum_fail}; note Foster sum is NOT fixed across stratum)")
        for E, u in worst_examples:
            print(f"     MAJ-FAIL: {E} c-vec={u}")

        # also: does the FULL 7-vec (include r=0) have a conserved sum (Foster)?
        sums = set(sum(conductance_vec(E)) for E in S)
        print(f"  distinct Foster total sums over stratum = {len(sums)} (e.g. {sorted(sums)[:3]} ...)")

        # ============ ATTACK 3: is measS7 monotone in the c-vector? ============
        # Test: order stratum by total conductance sum_{r=1..6} c_r; does measS7 track it?
        print("\n[ATTACK 3] Does measS7 increase with active-sector conductance C_act=sum_{r=1..6}c_r?")
        rows = []
        for E in S:
            cE = conductance_vec(E)
            Cact = sum(cE[1:])
            rows.append((Cact, measS7(E), E))
        rows.sort(key=lambda t: t[0], reverse=True)
        # is consec the max C_act?
        Cact_C = sum(cC[1:])
        max_Cact = rows[0][0]
        n_above = sum(1 for r in rows if r[0] > Cact_C)
        print(f"  C_act(consec)={Cact_C}={float(Cact_C):.4f}; max C_act in stratum={float(max_Cact):.4f}; #rivals with C_act>consec={n_above}")
        # correlation: does max-measS7 rival also have max C_act?
        best_meas = max(rows, key=lambda t: t[1])
        print(f"  argmax measS7 = {best_meas[2]} (measS7={float(best_meas[1]):.4f}), its C_act={float(best_meas[0]):.4f}")
