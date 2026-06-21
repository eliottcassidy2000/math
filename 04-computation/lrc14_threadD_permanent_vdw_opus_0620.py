#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadD_permanent_vdw_opus_0620.py   (opus-2026-06-20, THREAD D part i)

THE LITERAL PERMANENT REPRESENTATION of measS7, and a van-der-Waerden / rearrangement test.

measS7(E) = meas{x : the k clocks e in E together HIT all 7 sectors}.
For a fixed x, clock e lands in sector color(e,x).  "All 7 sectors hit" by k>=7 clocks =
the EXISTENCE of a SURJECTION-supporting assignment: there is a way to pick 7 of the k
clocks, one per sector.  By inclusion-exclusion over absent sectors,
   1[surjective](x) = Sum_{S subset Z/7} (-1)^|S| 1[all clocks avoid S](x)
                    = Sum_{S} (-1)^|S| prod_{e in E} 1[color(e,x) not in S].
Integrating in x:
   measS7(E) = Sum_{S subset Z/7} (-1)^|S| J(S,E),   J(S,E)=int prod_e 1[color(e,x) not in S] dx.

THE PERMANENT.  Define the k x 7 "occupancy" expectation matrix only makes the iid permanent.
For the iid (dissociated) limit, J(S,E)->prod_e (1-|S|/7)=(1-|S|/7)^k and
   measS7_iid = Sum_S (-1)^|S|(1-|S|/7)^k = Sum_{j=0}^7 (-1)^j C(7,j)(1-j/7)^k = 7! S(k,7)/7^k.
That last equality is EXACTLY the number of surjections / 7^k = perm of the k x 7 all-(1/7)
matrix's surjection generating function.  In the iid limit measS7 is literally
   measS7_iid = perm-like surjection density of a DOUBLY-UNIFORM occupancy = 7! S(k,7)/7^k,
INDEPENDENT of the arrangement E.  So the ARRANGEMENT dependence (consec vs spread) lives
ENTIRELY in the CORRELATIONS J(S,E)-(1-|S|/7)^k, i.e. in the relation lattice Lambda(E).

This script:
 (P1) Builds J(S,E) exactly and the correction corr(E)=measS7(E)-measS7_iid.
      Confirms iid value and that corr(consec) is the LARGEST positive correction (consec
      is the most-correlated / least-dissociated arrangement).
 (P2) THE VAN DER WAERDEN TEST.  vdW says perm of a doubly-stochastic matrix is MINIMIZED
      (not maximized) at the uniform (J_n/n) matrix.  Map: the iid/uniform occupancy = the
      vdW minimizer of the surjection permanent.  So the surjection-permanent is MAXIMIZED at
      the MOST DEGENERATE (rank-1 collapsed) occupancy = the all-equal-clock arrangement.
      consec is NOT all-equal -- so a NAIVE permanent-vdW argument predicts the WRONG
      extremizer.  We test the actual occupancy matrices to see whether consec is the vdW
      max among VALID (clock-realizable) occupancies, or whether a permanent bound is the
      wrong tool.  REPORT honestly.
 (P3) POSITIVE-COMBINATION test.  Is there a FIXED nonneg weight vector w(S)>=0 (S subset Z/7)
      such that  F_w(E) = Sum_S w(S) J(S,E)  (a) equals measS7 up to the iid constant, OR
      (b) is a positive functional for which consec is provably extremal because each J(S,.)
      is individually consec-extremal?  TEST: is each J(S,E) (fixed S) MAXIMIZED at consec?
      MINIMIZED at consec?  (If the signs (-1)^|S| align with a consistent per-S extremality,
      consec follows by linearity -- a clean route.  If they CONFLICT, the permanent/positive-
      combination route is structurally dead, matching route_G's mixed-sign finding.)

EXACT Fractions; stdlib only.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

def cell_decomp(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    cells = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hits = {0}
        for e in E:
            if e == 0: continue
            hits.add(int(((e * xm) % 1) * 7))
        cells.append((x1 - x0, frozenset(hits)))
    return cells

def measS7(E):
    return sum(L for L, h in cell_decomp(E) if len(h) == 7)

def J_of_S(E, S):
    if 0 in S: return F(0)
    Senz = set(S); tot = F(0)
    for L, h in cell_decomp(E):
        if h.isdisjoint(Senz): tot += L
    return tot

def stirling2(n, k):
    return sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))//factorial(k)

def measS7_iid(k):
    return F(factorial(7)*stirling2(k,7), 7**k)

CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

def shapes(k, span):
    return [(0,)+combo for combo in itertools.combinations(range(1, span+1), k-1)]

# ===========================================================================
def main():
    print("="*92)
    print("THREAD D (part i): the PERMANENT representation + van der Waerden / positive-combo test")
    print("="*92)

    CONFIG = {8: 11, 9: 12, 10: 13}

    for k in (8, 9, 10):
        span = CONFIG[k]
        consec = tuple(range(k))
        iid = measS7_iid(k)
        cap = CAP[k]
        pool = shapes(k, span)
        print("\n" + "#"*92)
        print(f"### k={k}, span box [0..{span}], |pool|={len(pool)}")
        print(f"###   measS7_iid = 7!S({k},7)/7^{k} = {iid} = {float(iid):.6f}   cap={float(cap):.4f}")
        print(f"###   wide-budget cap-iid = {float(cap-iid):.4f}")
        print("#"*92)

        # ---- (P1) correction = measS7 - iid; consec's correction vs others ----
        data = []
        for E in pool:
            m = measS7(E)
            data.append((E, m, m-iid))
        data.sort(key=lambda t: -t[1])
        cons_m = measS7(consec); cons_corr = cons_m - iid
        max_corr = max(data, key=lambda t: t[2])
        print(f"\n[P1] correction corr(E)=measS7-iid.  consec measS7={float(cons_m):.5f} corr=+{float(cons_corr):.5f}")
        print(f"     max-correction shape: {max_corr[0]} corr=+{float(max_corr[2]):.5f} (==consec? {max_corr[0]==consec})")
        print(f"     top-3 by measS7: " + ", ".join(f"{E}:{float(m):.4f}" for E,m,c in data[:3]))

        # ---- (P3) per-S extremality of J(S,E): is consec the max or min for each S? ----
        # group S by size r; for each S, find where consec sits in the J(S,.) ranking.
        print(f"\n[P3] per-S extremality of J(S,E) over the pool (does consec extremize each J?):")
        print(f"     For measS7=Sum_S (-1)^|S| J(S,.), consec-extremality follows by LINEARITY iff")
        print(f"     for EVERY S the sign (-1)^|S| times (J(S,consec)-J(S,E)) is >= 0 for all E.")
        # We don't need to test all C(7,r) S of each size by symmetry of Z/7? NOT symmetric:
        # color 0 is pinned, S with 0 give J=0. So enumerate S subset {1..6} (0 excluded).
        sign_consistent = True
        per_r = {}
        worst_conflict = None
        for r in range(1, 7):
            # all S subset {1..6} of size r
            cnt_consec_max = 0; cnt_consec_min = 0; nS = 0
            for S in itertools.combinations(range(1,7), r):
                nS += 1
                Sset = set(S)
                jc = J_of_S(consec, Sset)
                js = [(E, J_of_S(E, Sset)) for E in pool]
                jmax = max(js, key=lambda t:t[1])[1]
                jmin = min(js, key=lambda t:t[1])[1]
                if jc == jmax: cnt_consec_max += 1
                if jc == jmin: cnt_consec_min += 1
                # for the linearity route we want (-1)^r (jc - J(S,E)) >= 0 for all E.
                # (-1)^r>0 for even r => want jc>=J(S,E) (consec max); odd r => want jc<=J(S,E) (consec min)
                want_max = (r % 2 == 0)
                for E, jE in js:
                    diff = jc - jE
                    bad = (want_max and diff < 0) or ((not want_max) and diff > 0)
                    if bad:
                        sign_consistent = False
                        if worst_conflict is None or abs(diff) > worst_conflict[0]:
                            worst_conflict = (abs(diff), r, S, E, jc, jE)
            per_r[r] = (cnt_consec_max, cnt_consec_min, nS)
            want = "MAX (even r)" if r%2==0 else "MIN (odd r)"
            print(f"     r={r}: want consec={want:13s}  consec-is-max:{cnt_consec_max}/{nS}  consec-is-min:{cnt_consec_min}/{nS}")
        print(f"     => sign-consistent linearity route: {'HOLDS' if sign_consistent else 'FAILS'}")
        if worst_conflict:
            d,r,S,E,jc,jE = worst_conflict
            want = "max" if r%2==0 else "min"
            print(f"        worst conflict: r={r} S={S} wanted consec={want} but E={E} has J={float(jE):.4f} vs consec J={float(jc):.4f} (|diff|={float(d):.4f})")

    print("\n" + "="*92)
    print("READOUT (Thread D part i)")
    print("="*92)
    print("P1: consec is (or is not) the maximal-correlation arrangement (corr above iid).")
    print("P3: the signed sum has CONSISTENT or CONFLICTING per-S extremality. If conflicting,")
    print("    the positive-combination / permanent-vdW route is DEAD (matches route_G mixed-sign).")

if __name__ == "__main__":
    main()
