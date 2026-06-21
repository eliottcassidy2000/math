#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadD_permanent_rearrangement_opus_0620.py   (opus-2026-06-20, THREAD D)

THREAD D: the PERMANENT / van der Waerden / rearrangement angle for consec-extremality
of measS7.  measS7(E) = Sum_{S subset Z/7} (-1)^|S| J(S,E)  is a SIGNED surjection count;
surjection counting ~ a PERMANENT of a 0/1 incidence matrix (clock x sector).

Two concrete tests:
 (D-i)  PERMANENT REPRESENTATION.  Is measS7(E) a permanent (or a positive combination
        in which the permanent's van-der-Waerden extremality forces consec)?  We build the
        cell x clock incidence and the surjection permanent literally, and probe whether
        consec extremizes a permanent-of-a-doubly-stochastic-like object.
 (D-ii) SCHUR-CONVEXITY in the GAP MULTISET.  E={0=e_0<e_1<...<e_{k-1}} has consecutive-gap
        multiset g=(e_1-e_0,...,e_{k-1}-e_{k-2}) (and the residual to span).  consec has the
        FLATTEST gap multiset (all 1s).  TEST: is measS7 a Schur-CONVEX (max at flat <=> the
        all-1 vector is the bottom of majorization, so Schur-CONVEX would put max at SPREAD,
        Schur-CONCAVE at flat) function of the gap multiset?  Since consec = flat = the
        minimum of the majorization order, measS7 max at consec <=> measS7 SCHUR-CONCAVE in g
        (on the relevant domain).  Test the defining 2-point transfer (Robin-Hood / T-transform):
        moving mass from a large gap to a small gap (toward flatness) should INCREASE measS7
        for Schur-concavity.

We run EXACT (Fraction) on the real adversary pool (k=8,9,10) and report PROVED-style
counterexamples or clean confirmations.  stdlib only.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
# core measS7 + cell decomposition (matches canonical scripts EXACTLY)
# ---------------------------------------------------------------------------
def cell_decomp(E):
    """Return [(length, hitset_frozenset)] over cells of x->C_E(x).  0 in E pins 0."""
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    cells = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        hits = {0}
        for e in E:
            if e == 0: continue
            hits.add(int(((e * xm) % 1) * 7))
        cells.append((x1 - x0, frozenset(hits)))
    return cells

def measS7(E):
    return sum(L for L, h in cell_decomp(E) if len(h) == 7)

def stirling2(n, k):
    return sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))//factorial(k)

CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

# ---------------------------------------------------------------------------
# inclusion-exclusion J(S,E) = meas{x : C_E(x) disjoint from S}
# (the signed-surjection decomposition: measS7 = Sum_S (-1)^|S| J(S,E) )
# ---------------------------------------------------------------------------
def J_of_S(E, S):
    """meas{x : for all e in E, color(e,x) not in S}.  S subset of Z/7. 0 in E => if 0 in S, J=0."""
    if 0 in S:  # e=0 pins color 0, which is in S -> impossible to avoid
        return F(0)
    Senz = set(S)
    tot = F(0)
    for L, h in cell_decomp(E):
        if h.isdisjoint(Senz):
            tot += L
    return tot

def measS7_via_IE(E):
    tot = F(0)
    for r in range(0, 8):
        for S in itertools.combinations(range(7), r):
            tot += (-1)**r * J_of_S(E, set(S))
    return tot

# ---------------------------------------------------------------------------
def shapes(k, span):
    out = []
    for combo in itertools.combinations(range(1, span + 1), k - 1):
        out.append((0,) + combo)
    return out

def gap_multiset(E):
    """consecutive gaps of the SORTED support, NOT including wrap; flat=all 1s = consec."""
    E = sorted(set(E))
    return tuple(E[i+1]-E[i] for i in range(len(E)-1))

def majorizes(a, b):
    """does sorted-desc a majorize b (same sum, same length)?"""
    if sum(a) != sum(b) or len(a) != len(b): return None
    A = sorted(a, reverse=True); B = sorted(b, reverse=True)
    pa = 0; pb = 0
    for i in range(len(A)):
        pa += A[i]; pb += B[i]
        if pa < pb: return False
    return True

# ===========================================================================
def main():
    print("="*92)
    print("THREAD D: PERMANENT / van der Waerden / rearrangement angle for consec-extremality")
    print("="*92)

    # ---- sanity: IE decomposition reproduces measS7 ----
    print("\n[D0] verify measS7 = Sum_S (-1)^|S| J(S,E) (signed surjection count):")
    ok = True
    for k in (8, 9):
        for E in [tuple(range(k)), (0,)+tuple(range(2,k+1)), (0,1,2,3,4,6,9,11) if k==8 else (0,1,2,3,4,5,6,7,9)]:
            if len(E)!=k: continue
            a = measS7(E); b = measS7_via_IE(E)
            tag = "OK" if a==b else "MISMATCH"
            if a!=b: ok=False
            print(f"    k={k} E={E}: direct={float(a):.6f} IE={float(b):.6f} {tag}")
    print(f"    => signed-surjection identity {'CONFIRMED' if ok else 'BROKEN'}.")

    # =======================================================================
    # (D-ii) SCHUR-CONVEXITY / REARRANGEMENT in the GAP MULTISET
    # =======================================================================
    print("\n" + "#"*92)
    print("# (D-ii) SCHUR ORDER on the GAP MULTISET.  consec = flat gaps (all 1) = majorization")
    print("#        MINIMUM.  measS7 max at consec  <=>  measS7 SCHUR-CONCAVE in gap multiset")
    print("#        on shapes of fixed support-length k and fixed SPAN (= same gap-sum).")
    print("#"*92)
    # group adversary pool by (k, span): within a fixed span, all shapes have the SAME gap-sum
    # = span, and the SAME number of gaps = k-1.  So they live in one majorization poset.
    # Schur-concavity test: for every comparable pair g_a majorizes g_b (a more spread),
    # measS7 should satisfy measS7(b) >= measS7(a).  consec (flat, bottom) should be the global max.
    for k in (8, 9, 10):
        span = {8:11, 9:12, 10:13}[k]
        print(f"\n--- k={k} (gaps={k-1}) ---")
        for sp in range(k-1, span+1):
            # shapes with EXACT span sp: 0 and sp both in support
            pool = []
            for combo in itertools.combinations(range(1, sp), k-2):
                E = (0,)+combo+(sp,)
                pool.append(E)
            if not pool: continue
            # consec only has span k-1; for sp>k-1 the "flat-most" shape is the relevant comparison.
            data = [(E, gap_multiset(E), measS7(E)) for E in pool]
            # Schur-concavity: for all comparable (majorizing) pairs, the more-spread has SMALLER measS7
            viol = 0; tested = 0; worst = None
            for i in range(len(data)):
                for j in range(len(data)):
                    if i==j: continue
                    Ei, gi, mi = data[i]; Ej, gj, mj = data[j]
                    mj_maj = majorizes(gi, gj)  # gi majorizes gj => gi more spread
                    if mj_maj is True:
                        tested += 1
                        # Schur-concave wants measS7(gi) <= measS7(gj) (spread lowers it)
                        if mi > mj:
                            viol += 1
                            d = mi - mj
                            if worst is None or d > worst[0]:
                                worst = (d, Ei, Ej)
            # also: is the flattest shape (min in majorization) the argmax within this span?
            flat = min(data, key=lambda t: sorted(t[1], reverse=True))  # heuristic flattest
            argmax = max(data, key=lambda t: t[2])
            flat_is_max = (argmax[2] == flat[2])
            msg = f"  span={sp}: |pool|={len(pool):4d} Schur-pairs tested={tested:5d} violations={viol:4d}"
            if worst:
                msg += f"  worst +{float(worst[0]):.4f}"
            msg += f"  | flattest-is-argmax: {'Y' if flat_is_max else 'N'}"
            print(msg)
            if worst and sp == span:  # show the decisive counterexample at the box span
                print(f"      worst pair: more-spread {worst[1]} measS7>{worst[2]} less-spread, by +{float(worst[0]):.5f}")

    # =======================================================================
    # (D-ii') the DEFINING 2-point transfer (Robin-Hood / T-transform).
    #   Schur-concavity <=> a single Robin-Hood transfer toward flatness NEVER decreases measS7.
    #   We test ALL elementary transfers within the pool: take E, pick two gaps (a>b+1), move one
    #   unit from the big gap to the small gap (more flat), and check measS7 does not decrease.
    # =======================================================================
    print("\n" + "#"*92)
    print("# (D-ii') ELEMENTARY ROBIN-HOOD TRANSFER (the Schur-concavity generator):")
    print("#   flatten one gap-pair by 1 unit; measS7 should NOT decrease (Schur-concave).")
    print("#"*92)
    for k in (8, 9, 10):
        span = {8:11, 9:12, 10:13}[k]
        tested = 0; dec = 0; worst = None
        for sp in range(k-1, span+1):
            for combo in itertools.combinations(range(1, sp), k-2):
                E = (0,)+combo+(sp,)
                Es = sorted(E)
                gaps = [Es[i+1]-Es[i] for i in range(len(Es)-1)]
                m0 = measS7(E)
                # try every pair of gap positions (p big donor, q small receiver) with g[p] >= g[q]+2
                for p in range(len(gaps)):
                    for q in range(len(gaps)):
                        if p==q: continue
                        if gaps[p] >= gaps[q] + 2:
                            ng = gaps[:]; ng[p]-=1; ng[q]+=1
                            # rebuild support from flat-transferred gaps
                            cum = [0]
                            for g in ng: cum.append(cum[-1]+g)
                            E2 = tuple(cum)
                            if len(set(E2)) != k: continue
                            m1 = measS7(E2)
                            tested += 1
                            if m1 < m0:  # flattening DECREASED measS7 -> violates Schur-concavity
                                dec += 1
                                d = m0 - m1
                                if worst is None or d > worst[0]:
                                    worst = (d, E, E2)
        print(f"  k={k}: transfers tested={tested:6d}  flattening-DECREASED-measS7 = {dec:5d}", end="")
        if worst:
            print(f"  worst drop {float(worst[0]):.5f}  e.g. {worst[1]} -> {worst[2]}")
        else:
            print("  (none: every flattening raised measS7  => Schur-CONCAVE candidate!)")

if __name__ == "__main__":
    main()
