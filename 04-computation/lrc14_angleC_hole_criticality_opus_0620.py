#!/usr/bin/env python3
"""
lrc14_angleC_hole_criticality_opus_0620.py   (opus-2026-06-20)

ANGLE C -- HOLE-REMOVAL criticality (CA point-deletion).

SETUP. Fix span N. The "full board" is the consecutive run consec_{N+1} = {0,1,...,N}.
A bounded-span residual shape with |E|=k is  E = consec_{N+1} \ H  where H is the
set of HOLES (deleted clocks), |H| = N+1-k, with 0 not in H (we keep 0; if 0 in H,
scale/translate handles it -- but measS7 PINS residue 0 to the e=0 clock, so 0 must
stay; consec_{N+1} always contains 0).

The seven-sector cover measure:
    measS7(E) = meas{ x in [0,1) : { c(e,x) : e in E } = Z/7 },
    c(e,x) = floor(7 e x) mod 7   (e=0 always gives residue 0).

CA / coloring reframe: at a fixed x, each clock e "lights up" exactly one of 7 cells
(its residue color). measS7(E) is the measure of x where E's clocks light up ALL 7
cells. Removing a clock h (a hole) can only ever DECREASE the cover region:

    cover(E) subset cover(E + h)   for any extra clock h     (monotone in the clock set).

DEFINE the CRITICALITY (defect) of clock h relative to a clock-set F = E + {h}:
    Delta_h(F) := measS7(F) - measS7(F \ {h})
                = meas{ x : F covers Z/7  but  F\{h} does not }
                = meas{ x : h is the UNIQUE coverer of some residue at x, and
                          removing it breaks the full cover }.

This is exactly the "h is sole supporter of some color" event. It is a CA point-defect.

GOAL of Angle C: certify measS7(E) <= cap_k for all bounded-span residual shapes by a
HOLE-COUNTING argument:
    measS7(E) = measS7(consec_{N+1}) - [coverage LOST by removing the holes H].
The loss is  L(H) = measS7(consec_{N+1}) - measS7(consec_{N+1} \ H).  We want
    L(H) >= measS7(consec_{N+1}) - cap_k   for EVERY hole-set H with |H|=N+1-k (and 0 not in H),
equivalently measS7(E) <= cap_k.

A naive single-clock-criticality LOWER bound on L(H) would be additive:
    L_LB(H) = sum_{h in H} Delta_h(consec_{N+1})   ??? -- but removals INTERACT
(once one hole is removed, another may become critical or stop being critical).
By inclusion-exclusion / submodularity the additive single-hole sum is NOT in general
a valid bound on the joint loss in either direction.  This script TESTS:
  (1) the exact per-clock criticality Delta_h(consec_{N+1}) for each span N;
  (2) whether sum_{h in H} Delta_h >= required-drop holds (additive certificate);
  (3) the TRUE joint loss L(H) for all residual H (exhaustive), to see if
      measS7(E) <= cap is actually always satisfied with margin, and which hole-set
      is the worst (largest residual measS7).
  (4) the submodularity / supermodularity direction of L (does additive over- or
      under-estimate?), since that decides whether a clean hole-counting bound exists.

EXACT Fraction arithmetic throughout; stdlib only.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
# Canonical seven-sector cover measure (matches kps-S9 engine, verified).
# ---------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = set([F(0), F(1)])
    for e in Enz:
        for m in range(0, 7*e+1):
            bps.add(F(m, 7*e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        res = set(int(7*e*xm) % 7 for e in E)
        if len(res) == 7: total += x1 - x0
    return total

# cover MEASURE for the cell-decomposed engine, returns the cover region's cells
# (as list of (x0,x1)) so we can do exact set operations / criticality.
def cover_cells(E):
    """Return sorted list of (x0,x1, covered_bool) over a common refinement of all
    breakpoints of clocks in E.  covered_bool = full Z/7 cover on that cell."""
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = set([F(0), F(1)])
    for e in Enz:
        for m in range(0, 7*e+1):
            bps.add(F(m, 7*e))
    bps = sorted(bps)
    cells = []
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        res = frozenset(int(7*e*xm) % 7 for e in E)
        cells.append((x0, x1, res))
    return cells

# ---------------------------------------------------------------------------
# cap_k canonical (HYP-2603 / THM-535), == prompt's cap floats.
# cap_8=2243/5880=0.38146, cap_9=1979/4004=0.49426, cap_10=55/91=0.60440.
# ---------------------------------------------------------------------------
CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91),
       11: F(66,91), 12: F(6,7), 13: F(1)}

def measS7_AP(m):
    return measS7(tuple(range(m)))

def main():
    print("="*96)
    print("ANGLE C: HOLE-REMOVAL CRITICALITY  (opus-2026-06-20)")
    print("="*96)
    print("cap_k:", {k: f'{str(v)}={float(v):.5f}' for k,v in CAP.items()})
    consec = {m: measS7_AP(m) for m in range(1, 16)}
    print("consec measS7:", {m: float(consec[m]) for m in range(7,15)})
    print()

    # ----- per-clock criticality of consec_{N+1} -----
    print("="*96)
    print("(1) PER-CLOCK CRITICALITY  Delta_h(consec_{N+1}) = measS7(full) - measS7(full minus h)")
    print("    h ranges over 1..N (0 is pinned, removing it always destroys residue 0 -> trivial).")
    print("="*96)
    crit = {}  # crit[N][h]
    for N in range(7, 14):
        full = tuple(range(N+1))
        mfull = measS7(full)
        crit[N] = {}
        line = []
        for h in range(0, N+1):
            sub = tuple(e for e in full if e != h)
            d = mfull - measS7(sub)
            crit[N][h] = d
            line.append((h, d))
        print(f"  N={N:2d} (full=consec_{N+1}, measS7={float(mfull):.5f}):")
        # print criticality for h=0..N
        for h,d in line:
            tag = "  <-- e=0 pinned" if h==0 else ""
            print(f"      Delta_{h:<2} = {str(d):>14} = {float(d):.6f}{tag}")
        # sum of single criticalities (h=1..N), compare to (full - cap) thresholds
        s_all = sum(crit[N][h] for h in range(1,N+1))
        print(f"      sum_h Delta_h (h=1..N) = {float(s_all):.6f}   (vs measS7(full)={float(mfull):.6f})")
        print()

    # ----- EXHAUSTIVE residual: worst-case bounded-span shape per (k, span) -----
    print("="*96)
    print("(2) EXHAUSTIVE bounded-span residual: for each k and span N, the worst (max measS7)")
    print("    shape E = consec_{N+1} minus holes H (|H|=N+1-k, 0 not in H), and its margin to cap.")
    print("="*96)
    # residual regime per prompt: k=8 span 8..13, k=9 span 9..13, k=10 span 11..13.
    # We also include span==k-1 (=consec itself) as the baseline.
    regimes = {8: range(7, 14), 9: range(8, 14), 10: range(9, 14)}
    overall = {}
    for k in [8, 9, 10]:
        ck = CAP[k]
        print(f"  --- k={k}, cap_k={str(ck)}={float(ck):.5f} ---")
        gmax = F(0); gmaxE = None
        for N in regimes[k]:
            full = list(range(N+1))
            nholes = (N+1) - k
            if nholes < 0:
                continue
            # holes chosen from {1..N} (keep 0). number of shapes = C(N, nholes)
            best = F(-1); bestE = None
            cnt = 0; viol = 0
            for H in itertools.combinations(range(1, N+1), nholes):
                E = tuple(e for e in full if e not in set(H))
                # require span exactly N: max(E)=N means N not a hole
                if N in H:  # would reduce span
                    continue
                # primitivity not required for upper bound (scale-invariant; sub-shape ok)
                v = measS7(E)
                cnt += 1
                if v > best: best, bestE = v, E
                if v > ck: viol += 1
                if v > gmax: gmax, gmaxE = v, E
            if bestE is None:
                continue
            margin = ck - best
            flag = "OK" if best <= ck else "*** VIOLATION ***"
            print(f"    N={N:2d}: shapes={cnt:5d}  holes={nholes}  worst measS7={float(best):.5f} "
                  f"margin={float(margin):+.5f}  viol={viol}  {flag}")
            print(f"           worst E={bestE}")
        overall[k] = (gmax, gmaxE)
        print(f"    >>> k={k} overall worst over residual spans: measS7={float(gmax):.5f} at {gmaxE}, "
              f"margin to cap={float(ck-gmax):+.5f}")
        print()

    # ----- (3) hole-counting LOWER bound test: is sum of single-hole criticalities a
    #            valid lower bound on the JOINT loss?  (decides additive certificate) -----
    print("="*96)
    print("(3) ADDITIVE HOLE-COUNTING certificate test.")
    print("    Joint loss L(H) = measS7(full) - measS7(full minus H).")
    print("    Additive estimate L_add(H) = sum_{h in H} Delta_h(full).")
    print("    Need L(H) >= measS7(full) - cap_k  for ALL H to certify measS7(E)<=cap.")
    print("    Question: is L_add a valid LOWER bound on L (supermodular loss)?")
    print("="*96)
    for k in [8, 9, 10]:
        ck = CAP[k]
        print(f"  --- k={k}, cap_k={float(ck):.5f} ---")
        for N in regimes[k]:
            nholes = (N+1) - k
            if nholes < 1: continue
            full = list(range(N+1)); mfull = measS7(full)
            req_drop = mfull - ck          # need L(H) >= req_drop
            # additive lower bound is min over H of L_add; if even the min L_add >= req_drop,
            # AND L_add <= L (supermodular), the additive certificate closes it.
            min_Ladd = None; min_L = None; worst_gap_addvsL = None
            add_underest = True  # does L_add <= L hold everywhere?  (supermodular loss)
            worst_resid = F(-1)
            for H in itertools.combinations(range(1, N+1), nholes):
                if N in H: continue
                Ladd = sum(crit[N][h] for h in H)
                E = tuple(e for e in full if e not in set(H))
                L = mfull - measS7(E)
                if min_Ladd is None or Ladd < min_Ladd: min_Ladd = Ladd
                if min_L   is None or L    < min_L:    min_L = L
                if L < Ladd:  # additive OVER-estimates the loss -> NOT a valid lower bound
                    add_underest = False
                resid = mfull - L  # = measS7(E)
                if resid > worst_resid: worst_resid = resid
            cert_add = (min_Ladd is not None and min_Ladd >= req_drop)
            cert_true = (worst_resid <= ck)
            print(f"    N={N:2d} holes={nholes}: req_drop={float(req_drop):.5f}  "
                  f"min L_add={float(min_Ladd):.5f}  min L_true={float(min_L):.5f}")
            print(f"           L_add a valid LOWER bound on L (supermodular)? {add_underest}   "
                  f"additive cert closes? {cert_add}   true cert holds? {cert_true}")
    print()
    print("INTERPRETATION printed by the assembling agent.")
    print("DONE.")

if __name__ == "__main__":
    main()
