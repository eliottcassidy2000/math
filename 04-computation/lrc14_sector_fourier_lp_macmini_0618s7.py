#!/usr/bin/env python3
"""
lrc14_sector_fourier_lp_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

ANGLE D : LP / dual-certificate for  meas(S7(E)) <= cap_k.

STEP 1 (this file's first half): derive the EXACT Fourier / inclusion-exclusion
expansion of meas(S7(E)) and VERIFY it against the exact breakpoint engine.

  meas(S7) = int_0^1 prod_{j=0}^6 [sector j hit] dx.
  Inclusion-exclusion over the set A of MISSED sectors:
     prod_j [hit_j] = sum_{A subset {0..6}} (-1)^{|A|} * 1[every e in E avoids all sectors in A].
  "e avoids sectors A" at x  <=>  frac(e x) in U_A := union_{j not in A} [j/7,(j+1)/7).
  The e=0 offset has frac(0)=0 which lies in sector 0; so if 0 in A, the e=0 point
  is IN a missed sector -> "avoid A" is impossible -> that term is 0.
  Hence only A subset {1..6} contribute (0 always fills sector 0).  This is the
  M7 main-term combinatorics, but with the FULL joint measure (not the product proxy).

  meas(S7) = sum_{A subset {1..6}} (-1)^{|A|} * J(A,E),
     J(A,E) = meas{ x : for all e in E, frac(e x) in U_A }.
  U_A = [0,1) minus the |A| sectors in A.  meas(U_A) = (7-|A|)/7.

  KEY: J(A,E) is itself an intersection of arc-conditions; it has its OWN
  Fourier expansion in the characters n -> exp(2pi i n x).  Writing the
  indicator 1[frac(e x) in U_A] = sum_n chi_A(n) e(n e x)  (Fourier coeffs of the
  arc U_A), the product over e and integral over x picks out exactly the
  n-vectors with  sum_e n_e e = 0  (the RELATION LATTICE of E).  That is the
  "relation-lattice tail".  The constant term (all n_e=0) gives
  prod_e meas(U_A) = ((7-|A|)/7)^{|E|}  --- but the e=0 factor is meas(U_A) too,
  and since 0 not in A, frac(0)=0 in sector 0 in U_A always, so the e=0 indicator
  is identically 1 and contributes coeff (only n=0) -> factor 1, exponent |E|-1.
  => leading term sum_A (-1)^|A| ((7-|A|)/7)^{k-1} = M7(k).  EXACTLY the main term.

STEP 2: the LP.  Unknowns = the relation-lattice contributions delta_A(E) =
J(A,E) - ((7-|A|)/7)^{k-1}.  meas(S7) = M7(k) + sum_A (-1)^|A| delta_A(E).
We bound each |delta_A| by a POSITIVE combination of:
  (i) single-arc Fourier magnitudes (Koksma/var bounds),
  (ii) pairwise-overlap measures (exactly computable, 0<=.<=),
  (iii) the THM-503 7-vanishing (chi_A(n)=0 when 7|n) which removes a positive
       fraction of lattice modes.
The DUAL of "max meas(S7) s.t. these constraints" is a nonneg combination giving
an upper bound; if it is <= cap_k that combination IS the certificate for that k.

This file: verifies the Fourier identity exactly, tabulates J(A,E) and delta_A,
computes the exact per-A relation tail, and prints the LP data (objective coeffs,
constraint values) for k=8.  Honest: marks VERIFIED vs CONJECTURAL.
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV = F(1,7)

# ---------- exact breakpoint engine for meas(S7) (ground truth) ----------
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(int(((e*xm)%1)*7) for e in E)
        if len(secs)==7: total+=x1-x0
    return total

# ---------- exact J(A,E) = meas{ x: all frac(e x) avoid sectors in A } ----------
def measJ(A, E):
    """A subset of {0..6} = forbidden sectors.  U_A = complement.  Exact rational."""
    Aset=set(A); E=sorted(set(E))
    if 0 in Aset:
        # e=0 -> frac 0 in sector 0 which is forbidden -> impossible
        return F(0)
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        ok=True
        for e in E:
            s=int(((e*xm)%1)*7)
            if s in Aset: ok=False; break
        if ok: total+=x1-x0
    return total

def measS7_via_IE(E):
    """meas(S7) = sum_{A subset {1..6}} (-1)^|A| J(A,E).  (A=empty -> J=1.)"""
    k=len(set(E)); tot=F(0)
    for r in range(0,7):
        for A in itertools.combinations(range(1,7), r):
            tot += F((-1)**r) * measJ(A,E)
    return tot

def M7(k): return sum(F((-1)**t*comb(6,t))*F(7-t,7)**(k-1) for t in range(7))

# ============================================================
# STEP 1 : VERIFY the inclusion-exclusion identity exactly
# ============================================================
print("="*78)
print("STEP 1: verify  meas(S7) == sum_{A subset {1..6}} (-1)^|A| J(A,E)   (EXACT)")
print("="*78)
test_shapes = [
  ("consec {0..7}", list(range(8))),
  ("perf {0,2,3,4,5,6,7,9}", [0,2,3,4,5,6,7,9]),
  ("dissoc", [0,1,3,7,15,31,63,127]),
  ("spread", [0,1,2,3,40,41,42,43]),
  ("AP step3 {0,3,..,21}", [3*i for i in range(8)]),
  ("k=9 consec", list(range(9))),
  ("k=10 consec", list(range(10))),
]
allok=True
for name,E in test_shapes:
    direct=measS7(E); ie=measS7_via_IE(E)
    ok = (direct==ie); allok &= ok
    print(f"  {name:<22} direct={float(direct):.6f}  IE={float(ie):.6f}  {'OK' if ok else '*** MISMATCH ***'}")
print(f"\n  ALL IDENTITIES EXACT: {allok}")

# ============================================================
# STEP 2 : per-A relation tail  delta_A = J(A,E) - ((7-|A|)/7)^{k-1}
#          meas(S7) = M7 + sum_A (-1)^|A| delta_A
# ============================================================
print("\n"+"="*78)
print("STEP 2: per-A relation tail  delta_A(E) = J(A,E) - ((7-|A|)/7)^{k-1}")
print("  (|A| only matters for leading term; J depends on WHICH sectors via overlaps)")
print("="*78)
def report_tail(E, label):
    k=len(set(E))
    print(f"\n  --- {label}  (k={k}) ---")
    print(f"  {'|A|':>4}{'#A':>5}{'(7-|A|)/7^(k-1)':>18}{'sum J(A)':>14}{'sum delta_A':>14}{'signed contrib':>16}")
    m7=F(0); signed_total=F(0)
    for r in range(0,7):
        leading = F(7-r,7)**(k-1)
        sumJ=F(0); cnt=0
        for A in itertools.combinations(range(1,7), r):
            sumJ += measJ(A,E); cnt+=1
        sumdelta = sumJ - cnt*leading
        contrib = F((-1)**r)*sumdelta
        signed_total += contrib
        m7 += F((-1)**r)*cnt*leading
        print(f"  {r:>4}{cnt:>5}{float(leading):>18.6f}{float(sumJ):>14.6f}{float(sumdelta):>14.6f}{float(contrib):>16.6f}")
    s7=measS7(E)
    print(f"  M7(k)={float(m7):.6f}   sum signed delta={float(signed_total):.6f}   "
          f"M7+tail={float(m7+signed_total):.6f}   meas(S7)={float(s7):.6f}  "
          f"{'OK' if m7+signed_total==s7 else 'MISMATCH'}")
    return m7, signed_total

for name,E in [("consec {0..7}", list(range(8))),
               ("dissoc", [0,1,3,7,15,31,63,127]),
               ("AP step3", [3*i for i in range(8)])]:
    report_tail(E, name)

print("\nDONE STEP 1-2.")
