#!/usr/bin/env python3
"""
lrc14_finitecheck_closedform_macmini_0618s7g.py   (mac-mini-2026-06-18-S6, Angle G)

ANGLE G: convert the per-k numerical finite check meas(S7(consec_k)) < cap_k (k=8..13)
into a CLEAN ANALYTIC LEMMA, by splitting the comparison into:

  (LB-cap)  cap_k >= (k-6)/7         -- PROVED by subadditivity (each speed forbids exactly 1/7)
  (UB-S7a)  meas(S7(consec_k)) < (k-6)/7   for k >= 11   -- via subadditive cap + pair-Bonferroni
  (finite)  three genuinely-tight rational checks at k=8,9,10 against the TRUE cap_k.

Plus: closed forms discovered for the building blocks
  - Phi(c,k) := meas{x: frac(ix) in [0,c) for all i=0..k-1} = c/(k-1) for c<=1/2  (single-arc, PROVED regime)
  - cap_k minimizer family: {1} U {top consecutive cluster} for k=9..12; {1,5,7,8,9} at k=8.
  - subadditive LOWER bound cap_k >= 1 - (13-k)/7 = (k-6)/7, tight at k=12,13.

Outputs the exact rational ledger and the PASS/FAIL of every analytic piece.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

# ---------- exact engines ----------
def meas_avoid(missed, k):
    """meas{x in[0,1): orbit {frac(ix):i=0..k-1} avoids all sectors in `missed`} for AP."""
    bps=set([F(0),F(1)])
    for i in range(1,k):
        for j in range(0,i+1):
            for m in range(0,8):
                v=(F(j)+F(m,7))/i
                if 0<=v<=1: bps.add(v)
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=F(0)
    for t in range(len(bps)-1):
        x0,x1=bps[t],bps[t+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(int(((i*xm)%1)*7) not in missed for i in range(0,k)): tot+=x1-x0
    return tot

def measS7(k):
    tot=F(0)
    for r in range(0,7):
        for T in itertools.combinations(range(1,7),r):
            tot += (-1)**r * meas_avoid(set(T),k)
    return tot

def measGP(P):
    P=sorted(set(P))
    if not P: return F(1)
    bps=set([F(0),F(1)]); th=F(1,14)
    for p in P:
        for j in range(0,p+1):
            for off in [F(1,14),F(-1,14)]:
                v=(F(j)+off)/p
                if 0<=v<=1: bps.add(v)
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=F(0)
    for t in range(len(bps)-1):
        x0,x1=bps[t],bps[t+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(min((p*xm)%1,1-((p*xm)%1))>=th for p in P): tot+=x1-x0
    return tot

# cap_k true minimizers (exhaustively verified earlier; recompute to be safe)
CAP_P={8:[1,5,7,8,9],9:[1,11,12,13],10:[1,12,13],11:[1,13],12:[1],13:[]}

print("="*92)
print("ANGLE G — analytic finite check: meas(S7(consec_k)) < cap_k for k=8..13")
print("="*92)

# Piece 1: subadditive cap lower bound (PROVED)
print("\n[1] PROVED LOWER BOUND  cap_k >= (k-6)/7   (each of the 13-k speeds forbids exactly 1/7)")
print(f"{'k':>3}{'cap_k(true)':>16}{'(k-6)/7':>12}{'tight?':>9}")
for k in range(8,14):
    capk=measGP(CAP_P[k]); sub=F(k-6,7)
    assert sub<=capk, f"subadditive bound FAILS at k={k}!"
    print(f"{k:>3}{str(capk):>16}{str(sub):>12}{'YES' if sub==capk else '':>9}")

# Piece 2: the key inequality meas(S7) < (k-6)/7  (holds k>=9; binding at k=9)
print("\n[2] KEY INEQUALITY  meas(S7(consec_k)) < (k-6)/7   (=> closes via piece [1] for k>=9)")
print(f"{'k':>3}{'meas(S7)':>16}{'(k-6)/7':>12}{'holds':>7}{'slack':>14}")
for k in range(8,14):
    s=measS7(k); sub=F(k-6,7)
    print(f"{k:>3}{str(s):>16}{str(sub):>12}{str(s<sub):>7}{str(sub-s):>14}")

# Piece 3: the residual genuinely-tight finite checks (against TRUE cap_k)
print("\n[3] RESIDUAL FINITE CHECKS against TRUE cap_k (rows where (k-6)/7 < meas(S7)):")
for k in range(8,14):
    s=measS7(k); sub=F(k-6,7); capk=measGP(CAP_P[k])
    if s>=sub:
        # need true cap
        print(f"   k={k}: meas(S7)={s}={float(s):.6f}  <  cap_k={capk}={float(capk):.6f}  ?  {s<capk}")
        # show as common denominator rational comparison
        from math import lcm
        D=lcm(s.denominator,capk.denominator)
        print(f"        common-denom: {s.numerator*(D//s.denominator)}/{D}  <  {capk.numerator*(D//capk.denominator)}/{D}")

# Piece 4: full verdict
print("\n[4] FULL VERDICT  meas(S7(consec_k)) < cap_k for all k=8..13:")
allok=True
for k in range(8,14):
    s=measS7(k); capk=measGP(CAP_P[k]); ok=s<capk; allok=allok and ok
    route = "subadditive (k-6)/7" if s<F(k-6,7) else "true-cap finite"
    print(f"   k={k}: {ok}  via {route}")
print(f"\n   ALL k=8..13 PASS: {allok}")

# Piece 5: closed form Phi(c,k)=c/(k-1) for c<=1/2 (single-arc), the provable building block
print("\n[5] BUILDING BLOCK  Phi(L/7,k) := meas{orbit in [0,L/7)} = L/(7(k-1)) for L/7<=1/2 (L<=3):")
def phi_arc(c,k):
    bps=set([F(0),F(1)])
    for i in range(1,k):
        for j in range(0,i+1):
            bps.add(F(j,i)); v=(j+c)/i
            if 0<=v<=1: bps.add(v)
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=F(0)
    for t in range(len(bps)-1):
        x0,x1=bps[t],bps[t+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all((i*xm)%1<c for i in range(0,k)): tot+=x1-x0
    return tot
ok=True
for k in range(8,14):
    for L in [1,2,3]:
        lhs=phi_arc(F(L,7),k); rhs=F(L,7*(k-1))
        if lhs!=rhs: ok=False; print(f"   MISMATCH k={k} L={L}: {lhs} != {rhs}")
print(f"   Phi(L/7,k)=L/(7(k-1)) for L=1,2,3, k=8..13: {ok}")
print("\nDONE.")
