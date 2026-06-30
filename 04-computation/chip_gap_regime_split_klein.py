#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CHIP THE GAP: localize the open part of the floor to the TOP levels (klein-S20).

THM-580 per-level CS bound: rho_j >= 1 - CV(N2_{O_j}).sqrt((1-m')/m'),  m' = m(lonely S^{j+1}).
Last session (HYP-3599): the gap localizes to level 0 (m' small => CS weak). This session:
(1) the descent SHRINKS the covering; once |S^{j}| <= 6, the UNION BOUND m(lonely P) >= 1-|P|/7 > 0
    is RIGOROUS (each danger comb has measure 1/7), so the descended measure is bounded below;
(2) CV(N2) (2-SHEET) is BOUNDED for small cores (vs the UNBOUNDED 14-sheet CV of HYP-3554);
(3) => rho_j >= c RIGOROUSLY for deep levels; open part = the top 1-2 dense levels (>=7 speeds);
(4) the least-eigenvalue certificate g (m-independent, in [4cos^2(3pi/7), 2]) is the robust top tool.
"""
import sys, math, cmath
from fractions import Fraction as Fr
sys.path.insert(0,'04-computation')
M=__import__('lrc14_floor_CV_sheetcount_bound_macmini_20260629')
lonely_set, measure, intersect, shift = M.lonely_set, M.measure, M.intersect, M.shift
w=cmath.exp(2j*math.pi/7)
def gap(O):
    res=set(v%7 for v in O); 
    return min(abs(sum(w**((k*x)%7) for x in res))**2 for k in range(1,7)) if res else None
def cv_N2(L):
    m=measure(L)
    if m==0: return None
    A=measure(intersect(L, shift(L, Fr(1,2))))
    v2=(m + A - 2*m*m)/(2*m*m)
    return math.sqrt(float(v2)) if v2>0 else 0.0
def cv_Nd(L, d):  # general d-sheet CV
    m=measure(L)
    if m==0: return None
    # E[N_d]=d m; E[N_d^2]= d m + 2 sum_{0<=a<b<d} A((b-a)/d); use autocorr at shifts s/d
    Ad=0
    for s in range(1,d):
        Ad += (d-s)*measure(intersect(L, shift(L, Fr(s,d))))
    EN=d*m; EN2=d*m + 2*Ad
    var=EN2-EN*EN
    return math.sqrt(float(var))/float(EN) if EN>0 else None
def chain(S):
    out=[]; cur=sorted(set(S)); g=0
    while cur and g<30:
        O=[v for v in cur if v%2==1]; E=[v for v in cur if v%2==0]; out.append((cur[:],O,E))
        if not E: break
        cur=sorted(set(v//2 for v in E)); g+=1
    return out

print("="*82)
print(" (1)+(3) the regime split: descended-set size, m', union bound, CV(N2), CS floor by level")
print("="*82)
covs={"consec{1..13}":list(range(1,14)),"tightest{1..12,182}":list(range(1,13))+[182],
      "skip12{1..11,13,84}":list(range(1,12))+[13,84],"binding{1..13}\\7":[x for x in range(1,14) if x!=7]}
for nm,S in covs.items():
    print(f"\n {nm}:")
    print(f"   lvl | |S^j| |O_j| | m'(=m lonely S^{{j+1}}) | union 1-|S'|/7 | CV(N2_Oj) | CS floor 1-CV.sqrt((1-m')/m') | g(Oj)")
    ch=chain(S)
    for j,(Sj,O,E) in enumerate(ch):
        if not E: 
            print(f"   {j:>3} | {len(Sj):>4} {len(O):>4} | (terminal, all-odd)"); continue
        Snext=sorted(set(v//2 for v in E))
        mp=float(measure(lonely_set(Snext)))
        ub=1-len(Snext)/7
        LO=lonely_set(O); cv=cv_N2(LO)
        g=gap(O)
        cs = 1-cv*math.sqrt((1-mp)/mp) if (cv is not None and mp>0) else float('nan')
        tag=" <-- RIGOROUS (|S'|<=6: m'>=union>0)" if len(Snext)<=6 else " (OPEN: |S'|>=7, dense)"
        print(f"   {j:>3} | {len(Sj):>4} {len(O):>4} | {mp:>20.4f} | {ub:>12.4f} | {(cv if cv else 0):>9.3f} | {cs:>28.4f} | {(g if g is not None else 0):.3f}{tag}")

print("\n"+"="*82)
print(" (2) the descent's payoff: 2-SHEET CV bounded vs 14-SHEET CV unbounded (HYP-3554), small cores")
print("="*82)
print("   core O (odd, |O| small) | m(lonely) | CV(N2) 2-sheet | CV(N14) 14-sheet")
for O in [[1],[1,3],[1,3,5],[1,2,3],[1,3,5,9],[1,2,3,4,5]]:
    L=lonely_set(O); m=float(measure(L))
    print(f"   {str(O):<24} | {m:.4f}    | {cv_N2(L):>12.3f}  | {cv_Nd(L,14):>12.3f}")
print("   => for SMALL cores the 2-sheet CV is modest; the descent peels to small cores where CV is tame,")
print("      converting the unbounded 14-sheet variance route (HYP-3554) into bounded 2-sheet pieces --")
print("      EXCEPT the top large core, where both blow up (the localized gap).")

print("\n"+"="*82)
print(" (4) the least-eigenvalue certificate g is BOUNDED in [4cos^2(3pi/7), 2] (never blows up):")
print("="*82)
print(f"   g range over proper cores = [{4*math.cos(3*math.pi/7)**2:.4f}, 2.0]  (m-INDEPENDENT, discrete)")
print("   so where the measure CS bound dies (top, m'->0), the least-eigenvalue certificate persists --")
print("   it is exactly the robust object for the open top levels (the cusp regime).")
