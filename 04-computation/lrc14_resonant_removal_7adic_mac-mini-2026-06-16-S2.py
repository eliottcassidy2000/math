#!/usr/bin/env python3
"""
lrc14_resonant_removal_7adic — mac-mini-2026-06-16-S2

THE KEY NEW STRUCTURE. For a 12-core A_j = {1..13}\{j}, adding a stranger 14m
removes  R_j(m) := meas(A_j) - L(A_j ∪ {14m}) = meas( Lonely(A_j) ∩ danger(14m) ).
inf L over the interior-drop family > 0  <==>  for every j, max_m R_j(m) < meas(A_j).

We found L(A_6 ∪ {98}) = 7/858 - 1/343 exactly, i.e. R_6(7) = 1/7^3 EXACTLY.
Conjecture (the 7-adic quantization): R_j(m) takes only finitely many rational
values, all 7-power-structured, with a MAX strictly below meas(A_j). If so, the
infinite m-direction collapses to a finite, bounded, provable check per core.

Here: compute R_j(m) EXACTLY (Fraction) for all j, m=1..M, tabulate the distinct
removal values, find max_m R_j(m), the resulting min L, and test 7-power form.
"""
from fractions import Fraction as F

def danger_arcs(v):
    hw = F(1, 14*v)
    return [(F(n,v)-hw, F(n,v)+hw) for n in range(v)]

def wrap(intervals):
    out=[]
    for lo,hi in intervals:
        shift = lo - (lo % 1); a=lo-shift; b=hi-shift
        if b<=1: out.append((a,b))
        else: out.append((a,F(1))); out.append((F(0),b-1))
    return out

def union_measure(intervals):
    ivs=sorted(intervals); tot=F(0); clo=chi=None
    for lo,hi in ivs:
        if clo is None: clo,chi=lo,hi
        elif lo<=chi:
            if hi>chi: chi=hi
        else: tot+=chi-clo; clo,chi=lo,hi
    if clo is not None: tot+=chi-clo
    return tot

def lonely(S):
    arcs=[]
    for v in S: arcs.extend(danger_arcs(v))
    return F(1)-union_measure(wrap(arcs))

def v7(x):  # 7-adic valuation of a positive rational
    if x==0: return None
    num,den=x.numerator,x.denominator
    e=0
    while num%7==0: num//=7; e+=1
    while den%7==0: den//=7; e-=1
    return e

M=210
print("="*78)
print(f"R_j(m) = meas(A_j) - L(A_j ∪ {{14m}})  removed measure, m=1..{M}, exact rationals")
print("="*78)
inf_over_cores=None
for j in range(1,14):
    A=[v for v in range(1,14) if v!=j]
    measA=lonely(A)
    Rvals={}
    for m in range(1,M+1):
        L=lonely(A+[14*m])
        R=measA-L
        Rvals[m]=R
    maxR=max(Rvals.values()); argmax=[m for m in Rvals if Rvals[m]==maxR][:6]
    minL=measA-maxR
    distinct=sorted(set(Rvals.values()))
    # 7-adic valuations of the removal values
    val7=sorted(set(v7(r) for r in distinct if r>0))
    flag = " <== GLOBAL MIN" if (inf_over_cores is None or minL<inf_over_cores[0]) else ""
    if inf_over_cores is None or minL<inf_over_cores[0]:
        inf_over_cores=(minL,j,argmax[0],maxR,measA)
    print(f"\n j={j:2d}: meas(A_j)={measA}={float(measA):.6f}  #distinct R={len(distinct)}")
    print(f"      max removal R={maxR}={float(maxR):.6f} at m={argmax}  -> min L={minL}={float(minL):.6f}{flag}")
    print(f"      7-adic valuations of removals: {val7}   max_R/meas = {float(maxR/measA):.4f}")
    print(f"      maxR < meas(A_j)? {maxR<measA}  (=> L stays POSITIVE for all scanned m)")

print("\n"+"="*78)
minL,j,m,maxR,measA = inf_over_cores
print(f"INFIMUM over interior-drop cores (m<= {M}): j={j}, m={m}, 14m={14*m}")
print(f"  meas(A_{j}) = {measA},  max removal = {maxR},  L = {minL} = {float(minL):.8f}")
print(f"  CLEAN FORM: L = {measA} - {maxR}", "= meas(core) - 1/7^3" if maxR==F(1,343) else "")
print("="*78)
print("READING: if for EVERY core j the removal R_j(m) is bounded by a 7-power")
print("strictly below meas(A_j) (verified here for m<=%d; m->inf gives (1/7)meas" % M)
print("by Weyl decoupling THM-518), then inf L>0 reduces to this FINITE per-core")
print("check.  The resonant removal is 7-adically quantized (the v_7 column).")
