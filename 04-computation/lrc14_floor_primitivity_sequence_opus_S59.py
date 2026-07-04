#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE COVERING FLOOR IS A PRIMITIVE STATEMENT, and inf R' is the tight-locus rigidity (opus-S59).

Two facts behind the S59 global-min search:
 (1) R'(S) = meas(lonely S)/(m_R m_Q) = 0 EXACTLY on imprimitive tight families c*{1..13}: by
     scale-invariance meas(lonely(c V)) = meas(lonely V), and {1..13} is the AP with M=1/14 exactly
     (safe set = isolated points, measure 0). So the floor needs gcd=1 (kps HYP-4060, floor side).
 (2) A naive PRIMITIVE approach to the tight AP, S_c = {c+1, 2c, ..., 13c} (even c, primitive covering),
     does NOT drive R' -> 0: the coprime runner c+1 decorrelates (THM-611), and R' stabilizes ~0.53.
     So inf R' is NOT shown 0 by this route; whether the uniform floor is >0 is the open rigidity.
"""
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
BAND = Fr(1, 14)
def sa(v): return [((Fr(k)+BAND)/v, (Fr(k+1)-BAND)/v) for k in range(v)]
def inter(A, B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]);hi=min(A[i][1],B[j][1])
        if lo<hi:r.append((lo,hi))
        if A[i][1]<B[j][1]:i+=1
        else:j+=1
    return r
def meas(S):
    S=sorted(set(S)); a=sa(S[0])
    for v in S[1:]:
        a=inter(a,sa(v))
        if not a: return Fr(0)
    return sum(h-l for l,h in a)
def cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def prim(S): return reduce(gcd,S)==1
def Rp(S):
    R=[s for s in S if s%14]; Q=[s//14 for s in S if s%14==0]
    if not Q: return None
    mR,mQ=meas(R),meas(Q)
    if mR==0 or mQ==0: return (0.0, float(meas(S)))
    return float(meas(S)/(mR*mQ)), float(meas(S))

print("="*96)
print(" (1) R' = 0 on IMPRIMITIVE tight families c*{1..13}  (meas(lonely)=meas({1..13})=0)")
print("="*96)
for c in [1,2,6,14]:
    S=[c*j for j in range(1,14)]
    print(f"   c={c:>3}  S=c*[1..13]  covering={cov(S)}  primitive={prim(S)}  meas(lonely S)={float(meas(S)):.6f}")
print("   => the covering floor R'=meas/(m_R m_Q) -> 0 on these; the floor is a gcd=1 statement.")

print("\n"+"="*96)
print(" (2) primitive approach S_c = {c+1, 2c, ..., 13c}: does R' -> 0 as c grows?  (NO -- decorrelates)")
print("="*96)
print("   %4s %5s %5s %9s %10s %9s" % ('c','prim','cov',"R'",'meas(S)','maxspeed'))
for c in [2,4,6,8,10,12,16,20,30,40,60]:
    S=sorted([c+1]+[j*c for j in range(2,14)])
    r=Rp(S)
    print(f"   {c:>4} {str(prim(S)):>5} {str(cov(S)):>5} {r[0]:>9.5f} {r[1]:>10.6f} {max(S):>9}")
print("   => R' stabilizes ~0.53, NOT ->0: the coprime c+1 de-resonates (THM-611 single-coord decorrelation).")

print("\n"+"="*96)
print(" (3) the primitive global-min family (S59 descent) is near-imprimitive, small R', small magnitude")
print("="*96)
Smin=[2,4,5,6,8,10,14,16,18,20,22,24,26]
r=Rp(Smin)
print(f"   S={Smin}")
print(f"   = 2*({{1..13}}\\{{6}}) u {{5}};  covering={cov(Smin)} primitive={prim(Smin)}  R'={r[0]:.5f}  meas(S)={r[1]:.6f}")
print(f"   drop-7 {{1..6,8..13,14}} R'={Rp([1,2,3,4,5,6,8,9,10,11,12,13,14])[0]:.5f}  => drop-7 NOT global (0.076 < 0.315).")
print("   inf R' over primitive coverings: <= 0.076 (searched); >0 vs =0 is the tight-locus rigidity (HYP-4060). OPEN.")
print("DONE.")
