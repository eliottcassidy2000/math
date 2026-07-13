#!/usr/bin/env python3
"""
lrc14_Rct_cross_sector_klein_S274.py
====================================
klein-2026-07-13-S274  (owner: work on the cross-sector cancellation lemma)

The lemma (HYP-6285): |Error| = |Sum_s int_0^1 f_s(x) g_s(wx) dx| <= C0/w,
  f_s(x) = 1{E' misses exactly sector s},  g_s(y) = 1{y in [s/7,(s+1)/7)} - 1/7.

KEY REDUCTION (this session): the f_s are disjoint, and on each maximal interval I of
  R = {x : E' misses exactly one sector}, the MISSED SECTOR s_I is CONSTANT (moving to a
  different missed sector requires passing through miss-0 or miss-2, i.e. exiting R). Hence
     Error = Sum_{I in R} int_I g_{s_I}(w x) dx,   |int_I g_{s_I}(wx)dx| <= osc(G_s)/w = (6/49)/w,
  so   |Error| <= (6/49) * R_ct(E') / w,   R_ct(E') = #maximal intervals of R.
The whole lemma reduces to the GROWTH of R_ct. This script:
 (1) verifies "missed sector constant per interval" (the reduction's core claim);
 (2) measures R_ct(E') vs sum(e'), k=|E'|, diameter -- is it O(1), O(k), O(sqrt(sum)), O(sum)?
 (3) checks C0_emp = err*w against (6/49)*R_ct  (per-interval bound tight, or extra cross-interval
     cancellation?).
"""
import math
from math import gcd
from itertools import combinations

def sectors_occ(E, x):
    occ = 0
    for e in E: occ |= 1 << (int((e * x % 1.0) * 7.0) % 7)
    return occ

def miss_info(E, Ng):
    """Return (R_ct, list of (missed_sector changes within intervals -> nonconstant count)).
       R = x where exactly one of 7 sectors empty. Count maximal runs; track missed sector."""
    R_ct = 0
    nonconst = 0           # # intervals where the missed sector was NOT constant (should be 0)
    in_run = False
    run_missed = -1
    for k in range(1, Ng):
        x = k / Ng
        occ = sectors_occ(E, x)
        nfull = bin(occ).count("1")
        if nfull == 6:                       # exactly one sector empty
            missed = (~occ) & 0x7F
            ms = missed.bit_length() - 1     # the single empty sector index
            if not in_run:
                in_run = True; R_ct += 1; run_missed = ms
            elif ms != run_missed:
                nonconst += 1; run_missed = ms
        else:
            in_run = False; run_missed = -1
    return R_ct, nonconst

def moments2(E, Ng):
    s1 = s2 = 0
    for k in range(1, Ng):
        x = k / Ng; occ = sectors_occ(E, x)
        N = 7 - bin(occ).count("1"); s1 += N; s2 += N*(N-1)
    n = Ng-1; return s1/n, s2/n

def Phi_from_m(m1,m2,m3): return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def moments3(E,Ng):
    s1=s2=s3=0
    for k in range(1,Ng):
        x=k/Ng; occ=sectors_occ(E,x)
        N=7-bin(occ).count("1"); s1+=N; s2+=N*(N-1); s3+=N*(N-1)*(N-2)
    n=Ng-1; return s1/n,s2/n,s3/n
def Phi(E,Ng): m1,m2,m3=moments3(E,Ng); return Phi_from_m(m1,m2,m3)
def Phi_inf(C,Ng): m1,m2,m3=moments3(C,Ng); return Phi_from_m((6/7)*m1,(5/7)*m2,(4/7)*m3)

print("="*76)
print("(1)+(2) R_ct(E') and the 'missed sector constant per interval' check")
print("   families of varying sum(e'), k, diameter.  Ng scaled to resolve intervals.")
print("="*76)
def Rct(E):
    Ng = max(60000, 300*sum(E))
    return miss_info(E, Ng)
fams = [
  ("consec7 {0..6}",           [0,1,2,3,4,5,6]),
  ("consec8 {0..7}",           [0,1,2,3,4,5,6,7]),
  ("2*{0..6}",                 [0,2,4,6,8,10,12]),
  ("3*{0..6}",                 [0,3,6,9,12,15,18]),
  ("5*{0..6}",                 [0,5,10,15,20,25,30]),
  ("{0..6,20}",                [0,1,2,3,4,5,6,20]),
  ("{0..6,60}",                [0,1,2,3,4,5,6,60]),
  ("{0..6,200}",               [0,1,2,3,4,5,6,200]),
  ("dilated-AP7 span 70",      [0,10,20,30,40,50,60,70][:7]),
  ("spread7 {0,3,7,12,20,33,54}",[0,3,7,12,20,33,54]),
  ("consec9 {0..8}",           [0,1,2,3,4,5,6,7,8]),
  ("consec10 {0..9}",          list(range(10))),
]
print(f"  {'family':28s} {'k':>2} {'sum(e)':>6} {'diam':>5} {'R_ct':>5} {'nonconst':>8}  {'R_ct/sum':>8} {'R_ct/k':>6}")
for name,E in fams:
    rc,nc = Rct(E)
    s=sum(E); k=len(E); d=max(E)
    print(f"  {name:28s} {k:2d} {s:6d} {d:5d} {rc:5d} {nc:8d}  {rc/s:8.3f} {rc/k:6.2f}")

print()
print("="*76)
print("(3) R_ct GROWTH LAW: dilations c*{0..6} (fixed k=7, sum scales with c)")
print("="*76)
print(f"  {'c':>3} {'sum(e)':>7} {'R_ct':>5}  {'R_ct/sqrt(sum)':>14} {'R_ct/sum':>9}")
for c in [1,2,3,5,8,13,21]:
    E=[c*i for i in range(7)]
    rc,_=Rct(E); s=sum(E)
    print(f"  {c:3d} {s:7d} {rc:5d}  {rc/math.sqrt(s):14.3f} {rc/s:9.4f}")

print()
print("="*76)
print("(4) per-interval bound (6/49)*R_ct/w  vs  actual err*w  (cross-interval cancellation?)")
print("="*76)
for name,C in [("{0..6}",[0,1,2,3,4,5,6]),("2*{0..6}",[0,2,4,6,8,10,12]),("spread7",[0,3,7,12,20,33,54])]:
    Ngm=90000
    pinf=Phi_inf(C,Ngm)
    rc,_=Rct(C)
    bound_const=(6/49)*rc
    worst=0
    for w in [101,211,503,60,120,300]:
        Ng=max(90000,400*w)
        e=abs(Phi(C+[w],Ng)-Phi_inf(C,Ng)); worst=max(worst,e*w)
    print(f"  C={name:10s} R_ct={rc:3d}  (6/49)*R_ct = {bound_const:6.3f}   vs  worst err*w = {worst:.3f}   ratio bound/emp = {bound_const/max(worst,1e-9):.1f}")
print("\ndone.")
