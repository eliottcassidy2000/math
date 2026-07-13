#!/usr/bin/env python3
"""
lrc14_decay_bound_klein_S275.py
===============================
klein-2026-07-13-S275 (owner: work on the Sigma-e'-free decay bound |Error|<=C/w).

Error(E',w) = Phi(E' u {w}) - Phi_inf(E'), the two-scale far-element error (THM-725 object).
S274 showed the min-argument gives BOUNDEDNESS not decay. This session: is the decay bound TRUE?

PRIME grid Ng>>w (S274 lesson: Ng∝w aliases). Three tests:
 (A) DILATION counterexample: c*{0..6} at w=60c (=lcm). Error = const in c (dilation-inv) =>
     |Error| is NOT O(1/w) -- strict decay FALSE. (But dilations are IMPRIMITIVE => excluded.)
 (B) PRIMITIVE resonant growth: primitive 7-clusters of growing lcm, Error at w=lcm. O(1)/sqrt/linear?
 (C) THE ROW TARGET: |Error|<=0.097 (=> Phi<=Phi_inf+0.097<=cap9) on the moderate-w band 26<=w<7*diam,
     for primitive 7-clusters. Direct max-Phi check past the box.
"""
import math
from math import gcd
from functools import reduce
NG=999983  # large prime, >> all w and e'
def gcdall(xs): return reduce(gcd,[x for x in xs if x],0)
def lcm(xs):
    r=1
    for x in xs:
        if x: r=r*x//gcd(r,x)
    return r
def occ_of(E,x):
    o=0
    for e in E:o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def moments3(E):
    s1=s2=s3=0
    for k in range(1,NG):
        N=7-bin(occ_of(E,k/NG)).count("1"); s1+=N;s2+=N*(N-1);s3+=N*(N-1)*(N-2)
    n=NG-1;return s1/n,s2/n,s3/n
def Phi_from_m(m1,m2,m3):return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def Phi(E):m1,m2,m3=moments3(E);return Phi_from_m(m1,m2,m3)
def Phi_inf(C):m1,m2,m3=moments3(C);return Phi_from_m((6/7)*m1,(5/7)*m2,(4/7)*m3)
CAP9=1979/4004

print("="*74)
print("(A) DILATION (IMPRIMITIVE) counterexample: Error(c*{0..6}, w=60c) vs c")
print("    strict |Error|<=C/w would force Error->0; dilation-inv forces Error=const.")
print("="*74)
E0=[0,1,2,3,4,5,6]
for c in [1,2,4,8,16]:
    C=[c*e for e in E0]; w=60*c
    err=Phi(C+[w])-Phi_inf(C)
    print(f"  c={c:2d}: E'=c*{{0..6}} (sum={sum(C):3d}), w=60c={w:4d}, gcd(core)={gcdall(C+[w])}: Error={err:.5f}  Error*w={err*w:.3f}")
print("  => Error ~ const (dilation-inv), Error*w GROWS ~linearly => strict Sigma-free decay FALSE on dilations.")
print("     (all imprimitive gcd=c; EXCLUDED by the primitive reduction.)")

print()
print("="*74)
print("(B) PRIMITIVE resonant growth: primitive 7-clusters, Error at w=lcm(E'). O(1)/sqrt/linear?")
print("="*74)
prim7=[
  ("{0,1,2,3,4,5,6}",     [0,1,2,3,4,5,6]),        # lcm 60
  ("{0,1,2,3,4,5,7}",     [0,1,2,3,4,5,7]),        # lcm 420
  ("{0,1,2,4,7,10,13}",   [0,1,2,4,7,10,13]),
  ("{0,1,2,28,29,30,15}", [0,1,2,28,29,30,15]),    # 2blk, lcm 12180 (S274's ~3 case)
  ("{0,1,3,5,7,9,11}",    [0,1,3,5,7,9,11]),
  ("{0,1,2,3,5,7,11}",    [0,1,2,3,5,7,11]),        # lcm 2310
  ("{0,1,4,9,16,25,36}",  [0,1,4,9,16,25,36]),
]
print(f"  {'cluster':24s} {'sum':>5} {'lcm':>7} {'Error@lcm':>10} {'|Err|*lcm':>10} {'/sqrt(sum)':>10}")
for name,C in prim7:
    L=lcm(C)
    if L>=NG//30:
        print(f"  {name:24s} {sum(C):5d} {L:7d}  (lcm too big for grid)"); continue
    err=Phi(C+[L])-Phi_inf(C)
    print(f"  {name:24s} {sum(C):5d} {L:7d} {err:10.5f} {abs(err)*L:10.3f} {abs(err)*L/math.sqrt(sum(C)):10.3f}")

print()
print("="*74)
print("(C) ROW TARGET: max |Error| and max Phi over primitive 7-clusters, moderate-w band 26<=w<=8*diam")
print(f"    need |Error|<=cap9-Phi_inf; cap9={CAP9:.4f}. If max Phi<=cap9 with margin, row closes on the band.")
print("="*74)
band7=[
  [0,1,2,3,4,5,6],[0,1,2,3,4,5,7],[0,1,2,3,4,5,9],[0,1,2,3,4,5,13],
  [0,1,2,3,4,6,9],[0,1,2,4,6,8,10],[0,2,3,5,7,11,13],[0,1,2,3,7,8,9],
  [0,1,5,6,10,11,15],[0,1,2,3,4,8,12],
]
supPhi=0; supArg=None; supErr=0; supErrArg=None
for C in band7:
    if gcdall(C)!=1: continue
    d=max(C); pinf=Phi_inf(C)
    for w in [26,30,37,d+1,2*d,3*d,5*d,8*d]:
        if w<26 or w>=NG//30: continue
        E=C+[w]
        if gcdall(E)!=1: continue
        ph=Phi(E); err=ph-pinf
        if ph>supPhi: supPhi=ph; supArg=(C,w)
        if abs(err)>supErr: supErr=abs(err); supErrArg=(C,w,pinf)
print(f"  max Phi over band = {supPhi:.5f} at {supArg}   (cap9={CAP9:.5f}, margin +{CAP9-supPhi:.5f})")
print(f"  max |Error| over band = {supErr:.5f} at (C,w,Phi_inf)={supErrArg}")
print(f"  need |Error| <= cap9-Phi_inf; here Phi_inf(worst)~{supErrArg[2]:.3f} => budget {CAP9-supErrArg[2]:.3f}, used {supErr:.3f}")
print("\ndone.")
