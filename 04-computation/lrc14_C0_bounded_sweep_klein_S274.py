#!/usr/bin/env python3
"""
lrc14_C0_bounded_sweep_klein_S274.py
====================================
klein-2026-07-13-S274. DECISIVE test: is the two-scale constant C0 = err*w BOUNDED (uniform over
7-cluster shape AND diameter AND far-w), or does it grow?

Run-2 showed err*w ~constant within a shape-family as the scale grows (dilation ~0.25, {0..5,D} ~1.1).
This sweeps MANY shapes x growing diameter x clean+adversarial w to find sup C0. If sup C0 <= ~1.5,
the k=8 tail closes: d>25 => err <= C0/25 <= 0.06 => Phi <= Phi_inf_max(0.397)+0.06 = 0.457 < cap9.

Grid: Ng = max(120000, 150*w, 60*sum(C)) keeps frac(wx) resolved (no aliasing) and C resolved.
"""
import math
from math import gcd
def lcm(xs):
    r=1
    for x in xs:
        if x: r=r*x//gcd(r,x)
    return r
def sectors_occ(E,x):
    occ=0
    for e in E: occ|=1<<(int((e*x%1.0)*7.0)%7)
    return occ
def moments3(E,Ng):
    s1=s2=s3=0
    for k in range(1,Ng):
        x=k/Ng; occ=sectors_occ(E,x); N=7-bin(occ).count("1")
        s1+=N;s2+=N*(N-1);s3+=N*(N-1)*(N-2)
    n=Ng-1;return s1/n,s2/n,s3/n
def Phi_from_m(m1,m2,m3):return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def Phi(E,Ng):m1,m2,m3=moments3(E,Ng);return Phi_from_m(m1,m2,m3)
def Phi_inf(C,Ng):m1,m2,m3=moments3(C,Ng);return Phi_from_m((6/7)*m1,(5/7)*m2,(4/7)*m3)
CAP9=1979/4004

def shapes(D):
    """7-clusters of diameter D, varied shapes (all contain 0 and D)."""
    S=[]
    S.append([0,1,2,3,4,5,D])                     # compact-6 + far
    S.append([0,1,2,D-2,D-1,D,D//2])              # 2-block + mid
    S.append(sorted(set([round(i*D/6) for i in range(7)])))   # dilated-AP (equal spacing)
    S.append([0,1,2,3,D-1,D, (D+1)//2])           # 4 + 2 + 1
    S.append([0,3,7,12,20,D, (2*D)//3])           # irregular spread
    return [s for s in S if len(s)==7 and len(set(s))==7 and max(s)==D and min(s)==0]

print("="*72)
print("DECISIVE: sup C0 = err*w over shapes x diameter x (clean + adversarial w)")
print(f"cap9={CAP9:.5f}, Phi_inf_max~0.397, tail needs C0 <= 25*(cap9-0.397)={25*(CAP9-0.397):.3f}")
print("="*72)
supC0=0; suparg=None
for D in [30,60,120,240]:
    for C in shapes(D):
        L=lcm(C)
        # clean primes and adversarial lcm-multiples, all >> D, kept << grid
        base=[p for p in [1009,2003,4001] if p>20*D]
        adv=[m*L for m in [1,2] if 20*D < m*L < 200000]
        ws=(base+adv) or [4001]
        worst=0; wworst=None
        for w in ws:
            Ng=max(120000, 150*w, 60*sum(C))
            e=abs(Phi(C+[w],Ng)-Phi_inf(C,Ng))*w
            if e>worst: worst=e; wworst=w
        tag="COMPACT" if C[:6]==[0,1,2,3,4,5] else "spread"
        if worst>supC0: supC0=worst; suparg=(C,wworst)
        print(f"  D={D:4d} C={str(C):34s} sum={sum(C):4d}  sup err*w={worst:.3f}  (w={wworst})")
print("-"*72)
print(f"  SUP C0 over the whole sweep = {supC0:.3f}  at {suparg}")
need=25*(CAP9-0.397)
print(f"  tail-closure needs C0 <= {need:.3f}.  {'CLOSES (bounded, margin '+format(need-supC0,'.3f')+')' if supC0<need else 'does NOT close with this C0'}")
print()
print("  CLOSURE (if C0<=sup): any 8-core E, diam d=max(E)>25 => peel w=d, C=E\\{w} (7 elts),")
print(f"    err <= C0/d < {supC0:.2f}/25 = {supC0/25:.3f}, Phi(E) <= Phi_inf(C)+err <= 0.397+{supC0/25:.3f} = {0.397+supC0/25:.3f} < cap9={CAP9:.3f}.")
print("\ndone.")
