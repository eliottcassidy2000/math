#!/usr/bin/env python3
"""
lrc14_two_scale_constant_controlled_klein_S273.py
=================================================
klein-2026-07-12-S273  (owner: work on the explicit two-scale constant)

Follow-up to the first run, nailing three claims:
 (1) DILATION-INVARIANCE Phi(c*E)=Phi(E) EXACTLY -> the density rows reduce to PRIMITIVE cores;
     the large-diameter Phi=0.438 configs (4*{0..7}, 5*{0..7}) are IMPRIMITIVE (= consec-8 rescaled),
     NOT new tail cases.
 (2) CONTROLLED constant: measure err*w with Ng scaled to w (Ng=400*w) so frac(wx) is always
     resolved -> kills the w~Ng aliasing artifact of run 1. Get the TRUE uniform C_Phi over
     clean (prime) AND adversarial (lcm-multiple) w.
 (3) PRIMITIVE-band max: d=26..300 over PRIMITIVE structured 8-cores (gcd=1) only -> confirm the
     genuinely-spread primitive tail stays ~0.34 << cap9 (the 0.438 was imprimitive).
"""
import math
from math import gcd
from functools import reduce

CAP9=1979/4004
def gcdall(xs): return reduce(gcd,[x for x in xs if x],0)
def lcm(xs):
    r=1
    for x in xs:
        if x: r=r*x//gcd(r,x)
    return r
def moments(E,Ng=90000):
    s1=s2=s3=0
    for k in range(1,Ng):
        x=k/Ng; occ=0
        for e in E: occ|=1<<(int((e*x%1.0)*7.0)%7)
        N=7-bin(occ).count("1"); s1+=N; s2+=N*(N-1); s3+=N*(N-1)*(N-2)
    n=Ng-1; return s1/n,s2/n,s3/n
def Phi_from_m(m1,m2,m3): return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def Phi(E,Ng=90000): m1,m2,m3=moments(E,Ng); return Phi_from_m(m1,m2,m3)
def Phi_inf(C,Ng=90000): m1,m2,m3=moments(C,Ng); return Phi_from_m((6/7)*m1,(5/7)*m2,(4/7)*m3)

print("="*74)
print("(1) DILATION-INVARIANCE  Phi(c*E) = Phi(E)  (so reduce to PRIMITIVE cores)")
print("="*74)
E0=[0,1,2,3,4,5,6,7]
base=Phi(E0)
print(f"  Phi(consec-8 {{0..7}}) = {base:.6f}")
for c in [2,3,4,5,7]:
    Ec=[c*e for e in E0]
    print(f"  Phi({c}*{{0..7}}) (diam={7*c}, gcd={gcdall(Ec)}) = {Phi(Ec):.6f}   diff={abs(Phi(Ec)-base):.2e}")
print("  => the large-diameter 0.438 configs are DILATIONS of consec-8 (imprimitive), not new.")

print()
print("="*74)
print("(2) CONTROLLED constant err*w (Ng=400*w, so frac(wx) always resolved; no aliasing)")
print("="*74)
for C in [[0,1,2,3,4,5,6],[0,1,2,3,4,5,7],[0,2,4,6,8,10,12]]:
    L=lcm(C)
    Ngc=90000
    pinf=Phi_inf(C,Ngc)
    print(f"  C={C}  lcm={L}  Phi_inf={pinf:.5f}")
    worst=0
    # clean primes and adversarial lcm-multiples, all kept << grid via Ng=400*w
    tests=[("prime",p) for p in [37,101,211,503]]+[("lcm-mult",m*L) for m in [1,2,3,5]]
    for tag,w in tests:
        Ng=max(90000,400*w)
        e=abs(Phi(C+[w],Ng)-Phi_inf(C,Ng))
        worst=max(worst,e*w)
        print(f"    {tag:9s} w={w:5d} (Ng={Ng}): err={e:.5f}  err*w={e*w:.3f}")
    print(f"    -> worst err*w for this C = {worst:.3f}")

print()
print("="*74)
print("(3) PRIMITIVE-band max Phi, d=26..300 (gcd=1 structured cands ONLY; dilations excluded)")
print("="*74)
def prim_cands(d):
    C=[]
    C.append([0,1,2,3,4,5,6,d])
    C.append([0,1,2,3,4,5,d-1,d])
    C.append([0,1,2,3,4,d-2,d-1,d])
    C.append([0,1,2,3,d-3,d-2,d-1,d])
    C.append(sorted(set([0,1]+[round(i*d/6) for i in range(1,7)])))
    C.append(sorted(set([round(i*d/7) for i in range(8)])))              # near-dilated AP (prim if 7 nmid d)
    C.append(sorted(set([0,1]+[round(i*d/6) for i in range(1,6)]+[d])))
    C.append([0,1,2,3,4,5,6,d-1])                                        # slightly off
    out=[]
    for c in C:
        if len(c)==8 and len(set(c))==8 and max(c)==d and gcdall(c)==1:
            out.append(c)
    return out
bmax=-1; barg=None
for d in [26,27,29,31,37,41,50,61,80,101,151,201,300]:   # avoid multiples of small numbers to keep prim
    cs=prim_cands(d)
    if not cs:
        print(f"  d={d:4d}: (no primitive structured cand)"); continue
    mx=max(Phi(c) for c in cs); arg=max(cs,key=Phi)
    if mx>bmax: bmax=mx; barg=arg
    print(f"  d={d:4d}: max Phi over {len(cs)} PRIMITIVE cands = {mx:.5f}   arg={arg}")
print(f"  => PRIMITIVE-spread band max = {bmax:.5f} at {barg};  cap9={CAP9:.5f}, margin +{CAP9-bmax:.5f}")
print("     (contrast: run-1 band max 0.438 was 4*{0..7}/5*{0..7} = IMPRIMITIVE = consec-8 rescaled)")
print("\ndone.")
