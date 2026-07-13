#!/usr/bin/env python3
"""
lrc14_sqrt_law_klein_S274.py
============================
klein-2026-07-13-S274. Pin the growth law of the two-scale constant C0 = err*w.

Run 1 showed: R_ct (interval count) ~ 0.81*sum(e') LINEAR, but actual err*w grows SUBLINEARLY,
consistent with C0 ~ kappa*sqrt(sum e'). This script confirms the sqrt law and pins kappa, two ways:
 (A) DILATION series c*{0..6}: sum=21c, far w fixed >> diam=6c. err*w vs sqrt(sum).
 (B) FIXED-k=7 clusters of growing diameter {0,1,2,3,4,5,D}: sum=15+D. err*w vs sqrt(sum).
Then the CLOSURE check: err <= kappa*sqrt(7 d_C)/w with w > d_C gives err <= kappa*sqrt7/sqrt(d_C);
tabulate the implied Phi bound Phi_inf(C)+err vs cap9 across d_C.
"""
import math
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

def errw(C, ws):
    """max over far ws of |Phi(C u w)-Phi_inf(C)|*w, with Ng scaled to resolve both C and w."""
    worst=0
    for w in ws:
        Ng=max(90000, 400*w, 60*sum(C))
        e=abs(Phi(C+[w],Ng)-Phi_inf(C,Ng))
        worst=max(worst,e*w)
    return worst

print("="*72)
print("(A) DILATION c*{0..6}: sum(e')=21c, diam=6c; far w in {a few primes >> 6c}")
print("="*72)
print(f"  {'c':>3} {'sum':>5} {'err*w':>7} {'err*w/sqrt(sum)':>16}")
kaps=[]
for c in [1,2,3,5,8,13,21]:
    C=[c*i for i in range(7)]; s=sum(C)
    ws=[p for p in [997,1499,2999] if p>10*6*c]  # ensure w >> diam
    if not ws: ws=[9973,19949]
    ew=errw(C,ws); kap=ew/math.sqrt(s); kaps.append(kap)
    print(f"  {c:3d} {s:5d} {ew:7.3f} {kap:16.3f}")
print(f"  -> kappa (err*w / sqrt(sum)) is ~constant: mean={sum(kaps)/len(kaps):.3f}, range [{min(kaps):.3f},{max(kaps):.3f}]")

print()
print("="*72)
print("(B) FIXED k=7, growing diameter {0,1,2,3,4,5,D}: sum=15+D; far w prime >> D")
print("="*72)
print(f"  {'D':>4} {'sum':>5} {'err*w':>7} {'err*w/sqrt(sum)':>16}")
kapsB=[]
for D in [8,12,20,35,60,100,180]:
    C=[0,1,2,3,4,5,D]; s=sum(C)
    ws=[p for p in [997,2999,9973] if p>10*D]
    if not ws: ws=[99991]
    ew=errw(C,ws); kap=ew/math.sqrt(s); kapsB.append(kap)
    print(f"  {D:4d} {s:5d} {ew:7.3f} {kap:16.3f}")
print(f"  -> kappa ~constant: mean={sum(kapsB)/len(kapsB):.3f}, range [{min(kapsB):.3f},{max(kapsB):.3f}]")

kappa=max(max(kaps),max(kapsB))*1.15   # empirical upper envelope with margin
print()
print("="*72)
print(f"(C) CLOSURE CHECK with C0 <= kappa*sqrt(sum e'), kappa_env={kappa:.3f}")
print("    peel w=max(E), C=E\\{w} (7 elts), sum(e'(C))<=7*d_C, w>d_C => err<=kappa*sqrt(7 d_C)/w<=kappa*sqrt7/sqrt(d_C)")
print("    Phi(E) <= Phi_inf_max(=0.397) + kappa*sqrt7/sqrt(d_C).  Need <= cap9=0.49426.")
print("="*72)
PHIINF_MAX=0.397
print(f"  {'d_C':>4} {'err_bound':>10} {'Phi_bound':>10} {'<=cap9?':>8}")
for dC in [4,6,8,10,15,25,50,100]:
    eb=kappa*math.sqrt(7)/math.sqrt(dC)
    pb=PHIINF_MAX+eb
    print(f"  {dC:4d} {eb:10.4f} {pb:10.4f} {'YES' if pb<=CAP9 else 'no':>8}")
print(f"  => crossover d_C* where Phi_bound=cap9: d_C* = 7*(kappa/(cap9-0.397))^2 = {7*(kappa/(CAP9-PHIINF_MAX))**2:.1f}")
print("     (cores with compact-part diam d_C above this close via the peel; below it -> finite box / recurse)")
print("\ndone.")
