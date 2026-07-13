#!/usr/bin/env python3
"""
lrc14_delta_sep_klein_S280.py
=============================
klein-2026-07-13-S280. The large-sieve bound Q_s <= M_s(2 + (4/3)delta^-1), delta = min_{p!=p'} ||w(p-p')||
over R_s-endpoints. If delta = Omega(1) for non-resonant w, the sqrt-bound Q_s=O(M)=O(diam) is RIGOROUS.
Check delta empirically for clean w (and note behavior for resonant w).
Also report M_s (#endpoints) vs diam, and the implied rigorous Q_s bound vs the measured Q_s.
"""
import math
from math import gcd
NGRID=1500000
def sec(e,x): return int((e*x%1.0)*7.0)%7
def occ(E,x):
    o=0
    for e in E: o|=1<<sec(e,x)
    return o
def Rs_endpoints(E,s):
    pts=[]; inR=False
    for k in range(1,NGRID):
        x=k/NGRID; o=occ(E,x); miss=7-bin(o).count("1")
        cur=(miss==1) and not ((o>>s)&1)
        if cur!=inR: pts.append((k-0.5)/NGRID)
        inR=cur
    if inR: pts.append(1.0-0.5/NGRID)
    return pts
def frac(x): return x-math.floor(x)
def delta_min(pts,w):
    # min_{p!=p'} ||w(p-p')||; pts sorted
    P=sorted(pts); d=1.0
    for i in range(len(P)):
        for j in range(i+1,len(P)):
            v=min(frac(w*(P[j]-P[i])),1-frac(w*(P[j]-P[i])))
            if 0<v<d: d=v
    return d

print("delta = min_{p!=p'} ||w(p-p')|| over R_s-endpoints (clean w=997); is delta=Omega(1)?")
print("="*72)
print("  {:26s} {:>5} {:>6} {:>9} {:>18}".format("E'","diam","M_max","min_s delta","large-sieve Q<=M(2+4/3delta^-1)"))
w=997
for E in [[0,1,2,3,4,5,6],[0,1,2,3,4,5,25],[0,3,7,15,30,55,90],[0,10,27,55,99,150,199]]:
    dmin=1.0; Mmax=0; qb=0
    for s in range(7):
        pts=Rs_endpoints(E,s); M=len(pts)
        if M>=2:
            dl=delta_min(pts,w); dmin=min(dmin,dl); Mmax=max(Mmax,M)
    qb=Mmax*(2+ (4/3)/dmin) if dmin>0 else float('inf')
    print("  {:26s} {:5d} {:6d} {:9.4f} {:18.1f}".format(str(E),max(E),Mmax,dmin,qb))
print("-"*72)
print("  If delta bounded below (~O(1)): large sieve gives Q_s=O(M)=O(diam) RIGOROUSLY => sqrt-bound.")
print("  (compare large-sieve Q-bound to the measured maxQ_s ~1.5*diam from the 2nd-moment run.)")
print("\ndone.")
