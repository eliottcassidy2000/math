#!/usr/bin/env python3
"""
lrc14_second_moment_klein_S280.py
=================================
klein-2026-07-13-S280 (owner: prove the density sqrt-cancellation bound).

|S| <= (1/2pi^2) Sum_s Sum_ell |U_s(ell w)| |sin(pi ell/7)|/ell^2, U_s(N)=Sum_{R_s-endpoints} eps_p e(-N p).
Cauchy-Schwarz: Sum_ell |U_s| |sin|/ell^2 <= sqrt(Q_s)*sqrt(Sum sin^2/ell^2), Q_s = Sum_ell |U_s(ell w)|^2/ell^2.
So the sqrt-bound reduces to Q_s = O(diam) (=> |S|=O(sqrt(diam))). Large sieve predicts Q_s=O(M/delta),
M=#endpoints ~ diam, delta=min sep of {w p mod 1} ~ min||w/e'|| (O(1) for non-resonant w).

DECISIVE: compute Q_s vs diam for clean w. O(diam) => sqrt works (rigorous via large sieve).
O(diam^2) => fails (like covering averaging, mac-mini-S76). Also report |S| and the C-S bound.
"""
import math
from math import gcd
NGRID=1500000
LMAX=400
def sec(e,x): return int((e*x%1.0)*7.0)%7
def occ(E,x):
    o=0
    for e in E: o|=1<<sec(e,x)
    return o
def Rs_endpoints(E,s):
    """endpoints of R_s = {E misses exactly s}: list of (p, eps) eps=+1 arc-start, -1 arc-end."""
    eps=[]; inR=False
    for k in range(1,NGRID):
        x=k/NGRID; o=occ(E,x); miss=7-bin(o).count("1")
        cur = (miss==1) and not ((o>>s)&1)
        if cur and not inR: eps.append(((k-0.5)/NGRID,+1))
        if (not cur) and inR: eps.append(((k-0.5)/NGRID,-1))
        inR=cur
    if inR: eps.append((1.0-0.5/NGRID,-1))
    return eps
def U(endpts,N):
    re=im=0.0
    for p,e in endpts:
        a=-2*math.pi*N*p; re+=e*math.cos(a); im+=e*math.sin(a)
    return math.hypot(re,im)
def Qs(endpts,w):
    return sum(U(endpts,ell*w)**2/ell**2 for ell in range(1,LMAX+1))

SIN2=sum(math.sin(math.pi*ell/7)**2/ell**2 for ell in range(1,LMAX+1))
print("2nd moment Q_s = Sum_ell |U_s(ell w)|^2/ell^2 vs diameter (clean w=997); does Q_s ~ O(diam)?")
print("="*76)
print("  {:26s} {:>5} {:>6} {:>9} {:>9} {:>9}".format("cluster E'","diam","M_max","maxQ_s","Q/diam","|S|bound"))
w=997
fams=[
  [0,1,2,3,4,5,6],
  [0,1,2,3,4,5,12],
  [0,1,2,3,4,5,25],
  [0,1,2,4,8,16,32],
  [0,3,7,15,30,55,90],
  [0,5,13,28,54,88,140],
  [0,10,27,55,99,150,199],
]
for E in fams:
    d=max(E); maxQ=0; Mmax=0; Sbound=0.0
    for s in range(7):
        ep=Rs_endpoints(E,s); M=len(ep)
        q=Qs(ep,w); maxQ=max(maxQ,q); Mmax=max(Mmax,M)
        Sbound += math.sqrt(q)*math.sqrt(SIN2)
    Sbound/=(2*math.pi**2)
    print("  {:26s} {:5d} {:6d} {:9.2f} {:9.3f} {:9.3f}".format(str(E),d,Mmax,maxQ,maxQ/d,Sbound))
print("-"*76)
print("  Q/diam ~ const => Q_s=O(diam) => |S|=O(sqrt(diam)) [C-S], sqrt-cancellation HOLDS (rigorous via")
print("  large sieve, delta=min||w/e'||~O(1) for clean w). |S|bound = the C-S upper bound on |S|.")
print("\ndone.")
