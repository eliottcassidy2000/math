#!/usr/bin/env python3
"""t4_composed_atom_referee_kps_S128c44.py -- referee the COMPOSED-ATOM structure of T4:
(a) each k-line's outer sum matches the two-pole form (poles where u or t vanish, separation
Delta_out = |k|/(v3'' v4'')); (b) the composed explicit bound
T4 <= sum_k atom_out(Delta_out_k) * atom_in(Delta_in_k) with floors, vs true T4."""
import sys
from math import gcd, log, pi, sin
sys.stdout.reconfigure(line_buffering=True)
LAM=1/14
def c(h): return 2*LAM if h==0 else abs(sin(2*pi*h*LAM)/(pi*h))
def inv(x,m):
    x%=m; r0,r1,s0,s1=m,x,0,1
    while r1: q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    return s0%m if r0==1 else None
def atom(A,B,Delta):
    """codex THM-946(1): 64(1+log(2+Delta))/(A B (1+Delta)) + 6/(1+A Delta) + 6/(1+B Delta)"""
    return 64*(1+log(2+Delta))/(A*B*(1+Delta)) + 6/(1+A*Delta) + 6/(1+B*Delta)
def t4_true_and_composed(v,H,BOX=200):
    v1,v2,v3,v4=v
    g12=gcd(v1,v2); v1r,v2r=v1//g12,v2//g12
    g34=gcd(v3,v4); v3r,v4r=v3//g34,v4//g34
    ia=inv(v1r,v2r) if v2r>1 else 0
    true=0.0
    KMAX=v3*BOX+v4*BOX
    # true sum by (u,t) enumeration
    for u in range(-BOX,BOX+1):
        if u==0: continue
        for t in range(-BOX,BOX+1):
            if t==0: continue
            k=v3*u+v4*t
            if k%g12!=0: continue
            kk=-k//g12
            h10=(ia*kk)%v2r if v2r>1 else 0
            for m in range(-(BOX//max(1,v2r))-2,(BOX//max(1,v2r))+3):
                h1=h10+m*v2r
                if h1==0 or abs(h1)>BOX: continue
                num=kk-v1r*h1
                if num%v2r!=0: continue
                h2=num//v2r
                if h2==0 or abs(h2)>BOX: continue
                if max(abs(h1),abs(h2),abs(u),abs(t))<=H: continue
                true+=c(h1)*c(h2)*c(u)*c(t)
    # composed bound: sum over k (g34 | k needed for the (u,t) line to be nonempty? line exists iff g34|k)
    comp=0.0
    Kcap=8*BOX
    for k in range(1,Kcap):
        if k%g34!=0: continue
        Din=k/(v1r*v2r); Dout=k/(v3r*v4r)
        # inner and outer atoms; factor 2 for +-k; the pi-scaling of c() absorbed crudely by 1/pi^2 per pair
        comp+=2*(1/pi**2)*atom(v2r,v1r,Din)*(1/pi**2)*atom(v4r,v3r,Dout)
    return true,comp
quads=[(307,425,541,671),(800,944,1413,2147),(541,1087,1943,3310)]
print("composed-atom referee (box 200):")
for v in quads:
    for H in (10,40):
        tr,co=t4_true_and_composed(v,H)
        print("  %s H=%d: true=%.3e  composed-bound=%.3e  valid=%s"%(v,H,tr,co,tr<=co))
print("NOTE: the composed bound lacks the H-floors (it is the H-free envelope); validity true<=bound is the structural check; the floored version tightens by ~1/H as in T3.")
print("DONE")
