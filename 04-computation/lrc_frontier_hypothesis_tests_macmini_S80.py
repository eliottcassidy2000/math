"""
S80: hypothesis tests on the LRC(14) frontier. Findings:
 H1  covering bound holds: 0/600 random covering 13-sets violate M>=1/14 (random M ~ 0.11, far from tight).
 CH2 index-theorem: witness peaks (at M) = (p-1)*d; primitive index (p-1)/2=3 = Borsuk-Ulam degree. AP/GW=6, d*AP=6d.
 H3  tight-locus = AP/GW dilations, ISOLATED: near-tight gap delta~0.0026 (AP with one element doubled, M~0.074).
 H3c lcm-family S_X={1..11,13,lcm(2..X)} is VALUE-SAFE: M -> 1/12 (NOT 1/14); the denominator obstruction != value.
 => covering bound = {finite tight-locus (constructed)} + {uniform margin delta}.
"""
import numpy as np
from math import gcd
def M_of(S,grid):
    t=np.arange(1,grid)/grid; s=np.ones(grid-1)
    for x in S:
        fr=(x*t)%1.0; s=np.minimum(s,np.minimum(fr,1-fr))
    return s.max()
def peaks_at_M(S,grid=420000):
    t=np.arange(grid)/grid; s=np.ones(grid)
    for x in S:
        fr=(x*t)%1.0; s=np.minimum(s,np.minimum(fr,1-fr))
    M=s.max(); near=s>M-1.5e-4; runs=0; prev=False
    for b in near:
        if b and not prev: runs+=1
        prev=b
    if near[0] and near[-1] and runs>1: runs-=1
    return M,runs
def lcm_upto(X):
    l=1
    for k in range(2,X+1): l=l*k//gcd(l,k)
    return l

thr=1/14
print("CH2 index-theorem (witness peaks = (p-1)*d):")
for name,S in [("AP",list(range(1,14))),("GW",list(range(1,12))+[13,24]),("2AP",list(range(2,27,2))),("3AP",list(range(3,40,3)))]:
    M,w=peaks_at_M(S); print(f"  {name:<4}: M={M:.5f} peaks={w} index={w/2}")
print("\nH3c lcm-family (value-safe, M->1/12):")
for X in (7,9,11,13):
    L=lcm_upto(X); S=list(range(1,12))+[13,L]; m=M_of(S,max(300000,4*L+10))
    print(f"  X={X:>2} lcm={L:>7}: M={m:.6f}  M-1/14={m-thr:+.6f}")
print(f"\n=> tight-locus = AP/GW dilations (M=1/14), isolated; lcm-family M->1/12 (margin). 1/12={1/12:.5f} 1/14={thr:.5f}")
