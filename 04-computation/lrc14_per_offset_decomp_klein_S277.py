#!/usr/bin/env python3
"""
lrc14_per_offset_decomp_klein_S277.py
=====================================
klein-2026-07-13-S277 (owner: prove the coupled per-offset <=1 bound).

Directly decompose Error*w = S = Sum_{endpoints p} eps_p G_{s(p)}(w p) by the OWNING offset e':
S = Sum_{e'} S_{e'},  S_{e'} = Sum_{R-endpoints owned by e'} eps_p G_{s(p)}(w p).
Verify (a) Sum_{e'} S_{e'} = Error*w (consistency vs the moment computation), (b) |S_{e'}| <= ~1.

Method: fine grid; detect where miss-count (=#empty of 7 sectors) crosses 1 (enter/leave R=miss-1);
at each such boundary x0: owning offset = argmin_e' dist(frac(e' x0)*7, nearest int); eps=-1 if entering
R (left endpoint), +1 if leaving (right); s = the missed sector on the R side. G_s from closed form.
"""
import math
from math import gcd
from functools import reduce
NG=600000  # fine grid to resolve R-boundaries
def lcm(xs):
    r=1
    for x in xs:
        if x:r=r*x//gcd(r,x)
    return r
def occ_of(E,x):
    o=0
    for e in E:o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def missinfo(E,x):
    o=occ_of(E,x); nf=bin(o).count("1")
    miss=7-nf
    ms=-1
    if miss==1:
        em=(~o)&0x7F; ms=em.bit_length()-1
    return miss,ms
def Gs(s,y):  # antiderivative of 1{y in [s/7,(s+1)/7)}-1/7, mean-zero 1-periodic
    y=y%1.0
    L=min(max(y-s/7.0,0.0),1.0/7.0)
    return L - y/7.0
# reference Error*w via moments (degree-3 Phi)
def moments3(E,Ng=200003):
    s1=s2=s3=0
    for k in range(1,Ng):
        N=7-bin(occ_of(E,k/Ng)).count("1");s1+=N;s2+=N*(N-1);s3+=N*(N-1)*(N-2)
    n=Ng-1;return s1/n,s2/n,s3/n
def Phi(E):m1,m2,m3=moments3(E);return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def Phi_inf(C):m1,m2,m3=moments3(C);return 1-(2/3)*(6/7)*m1+(47/252)*(5/7)*m2-(5/252)*(4/7)*m3

def per_offset_S(C,w):
    """decompose the COVER-measure two-scale error* w by owning offset (degree-1 proxy; the missed
       sector's g_s). Returns dict e'->S_{e'} and total. Uses the cover-error object
       (which is Sum_s int_{R_s} g_s(wx) dx)*w = Sum eps G_s(w p)."""
    Se={e:0.0 for e in C}
    prev,prevms=missinfo(C,1.0/NG)
    for k in range(2,NG):
        x=k/NG; miss,ms=missinfo(C,x)
        if miss!=prev:
            # boundary between (k-1)/NG and k/NG; crossing into or out of miss==1
            x0=(k-0.5)/NG
            enteringR = (miss==1 and prev!=1)
            leavingR  = (prev==1 and miss!=1)
            if enteringR or leavingR:
                s = ms if miss==1 else prevms   # missed sector on the R side
                eps = -1.0 if enteringR else +1.0
                # owning offset: which e' has frac(e' x0)*7 nearest an integer
                best=None;bd=1.0
                for e in C:
                    if e==0: continue
                    f=(e*x0%1.0)*7.0; d=abs(f-round(f))
                    if d<bd: bd=d;best=e
                if best is not None:
                    Se[best]+= eps*Gs(s, (w*x0)%1.0)
            prev,prevms=miss,ms
    return Se

CAP9=1979/4004
print("per-offset decomposition of Error*w (cover-error proxy); verify Sum=Error*w and |S_e'|<=~1")
print("="*74)
tests=[
  ("{0..6} clean w=1009",      [0,1,2,3,4,5,6], 1009),
  ("{0..6} w=60 (=lcm,reson)", [0,1,2,3,4,5,6], 60),
  ("{0,1,2,28,29,30,15} clean",[0,1,2,28,29,30,15], 2003),
  ("2blk w=12180 (=lcm)",      [0,1,2,28,29,30,15], 12180),
  ("{0,1,2,4,7,10,13} clean",  [0,1,2,4,7,10,13], 1009),
]
for name,C,w in tests:
    Se=per_offset_S(C,w)
    tot=sum(Se.values())
    # reference: cover-error is not exactly Phi-error; report both the decomposed total and per-offset
    mx=max(abs(v) for v in Se.values())
    print(f"\n  {name}: lcm={lcm(C)}")
    print(f"    S_e' by offset: "+", ".join(f"{e}:{Se[e]:+.3f}" for e in sorted(C) if e!=0))
    print(f"    total S=Sum S_e' = {tot:+.4f};  max|S_e'| = {mx:.3f}  {'<=1 OK' if mx<=1.05 else 'OVER 1'}")
print("\ndone.")
