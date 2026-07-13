#!/usr/bin/env python3
"""
lrc14_per_offset_C_klein_S277.py
================================
klein-2026-07-13-S277. The S276 "per-offset <=1" is FALSE (per-offset reaches ~1.7). Pin the TRUE
per-offset absolute constant C = max_{e'} |S_{e'}| over a broad sweep, and reconfirm the ROBUST total
bound |S_cover| <= c*R (the real result). If per-offset C is bounded (~2), |S|<=C*k=O(k) still closes.

S_cover = (p0(E'u w)-Plat(E'))*w  [cover-measure two-scale error * w; robust, no attribution].
R(E',w) = sum_{e'>=2} min(1, 1/(e'||w/e'||)).  Per-offset S_{e'} via endpoint attribution.
"""
import math
from math import gcd
from functools import reduce
def lcm(xs):
    r=1
    for x in xs:
        if x:r=r*x//gcd(r,x)
    return r
def occ(E,x):
    o=0
    for e in E:o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def Gs(s,y):
    y=y%1.0; L=min(max(y-s/7.0,0.0),1.0/7.0); return L-y/7.0
def ressum(C,w):
    S=0.0
    for e in C:
        if e>=2:
            f=(w/e)%1.0;nm=min(f,1-f); S+=min(1.0,1.0/(e*nm)) if nm>0 else 1.0
    return S
def scover_and_peroffset(C,w,NG):
    # p0(E' u w), p0(E'), p1(E') on the grid; and per-offset endpoint decomposition
    c0p=c0=c1=0
    prev=None; prevms=-1; Se={e:0.0 for e in C if e!=0}
    for k in range(1,NG):
        x=k/NG; o=occ(C,x); nf=bin(o).count("1"); miss=7-nf
        ow=o|(1<<(int((w*x%1.0)*7.0)%7))
        if ow==0x7F: c0p+=1
        if nf==7: c0+=1
        elif nf==6: c1+=1
        ms=((~o)&0x7F).bit_length()-1 if miss==1 else -1
        if prev is not None and miss!=prev:
            enter=(miss==1 and prev!=1); leave=(prev==1 and miss!=1)
            if enter or leave:
                x0=(k-0.5)/NG; s=ms if miss==1 else prevms; eps=-1.0 if enter else 1.0
                best=None;bd=1.0
                for e in C:
                    if e==0:continue
                    f=(e*x0%1.0)*7.0;d=abs(f-round(f))
                    if d<bd:bd=d;best=e
                if best is not None: Se[best]+=eps*Gs(s,(w*x0)%1.0)
        prev,prevms=miss,ms
    n=NG-1; plat=c0/n+(1/7)*(c1/n); scov=(c0p/n-plat)*w
    return scov, Se

print("broad sweep: total S_cover vs R, and per-offset max C = max|S_e'|")
print("="*76)
clusters=[
  [0,1,2,3,4,5,6],[0,1,2,3,4,5,7],[0,1,2,4,7,10,13],[0,1,2,28,29,30,15],
  [0,1,3,5,7,9,11],[0,1,2,3,5,7,11],[0,2,4,6,8,10,12],[0,1,2,3,4,10,20],
  [0,1,5,6,10,11,15],[0,3,6,9,12,15,18],
]
supC=0; supArg=None; supRatio=0
for C in clusters:
    L=lcm(C); NG=max(200000, 800*max(C))
    for w in [1009, 2003, L if L<40000 else 2*max(C)+1]:
        if w>=NG//8: continue
        scov,Se=scover_and_peroffset(C,w,NG)
        R=ressum(C,w); mx=max(abs(v) for v in Se.values()); arg=max(Se,key=lambda e:abs(Se[e]))
        ratio=abs(scov)/max(R,1e-9)
        supRatio=max(supRatio,ratio)
        if mx>supC: supC=mx; supArg=(C,w,arg)
        tag="reson" if w==L else "clean"
        print(f"  C={str(C):24s} w={w:6d}({tag}) |S_cov|={abs(scov):.3f} R={R:.2f} |S|/R={ratio:.3f} maxС_e'={mx:.3f}@e'={arg}")
print("-"*76)
print(f"  PER-OFFSET max C = {supC:.3f} at (C,w,e')={supArg}")
print(f"  TOTAL |S_cover|/R max ratio = {supRatio:.3f}  (S276 claimed <=0.61 -- check)")
print(f"  => per-offset bounded by ~{supC:.1f} (NOT <=1); |S|<={supC:.1f}*k=O(k) closes row (bigger box).")
print("\ndone.")
