#!/usr/bin/env python3
"""Verify the SHEET-COUNTING argument closing the comparable-V boundary core:
12 speeds are multiples of b; per sheet n, V's safe set has >= 6b/7-1 good sheets.
kind-mendel-2026-06-22-S3."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
def gl(xs): return reduce(lambda a,b:a*b//gcd(a,b),[x for x in xs if x],1)
def nrm(y): f=y-floor(y); return min(f,1-f)

def G12_arcs(thr=F(1,14)):
    S=list(range(1,13)); bps=set([F(0),F(1)])
    for s in S:
        for n in range(s+1):
            for sg in (F(-1),F(1)):
                x=F(n,s)+sg*thr/s
                if 0<=x<=1: bps.add(x)
    bps=sorted(bps); arcs=[]; cur=None
    for a,b in zip(bps,bps[1:]):
        ok=all(nrm(s*(a+b)/2)>=thr for s in S)
        if ok:
            if cur and cur[1]==a: cur=(cur[0],b)
            else:
                if cur: arcs.append(cur)
                cur=(a,b)
        else:
            if cur: arcs.append(cur); cur=None
    if cur: arcs.append(cur)
    return arcs

arcs=G12_arcs()
arc=max(arcs,key=lambda t:t[1]-t[0])    # widest arc of G_12
ulo,uhi=arc; w=uhi-ulo
print(f"G_12 has {len(arcs)} arcs; widest = [{float(ulo):.5f},{float(uhi):.5f}] width w={float(w):.6f}")
print("near 1/13? 1/13=%.5f\n"%(1/13))

def good_sheets(b,V,thr=F(1,14)):
    "count n in 0..b-1 s.t. interval {V(n+u)/b: u in arc} (mod1) meets [thr,1-thr]"
    good=0
    for n in range(b):
        lo=V*(n+ulo); hi=V*(n+uhi)        # before /b; image length = V*w
        L=F(V*w, b)
        if L>=F(1):            # wraps fully -> meets safe
            good+=1; continue
        p=(F(V*(n)+V*ulo, b))%1            # left endpoint mod1
        q=(p+L)                            # right endpoint (may exceed 1)
        # danger arc = (-thr,thr) mod 1 = [0,thr) U (1-thr,1). interval [p,q] entirely in danger?
        # safe = [thr, 1-thr]; interval meets safe iff NOT contained in danger.
        # contained in danger: [p,q] subset of (1-thr,1+thr) (i.e. around 0)
        # normalize: shift so danger centered. Check if [p,q] (mod1, length L<1) avoids meeting [thr,1-thr].
        seg_lo=p; seg_hi=p+L
        meets_safe=False
        # sample endpoints + check overlap with [thr,1-thr] over the (possibly wrapped) segment
        pts=[seg_lo, seg_hi, F(thr), F(1)-thr]
        # robust: does [seg_lo,seg_hi] (real interval) intersect [thr,1-thr] + Z ?
        for k in (-1,0,1):
            a=max(float(seg_lo), float(thr)+k); bb=min(float(seg_hi), float(1-thr)+k)
            if bb>a: meets_safe=True; break
        if meets_safe: good+=1
    return good

print(f"{'b':>3} {'V':>4} {'gcd':>3} {'good_sheets':>11} {'6b/7-1':>7} {'>=bound?':>8} {'>=1?':>5}")
import math
for b in [2,3,5,7,10,13,20,30,49]:
    for V in [14,28,98,14*b+14]:
        if gcd(b,V)!=1: continue
        gs=good_sheets(b,V); bnd=6*b/7-1
        print(f"{b:>3} {V:>4} {gcd(b,V):>3} {gs:>11} {bnd:>7.2f} {str(gs>=math.floor(bnd) if bnd>0 else True):>8} {str(gs>=1):>5}")
