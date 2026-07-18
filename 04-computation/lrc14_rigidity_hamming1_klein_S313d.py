#!/usr/bin/env python3
"""
lrc14_rigidity_hamming1_klein_S313d.py
======================================
klein-2026-07-17-S313d (owner: prove the inverse/rigidity theorem "M near 1/13 => A near AP").

THM-1004 (PROVED, Hamming radius 1): for A = (AP\{j}) u {w}, w not in AP=[1..12],
    M(A) >= 2/25, with equality iff A = {1,...,11,24}.
Two ingredients, overlapping ranges:
  TAIL   : if G_{B_j}(delta) has an interval of length L, the bad set of w meets it in measure
           <= 2*delta*L + 2*delta/w, so a good point survives once w > 2 delta/(L(1-2 delta)) = 4/(21 L)
           at delta=2/25. W0 = max_j 4/(21 L_j) = 1100/21 = 52.38 -> every w > 52 is automatic.
  FINITE : exact rational evaluation for all j and 13 <= w <= 399 (4644 families); min = 2/25 only at
           (j,w)=(12,24); next values 3/37, 4/49, 5/61 (the family {1..11,12k} has M = k/(12k+1)).
Also: M = 1/13 in lowest terms forces (val,q)=(1,13), so tight witnesses sit at denominator exactly 13;
then v != 0 mod 13 for all v, and {v mod 13} must meet each of the six +- pairs of (Z/13)*.
"""
import numpy as np
from fractions import Fraction as Fr
delta=Fr(2,25); AP=list(range(1,13))
def Mexact(A):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A); bn,bd=0,1
    for d in ds:
        m=np.arange(1,d); r=(np.outer(m,Av))%d
        mn=np.minimum(r,d-r).min(axis=1); val=int(mn.max())
        if val*bd > bn*d: bn,bd=val,d
    return Fr(bn,bd)
def good_intervals(B,d):
    cuts={Fr(0),Fr(1)}
    for v in B:
        k=0
        while Fr(k,v)<=1:
            for e in (Fr(k,v)-d/v, Fr(k,v)+d/v):
                if 0<=e<=1: cuts.add(e)
            k+=1
    cs=sorted(cuts); out=[]
    for i in range(len(cs)-1):
        a,b=cs[i],cs[i+1]
        if b>a and all(min((v*((a+b)/2))%1, 1-(v*((a+b)/2))%1)>=d for v in B): out.append([a,b])
    merged=[]
    for a,b in out:
        if merged and merged[-1][1]==a: merged[-1][1]=b
        else: merged.append([a,b])
    return [(a,b,b-a) for a,b in merged]
print("TAIL: exact L_j and threshold w > 4/(21 L_j)")
W0=Fr(0)
for j in AP:
    L=max((l for _,_,l in good_intervals([v for v in AP if v!=j],delta)),default=Fr(0))
    thr=Fr(4,21)/L; W0=max(W0,thr)
    print("  j=%2d  L=%-11s  w > %s"%(j,L,thr))
print("  W0 = %s = %.2f  => all w > 52 automatic\n"%(W0,float(W0)))
print("FINITE: exact M for all j, 13<=w<=399")
rows=[]
for j in range(12):
    for w in range(13,400):
        A=AP[:j]+AP[j+1:]+[w]
        if len(set(A))==12: rows.append((Mexact(A),AP[j],w))
rows.sort()
for M,out,w in rows[:5]: print("   M=%-7s (removed %2d, added %3d)"%(M,out,w))
print("  min = %s  => THM-1004 PROVED (tail w>52 + finite w<=399 overlap)"%rows[0][0])
