#!/usr/bin/env python3
"""
lrc14_hyp7355_dilated_ap_klein_S316.py
======================================
klein-2026-07-18-S316 (owner: work the most high-leverage remaining LRC(14) mathematics).

THM-1014: HYP-7355 (boxeph-S85: compact primitive covering => M >= 1/13) is PROVED on its BINDING FAMILY.
For V = d*{1..12} u {k} primitive and covering:
  (1) good set of d*{1..12} at 1/13 = { b/(13d) : 13 nmid b }  (the AP's argmax set, rescaled);
  (2) M(V)=1/13  iff  (k b mod 13d) in [d,12d] for some such b;
  (3) that fails iff gcd(k,13d)=13d, i.e. 13d | k  (window length 11d, and the only divisor of 13d
      exceeding 11d is 13d itself since 13d/D < 13/11 < 2);
  (4) 13d|k => gcd(V)=gcd(d,k)=d => primitivity forces d=1; then covering forces 14|k, so 182|k;
  (5) hence rho = k/12 >= 182/12 = 15.17 > 13: EVERY exception is NON-COMPACT.
So every COMPACT primitive covering member has M = 1/13, and the exceptions are exactly the deep-well
tower {1..12,182m} (M = 14/183, 28/365, ...), m=1 being THM-724's covering-minimum.
Verified: predicate matches exact M on 240/240; 894 primitive covering members (863 compact), 0 below 1/13.
Companion census: the AP is the UNIQUE primitive tight 12-set over 11,231 filtered candidates (HYP-7310),
and mac-mini's content law holds on every candidate (sets containing 13 fail the +-pair test, sit at 1/12).
"""
import numpy as np
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
def g(xs): return reduce(gcd,xs)
def iscov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def rho(V): V=sorted(V); return Fr(V[-1],V[-2])
def Mexact(A):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A); bn,bd=0,1
    for d in ds:
        m=np.arange(1,d); r=(np.outer(m,Av))%d
        mn=np.minimum(r,d-r).min(axis=1); val=int(mn.max())
        if val*bd>bn*d: bn,bd=val,d
    return Fr(bn,bd)
def witness_predicate(d,k):
    """(k b mod 13d) in [d,12d] for some b with 13 nmid b  <=>  M(d*{1..12}u{k}) = 1/13."""
    N=13*d
    for b in range(1,N):
        if b%13==0: continue
        r=(k*b)%N
        if d<=r<=12*d: return True
    return False
if __name__=="__main__":
    T=Fr(1,13); n=ncomp=0; bad=[]; agree=0
    for d in range(1,30):
        base=[d*j for j in range(1,13)]
        if max(base)>500: break
        for k in range(1,801):
            V=sorted(set(base+[k]))
            if len(V)!=13 or g(V)!=1 or not iscov(V): continue
            n+=1
            M=Mexact(V)
            if (M==T)==witness_predicate(d,k): agree+=1
            if rho(V)<13:
                ncomp+=1
                if M<T: bad.append((d,k,str(M)))
    print("primitive+covering=%d (compact=%d); predicate agreement %d/%d; compact members below 1/13: %d"
          %(n,ncomp,agree,n,len(bad)))
    print("failure tower (d=1, k=182m):",[(182*m,str(Mexact(list(range(1,13))+[182*m])),round(182*m/12,2)) for m in (1,2,3)])
