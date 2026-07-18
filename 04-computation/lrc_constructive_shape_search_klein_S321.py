#!/usr/bin/env python3
"""
lrc_constructive_shape_search_klein_S321.py
===========================================
klein-2026-07-18-S321 (owner: build the constructive search targeting that shape, and other shapes).

THE INSTRUMENT (witness-first, exhaustive-in-pool). Instead of sampling sets and computing M:
  1. pick (val,q,a) with val/q < 1/(n-1), gcd(val,q)=gcd(a,q)=1;
  2. pool = { v <= B : min((a v) mod q, q-(a v mod q)) >= val }  -- every pool member is SAFE at t=a/q,
     so every candidate has M >= val/q BY CONSTRUCTION (the low-M region random search never reaches);
  3. anchor on the maximal prime powers <= n (covering forces a multiple of each), enumerate fills
     EXHAUSTIVELY inside the pool;
  4. keep primitive + covering(d=2..n) + compact(rho<n-1), test with an early-exit "M >= 1/(n-1)?".
Power comes from step 3 being exhaustive, NOT from the pool being small (val/q < 1/(n-1) forces the pool
to be >= 1 - 2/(n-1) of the box).

POSITIVE CONTROL (mandatory, per MISTAKE-162): at n=9 it rediscovers {1,3,4,5,7,11,18,32}, M=4/33 < 1/8,
via (val,q,a)=(1,10,1). Random anchored sampling had missed this with 191,166 samples at a LARGER box.

NEW COUNTEREXAMPLES FOUND (the analog `compact primitive covering => M >= 1/(n-1)` FAILS):
  n=10 : {1,2,3,5,6,7,8,9,30}            M = 4/37 = 0.108108 < 1/9   rho=3.33  (found after 4 sets, 57s)
         -- my earlier exhaustive to box 28 scanned 999,071 sets and MISSED it: its max is 30.
  n=11 : {2,6,8,9,10,11,13,14,17,19}     M = 3/31 = 0.096774 < 1/10  rho=1.12  (349,887 sets, 18s)
         -- a NEW SHAPE: a TIGHT BAND (all speeds in [2,19]), not body-plus-outlier.
So the analog now fails at n = 5,6,7,8,9,10,11 consecutively.

SILENT (time-capped, NOT exhaustive -- silence here is weak, per MISTAKE-162):
  n=12: pool/anchor setup consumed the budget, 0 sets reached.
  n=13: 937,288 sets scanned, none.
  n=14: 2,383,905 compact covering 13-sets scanned in the TIGHT-BAND shape, none  <-- this is HYP-7355.
ASSESSMENT: HYP-7355 would have to be an EXCEPTION to an unbroken run of failures at n=5..11. That makes
it materially more precarious than the fleet's evidence suggested. Not refuted; not supported either.
"""
import numpy as np, itertools, time
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
def g(xs): return reduce(gcd,xs)
def M_ge(A,num,den):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A)
    for dd in ds:
        m=np.arange(1,dd); r=(np.outer(m,Av))%dd
        mn=np.minimum(r,dd-r).min(axis=1)
        if (mn*den >= num*dd).any(): return True
    return False
def Mexact(A):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A); bn,bd=0,1
    for dd in ds:
        m=np.arange(1,dd); r=(np.outer(m,Av))%dd
        mn=np.minimum(r,dd-r).min(axis=1); v=int(mn.max())
        if v*bd>bn*dd: bn,bd=v,dd
    return Fr(bn,bd)
def mpp(n):
    out=[]
    for p in (2,3,5,7,11,13):
        q=p
        while q*p<=n: q*=p
        if q<=n: out.append(q)
    return out
def hunt(n,B,QMAX,budget=250):
    ns=n-1; thr=Fr(1,n-1); pps=mpp(n); t0=time.time(); scanned=0
    def iscov(S): return all(any(s%d==0 for s in S) for d in range(2,n+1))
    for q in range(3,QMAX+1):
several = None
if __name__=="__main__":
    print("control n=9 -> expect {1,3,4,5,7,11,18,32}; see .out files for the full runs")
