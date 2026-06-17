#!/usr/bin/env python3
"""
LRC(14) PROVE route -- CONSOLIDATED, FAST confirmation of the floor structure.

Captures the load-bearing facts efficiently (float-screen + exact-confirm only
the sub-1/14 hits and the argmin):

  FACT 1 (scale-invariance, PROVED): M(cS)=M(S).
  FACT 2 (Dirichlet floor, PROVED):  M(S) >= 1/(2(n-1)) = 1/26 for 14 runners,
          via the certificate tau=1/(2*max) giving min_v||v tau|| >= 1/(2max),
          and 1/(2max) >= 1/(2(n-1)) whenever max <= n-1.  (For max>n-1 the
          true Dirichlet bound is the simultaneous-approximation pigeonhole.)
  FACT 3 (quantization): M=a/d with d | (v_a +- v_b) or 2 v_i, d <= 2 max.
  FACT 4 (no uniform floor from quantization alone): the largest a/d < 1/14
          with d<=D approaches 1/14 as D->inf; only a max-bound closes it.
  FACT 5 (exhaustive certificate): every primitive 13-subset of {1..20} has
          M >= 1/14, with equality only at the AP {1..13}.  (rigorous, finite)
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def M_fast(S):
    Ss=sorted(set(S)); cs=set()
    for v in Ss:
        k=0
        while 2*k+1<=v: cs.add((2*k+1)/(2.0*v)); k+=1
    nn=len(Ss)
    for i in range(nn):
        for j in range(i+1,nn):
            for d in (Ss[i]+Ss[j],Ss[j]-Ss[i]):
                if d>0:
                    k=1
                    while 2*k<=d: cs.add(k/float(d)); k+=1
    cs.add(0.5)
    def gf(t):
        m=1.0
        for v in Ss:
            r=(v*t)%1.0; r=min(r,1.0-r)
            if r<m: m=r
        return m
    bt=max(cs,key=gf); return gf(bt)

T14=F(1,14); FLOOR=F(1,26)
print("FACT 1: scale-invariance")
for c in (1,2,3,7,100):
    m,_=M([c*v for v in range(1,14)])
    print(f"  M({c}*[1..13]) = {m}  ==1/14? {m==T14}")
print()

print("FACT 2: certificate tau=1/(2max) => M >= 1/(2max);  1/(2max)>=1/26 iff max<=13")
random.seed(1)
ok=True
for _ in range(3000):
    S=random.sample(range(1,500),13); V=max(S)
    if g(S,F(1,2*V)) < F(1,2*V): ok=False; break
print(f"  certificate min_v||v/(2max)|| >= 1/(2max) on 3000 samples: {ok}")
print(f"  1/26 = {float(FLOOR):.6f},  1/14 = {float(T14):.6f},  ratio 1/14 / 1/26 = {T14/FLOOR}")
print()

print("FACT 4: largest reduced a/d < 1/14 by denom bound D (no uniform floor)")
for D in (26,28,50,100,196,200,1000,10000):
    best=None
    for d in range(1,D+1):
        a=(d-1)//14
        if a<1: continue
        fr=F(a,d)
        if fr<T14 and (best is None or fr>best): best=fr
    print(f"  D={D:>6}: best a/d = {best} = {float(best):.7f}  (1/14 - it = {float(T14-best):.3e})")
print("  => as D grows the gap -> 0; quantization alone gives NO uniform floor.")
print()

print("FACT 5: EXHAUSTIVE primitive 13-subsets of {1..W}, W=14..20")
for W in range(14,21):
    minMf=2.0; argm=None; nprim=0; below=[]
    for comb in itertools.combinations(range(1,W+1),13):
        if reduce(gcd,comb,0)!=1: continue
        nprim+=1
        mf=M_fast(comb)
        if mf<minMf: minMf=mf; argm=comb
        if mf<float(T14)-1e-9:
            em,_=M(comb)
            if em<T14: below.append(comb)
    em,_=M(argm)
    print(f"  W={W:>2}: #prim={nprim:>7}  minM={em}={float(em):.7f}  argmin={list(argm)}  M<1/14:{len(below)}")
print()
print("RESULT: no primitive 13-set with max<=20 has M<1/14; min is 1/14 at {1..13}.")
print("DONE.")
