#!/usr/bin/env python3
"""
LRC(14) DISPROVE - MINIMAL FAST breadth (disprove-relevant families, small ranges).
Avoids the high-range float-screen blowup. Covers: random sweeps in ranges<=40,
near-multiples-of-14, small-lcm subsets of {1..24}, multiples-of-base+extra,
geometric/interleaved. All exact-confirmed if float M < 1/14+eps.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def gE(S,t): return min(nrm(v*t) for v in S)
def candE(S):
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
def Mexact(S):
    b=F(0); at=None
    for t in candE(S):
        v=gE(S,t)
        if v>b: b=v; at=t
    return b,at
def Mfloat(S):
    Ss=sorted(set(S)); Cs=set([0.5])
    for v in Ss:
        k=0
        while 2*k+1<=v: Cs.add((2*k+1)/(2.0*v)); k+=1
    for i in range(len(Ss)):
        for j in range(i+1,len(Ss)):
            for d in (Ss[i]+Ss[j],Ss[j]-Ss[i]):
                if d>0:
                    k=1
                    while 2*k<=d: Cs.add(k/float(d)); k+=1
    best=0.0
    for tf in Cs:
        m=1.0
        for v in Ss:
            x=v*tf; r=x-int(x)
            if r<0: r+=1
            dd=r if r<=0.5 else 1-r
            if dd<m: m=dd
            if m<=best: break
        if m>best: best=m
    return best
def is_prim(S): return reduce(gcd,S)==1
THRf=1.0/14.0; THRESH=F(1,14)
bestf=[1.0,None]; ce=[]; seen=set()
def consider(S):
    S=tuple(sorted(set(int(x) for x in S)))
    if len(S)!=13 or any(v<=0 for v in S): return
    if not is_prim(S): return
    if S in seen: return
    seen.add(S)
    mf=Mfloat(S)
    if mf<bestf[0]: bestf[0]=mf; bestf[1]=S
    if mf<THRf+1e-9:
        m,at=Mexact(S)
        if m<THRESH: ce.append((m,S,at))

random.seed(11)
print("(f) random ranges<=40",flush=True)
for hi in [16,20,25,30,40]:
    for _ in range(80000): consider(random.sample(range(1,hi+1),13))
print(f"  bestf={bestf[0]:.6f} ce={len(ce)}",flush=True)
print("(e) near-14k sweep",flush=True)
pool=sorted(set([14*a+b for a in range(0,6) for b in range(-5,6) if 14*a+b>0]+list(range(1,28))))
for _ in range(150000): consider(random.sample(pool,13))
print(f"  bestf={bestf[0]:.6f} ce={len(ce)}",flush=True)
print("(c') small-lcm subsets of {1..24} + multiples+extra",flush=True)
for _ in range(120000): consider(random.sample(range(1,25),13))
for bm in range(1,7):
    for extra in range(1,80): consider([bm*k for k in range(1,13)]+[extra])
print(f"  bestf={bestf[0]:.6f} ce={len(ce)}",flush=True)
print("(g) geometric & interleaved APs",flush=True)
for r in range(2,6):
    geo=[r**k for k in range(0,13)]; consider(geo)
    for sh in range(0,12): consider([x+sh for x in geo[:13]])
for _ in range(40000):
    s1=random.randint(1,10); st1=random.randint(1,5); l1=random.randint(4,9)
    s2=random.randint(1,20); st2=random.randint(1,5)
    A=[s1+st1*k for k in range(l1)]; B=[s2+st2*k for k in range(13-l1)]
    consider(set(A)|set(B))
print(f"  bestf={bestf[0]:.6f} ce={len(ce)}",flush=True)
print("="*70,flush=True)
print(f"Counterexamples (exact M<1/14): {len(ce)}",flush=True)
ce.sort()
for m,S,at in ce[:25]: print(f"  M={m}={float(m):.6f} S={S} at {at}",flush=True)
print(f"\nFloat champion: {bestf[0]:.8f} (1/14={THRf:.8f})  S={bestf[1]}",flush=True)
if bestf[1]:
    mm,att=Mexact(bestf[1])
    print(f"  EXACT M of champion: {mm}={float(mm):.8f} at {att}; >=1/14? {mm>=THRESH}",flush=True)
print(f"\nVERDICT: exact M<1/14 found? {'YES - COUNTEREXAMPLE' if ce else 'NO'}",flush=True)
