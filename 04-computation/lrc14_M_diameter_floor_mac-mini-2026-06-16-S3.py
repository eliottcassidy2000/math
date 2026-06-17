#!/usr/bin/env python3
"""
LRC(14) PROVE route, part 3: hunting a SCALE-COMPATIBLE forced floor on M.

The blocker: M(S) >= 1/(2 max(S)) is real but max(S) is unbounded among
primitive reps.  We need a floor whose RHS is controlled by an invariant that
IS effectively bounded for a minimizer.  Candidates to test for a forced floor:

  (F1) M(S) >= 1/(2 max(S))                      [confirmed part1]
  (F2) M(S) >= 1/(2 (max(S)-min(S)))? (diameter) -- if min>0, diam<max
  (F3) M(S) >= 1/(2 * (#distinct pairwise sums/diffs))? (combinatorial)
  (F4) M(S) >= c / (number of envelope vertices in (0,1/2])?
  (F5) two-runner local floor near each tight rational.

Also: the DISPROVE-seed reverse-engineered.  A counterexample M=a/d<1/14 forces
d>14 and d | (v_a+-v_b) or 2 v_i.  We enumerate ALL (a,d) with a/d in (0,1/14),
d<=2*MAXBOUND, and check Stern-Brocot structure: the fractions in (0,1/14) with
smallest denominators are 1/15,1/16,...,2/29,... -- list them and see what
denominators are 'cheap'.  A would-be counterexample lives at one of these and
must be realizable as min_v ||v tau*|| at an envelope vertex with that denom.

Finally: confirm the EXACT tight extremizers and the value 13/7 = M*2max bound.
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

print("="*70)
print("PART A: test candidate forced floors over wide random + structured sets")
print("="*70)
random.seed(2026)
samples=[list(range(1,14)),[1,2,3,4,5,6,7,8,9,10,11,13,24],[1,2,3,4,5,6,7,8,9,10,11,13,36]]
for _ in range(300):
    samples.append(sorted(random.sample(range(1,60),13)))
for _ in range(150):
    samples.append(sorted(random.sample(range(1,500),13)))

worst={'F1':(None,None),'F2':(None,None)}
for S in samples:
    m,_=M(S)
    mx=max(S); mn=min(S); diam=mx-mn
    r1=m*2*mx
    r2=m*2*diam if diam>0 else None
    if worst['F1'][0] is None or r1<worst['F1'][0]: worst['F1']=(r1,S)
    if r2 is not None and (worst['F2'][0] is None or r2<worst['F2'][0]): worst['F2']=(r2,S)
print(f"samples: {len(samples)}")
print(f"F1: min over samples of M*2max   = {worst['F1'][0]} = {float(worst['F1'][0]):.5f}")
print(f"     at {worst['F1'][1]}")
print(f"     -> floor M>=1/(2max) {'HOLDS' if worst['F1'][0]>=1 else 'FAILS'} (ratio>=1 means holds)")
print(f"F2: min over samples of M*2diam  = {worst['F2'][0]} = {float(worst['F2'][0]):.5f}")
print(f"     at {worst['F2'][1]}")
print(f"     -> floor M>=1/(2diam) {'HOLDS' if worst['F2'][0]>=1 else 'FAILS'}")
print()
print("Interpretation: F1 (max) holds with extremal ratio 13/7 at the tight AP.")
print("F2 (diameter) is STRICTER (diam<=max); if it also holds it's a better")
print("floor, but diam is still unbounded among primitive reps -- same blocker.")
print()

print("="*70)
print("PART B: Stern-Brocot fractions in (0,1/14) by denominator (disprove map)")
print("="*70)
print("""
A counterexample value M=a/d lies in (0,1/14).  List the fractions in (0,1/14)
with the SMALLEST denominators -- these are the 'cheapest' targets a disproof
could aim at.  For each, denom d must divide some (v_a+-v_b) or 2 v_i.
""")
THR=F(1,14)
byd=[]
for d in range(15, 60):
    # all a with 0<a/d<1/14 i.e. 1<=a<=floor((d-1)/14)
    amax=(d-1)//14
    for a in range(1, amax+1):
        fr=F(a,d)
        if fr<THR and gcd(a,d)==1:
            byd.append((fr, a, d))
byd.sort()
print(f"reduced fractions a/d in (0,1/14) with d<=59:  {len(byd)} of them")
print("the 15 LARGEST (closest to 1/14 from below):")
for fr,a,d in sorted(byd)[-15:]:
    print(f"   {a}/{d} = {float(fr):.6f}   (1/14-it = {float(THR-fr):.6f})")
print()
print("the 10 with SMALLEST denominator:")
for fr,a,d in sorted(byd,key=lambda x:(x[2],x[0]))[:10]:
    print(f"   {a}/{d} = {float(fr):.6f}")
print()

print("="*70)
print("PART C: can M actually EQUAL one of these (0,1/14) values? exhaustive max<=24")
print("="*70)
print("""
Extend the exhaustive sweep to max<=24 (includes the sporadic tight {1..11,13,24}).
Report min M and any M<1/14.  C(24,13)=2496144 -- use float screen, exact-confirm
only sub-1/14 hits and the global argmin.
""")
def M_fast(S):
    Ss=sorted(set(S)); cs=set()
    for v in Ss:
        k=0
        while 2*k+1<=v: cs.add((2*k+1)/(2.0*v)); k+=1
    n=len(Ss)
    for i in range(n):
        for j in range(i+1,n):
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

THRf=float(THR)
for W in (22,23,24):
    minM=2.0; argm=None; nprim=0; below=[]
    for comb in itertools.combinations(range(1,W+1),13):
        if reduce(gcd,comb,0)!=1: continue
        nprim+=1
        mf=M_fast(comb)
        if mf<minM: minM=mf; argm=comb
        if mf<THRf-1e-9:
            em,_=M(comb)
            if em<THR: below.append((comb,em))
    em,eat=M(argm)
    print(f"W={W}: #prim={nprim}  min M={em}={float(em):.7f} at {list(argm)}  (<1/14? {em<THR})")
    if below:
        print(f"   *** {len(below)} sub-1/14 EXACT! ***")
        for c,e in below[:10]: print(f"      {list(c)} M={e}")
print()
print("DONE part 3.")
