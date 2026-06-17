#!/usr/bin/env python3
"""
LRC(14) PROVE route, part 5: the Dirichlet floor M >= 1/(2(n-1)) = 1/26 and the
factor-of-2 gap to the conjecture 1/14.

CLASSICAL TRIVIAL BOUND (Dirichlet pigeonhole): for n runners (n-1 nonzero
speeds), there is tau with all ||v tau|| >= 1/(2(n-1)).  For n=14: M(S) >= 1/26.
LRC(14) asks to improve 1/26 -> 1/14, exactly a factor 13/7 ~ 1.857.

Note: the EXTREMAL ratio M*2max = 13/7 found in part1 (at the tight AP) is EXACTLY
this factor.  Coincidence?  For the AP, max=13=n-1, so 1/(2max)=1/(2(n-1))=1/26
and M=1/14, ratio (1/14)/(1/26)=26/14=13/7.  So the AP is the config where the
Dirichlet floor is *tightest* (closest to being achieved), yet still beaten by
the true M.  The conjecture = "even the worst config beats Dirichlet by 13/7".

This script:
 1. Verifies M(S) >= 1/(2(n-1)) = 1/26 holds (the proven floor) on a big sample,
    and finds the config(s) closest to the floor (smallest M*26).
 2. Establishes which configs sit between 1/26 and 1/14 (the 'hard middle band')
    -- these are the configs the conjecture's proof must handle, and the only
    place a counterexample could (cannot, if LRC true) live.
 3. Quantifies the gap: smallest M observed and its ratio to 1/26 and to 1/14.

DISPROVE-seed: a counterexample has M in (1/26, 1/14) ... no -- BELOW 1/14.  Since
M>=1/26 is PROVEN, any counterexample has M in [1/26, 1/14).  We enumerate the
RATIONALS a/d in [1/26,1/14) with bounded denom (the only possible counterexample
M-values) and check exhaustively that none is realized for small max.
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

n=14
FLOOR=F(1,2*(n-1))   # 1/26
TARGET=F(1,n)        # 1/14
print(f"n={n}: Dirichlet floor 1/(2(n-1)) = {FLOOR} = {float(FLOOR):.6f}")
print(f"        conjecture target 1/n     = {TARGET} = {float(TARGET):.6f}")
print(f"        gap factor TARGET/FLOOR    = {TARGET/FLOOR} = {float(TARGET/FLOOR):.6f}")
print()

print("="*70)
print("PART 1: verify M >= 1/26 and find configs closest to the floor")
print("="*70)
random.seed(101)
samples=[list(range(1,14))]
for _ in range(2000):
    samples.append(sorted(random.sample(range(1,80),13)))
for _ in range(1000):
    samples.append(sorted(random.sample(range(1,1000),13)))
worst=None; worstS=None; below_floor=0; below_target=[]
for S in samples:
    m,_=M(S)
    r=m/FLOOR
    if worst is None or r<worst: worst=r; worstS=S[:]
    if m<FLOOR: below_floor+=1
    if m<TARGET: below_target.append((S[:],m))
print(f"samples: {len(samples)}")
print(f"M < 1/26 (floor violations): {below_floor}  (should be 0)")
print(f"smallest M*26 = {worst} = {float(worst):.5f}  at {worstS}")
print(f"M < 1/14 (would-be counterexamples): {len(below_target)}")
for S,m in below_target[:10]:
    print(f"   {S} M={m}={float(m):.6f}")
print()

print("="*70)
print("PART 2: the hard middle band [1/26, 1/14) -- rationals a counterexample")
print("        could equal, by denominator")
print("="*70)
print("M is proven >= 1/26.  A counterexample has M in [1/26, 1/14).")
print("List reduced a/d in [1/26,1/14) with smallest denominators:")
inband=[]
for d in range(14,80):
    for a in range(1,d):
        fr=F(a,d)
        if FLOOR<=fr<TARGET and gcd(a,d)==1:
            inband.append((d,fr))
inband.sort()
print(f"  reduced fractions in [1/26,1/14) with d<=79: {len(inband)}")
for d,fr in inband[:20]:
    print(f"   {fr} = {float(fr):.6f}  (denom {d})")
print()
print("Each such M-value needs denom d | (v_a +- v_b) or 2 v_i.  The smallest")
print("denominators (15,16,...) are easy to realize as denominators of envelope")
print("vertices; the question is whether the MIN over runners ever lands there.")
print()

print("="*70)
print("PART 3: exhaustive max<=22 -- min M and full middle-band occupancy")
print("="*70)
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

THRf=float(TARGET)
distinctM={}
for W in (18,20,22):
    minM=2.0; argm=None; nprim=0; below=0
    for comb in itertools.combinations(range(1,W+1),13):
        if reduce(gcd,comb,0)!=1: continue
        nprim+=1
        mf=M_fast(comb)
        if mf<minM: minM=mf; argm=comb
        if mf<THRf-1e-9:
            em,_=M(comb)
            if em<TARGET: below+=1; distinctM[em]=distinctM.get(em,0)+1
    em,_=M(argm)
    print(f"W={W}: #prim={nprim}  min M={em}={float(em):.7f}  M<1/14 count={below}")
print()
if distinctM:
    print("distinct M-values below 1/14:", {str(k):v for k,v in distinctM.items()})
else:
    print("NO primitive 13-set with max<=22 has M<1/14.  Min is 1/14 (tight AP).")
print()
print("DONE part 5.")
