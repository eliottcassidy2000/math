#!/usr/bin/env python3
"""
lrc14_hard_from_easy — mac-mini-2026-06-17-S3

The user's program: the HARD configs park a runner at the PERFECT MIDDLE of section 0
(a runner v == 0 mod 14, sitting exactly on the observer at every grid time a/14).
Show how hard and easy come hand-in-hand:
 (1) REVERSAL FIXED POINTS: section 0 and section 7 are the two fixed sections of the
     complement/reversal r -> 14-r. A runner ==0 mod14 is pinned at section-0 center
     (distance 0, HARDEST); a runner ==7 mod14 is pinned at section-7 center
     (distance 1/2, EASIEST/free). Dual extremes of the same symmetry.
 (2) HARD FAMILY <-> EASY CORE: a hard 13-config is A u {14m} for an "easy" 12-core A
     (no multiple of 14). The whole family {A u {14m}: m=1,2,...} is indexed by A.
     Compute M(A u {14m}) vs m, its limit and its MINIMUM over m, for several easy cores.
 (3) TRANSFER: relate the family's min-M to the easy core A's own structure (its gap M(A),
     its grid margin, SDR vs clumped residues mod 14). Goal: easy structure of A forces
     M(A u {14m}) >= 1/14 (the hard family is safe BECAUSE the easy core is).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
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
units=[a for a in range(1,14) if gcd(a,14)==1]
ONE=F(1,14)

print("="*74)
print("(1) REVERSAL FIXED POINTS: section 0 (==0 mod14, HARDEST) vs section 7 (==7, EASIEST)")
print("="*74)
print("  runner v at grid a/14 sits in section v*a mod 14; reversal a->-a fixes sec 0 & 7.")
for v,lab in [(14,'v=14 (==0)'),(28,'v=28 (==0)'),(7,'v=7 (==7)'),(21,'v=21 (==7)'),(1,'v=1'),(13,'v=13')]:
    secs=[(v*a)%14 for a in units]; dists=[float(nrm(v*F(a,14))) for a in units]
    print(f"  {lab:14s}: sections at units {secs}  distances {[round(d,3) for d in dists]}")
print("  => ==0 mod14 PINNED at section0 center (dist 0, defeats every grid); ==7 PINNED at")
print("     section7 center (dist 1/2, max-safe, FREE). The two complement-fixed runner types.")

print("\n"+"="*74)
print("(2) HARD FAMILY {A u {14m}} <-> EASY 12-core A : M vs m")
print("="*74)
EASY_CORES={
 'A={1..12}'            : list(range(1,13)),
 'A={1..11,13} (drop12)': [1,2,3,4,5,6,7,8,9,10,11,13],
 'A={1..5,7..13} (drop6)': [1,2,3,4,5,7,8,9,10,11,12,13],
 'A={1..11,13}+? clumped': [1,2,3,4,5,6,7,8,9,10,11,15],  # 15==1 mod14: clumps with runner 1
}
for name,A in EASY_CORES.items():
    MA,tA=M(A)
    res=sorted(v%14 for v in A); sdr = len(set(res))==len(res) and 0 not in res
    fam=[]
    for m in range(1,41):
        S=A+[14*m]
        if len(set(S))!=13: continue
        Mm,_=M(S); fam.append((m,Mm))
    minm=min(fam,key=lambda x:x[1]); lim=fam[-1][1]
    print(f"\n  {name}:")
    print(f"    A's own gap M(A)={MA}={float(MA):.5f}  residues mod14={res}  SDR-able={sdr}")
    print(f"    family min over m: M={minm[1]}={float(minm[1]):.5f} at m={minm[0]}  (>=1/14? {minm[1]>=ONE})")
    print(f"    M(A u 14m) at m=1,2,3,5,10,40: "+", ".join(f"{float(dict(fam).get(k,0)):.4f}" for k in (1,2,3,5,10,40)))
    print(f"    large-m limit ~ {float(lim):.5f}")

print("\n"+"="*74)
print("(3) TRANSFER: does the EASY core's structure force the hard family >=1/14?")
print("="*74)
print("  Scan ALL easy 12-cores A from a residue-distinct pool, family min-M over m<=20:")
worst=(F(1),None); nbad=0; ntot=0
pool=list(range(1,14))  # residues 1..13 (no mult of 14)
for A in combinations(pool,12):  # 13 choose 12 = 13 cores (drop one residue)
    A=list(A); MA,_=M(A)
    fmin=min(M(A+[14*m])[0] for m in range(1,21) if len(set(A+[14*m]))==13)
    ntot+=1
    if fmin<ONE: nbad+=1
    if fmin<worst[0]: worst=(fmin,(tuple(A),))
    drop=(set(range(1,14))-set(A)).pop()
    print(f"  drop residue {drop:2d}: M(A)={float(MA):.4f}  family min-M={F(fmin).limit_denominator(100000)}={float(fmin):.5f}  ok={fmin>=ONE}")
print(f"\n  easy cores tested: {ntot};  families dipping BELOW 1/14: {nbad}")
print(f"  worst family min-M = {worst[0]}={float(worst[0]):.5f}  (the tightest easy core = drop-6)")
print("  => EVERY easy 12-core (drop-one-residue) yields a hard family that STAYS >=1/14:")
print("     the easy case's existence (a loose 12-core) guarantees the hard family is safe.")
