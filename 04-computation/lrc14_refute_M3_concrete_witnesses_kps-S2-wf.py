#!/usr/bin/env python3
"""
Confirm the M3 refutation with CONCRETE PRIMITIVE INTEGER speed sets (real LRC inputs),
running the ORIGINAL build_tournament_M3 (slowest = reference, exact Fractions).

Residue-witnesses found in Phase 3 (ordered by speed, pos0 = slowest = reference):
  n=4: residues (2,2,4,6) -> realizes forbidden (H=5,#3cyc=2,score(1,1,2,2)) [regular near-3-cycle]
  n=5: residues (2,2,2,4,6) -> realizes forbidden (H=9,#3cyc=3,score(1,1,2,3,3))
       residues (2,2,4,4,6) -> realizes forbidden (H=11,#3cyc=4,score(1,2,2,2,3))
We pick actual integer speeds with those residues mod 14, ASCENDING (so reference=slowest is right residue),
primitive (gcd=1), and DISTINCT positive integers.
"""
from fractions import Fraction as F
from math import gcd
from itertools import permutations
from functools import reduce

UNITS=[a for a in range(1,14) if gcd(a,14)==1]
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def section(v,a): return (v*a)%14

def build_tournament_M3(speeds):
    S=sorted(set(speeds)); n=len(S); v0=S[0]
    w={a:(1 if section(v0,a) in (1,13) else -1) for a in UNITS}
    fr={v:{a:nrm(F(v*a,14)) for a in UNITS} for v in S}
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            x=S[i]; y=S[j]; D=0
            for a in UNITS:
                dx=fr[x][a]; dy=fr[y][a]
                s=1 if dx>dy else (-1 if dx<dy else 0)
                D+=w[a]*s
            if D>0: adj[i][j]=1
            elif D<0: adj[j][i]=1
            else:
                adj[i][j]=1 if x>y else 0
                if x<=y: adj[j][i]=1
    return S,adj

def scores(adj,n): return [sum(adj[i][j] for j in range(n)) for i in range(n)]
def c3(adj,n,sc): return n*(n-1)*(n-2)//6 - sum(s*(s-1)//2 for s in sc)
def hpaths(adj,n):
    cnt=0
    for p in permutations(range(n)):
        if all(adj[p[k]][p[k+1]] for k in range(n-1)): cnt+=1
    return cnt
def primitive(S): return reduce(gcd,S)==1
def desc(S):
    _,adj=build_tournament_M3(S); n=len(S); sc=scores(adj,n)
    return (hpaths(adj,n), c3(adj,n,sc), tuple(sorted(sc)), [section(v,1)%14 for v in sorted(S)], adj)

def find_speeds_for_residues(res_ascending):
    """Find DISTINCT positive integers, ascending, with given residues mod 14 (in ascending order), primitive.
       Greedy: smallest strictly-increasing distinct speeds matching residues, then bump first speed
       by +14 across a few options to fix primitivity without breaking ascending residue match."""
    n=len(res_ascending)
    def greedy(start_bump):
        speeds=[]; cur=0
        for idx,r in enumerate(res_ascending):
            c=(r if r>0 else 14)
            while c<=cur:
                c+=14
            if idx==0:
                c+=14*start_bump
            speeds.append(c); cur=c
        return speeds
    for bump in range(0,60):
        speeds=greedy(bump)
        if len(set(speeds))==n and all(s>0 for s in speeds) and primitive(speeds) \
           and [s%14 for s in speeds]==list(res_ascending) and speeds==sorted(speeds):
            return speeds
    # fallback: bump the LAST speed instead
    for bump in range(1,60):
        speeds=greedy(0); speeds[-1]+=14*bump; speeds=sorted(speeds)
        if len(set(speeds))==n and primitive(speeds) and [s%14 for s in speeds]==list(res_ascending):
            return speeds
    return None

def covering_check(S):
    """Does S contain a multiple of every q in 2..14? (HARD/covering test)"""
    return {q:any(v%q==0 for v in S) for q in range(2,15)}

print("UNITS:",UNITS)
print("="*70)

targets=[
  ("n=4 forbidden (5,2,(1,1,2,2)) regular near-3-cycle", (2,2,4,6)),
  ("n=5 forbidden (9,3,(1,1,2,3,3))", (2,2,2,4,6)),
  ("n=5 forbidden (11,4,(1,2,2,2,3))", (2,2,4,4,6)),
]

for name,res in targets:
    res=tuple(sorted(res))
    speeds=find_speeds_for_residues(res)
    print(f"\nTARGET: {name}")
    print(f"  residues (ascending): {res}")
    if speeds is None:
        print("  could not realize as distinct primitive integers (unexpected)"); continue
    H,cc,sc,sects,adj=desc(speeds)
    print(f"  CONCRETE primitive speed set S = {speeds}  (gcd={reduce(gcd,speeds)})")
    print(f"  residues mod14 = {[s%14 for s in speeds]}")
    print(f"  M3 tournament -> H={H}, #3cyc={cc}, score={sc}")
    print(f"  adjacency rows: {adj}")
    cov=covering_check(speeds)
    missing=[q for q in range(2,15) if not cov[q]]
    print(f"  covering 2..14? missing divisors: {missing if missing else 'NONE (full covering set!)'}")

print("\n" + "="*70)
print("Now also try to realize forbidden classes with a GENUINE COVERING / HARD set (contains 84-type park).")
# build a covering n=5 set whose M3 class is high-3cyc. Search small.
from itertools import combinations
best=None
# we want a primitive 5-set, ideally covering many q, that realizes H>=9
import random
rng=random.Random(7)
hits=[]
# brute over residue patterns we know give H>=9: must have multiple residue-2 (or 12) runners etc.
# realize concretely with speeds up to ~120, search for covering ones
pool=list(range(1,400))
random_found=0
for _ in range(200000):
    S=sorted(rng.sample(pool,5))
    if reduce(gcd,S)!=1: continue
    H,cc,sc,sects,adj=desc(S)
    if H>=9:
        random_found+=1
        cov=covering_check(S); miss=[q for q in range(2,15) if not cov[q]]
        if random_found<=8 or len(miss)<=3:
            hits.append((S,H,cc,sc,len(miss),miss))
print(f"random 5-sets (speeds<=139) realizing H>=9 found in 40000 samples: {random_found}")
for S,H,cc,sc,nm,miss in hits[:12]:
    print(f"   S={S} H={H} c3={cc} score={sc} missing-div={miss}")
