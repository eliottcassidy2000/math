#!/usr/bin/env python3
"""
VERIFY the residue-reduction is EXACT, then push the search to n=6,7 (sampled over real primitive sets).

Step A: For random primitive distinct-speed sets, confirm build_tournament_M3(speeds) (slowest=ref, exact)
        EQUALS build_adj(residues-in-ascending-speed-order). If they always match, the residue enumeration
        in the primitive_exhaustive script is a faithful (complete) model of ALL primitive LRC inputs.

Step B: Directly search n=6 and n=7 over real primitive integer speed sets (random + structured covering/tight)
        for the high-3-cycle / regular end (H large, score near-regular). Report the max H and #3cyc reachable.
"""
from fractions import Fraction as F
from math import gcd
from itertools import permutations, combinations
from functools import reduce
import random, sys

UNITS=[a for a in range(1,14) if gcd(a,14)==1]
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def section(v,a): return (v*a)%14
FRAC={r: tuple(nrm(F(r*a,14)) for a in UNITS) for r in range(14)}
WPAT={r: tuple(1 if (r*a)%14 in (1,13) else -1 for a in UNITS) for r in range(14)}

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
                if x>y: adj[i][j]=1
                else: adj[j][i]=1
    return S,adj

def build_adj_res(res):
    n=len(res); v0=res[0]; w=WPAT[v0]; fr=[FRAC[r] for r in res]
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            D=0
            for k in range(6):
                if fr[i][k]>fr[j][k]: D+=w[k]
                elif fr[i][k]<fr[j][k]: D-=w[k]
            if D>0: adj[i][j]=1
            elif D<0: adj[j][i]=1
            else: adj[j][i]=1
    return adj

def scores(adj,n): return [sum(adj[i][j] for j in range(n)) for i in range(n)]
def c3(adj,n,sc): return n*(n-1)*(n-2)//6 - sum(s*(s-1)//2 for s in sc)
def hpaths(adj,n):
    cnt=0
    for p in permutations(range(n)):
        if all(adj[p[k]][p[k+1]] for k in range(n-1)): cnt+=1
    return cnt
def desc(adj,n):
    sc=scores(adj,n); return (hpaths(adj,n), c3(adj,n,sc), tuple(sorted(sc)))
def prim(S): return reduce(gcd,S)==1

# ---------- Step A: verify reduction exact ----------
print("STEP A: verify build_tournament_M3 == build_adj_res(residues asc) on random primitive sets")
rng=random.Random(123); mism=0; checked=0
for n in [3,4,5,6,7,8]:
    for _ in range(3000):
        S=sorted(rng.sample(range(1,300),n))
        if not prim(S): continue
        Ssrt,adj1=build_tournament_M3(S)
        res=[s%14 for s in Ssrt]      # ascending speed order
        adj2=build_adj_res(res)
        checked+=1
        if adj1!=adj2:
            mism+=1
            if mism<=5: print("   MISMATCH:",S,res)
print(f"   checked {checked} primitive sets, mismatches={mism}")
print(f"   => reduction is {'EXACT (residue enumeration is complete)' if mism==0 else 'NOT exact (!!)'}")
print()

# ---------- Step B: n=6,7 direct search for high-3cyc / regular end ----------
print("STEP B: search n=6,7 real primitive sets for high-3cyc / near-regular end")
for n in [6,7]:
    rng=random.Random(999+n)
    maxH=0; maxc3=0; seen=set(); best=None
    # random
    for _ in range(120000):
        S=sorted(rng.sample(range(1,400),n))
        if not prim(S): continue
        _,adj=build_tournament_M3(S)
        d=desc(adj,n); seen.add(d)
        if d[0]>maxH: maxH=d[0]; best=(S,d)
        if d[1]>maxc3: maxc3=d[1]
    # structured: covering sets {1..n-1, big multiples}, tight AP {1..n}, GW-type
    structured=[
        list(range(1,n+1)),                       # AP {1..n}
        list(range(1,n))+[14],                     # park one at residue 0
        list(range(1,n))+[84],                     # 84-park
        [1,2,3,4,5,6,7][:n],
        [1,3,5,7,9,11,13][:n],                     # odd AP
        [1,2,4,8,16,32,64][:n],                    # geometric
        [2,3,5,7,11,13,17][:n],                    # primes
    ]
    for S in structured:
        S=sorted(set(S))
        if len(S)!=n or not prim(S): continue
        _,adj=build_tournament_M3(S); d=desc(adj,n); seen.add(d)
        if d[0]>maxH: maxH=d[0]; best=(S,d)
        if d[1]>maxc3: maxc3=d[1]
    print(f"  n={n}: distinct (H,c3,score) classes seen = {len(seen)}; MAX H={maxH}, MAX c3={maxc3}")
    print(f"        argmax-H example: {best}")
    # is the regular tournament reachable? regular n=7 has score all 3, H=189; n=6 has no regular
    reg_scores={6:None, 7:(3,3,3,3,3,3,3)}
    if n in reg_scores and reg_scores[n] is not None:
        reg_seen=any(d[2]==reg_scores[n] for d in seen)
        print(f"        regular score {reg_scores[n]} reachable? {reg_seen}")
sys.stdout.flush()
