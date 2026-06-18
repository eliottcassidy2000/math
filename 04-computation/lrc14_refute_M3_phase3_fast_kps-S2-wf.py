#!/usr/bin/env python3
"""
FAST canonical exhaustive Phase-3 verifier for the M3 Walsh-weighted-lead forbidden-class claim.

Key fact (proved by direct check): M3 sign-votes depend ONLY on speeds mod 14.
The only actual-speed dependence is (a) reference = slowest speed's residue, (b) tie-break = actual speed order.
So EVERY achievable M3 tournament for ANY primitive LRC input (any size, any covering/tight/sporadic case)
arises from some ordered residue-tuple (position 0 = slowest = reference; later positions = faster speeds,
used for tie-break). Enumerate ALL 14^n ordered residue-tuples -> superset of all real inputs.

Drops the expensive canonical-signature; forbidden check needs only (H, c3, score).
Unbuffered output (run with python3 -u or rely on prints + flush).
"""
from fractions import Fraction as F
from math import gcd
from itertools import product, permutations
import sys

UNITS=[a for a in range(1,14) if gcd(a,14)==1]

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r

# Precompute frac rows per residue 0..13 (as tuples of Fractions over UNITS)
FRAC={r: tuple(nrm(F(r*a,14)) for a in UNITS) for r in range(14)}
# Precompute weight pattern per reference residue
WPAT={r: tuple(1 if (r*a)%14 in (1,13) else -1 for a in UNITS) for r in range(14)}

def build_adj(res):
    n=len(res); v0=res[0]; w=WPAT[v0]
    fr=[FRAC[r] for r in res]
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        fri=fr[i]
        for j in range(i+1,n):
            frj=fr[j]
            D=0
            for k in range(len(UNITS)):
                dx=fri[k]; dy=frj[k]
                if dx>dy: D+=w[k]
                elif dx<dy: D-=w[k]
            if D>0: adj[i][j]=1
            elif D<0: adj[j][i]=1
            else: adj[j][i]=1   # tie-break x>y false (speed_i<speed_j) -> y->x
    return adj

def scores(adj,n):
    return [sum(adj[i][j] for j in range(n)) for i in range(n)]

def c3(adj,n,sc):
    total=n*(n-1)*(n-2)//6
    return total - sum(s*(s-1)//2 for s in sc)

def hpaths(adj,n):
    cnt=0
    for perm in permutations(range(n)):
        ok=True
        for k in range(n-1):
            if not adj[perm[k]][perm[k+1]]: ok=False;break
        if ok: cnt+=1
    return cnt

def desc(adj,n):
    sc=scores(adj,n)
    return (hpaths(adj,n), c3(adj,n,sc), tuple(sorted(sc)))

FORB={
4:{(5,2,(1,1,2,2))},
5:{(9,3,(1,1,2,3,3)),(11,4,(1,2,2,2,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))},
6:None,  # we'll just report all classes seen, flag the regular/high-3cyc end
}

def run(n):
    forb=FORB.get(n)
    classes={}
    found=[]
    total=0
    for res in product(range(14),repeat=n):
        total+=1
        adj=build_adj(list(res))
        d=desc(adj,n)
        classes[d]=classes.get(d,0)+1
        if forb and d in forb:
            found.append((res,d))
    print(f"=== n={n}: enumerated {total} = 14^{n} ordered residue-tuples ===")
    print(f"distinct (H,c3,score) classes realized: {len(classes)}")
    # sort by H then c3
    for d in sorted(classes):
        print(f"   {d}: {classes[d]}")
    if forb is not None:
        print(f"FORBIDDEN targets: {sorted(forb)}")
        print(f"FORBIDDEN realized? {len(found)>0}")
        for res,d in found[:20]:
            print(f"   WITNESS residues={res} -> {d}")
    # report max H achieved and whether any 3-cycle-dense class appears
    maxH=max(d[0] for d in classes)
    maxc3=max(d[1] for d in classes)
    print(f"MAX H realized = {maxH}; MAX #3cyc realized = {maxc3}")
    sys.stdout.flush()
    return classes,found

if __name__=="__main__":
    for n in [4,5]:
        run(n)
    print("\n[Optional] n=6 class census (no fixed forbidden list; show the realized spectrum):")
    run(6)
