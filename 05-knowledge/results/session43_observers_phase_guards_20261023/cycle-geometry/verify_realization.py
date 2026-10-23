#!/usr/bin/env python3
"""Exhaustive lower-envelope criterion vs independent exponent-box search."""
from itertools import combinations,product
import json

def check(x,label):
    if not x: raise RuntimeError(label)
def factor(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def lower(rows):
    result={}
    for b,k in rows:
        for p,f in factor(b).items():result[p]=max(result.get(p,0),k*f)
    return result

def height(e,b):
    return min(e.get(p,0)//f for p,f in factor(b).items())
def criterion(rows):
    e=lower(rows)
    return all(height(e,b)==k for b,k in rows),e

bases=[2,3,4,6,8,9,12]
checked=0;possible=0;impossible=0
# Any feasible triple with k<=3 has minimal exponents <=9 at2, <=6 at3.
# Exhaustively test a larger fixed box, not just the proposed minimal vector.
box=[({2:e,3:f},tuple(height({2:e,3:f},b) for b in bases)) for e in range(12) for f in range(9)]
for bs in combinations(bases,3):
    inds=[bases.index(b) for b in bs]
    attained={tuple(h[i] for i in inds) for e,h in box}
    for ks in product(range(4),repeat=3):
        rows=list(zip(bs,ks));yes,e=criterion(rows)
        check(yes==(ks in attained),('criterion mismatch',rows))
        if yes:
            possible+=1
            for f,hs in box:
                if tuple(hs[i] for i in inds)==ks:
                    check(all(f[p]>=v for p,v in e.items()),'least property')
        else: impossible+=1
        checked+=1
cycle_rows=[(4,2),(2,5),(5,3),(3,7),(7,4)]
yes,e=criterion(cycle_rows);check(yes,'cycle feasible')
N0=1
for p,k in e.items():N0*=p**k
check(N0==39046339332000//(11*169),'least cycle realization')
# Pairwise individually possible but shared-prime incompatible.
hostile=[(2,3),(4,2)]
check(not criterion(hostile)[0],'shared prime hostile')
# A nonminimal shared base can be repaired at a different prime.
repair=[(2,3),(6,2)]
check(criterion(repair)[0] and lower(repair)=={2:3,3:2},'two-prime repair')
# Zero heights are legitimate at N=1; do not claim a least N>1.
check(criterion([(2,0),(3,0)])[0],'zero targets')
check(all(v==0 for v in lower([(2,0),(3,0)]).values()),'zero minimum')
print(json.dumps({'status':'PASS','constraint_triples':checked,'realizable':possible,'infeasible':impossible,'independent_box_exponent_bounds':{'2':[0,11],'3':[0,8]},'least_cycle_N':N0,'least_cycle_exponents':e,'shared_prime_hostile':hostile,'repair':repair},indent=2,sort_keys=True))
