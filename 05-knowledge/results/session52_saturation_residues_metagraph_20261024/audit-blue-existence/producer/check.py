#!/usr/bin/env python3
"""Small exact controls only; no all-order assertion is inferred."""
from itertools import combinations,permutations,product
import json

def check(n):
    q=lambda v:n-1-v
    pairs=list(combinations(range(n),2)); unseen=set(pairs); orbits=[]
    while unseen:
        e=min(unseen); f=tuple(sorted((q(e[0]),q(e[1]))))
        orb=sorted({e,f}); unseen.difference_update(orb); orbits.append(orb)
    total=balanced=0; minimum=None; centers=set()
    for bits in product(range(2),repeat=len(orbits)):
        a=[[0]*n for _ in range(n)]
        for orb,b in zip(orbits,bits):
            for u,v in orb:a[u][v]=b;a[v][u]=1-b
        total+=1
        d=list(map(sum,a));p=n//2
        if sorted(d)!=[p-1]*p+[p]*p:continue
        balanced+=1; counts=[0]*p
        for order in permutations(range(p)):
            for signs in product(range(2),repeat=p):
                half=[q(v) if s else v for v,s in zip(order,signs)]
                path=half+[q(v) for v in reversed(half)]
                if d[path[0]]!=p:continue
                if all(a[u][v] for u,v in zip(path,path[1:])):
                    counts[min(path[p-1],q(path[p-1]))]+=1
        if min(counts)<1:raise RuntimeError(('missing central arc',n,bits,counts))
        minimum=min(sum(counts),minimum if minimum is not None else sum(counts))
    return dict(n=n,anti_invariant_matrices=total,balanced_matrices=balanced,minimum_symmetric_downhill_paths=minimum,every_paired_arc_realized=True)
if __name__=='__main__':print(json.dumps({'status':'FINITE-EXACT','rows':[check(n) for n in (2,4,6)]},sort_keys=True,indent=2))
