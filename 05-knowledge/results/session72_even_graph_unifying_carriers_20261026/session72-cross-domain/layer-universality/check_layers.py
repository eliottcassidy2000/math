#!/usr/bin/env python3
from itertools import combinations
from math import prod
from pathlib import Path
import hashlib,json,sys

def need(ok,why):
    if not ok:raise RuntimeError(why)
if '--inject-failure' in sys.argv:need(False,'injected layer failure')
parent=Path(__file__).resolve().parent.parent/'THEOREM-CYCLIC-ALIGNMENT.md'
need(hashlib.sha256(parent.read_bytes()).hexdigest()=='c9cdb86dfd94929b7bdf5a3a9145bf42fc95be4e6f10f2305197cb39d96650cf','parent CA2/CA3 source pin')
primes=[3,5,7,11,13,17,19,23]
need(all(all(p%d for d in range(2,p)) for p in primes),'distinct odd prime choices')

def connected(S,E):
    if not S:return True
    seen={min(S)}
    while True:
        new=seen|{j for i,j in E if i in seen and j in S}|{i for i,j in E if j in seen and i in S}
        if new==seen:return seen==S
        seen=new
families=cases=0
for r in range(1,4):
    subsets=[{i for i in range(r) if mask>>i&1} for mask in range(1<<r)]
    edges=list(combinations(range(r),2))
    for fm in range(1<<len(subsets)):
        sets=[s for j,s in enumerate(subsets) if fm>>j&1]
        ps=primes[:len(sets)]
        ds=[prod(p for p,S in zip(ps,sets) if i in S) for i in range(r)]
        m=prod(ps) if ps else 3
        need(m>=3 and m%2==1 and all(m%d==0 for d in ds),'odd-carousel period scope')
        layers=[{i for i,d in enumerate(ds) if d%p==0} for p in ps]
        need(layers==sets,'exact prescribed prime layers')
        need(all(all(d%(p*p)!=0 for d in ds) for p in ps),'no higher prime layers')
        for em in range(1<<len(edges)):
            E=[e for j,e in enumerate(edges) if em>>j&1]
            need(all(connected(s,E) for s in sets)==all(connected(s,E) for s in layers),'simultaneous induced connectivity representation')
            cases+=1
        families+=1
sets=[{0,2},{0,2},set(),{1}];ps=primes[:4];ds=[prod(p for p,S in zip(ps,sets) if i in S) for i in range(3)]
need([{i for i,d in enumerate(ds) if d%p==0} for p in ps]==sets,'duplicate and empty indexed supports')
print(json.dumps({'status':'PASS','families':families,'graph_family_cases':cases,'index_sizes':[1,2,3],'duplicates_and_empty_supports':True,'parent_theorem_sha256':'c9cdb86dfd94929b7bdf5a3a9145bf42fc95be4e6f10f2305197cb39d96650cf','scope':'Representation universality only; no complexity or huge native graph census claim'},sort_keys=True,indent=2))
