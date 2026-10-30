"""Independent incidence-only enumeration for the constant-antichain lemma."""
from itertools import combinations,product
import json

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
checks=0;topologies=0;cases=[]
for n in range(1,6):
    choices=[list(combinations(range(i+2),2)) for i in range(n)]
    kept=tests=0
    for parents in product(*choices):
        ancestors=[]
        for i,ps in enumerate(parents):
            ancestors.append(set(ps)|set().union(*(ancestors[p-2] for p in ps if p>=2)))
        if not set(range(2,n+1)).issubset(ancestors[-1]):continue
        kept+=1
        for mask in range(1,1<<n):
            seeds={i+2 for i in range(n) if mask>>i&1}
            if any(u in ancestors[v-2] for u in seeds for v in seeds if u!=v):continue
            C=set(seeds)
            for i,ps in enumerate(parents):
                if set(ps)<=C:C.add(i+2)
            if n+1 in C:continue
            B={i+2 for i,ps in enumerate(parents) if i+2 not in C and set(ps)&C}
            require(len(C)+len(B)>=2*len(seeds),(parents,seeds,C,B))
            tests+=1
    checks+=tests;topologies+=kept;cases.append({'gates':n,'live_topologies':kept,'eligible_antichains':tests})
# Comparable seeds invalidate the proposed factor two: a*b -> (a*b)*c -> XOR x.
hostile_parents={4:(0,1),5:(4,2),6:(5,3)}
seeds={4,5};C=set(seeds)
for node,ps in hostile_parents.items():
    if set(ps)<=C:C.add(node)
B={node for node,ps in hostile_parents.items() if node not in C and set(ps)&C}
seed_count=len(seeds);actual_deletions=len(C|B)
require(actual_deletions<2*seed_count,'hostile must disprove missing-antichain variant')
for a,b,c,x in product([0,1],repeat=4):
    if a==0:require((((a&b)&c)^x)==x,'hostile output stays nonconstant')
print(json.dumps({'status':'PASS','scope':'all binary distinct-pin gate topologies through5gates over2 external inputs, one output, live DAGs; symbolic seed constancy only','topologies':topologies,'eligible_antichain_cases':checks,'cases':cases,'comparable_seed_hostile':{'seed_count':seed_count,'deletions':actual_deletions}},indent=2,sort_keys=True))
