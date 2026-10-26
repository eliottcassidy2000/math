"""Independent literal class orbits, pair-free survivors, and naive transform hostile."""
from itertools import combinations, permutations
from math import comb
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def census(n):
    pairs=list(combinations(range(n),2));index={p:i for i,p in enumerate(pairs)}
    maps=[]
    for perm in permutations(range(n)):
        relabel=[]
        for a,b in pairs:
            x,y=perm[a],perm[b]
            relabel.append((index[min(x,y),max(x,y)],int(x>y)))
        maps.append(relabel)
    seen=set();ordinary=core=0
    for mask in range(2**len(pairs)):
        if mask in seen:continue
        ordinary+=1
        def arc(a,b):
            return ((mask>>index[min(a,b),max(a,b)])&1) ^ int(a>b)
        pairfree=not any(all(arc(a,c)==arc(b,c) for c in range(n) if c not in (a,b)) for a,b in pairs)
        core+=pairfree
        orbit={sum((((mask>>old)&1)^flip)<<j for j,(old,flip) in enumerate(mapping)) for mapping in maps}
        need(orbit.isdisjoint(seen),'literal relabeling orbit disjointness')
        seen.update(orbit)
    need(len(seen)==2**len(pairs),'complete labelled universe')
    return ordinary,core,len(seen)

rows=[]
for n in range(7):
    ordinary,core,labelled=census(n)
    rows.append({'n':n,'ordinary':ordinary,'pair_free_core':core,'labelled':labelled})
need([r['ordinary'] for r in rows]==[1,1,1,2,4,12,56],'ordinary literal count')
need([r['pair_free_core'] for r in rows]==[1,1,0,1,2,5,26],'core literal count')
naive={n:sum((-1)**(n-k)*comb(n-1,k-1)*rows[k]['ordinary'] for k in range(1,n+1)) for n in range(1,7)}
need(next(n for n in naive if naive[n]!=rows[n]['pair_free_core'])==4,'first scalar substitution failure')
need(naive[4]==0 and naive[6]==20,'scalar binomial substitution misses cycle inventory')
print(json.dumps({'status':'PASS','orders':[0,6],'rows':rows,'naive_scalar_substitution_hostiles':{'order4':{'naive':naive[4],'actual':2},'order6':{'naive':naive[6],'actual':26}},'scope':'Literal tournament relabeling orbits and homogeneous-pair test; no inherited counting or clone routines imported.'},indent=2,sort_keys=True))
