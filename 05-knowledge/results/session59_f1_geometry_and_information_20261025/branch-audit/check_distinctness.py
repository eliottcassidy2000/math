#!/usr/bin/env python3
"""Independent minimal kernel-union audit for nine-entry injectivity."""
from itertools import product,combinations
from collections import Counter
from math import gcd,prod
from functools import reduce
import json
from pathlib import Path
from check_independent import invfac,need

def canonical(row):
    first=next(x for x in row if x)
    return row if first>0 else tuple(-x for x in row)

points=list(product(range(-1,2),repeat=2))
all_differences=sorted({canonical((x[0]-y[0],x[1]-y[1])) for x,y in combinations(points,2)})
rows=[(2,0),(0,2),(2,2),(2,-2),(2,1),(2,-1),(1,2),(1,-2)]
need(len(all_differences),12,'twelve canonical pair differences')
need(set(rows).issubset(set(all_differences)),True,'retained rows are actual differences')
discarded=[r for r in all_differences if r not in rows]
need(set(discarded),{(1,0),(0,1),(1,1),(1,-1)},'four redundant primitive rows')
for r in discarded:
    need(tuple(2*x for x in r) in rows,True,('kernel inclusion',r))
# Since ker(r) lies in ker(2r) for every abelian group, both unions coincide.
agg=Counter();records=[]
for mask in range(1<<len(rows)):
    sub=[r for i,r in enumerate(rows) if (mask>>i)&1]
    ds=invfac(sub,2);agg[ds]+=(-1)**len(sub)
    records.append({'mask':mask,'invariants':ds})
expected={():1,(1,):-4,(2,):-4,(1,1):-2,(1,2):-4,(1,4):4,(1,5):2,(1,6):4,(2,2):2,(2,4):1}
need({k:v for k,v in agg.items() if v},expected,'signed eight-row census')

def h(d,mods):return prod(gcd(d,m) for m in mods)
def formula(N,mods):
    h2,h4,h5,h6=(h(d,mods) for d in [2,4,5,6])
    return N*N-4*(1+h2)*N-2-4*h2+4*h4+2*h5+4*h6+2*h2*h2+h2*h4

def count(mods):
    els=list(product(*(range(m) for m in mods)));N=len(els)
    full=distinct=nonidentity=0
    for x,z in product(els,repeat=2):
        off={tuple((i*a+j*b)%m for a,b,m in zip(x,z,mods)) for i,j in points}
        if len(off)!=9:continue
        full+=1
        for e in els:
            vals={tuple((a+b)%m for a,b,m in zip(e,v,mods)) for v in off}
            need(len(vals),9,('translation injectivity',mods))
            distinct+=1
            nonidentity+=((0,)*len(mods) not in vals)
    need(full,formula(N,mods),('pair formula',mods))
    need(full,sum(c*N**(2-len(ds))*prod(h(d,mods) for d in ds) for ds,c in agg.items()),('IE formula',mods))
    need(distinct,N*full,('all centers',mods))
    need(nonidentity,(N-9)*full,('nonidentity centers',mods))
    return {'factors':mods,'N':N,'injective_pairs':full,'distinct_arrays':distinct,'distinct_nonidentity_arrays':nonidentity}

modslist=[(n,) for n in range(1,41)]+[(2,2),(2,4),(4,4),(2,2,2),(2,2,4),(2,3,4),(2,4,4),(3,3),(3,9),(3,6),(5,5),(2,2,2,2)]
groups=[count(mods) for mods in modslist]
field_tests=[]
for p,k in [(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(5,1),(5,2),(7,1),(7,2),(11,1)]:
    N=p**k;actual=formula(N,(p,)*k)
    expected=0 if p==2 else (N-1)*(N-(3 if p==3 else 5 if p==5 else 7))
    need(actual,expected,('field specialization',p,k))
    field_tests.append({'p':p,'dimension':k,'N':N,'I':actual})
print(json.dumps({'status':'PASS','canonical_differences':all_differences,'retained_rows':rows,'discarded_rows':discarded,'signed_SNF':{str(k):v for k,v in sorted(agg.items())},'records':records,'groups':groups,'field_tests':field_tests},sort_keys=True,indent=2))
