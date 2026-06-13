"""C=39=3*13 confirmation: predict deficiency=2^{(13-1)/2}=64, flips {13} and {3,6,9,12,15,18}."""
from itertools import product
import math
def folded(V,eps,C):
    out=[]
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            f=(eps[i]*V[i]-eps[j]*V[j])%C; out.append(min(f,C-f))
    return tuple(sorted(out))
n=20;C=39;V=list(range(1,n))
fm={}
for bits in product([0,1],repeat=n-2):
    eps=tuple([1]+[1 if b else -1 for b in bits])
    fm.setdefault(folded(V,list(eps),C),[]).append(eps)
coll=[g for g in fm.values() if len(g)>1]
defic=sum(len(g)-1 for g in coll)
stats={}
for g in coll:
    b=g[0]
    for o in g[1:]:
        D=[V[i] for i in range(len(V)) if b[i]!=o[i]]
        D=min([D,[V[i] for i in range(len(V)) if b[i]==o[i]]],key=len)
        gc=tuple(sorted(set(math.gcd(d,C) for d in D)))
        stats[(len(D),gc)]=stats.get((len(D),gc),0)+1
print(f"C=39 n=20: deficiency={defic} (predict 2^6=64); class-sizes={sorted((len(g) for g in coll),reverse=True)[:5]}...")
print(f"Diff-stats={stats}")
