#!/usr/bin/env python3
from itertools import combinations,product
from fractions import Fraction
from pathlib import Path
import json
q=3;n=6;r=3

def Q(v):return sum(v[i]*v[i+r] for i in range(r))%q

def bil(v,w):return sum(v[i]*w[i+r]+v[i+r]*w[i] for i in range(r))%q

def contains(w,v):
    v=list(v)
    for row in w:
        pivot=next(i for i,x in enumerate(row) if x)
        z=v[pivot]
        v=[(a-z*b)%q for a,b in zip(v,row)]
    return not any(v)

spaces=0;maximal=[]
for pivots in combinations(range(n),r):
    free=[(i,j) for i,p in enumerate(pivots) for j in range(p+1,n) if j not in pivots]
    for vals in product(range(q),repeat=len(free)):
        w=[[int(j==p) for j in range(n)] for p in pivots]
        for (i,j),v in zip(free,vals):w[i][j]=v
        spaces+=1
        if all(Q(v)==0 for v in w) and all(bil(v,u)==0 for v,u in combinations(w,2)):
            maximal.append(w)
e1=(1,0,0,0,0,0);e2=(0,1,0,0,0,0);e3=(0,0,1,0,0,0);e12=(1,1,0,0,0,0)
dep=[e1,e2,e12];ind=[e1,e2,e3]
a=sum(not any(contains(w,v) for v in dep) for w in maximal)
b=sum(not any(contains(w,v) for v in ind) for w in maximal)
assert all(bil(v,u)==0 for B in [dep,ind] for v,u in combinations(B,2))
assert spaces==33880 and len(maximal)==80 and a==60 and b==61
out={'all_three_subspaces':spaces,'maximal_isotropic_three_spaces':len(maximal),'dependent_K3_avoiding':a,'independent_K3_avoiding':b,'dependent_probability':str(Fraction(a,len(maximal))),'independent_probability':str(Fraction(b,len(maximal))),'status':'all assertions passed'}
print(json.dumps(out,indent=2))
Path(__file__).with_name('graph_counterexample_results.json').write_text(json.dumps(out,indent=2)+'\n')
