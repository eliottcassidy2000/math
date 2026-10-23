#!/usr/bin/env python3
"""Exact finite E8 root-grading and arithmetic checks; standard library only."""
from itertools import combinations, product
from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path
roots=[] # twice the standard root coordinates
for i,j in combinations(range(8),2):
 for s,t in product((-2,2),repeat=2):
  r=[0]*8;r[i]=s;r[j]=t;roots.append(tuple(r))
for r in product((-1,1),repeat=8):
 if sum(x<0 for x in r)%2==0: roots.append(r)
R=set(roots);assert len(R)==240
assert all(sum(x*x for x in r)==8 for r in R)
grade=lambda r:(r[0]+r[1])//2
counts=Counter(map(grade,R));assert dict(counts)=={-2:1,-1:56,0:126,1:56,2:1}
theta=(2,2,0,0,0,0,0,0)
positive=[r for r in R if grade(r)==1]
pairs={tuple(sorted((r,tuple(a-b for a,b in zip(theta,r))))) for r in positive}
assert len(pairs)==28
assert all(t in R and grade(t)==1 for pair in pairs for t in pair)
# All root sums of two degree-one roots that are roots equal theta.
sums={tuple(a+b for a,b in zip(r,s)) for r in positive for s in positive}
assert sums & R == {theta}
# Root spaces plus the appropriate Cartan subspaces.
dims={'e7':counts[0]+7,'heisenberg':counts[1]+counts[2], 'intermediate':counts[0]+7+counts[1]+counts[2], 'parabolic':counts[0]+8+counts[1]+counts[2], 'e8':len(R)+8}
assert dims=={'e7':133,'heisenberg':57,'intermediate':190,'parabolic':191,'e8':248}
assert 57==3*19==2**5+5**2
assert all(p%4==3 and all(p%d for d in range(2,int(p**.5)+1)) for p in (3,19))
# Minimal-model degree-zero dimension N=57 has all coprime p<q realizations.
models=[]
from math import gcd
for p in range(2,116):
 for q in range(p+1,116):
  if (p-1)*(q-1)==114 and gcd(p,q)==1:
   c=1-F(6*(p-q)**2,p*q);N=57;cb=1-F(6*(2-115)**2,2*115)
   A=F(N*(2*N+1),8)*(c-cb)
   models.append({'p':p,'q':q,'c':str(c),'N':N,'first_transfer_coefficient':str(A)})
assert [(m['p'],m['q']) for m in models]==[(2,115),(3,58),(4,39),(7,20)]
result={'root_grades':dict(sorted(counts.items())), 'symplectic_pairs':len(pairs), 'dimensions':dims,'arithmetic_57':{'Blum_primes':[3,19],'Leyland_pair':[2,5]},'formal_level_one_c':str(F(190,25)),'minimal_models_with_C2_dimension_57':models}
Path(__file__).with_name('exact_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
