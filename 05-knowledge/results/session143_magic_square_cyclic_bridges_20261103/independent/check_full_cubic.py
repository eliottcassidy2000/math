#!/usr/bin/env python3
"""Independent exact central-symmetry candidate audit for all cubic assignments.
This eliminates B!=0; check_independent.py handles the remaining odd plane.
No producer code imported; standard library; explicit exceptions survive -O.
"""
from fractions import Fraction as F
from itertools import permutations,product
from collections import Counter
import json,sys

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def pairings(xs):
 if not xs:yield ();return
 a=xs[0]
 for j,b in enumerate(xs[1:],1):
  for tail in pairings(xs[1:j]+xs[j+1:]):yield ((a,b),)+tail
def rref(rows):
 rows=[list(map(F,row)) for row in rows];pivots=[];r=0
 for c in range(3):
  piv=next((i for i in range(r,len(rows)) if rows[i][c]),None)
  if piv is None:continue
  rows[r],rows[piv]=rows[piv],rows[r]
  den=rows[r][c];rows[r]=[x/den for x in rows[r]]
  for i in range(len(rows)):
   if i!=r:
    factor=rows[i][c];rows[i]=[x-factor*y for x,y in zip(rows[i],rows[r])]
  pivots.append(c);r+=1
  if r==len(rows):break
 return rows[:r],pivots
def nullspace(rows):
 R,pivots=rref(rows);basis=[]
 for c in range(3):
  if c in pivots:continue
  v=[F(0)]*3;v[c]=1
  for row,p in zip(R,pivots):v[p]=-row[c]
  basis.append(tuple(v))
 return tuple(basis)
def canonical_ray(v):
 first=next(x for x in v if x)
 return tuple(x/first for x in v)
spaces=Counter();candidates={};count=0
for center in range(-4,5):
 for pairs in pairings(tuple(k for k in range(-4,5) if k!=center)):
  rows=[(center**3,F(center*center)-F(20,3),center)]
  for a,b in pairs:rows.append((a**3+b**3,F(a*a+b*b)-F(40,3),a+b))
  basis=nullspace(rows);spaces[len(basis)]+=1;count+=1
  if any(v[1] for v in basis):
   require(len(basis)==1,'nonodd candidate must be one ray')
   ray=canonical_ray(basis[0]);candidates.setdefault(ray,[]).append({'center':center,'pairs':pairs})
require(count==945,'all9*105center-pair systems')
require(dict(spaces)=={0:914,1:30,2:1},'complete nullity census')
expected={(F(1),F(3),F(-18)),(F(1),F(-3),F(-18)),(F(1),F(3,2),F(-39,2)),(F(1),F(-3,2),F(-39,2))}
if '--mutate-ray' in sys.argv:expected.remove((F(1),F(3),F(-18)))
require(set(candidates)==expected,'four nonodd central-symmetry rays')
rejections=[]
for ray,examples in sorted(candidates.items()):
 A,B,C=ray; values=[A*k**3+B*k*k+C*k-F(20,3)*B for k in range(-4,5)]
 require(sum(values)==0,'centered values')
 for example in examples:
  quartet=[values[a+4] for a,b in example['pairs']]
  survivors=[]
  for ps in permutations(quartet):
   for signs in product((-1,1),repeat=4):
    u,v,w,z=(s*x for s,x in zip(signs,ps))
    if u+v==w and u-v==z:survivors.append((u,v,w,z))
  require(not survivors,'nonodd central-symmetry ray cannot be magic')
 # Any magic quartet abs-values has max=abs(u)+abs(v), so max is at least
 # the sum of the two smallest members of the quartet (triangle constraint).
 absq=sorted(abs(values[a+4]) for a,b in examples[0]['pairs'])
 require(absq[-1]<absq[0]+absq[1],'strict triangle hostile obstruction')
 rejections.append({'ray':list(map(str,ray)),'absolute_quartet':list(map(str,absq)),'center_pairings':len(examples)})
print(json.dumps({'status':'PASS','systems':count,'nullity_census':dict(sorted(spaces.items())),
 'nonodd_candidates_rejected':rejections,'conclusion':'B must vanish for any nonconstant magic cubic. Combined with independently checked odd-plane classification: A=B=0, or B=0 and C=-15A.',
 'scope':'All real cubic numerical assignments on indices -4..4; not a squarehood or arbitrary MSOS theorem.'},indent=2,sort_keys=True))
