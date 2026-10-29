#!/usr/bin/env python3
"""Independent primitive mod-eight exclusions for E:y^2=x^3-x."""
from itertools import product
from math import isqrt
import json
N=0

def ck(v,msg):
 global N
 N+=1
 if not v:raise RuntimeError(msg)

S=(-2,-1,1,2)
candidates=[d for d in product(S,repeat=3) if d[0]*d[1]*d[2]>0 and isqrt(d[0]*d[1]*d[2])**2==d[0]*d[1]*d[2]]
real=[d for d in candidates if d[2]>0 and d[0]*d[1]>0]
ck(len(candidates)==16,'candidate cardinality');ck(len(real)==8,'real survivors')
counts={}
for d in real:
 count=0
 for z1,z2,z3,t in product(range(8),repeat=4):
  if not any(v%2 for v in (z1,z2,z3,t)):continue
  count+=int((d[0]*z1*z1-d[1]*z2*z2-t*t)%8==0 and (d[2]*z3*z3-d[0]*z1*z1-t*t)%8==0)
 counts[d]=count
 ck((count>0)==(abs(d[0])==1),('local exclusion',d))
expected={(-1,-2,2):256,(-1,-1,1):128,(1,1,1):128,(1,2,2):256}
ck({d:c for d,c in counts.items() if c}==expected,'literal modular counts')
# Kummer values at branch points, replacing the zero factor by product of differences.
roots=(0,1,-1)
images={(1,1,1)}
def sf(x):
 sign=-1 if x<0 else 1;x=abs(x);ans=sign;p=2
 while p*p<=x:
  power=0
  while x%p==0:x//=p;power+=1
  if power%2:ans*=p
  p+=1
 return ans*x
for a in roots:
 vals=[]
 for b in roots:
  vals.append(sf(a-b) if a!=b else sf(__import__('functools').reduce(lambda x,y:x*y,(a-c for c in roots if c!=a),1)))
 images.add(tuple(vals))
ck(images==set(expected),'actual rational Kummer images')
for a,b in product(images,repeat=2):
 ck(tuple(sf(x*y) for x,y in zip(a,b)) in images,'Kummer subgroup closure')
print(json.dumps({'status':'PASS','checks':N,'primitive_mod8_counts':{str(d):c for d,c in counts.items()},'scope':'finite local exclusion and four actual rational images; full descent identification cited separately'},sort_keys=True))
