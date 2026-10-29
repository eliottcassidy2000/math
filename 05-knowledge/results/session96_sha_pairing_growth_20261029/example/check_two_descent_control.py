#!/usr/bin/env python3
from itertools import product
import json
checks=0

def check(ok,msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)

def squareclass(a):
 sign=1 if a>0 else -1;a=abs(a);v=0
 while a%2==0:a//=2;v+=1
 check(a==1,'only 2-supported inputs')
 return sign*(2 if v%2 else 1)

classes=(1,2,-1,-2)
triples=[(a,b,squareclass(a*b)) for a,b in product(classes,repeat=2)]
check(len(set(triples))==16,'full source bank')
real=[d for d in triples if d[2]>0 and d[0]*d[1]>0]
check(len(real)==8,'real sign bank')
rows=[];survive=[];tested=0
for d1,d2,d3 in real:
 count=0
 for z1,z2,z3,t in product(range(8),repeat=4):
  tested+=1
  primitive=any(v%2 for v in (z1,z2,z3,t))
  equations=(d1*z1*z1-d2*z2*z2-t*t)%8==0 and (d3*z3*z3-d1*z1*z1-t*t)%8==0
  if primitive and equations:count+=1
  if not primitive and not any((z1,z2,z3,t)):check(equations,'zero-vector hostile')
 check((count>0)==(abs(d1)==1),('local exclusion',(d1,d2,d3)))
 if count:survive.append((d1,d2,d3))
 rows.append({'triple':[d1,d2,d3],'primitive_mod8_points':count})
witnesses={(1,1,1):(1,1,1,0),(-1,-1,1):(0,1,1,1),(1,2,2):(1,0,1,1),(-1,-2,2):(1,1,0,1)}
check(set(survive)==set(witnesses),'exact survivors')
for (d1,d2,d3),(z1,z2,z3,t) in witnesses.items():
 check(d1*z1*z1-d2*z2*z2==t*t,'first rational equation')
 check(d3*z3*z3-d1*z1*z1==t*t,'second rational equation')
 check(any((z1,z2,z3,t)),'projective nonzero')
for a,b in product(survive,repeat=2):
 check(tuple(squareclass(x*y) for x,y in zip(a,b)) in survive,'Kummer subgroup')
# Derive all root images from root differences rather than stored signs.
roots=(0,1,-1);images=[(1,1,1)]
for root in roots:
 values=[]
 for i,e in enumerate(roots):
  if root!=e:values.append(squareclass(root-e))
  else:
   val=1
   for other in roots:
    if other!=e:val*=e-other
   values.append(squareclass(val))
 images.append(tuple(values))
check(set(images)==set(survive),'literal root Kummer images')
print(json.dumps({'status':'PASS','checks':checks,'projective_residue_vectors_tested':tested,'candidate_triples':16,'real_survivors':8,'selmer_classes':rows,'rational_images':images,'conclusion':'Sha[2^infinity]=0 for y^2=x^3-x using cited full 2-descent; no odd-primary conclusion'},sort_keys=True))
