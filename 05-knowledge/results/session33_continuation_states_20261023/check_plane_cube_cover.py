#!/usr/bin/env python3
import itertools,json
from fractions import Fraction as F
checks=0
def ck(x,m):
 global checks
 checks+=1
 if not x:raise ValueError(m)
v=[(F(1-t*t,1+t*t),F(2*t,1+t*t)) for t in range(1,8)]
p=[tuple(sum((v[i][j] for i in range(7) if m>>i&1),F(0)) for j in range(2)) for m in range(128)]
ck(len(set(p))==128,'distinct positions');edges=0
for a,b in itertools.combinations(range(128),2):
 unit=sum((p[a][j]-p[b][j])**2 for j in range(2))==1
 ck(unit==((a^b).bit_count()==1),'induced unit-distance Q7')
 edges+=unit
normals=[4,1,3,7,6,5,2]
def s(m):
 o=0
 for i,a in enumerate(normals):
  if m>>i&1:o^=a
 return o
for m in range(128):
 ck({s(m^(1<<i)) for i in range(7)}==set(range(8))-{s(m)},'K8 regular cover')
 ck(all(m.bit_count()%2!=(m^(1<<i)).bit_count()%2 for i in range(7)),'proper two-coloring')
# Triangle 0 -> normal0 -> normal0 XOR normal1 -> 0 lifts open.
a,b=normals[:2];c=a^b;i=normals.index(c);end=1|2|(1<<i)
ck(s(end)==0 and end!=0 and end.bit_count()==3,'triangle holonomy')
print(json.dumps({'status':'PASS','checks':checks,'vertices':128,'edges':edges,'directions':[[str(x),str(y)] for x,y in v],'source_chromatic':2,'quotient_chromatic':8,'triangle_lift_endpoint':end,'scope':'exact induced unit-distance realization of Q7; local bijection to K8; no HN lower bound'},indent=2,sort_keys=True))
