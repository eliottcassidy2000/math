#!/usr/bin/env python3
import json
from pathlib import Path

def need(x):
 if not x:raise RuntimeError('audit failure')
Q=[];P=[]
for r in range(100):
 # Build shells from ordered lists rather than rank/isqrt formulas.
 a=[(r,j) for j in range(r+1)]+[(i,r) for i in range(r-1,-1,-1)]
 if r%2:a=[(y,x) for x,y in a]
 Q+=a
 m=r//2
 if r%2==0:
  P += [(-m,y) for y in range(m,-m-1,-1)]+[(x,-m) for x in range(-m+1,m+1)]
 else:
  P += [(m+1,y) for y in range(-m,m+2)]+[(x,m+1) for x in range(m,-m-1,-1)]
 need(len(set(Q))==len(Q)==(r+1)**2)
 need(set(P)=={(x,y) for x in range(-m,(r+1)//2+1) for y in range(-m,(r+1)//2+1)})
for path in [Q,P]:need(all(abs(a-c)+abs(b-d)==1 for (a,b),(c,d) in zip(path,path[1:])))
mapQP=dict(zip(Q,P));bad=None
for q,p in mapQP.items():
 q2=(q[0],q[1]+1)
 if q2 in mapQP:
  p2=mapQP[q2]
  if sum(abs(a-b) for a,b in zip(p,p2))!=1:
   bad={'quadrant_doubling_edge':[q,q2],'plane_positions':[p,p2]};break
need(bad is not None)
labels=[(2*j+1)*2**k for j,k in Q];need(len(set(labels))==10000)
cycle=[1,3,6,4,2,1]
need(all(abs(a-b)==2 or a==2*b or b==2*a for a,b in zip(cycle,cycle[1:])))
out={'status':'PASS','complete_shells':100,'path_vertices_each':10000,'all_dyadic_edges_not_preserved':bad,'literal_plus2_doubling_odd_cycle':cycle}
Path(__file__).with_name('PLANE_RESULT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
