#!/usr/bin/env python3
from itertools import combinations
from math import isqrt
import json
checks=0

def req(v,msg):
 global checks
 checks+=1
 if not v:raise RuntimeError(msg)

rows=[]
for p in [3,5,11,13,19,29]:
 req(2 not in {z*z%p for z in range(p)},'nonresidue')
 for n in range(1,25):
  L=p
  while L<=4*n+2:L*=p*p
  qs=[p*i+1 for i in range(n)]
  labels=[4*L*q for q in qs]
  target=set()
  for i,j in combinations(range(n),2):
   x=(L*(qs[j]-qs[i])-1)**2-labels[i]
   k=len(labels);labels.append(x);target.update([(i,k),(j,k)])
  req(len(labels)==len(set(labels)) and min(labels)>0,'collision/positivity')
  req(max(labels)<=36*p**6*n**4,'height')
  for i,j in combinations(range(len(labels)),2):
   z=labels[i]+labels[j]
   req((isqrt(z)**2==z)==((i,j) in target),('adjacency',p,n,i,j))
  rows.append({'prime':p,'n':n,'vertices':len(labels),'height':max(labels)})
print(json.dumps({'status':'PASS','always_active_checks':checks,'complete_graph_rows':rows},sort_keys=True,indent=2))
