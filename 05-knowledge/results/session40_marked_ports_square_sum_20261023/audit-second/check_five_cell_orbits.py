#!/usr/bin/env python3
import itertools,json
from pathlib import Path
P=list(itertools.permutations(range(5)))
def gauge(p,a,b):return tuple((p[(i+a)%5]+b)%5 for i in range(5))
def canon(p):return min(gauge(p,a,b) for a in range(5) for b in range(5))
def converse(p):
 inv=[p.index(i) for i in range(5)];return tuple((-inv[(-i)%5])%5 for i in range(5))
def parity(p):return sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2
orbits={}
for p in P:orbits.setdefault(canon(p),[]).append(p)
if sorted(map(len,orbits.values()))!=[5]*4+[25]*4:raise RuntimeError('orbit sizes')
rows=[]
for p,o in sorted(orbits.items()):
 if len({parity(q) for q in o})!=1:raise RuntimeError('parity gauge')
 q=canon(converse(p));J=sum(1 if (p[j]-p[i])%5 in(1,2) else -1 for i in range(5) for j in range(5) if(j-i)%5 in(1,2));rows.append({'representative':p,'size':len(o),'parity':parity(p),'converse_representative':q,'J':J,'J_mod8':J%8,'J_mod16':J%16})
merged={min(p,canon(converse(p))) for p in orbits}
for row in rows:
 p=row['representative'];A=[set() for _ in range(10)]
 for i in range(5):
  for j in range(5):
   if (j-i)%5 in (1,2):A[i].add(j);A[5+i].add(5+j)
   if j==p[i]:A[i].add(5+j)
   else:A[5+j].add(i)
 def count(v,mask):
  if mask==1023:return 1
  return sum(count(w,mask|(1<<w)) for w in A[v] if not(mask>>w&1))
 row['H']=sum(count(v,1<<v) for v in range(10))
if [row['H'] for row in rows]!=[8265,8245,8221,8201,8193,8175,8175,8245]:raise RuntimeError('path counts')
require = lambda b: None if b else (_ for _ in ()).throw(RuntimeError('J classification'))
for a,b in itertools.combinations(rows,2):
 same=min(tuple(a['representative']),tuple(a['converse_representative']))==min(tuple(b['representative']),tuple(b['converse_representative']))
 require((a['J']==b['J'])==same);require((a['J_mod16']==b['J_mod16'])==same)
require(len({a['J_mod8'] for a in rows})==4)
out={'status':'PASS','ordinary_types':len(orbits),'merged_types':len(merged),'rows':rows}
Path(__file__).with_name('five-cell-orbits-audit-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
