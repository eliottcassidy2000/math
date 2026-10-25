#!/usr/bin/env python3
from itertools import product
import json

def need(v,m):
 if not v:raise RuntimeError(m)
def rank(cs):
 b={}
 for c in cs:
  while c:
   h=c.bit_length()-1
   if h not in b:b[h]=c;break
   c^=b[h]
 return len(b)
records=[]
for k in (2,3,4):
 points=list(product(range(3),repeat=k));idx={x:i for i,x in enumerate(points)};n=len(points)
 translations=[];cycles=[]
 for a in points:
  perm=[idx[tuple((x[j]+a[j])%3 for j in range(k))] for x in points]
  translations.append(perm)
  unseen=set(range(n));c=0
  while unseen:
   x=min(unseen);c+=1
   while x in unseen:unseen.remove(x);x=perm[x]
  cycles.append(c)
 gens=[]
 for j in range(k):
  a=tuple(int(i==j) for i in range(k));gens.append(translations[idx[a]])
 # Bouquet coboundary V -> V^k with blocks(I+g).
 cols=[sum(((1<<i)^(1<<g[i]))<<(j*n) for j,g in enumerate(gens)) for i in range(n)]
 need(rank(cols)==n-1,'one common invariant')
 hv=k*n-rank(cols)
 # E basis e_i+e_last; output E coordinates are first n-1 bits per block.
 ecols=[]
 for i in range(n-1):
  c=0
  for j,g in enumerate(gens):
   y=(1<<i)^(1<<(n-1))^(1<<g[i])^(1<<g[n-1])
   need(y.bit_count()%2==0,'augmentation preserved')
   c|=(y&((1<<(n-1))-1))<<(j*(n-1))
  ecols.append(c)
 he=k*(n-1)-rank(ecols)
 need(min(cycles)==3**(k-1),'minimum individual fixed space')
 need(hv==n*(k-1)+1 and he==(n-1)*(k-1),'gluing dimensions')
 records.append({'k':k,'gap_dimension':n,'common_fixed_dimension':1,'cancellation_order':min(cycles),'full_H1':hv,'forgotten_H1':he})
print(json.dumps({'status':'PASS','universe':'regular C3^k actions for k2,3,4; all group elements; direct bouquet coboundaries overF2','records':records},indent=2))
