#!/usr/bin/env python3
"""Exact integer-kernel + modular-rank audit of guarded cube leakage."""
from itertools import combinations,product
from random import Random
import json
CHECKS=0
def ck(x,m):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(m)
def components(n,edges):
 adj=[set() for _ in range(n)]
 for u,v in edges:adj[u].add(v);adj[v].add(u)
 out=[];unseen=set(range(n))
 while unseen:
  v=unseen.pop();color={v:1};q=[v];ok=True
  while q:
   u=q.pop()
   for w in adj[u]:
    if w not in color:color[w]=-color[u];unseen.discard(w);q.append(w)
    elif color[w]==color[u]:ok=False
  out.append((color,ok))
 return out

def rank_mod(rows,p=1000003):
 piv={}
 for row in rows:
  v=[x%p for x in row]
  while any(v):
   j=next(i for i,x in enumerate(v) if x)
   if j not in piv:
    z=pow(v[j],-1,p);piv[j]=[(x*z)%p for x in v];break
   z=v[j];v=[(a-z*b)%p for a,b in zip(v,piv[j])]
 return len(piv)

def inspect(n,edges):
 legal=[s for s in range(1<<n) if all(not(s>>u&1 and s>>v&1) for u,v in edges)];S=set(legal)
 bad=sorted({s^(1<<v) for s in legal for v in range(n)}-S)
 C=[[int((s^z).bit_count()==1) for z in bad] for s in legal]
 cols=list(zip(*C)) if bad else []
 ck(all(sum(c) in (1,2) for c in cols),'one/two predecessor column')
 comps=components(n,edges);bip=[c for c,ok in comps if ok];b=len(bip)
 # Every basis vector is checked against literal incidence C, over integers.
 basis=[]
 for mask in range(1<<b):
  chosen=[bip[i] for i in range(b) if mask>>i&1];v=[]
  for s in legal:
   total=1;seen=0
   for c in chosen:
    hit=[u for u in c if s>>u&1]
    if len(hit)!=1:total=0;break
    total*=c[hit[0]];seen|=1<<hit[0]
   v.append(total if seen==s else 0)
  ck(any(v),'nonzero kernel basis vector')
  for col in cols:ck(sum(a*x for a,x in zip(v,col))==0,'exact left kernel')
  basis.append(v)
 ck(all(not any(x and y for x,y in zip(u,v)) for u,v in combinations(basis,2)),'disjoint kernel supports')
 rank=rank_mod(C)
 ck(rank==len(legal)-(1<<b),'rank matches integer nullity bound')
 e=sum((s^(1<<v)) in S for s in legal for v in range(n))//2
 trace=sum(sum(row) for row in C)
 ck(trace==n*len(legal)-2*e,'trace exit degree identity')
 return dict(states=len(legal),rank=rank,nullity=1<<b,trace=trace,edges=e)

def main():
 total=0;by_n={}
 for n in range(6):
  pairs=list(combinations(range(n),2));count=0
  for mask in range(1<<len(pairs)):
   inspect(n,[e for i,e in enumerate(pairs) if mask>>i&1]);count+=1
  by_n[n]=count;total+=count
 rng=Random(913)
 for n in range(6,10):
  for _ in range(30):inspect(n,[e for e in combinations(range(n),2) if rng.randrange(4)==0]);total+=1
 cycles=[]
 for n in range(3,11):
  x=inspect(n,sorted({tuple(sorted((i,(i+1)%n))) for i in range(n)}));x['n']=n
  a,b=0,1
  for _ in range(n-1):a,b=b,a+b
  ck(x['edges']==n*a,'cycle edge formula')
  cycles.append(x)
 print(json.dumps(dict(status='PASS',checks=CHECKS,graph_cases=total,exhaustive=by_n,cycles=cycles),indent=2))
if __name__=='__main__':main()
