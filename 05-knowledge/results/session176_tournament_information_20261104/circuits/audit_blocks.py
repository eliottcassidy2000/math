#!/usr/bin/env python3
import itertools,json
from explore import kernel

def req(p,m):
 if not p:raise RuntimeError(m)

def block_tournament(k,mask):
 n=3*k;a=[0]*n
 for u,v in itertools.combinations(range(n),2):
  if u//3<v//3: a[u]|=1<<v
  else:
   edge=(v-u==1)
   if mask>>(u//3)&1:edge=not edge
   if edge:a[u]|=1<<v
   else:a[v]|=1<<u
 return a

def literal_K(a):
 n=len(a);K=[[0]*n for _ in a]
 for p in itertools.permutations(range(n)):
  for s in range(1,n):
   if all(a[p[i]]>>p[i+1]&1 for i in range(n-1) if i!=s-1):K[p[s-1]][p[s]]+=1
 return K

rows=[]
for k in range(1,4):
 signatures=set()
 for mask in range(1<<k):
  a=block_tournament(k,mask);H,K=kernel(a)
  req(H==3**k,'H')
  if k<=2:req(K==literal_K(a),'literal cover')
  n=len(a);W={}
  for u,v in itertools.combinations(range(n),2):
   w=(K[u][v]+K[v][u])//2;W[u,v]=w
   if u//3==v//3:req(w==3**(k-1),'internal')
   else:
    t=v//3-u//3-1
    req(K[u][v]==(3**(k-2) if t==0 else 0),'cross forward')
    req(K[v][u]==3**k*4**t,'cross reverse')
  deg=[sum(W[min(u,v),max(u,v)] for v in range(n) if a[u]>>v&1) for u in range(n)]
  signatures.add((H,tuple(W.values()),tuple(deg)))
 req(len(signatures)==1,'response collision')
 rows.append(dict(blocks=k,orientations=1<<k,distinct_signatures=len(signatures)))
print(json.dumps(dict(status='PASS',families=rows,independent_literal_covers_through_blocks=2),indent=2))
