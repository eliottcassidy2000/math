#!/usr/bin/env python3
"""Independent direct sparse-pair fibre census; producer module tree reused."""
from itertools import combinations,permutations
from collections import Counter
from pathlib import Path
import importlib.util,json
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('f',root/'triples/colored_fibers.py')
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
def check(b,m):
 if not b:raise RuntimeError(m)
def signature(a,n):
 return tuple(all(sum((a[v]>>w)&1 for w in s if w!=v)==1 for v in s) for s in combinations(range(n),3))
def predicted(t,arcs):
 kind,mask,kids=t
 if not kids:return 1
 local=[];below=[[] for _ in kids]
 for u,v in arcs:
  i=next(i for i,k in enumerate(kids) if k[1]>>u&1)
  j=next(i for i,k in enumerate(kids) if k[1]>>v&1)
  if i==j:below[i].append((u,v))
  else:local.append((i,j))
 val=1
 for k,edges in zip(kids,below):val*=predicted(k,edges)
 if kind=='prime':return val*(1 if local else 2)
 le=0
 for perm in permutations(range(len(kids))):
  pos={v:i for i,v in enumerate(perm)}
  le+=all(pos[i]<pos[j] for i,j in local)
 return val*le
def certificate(t,cols):
 kind,mask,kids=t
 edges=[]
 for child in kids:edges+=certificate(child,cols)
 candidates=list(combinations(range(len(kids)),2)) if kind=='prime' else [(i,i+1) for i in range(len(kids)-1)]
 if kind=='leaf':return edges
 for i,j in candidates:
  uv=next(((u,v) for u in range(len(cols)) for v in range(len(cols)) if kids[i][1]>>u&1 and kids[j][1]>>v&1 and cols[u]!=cols[v]),None)
  if uv:
   edges.append(tuple(sorted(uv)))
   if kind=='prime':break
  elif kind=='linear':raise RuntimeError('missing split')
 return edges
counts=[];designed=0
for n in range(1,6):
 pairs=list(combinations(range(n),2));data=[]
 for code in range(1<<len(pairs)):
  a=f.rows(n,code);data.append((code,a,signature(a,n),f.tree(a,(1<<n)-1)))
 masks=range(1<<len(pairs)) if n<=4 else sorted(set([0,(1<<len(pairs))-1]+[1<<j for j in range(len(pairs))]+[(i*137)%1024 for i in range(40)]))
 tested=0
 for mask in masks:
  direct=Counter((sig,code&mask) for code,a,sig,t in data)
  for code,a,sig,t in data:
   arcs=[(u,v) if a[u]>>v&1 else (v,u) for j,(u,v) in enumerate(pairs) if mask>>j&1]
   check(predicted(t,arcs)==direct[sig,code&mask],('product',n,code,mask));tested+=1
 for code,a,sig,t in data:
  cols=[-1]*n;f.engineered(t,cols);edges=certificate(t,cols)
  check(len(edges)<=n-1,'budget')
  mask=sum(1<<pairs.index(edge) for edge in edges)
  check(sum(sig==other and code&mask==c&mask for c,_,other,_ in data)==1,('certificate',n,code));designed+=1
 counts.append(dict(n=n,masks=len(masks),cases=tested))
print(json.dumps(dict(status='PASS',sparse_counts=counts,designed_certificates=designed,scope='All pair constraint subsets through n4; deterministic 51 masks at n5. Literal C3 grouping and permutation linear-extension count; producer tree/color routines reused.'),indent=2))
