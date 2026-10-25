#!/usr/bin/env python3
"""Exact local-coefficient audit; finite tests support, not replace, the proof."""
from itertools import combinations, permutations, product
from collections import Counter
import json
P=tuple(permutations(range(3))); ID=(0,1,2)
E=(0,3,5,6); ecoord={0:0,3:1,5:2,6:3}
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def compose(p,q):return tuple(p[q[i]] for i in range(3))
def inv(p):return tuple(p.index(i) for i in range(3))
def act(p,b):return sum(((b>>i)&1)<<p[i] for i in range(3))
def rank(columns):
 basis={}
 for c in columns:
  while c:
   i=c.bit_length()-1
   if i not in basis:basis[i]=c;break
   c^=basis[i]
 return len(basis)
def mat(n,es,ps,kind):
 dim={'V':3,'E':2,'scalar':1}[kind]
 bs={'V':(1,2,4),'E':(3,5),'scalar':(1,)}[kind]
 cols=[]
 for v in range(n):
  for b in bs:
   col=0
   for j,((u,w),p) in enumerate(zip(es,ps)):
    val=(b if w==v else 0)^(act(p,b) if u==v and kind!='scalar' else (b if u==v else 0))
    if kind=='E':val=ecoord[val]
    col|=val<<(dim*j)
   cols.append(col)
 return cols
def holonomy(n,es,ps):
 adj=[[] for _ in range(n)]
 for (u,v),p in zip(es,ps):adj[u].append((v,p));adj[v].append((u,inv(p)))
 T={0:ID};todo=[0]
 while todo:
  u=todo.pop()
  for v,p in adj[u]:
   if v not in T:T[v]=compose(p,T[u]);todo.append(v)
 if len(T)<n:return None
 gens=[compose(inv(T[v]),compose(p,T[u])) for (u,v),p in zip(es,ps)]
 H={ID};todo=[ID]
 while todo:
  h=todo.pop()
  for g in gens:
   t=compose(g,h)
   if t not in H:H.add(t);todo.append(t)
 return H
cases=0;hist=Counter()
for n in range(1,5):
 possible=list(combinations(range(n),2))
 for em in range(1<<len(possible)):
  es=[e for i,e in enumerate(possible) if em>>i&1]
  if holonomy(n,es,[ID]*len(es)) is None:continue
  for ps in product(P,repeat=len(es)):
   H=holonomy(n,es,ps)
   fix=[b for b in E if all(act(p,b)==b for p in H)]
   f=len(fix).bit_length()-1
   beta=len(es)-n+1
   re=rank(mat(n,es,ps,'E'));rv=rank(mat(n,es,ps,'V'));rs=rank(mat(n,es,ps,'scalar'))
   he=2*len(es)-re;hv=3*len(es)-rv
   need(he==2*(beta-1)+f,'twisted Euler formula')
   need(hv==he+beta,'equivariant cohomology split')
   need(rv==re+rs,'direct sum differential')
   need(3*n-rv==1+f,'affine solution dimension')
   need(f=={1:2,2:1,3:0,6:0}[len(H)],'complete subgroup fixed-space table')
   hist[(beta,len(H),he)]+=1;cases+=1
# All affine transports on a cycle after tree gauge; brute root frame counts.
cycle_records=[]
for sigma in P:
 for w in range(8):
  roots=[b for b in range(8) if b^act(sigma,b)==w]
  projected=(w.bit_count()%2==0)
  need((len(roots)>0)==(rank([b^act(sigma,b) for b in (1,2,4)]+[w])==rank([b^act(sigma,b) for b in (1,2,4)])),'brute affine solvability')
  if roots:
   need(projected,'full section implies parity section')
   need(sum(b.bit_count()%2==0 for b in roots)==len(roots)//2,'equivariant odd split balances frames')
  if sigma in ((0,1,2),(1,0,2),(1,2,0)):
   cycle_records.append({'sigma':sigma,'gain':w,'scalar_compatible':projected,'frame_count':len(roots)})
need(not [b for b in range(8) if b^act(ID,b)==3],'identity false positive')
trans=(1,0,2);rot=(1,2,0)
need(not [b for b in range(8) if b^act(trans,b)==5],'transposition false positive')
need(all(any(b^act(rot,b)==w for b in range(8)) for w in E),'three-cycle repairs every even gain')
# Two cycles: theta graph with edges12,13,14,23,24 after relabeling0..3.
es=[(0,1),(0,2),(0,3),(1,2),(1,3)];ps=[ID,ID,ID,rot,rot]
need(2*len(es)-rank(mat(4,es,ps,'E'))==2,'two cycles retain two obstructions')
# All48 affine cube maps preserve adjacency and act on Fano halves by parity.
for sigma in P:
 for omega in range(8):
  f=lambda x: omega^act(sigma,x)
  need(len({f(x) for x in range(8)})==8,'affine cube bijection')
  need(all(((f(x)^f(y)).bit_count()==1)==((x^y).bit_count()==1) for x,y in combinations(range(8),2)),'cube adjacency preserved')
  need(all(f(x).bit_count()%2==(x.bit_count()+omega.bit_count())%2 for x in range(8)),'Fano halves transform by parity')
for omega in range(1,8):
 need(not any(x^omega==x for x in range(8)),'translation has no fixed cube state')
print(json.dumps({'status':'PASS','checks':checks,'connected_graph_transport_cases':cases,'universe':'all connected labeled simple graphs on1..4 vertices; every S3 edge transport; all48 affine one-cycle transports','dimension_census':[{'beta':b,'holonomy_order':h,'lost_dimension':d,'cases':v} for (b,h,d),v in sorted(hist.items())],'cycle_hostiles':cycle_records},indent=2))
