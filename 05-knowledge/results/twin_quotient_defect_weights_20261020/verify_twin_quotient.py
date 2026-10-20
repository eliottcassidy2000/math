#!/usr/bin/env python3
"""Exact weighted-quotient controls; no repository imports or mutations."""
from itertools import combinations,permutations,product
from fractions import Fraction
from pathlib import Path
import json
checks=0

def check(ok,msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)

def graph_edges(q,mask):
 return {e for i,e in enumerate(combinations(range(q),2)) if mask>>i&1}

def autos(q,E,c):
 return [p for p in permutations(range(q)) if all(c[i]==c[p[i]] for i in range(q)) and {tuple(sorted((p[i],p[j]))) for i,j in E}==E]

def sign(E,c,p):
 return -1 if sum(c[i]*c[j] for i,j in E if p[i]>p[j])%2 else 1

def transform(a,p):
 b=[0]*len(a)
 for i,j in enumerate(p): b[j]=a[i]
 return tuple(b)

def calc(E,c):
 L=autos(len(c),E,c); ks=tuple(x//3 for x in c)
 chi={p:sign(E,c,p) for p in L}
 even=all(v==1 for v in chi.values())
 seen=set(); b=0; orbit_sum=0
 for a in product(*(range(k+1) for k in ks)):
  if a in seen:continue
  seen.update(transform(a,p) for p in L)
  stabilizer=[p for p in L if transform(a,p)==a]
  if all(chi[p]==1 for p in stabilizer):
   orbit_sum+=(-1)**sum(a)
   if sum(a):b+=(-1)**(sum(a)+1)
 total=0
 for p in L:
  used=set(); value=chi[p]
  for i in range(len(c)):
   if i in used:continue
   cyc=[];j=i
   while j not in used:used.add(j);cyc.append(j);j=p[j]
   k=ks[i];value*= k+1 if len(cyc)%2==0 else int(k%2==0)
  total+=value
 D=Fraction(total,len(L))
 check(D.denominator==1,'signed average integral')
 check(D==orbit_sum,'Burnside equals direct packing-vector orbit sum')
 check(b==int(even)-D,'closed b formula')
 if not any(ks):check(b==0,'no triple boundary')
 return dict(b=b,e=int(even and any(ks)),D=int(D),quotient_aut=len(L))

cases=0
# All false-twin-free quotient graphs up to four vertices, selected weights.
for q in range(1,5):
 for mask in range(1<<(q*(q-1)//2)):
  E=graph_edges(q,mask)
  nb=[frozenset(j if i==v else i for i,j in E if v in (i,j)) for v in range(q)]
  if len(set(nb))<q:continue
  for c in product((1,3,4,6),repeat=q):
   calc(E,c);cases+=1
# Full original and contracted graph sign computed on explicit vertex permutations.
contractions=0
for c in ((3,3),(6,6),(3,3,6),(3,3,6,12),(1,3,4),(3,4,6)):
 q=len(c);E=set(combinations(range(q),2)) if q!=4 else {(0,1),(0,3),(1,3),(2,3)}
 offsets=[];v=0
 for size in c:offsets.append(v);v+=size
 edges=[(offsets[i]+u,offsets[j]+w) for i,j in E for u in range(c[i]) for w in range(c[j])]
 for a in product(*(range(x//3+1) for x in c)):
  blocks=[]
  for i,size in enumerate(c):
   blocks.extend(tuple(offsets[i]+3*t+j for j in range(3)) for t in range(a[i]))
   blocks.extend((offsets[i]+t,) for t in range(3*a[i],size))
  owner={v:i for i,B in enumerate(blocks) for v in B}
  CE={tuple(sorted((owner[u],owner[w]))) for u,w in edges}
  for p in autos(q,E,c):
   if transform(a,p)!=a:continue
   g={offsets[i]+u:offsets[p[i]]+u for i in range(q) for u in range(c[i])}
   alpha={i:owner[g[B[0]]] for i,B in enumerate(blocks)}
   old=(-1)**sum(g[u]>g[w] for u,w in edges)
   new=(-1)**sum(alpha[u]>alpha[w] for u,w in CE)
   check(old==new,'explicit odd-module contraction sign');contractions+=1
balanced=[]
for c in range(1,61):
 z=calc({(0,1)},(c,c));k=c//3;s=int(k%2==0)
 check(z['D']==Fraction(s+(-1)**c*(k+1),2),'balanced complete bipartite formula')
 balanced.append(dict(c=c,**z))
family=[]
for r in range(1,21):
 X=calc({(0,1)},(6*r,6*r))
 if r==1:Y=calc({(0,1),(0,2),(1,2)},(3,3,6))
 else:Y=calc({(0,1),(0,3),(1,3),(2,3)},(3,3,6*r-6,6*r))
 check(X['b']-X['e']==-(r+1),'unbounded negative defect')
 check(Y['b']-Y['e']==1,'split join positive defect remains one')
 family.append(dict(r=r,source=X,target=Y))
result=dict(status='PASS',checks=checks,weighted_quotient_cases=cases,explicit_contractions=contractions,balanced_family=balanced,split_join_family=family)
out=json.dumps(result,indent=2,sort_keys=True)+'\n'
Path(__file__).with_name('verification.json').write_text(out)
print(json.dumps(dict(status='PASS',checks=checks,weighted_quotient_cases=cases,explicit_contractions=contractions),sort_keys=True))
