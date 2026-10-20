#!/usr/bin/env python3
"""Exact tournament Hopf / connected sign-even graph primitive bridge.
No third-party dependencies. Default finite matrix census through order five.
"""
from __future__ import annotations
import argparse,itertools,json
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict

CHECKS=0
def check(p,msg):
 global CHECKS
 CHECKS+=1
 if not p: raise AssertionError(msg)

@lru_cache(None)
def pairs(n): return tuple(itertools.combinations(range(n),2))

@lru_cache(None)
def actions(n):
 ps=pairs(n); ix={p:i for i,p in enumerate(ps)}; ans=[]
 for perm in itertools.permutations(range(n)):
  dest=[]; rev=0
  for k,(i,j) in enumerate(ps):
   u,v=perm[i],perm[j]
   dest.append(1<<ix[tuple(sorted((u,v)))])
   if u>v: rev|=1<<k
  revout=sum(dest[k] for k in range(len(ps)) if rev>>k&1)
  ans.append((tuple(dest),rev,revout))
 return tuple(ans)

def transform(mask,dest):
 out=0
 while mask:
  low=mask&-mask; out|=dest[low.bit_length()-1]; mask-=low
 return out

@lru_cache(None)
def catalog(n,tournament):
 N=1<<len(pairs(n)); canon={}; reps=[]; evens=[]
 for mask in range(N):
  if mask in canon: continue
  reps.append(mask); even=True
  for dest,rev,revout in actions(n):
   new=transform(mask,dest)^(revout if tournament else 0)
   canon[new]=mask
   if not tournament and new==mask and (mask&rev).bit_count()%2: even=False
  if even: evens.append(mask)
 return tuple(reps),canon,tuple(evens)

@lru_cache(None)
def restrict(n,mask,vertices,tournament=True):
 ix={v:i for i,v in enumerate(vertices)}; pi={p:i for i,p in enumerate(pairs(len(vertices)))}; out=0
 for k,(i,j) in enumerate(pairs(n)):
  if i in ix and j in ix and mask>>k&1: out|=1<<pi[(ix[i],ix[j])]
 return catalog(len(vertices),tournament)[1][out]

def connected(n,mask):
 if n==0:return False
 adj=[0]*n
 for k,(i,j) in enumerate(pairs(n)):
  if mask>>k&1:adj[i]|=1<<j;adj[j]|=1<<i
 seen=1;old=0
 while old!=seen:
  old=seen
  for i in range(n):
   if seen>>i&1:seen|=adj[i]
 return seen==(1<<n)-1

def strong(n,mask):
 if n==0:return False
 reach=[1<<i for i in range(n)]
 for k,(i,j) in enumerate(pairs(n)):
  if mask>>k&1:reach[i]|=1<<j
  else:reach[j]|=1<<i
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 return all(r==(1<<n)-1 for r in reach)

@lru_cache(None)
def join(a,b):
 n,x=a;m,y=b; ix={p:i for i,p in enumerate(pairs(n+m))}; out=0
 for k,p in enumerate(pairs(n)):
  if x>>k&1:out|=1<<ix[p]
 for k,(i,j) in enumerate(pairs(m)):
  if y>>k&1:out|=1<<ix[(i+n,j+n)]
 for i in range(n):
  for j in range(n,n+m):out|=1<<ix[(i,j)]
 return n+m,catalog(n+m,True)[1][out]

def plus(a,b,scale=F(1)):
 c=defaultdict(F,a)
 for k,v in b.items(): c[k]+=scale*v
 return {k:v for k,v in c.items() if v}
def mul(a,b):
 c=defaultdict(F)
 for x,u in a.items():
  for y,v in b.items():c[join(x,y)]+=u*v
 return {k:v for k,v in c.items() if v}

@lru_cache(None)
def ordered_partitions(n):
 out=[()]
 for v in range(n):
  nex=[]
  for p in out:
   for i in range(len(p)):
    nex.append(p[:i]+(p[i]+(v,),)+p[i+1:])
   for i in range(len(p)+1):nex.append(p[:i]+((v,),)+p[i:])
  out=nex
 return tuple(out)

@lru_cache(None)
def eulerian(letter):
 n,mask=letter; out=defaultdict(F)
 for blocks in ordered_partitions(n):
  k=len(blocks); term=(0,0)
  for b in blocks:term=join(term,(len(b),restrict(n,mask,b)))
  out[term]+=F((-1)**(k-1),k)
 return {k:v for k,v in out.items() if v}

def coproduct(poly,reduced=False):
 out=defaultdict(F)
 for (n,mask),c in poly.items():
  for bits in range(1<<n):
   if reduced and bits in (0,(1<<n)-1):continue
   A=tuple(i for i in range(n) if bits>>i&1); B=tuple(i for i in range(n) if not bits>>i&1)
   a=(len(A),restrict(n,mask,A));b=(len(B),restrict(n,mask,B))
   out[(a,b)]+=c
 return {k:v for k,v in out.items() if v}

def walsh(n,g):
 out=defaultdict(F); den=1<<len(pairs(n))
 for t,rep in catalog(n,True)[1].items():out[(n,rep)]+=F((-1)**((g&~t).bit_count()),den)
 return {k:v for k,v in out.items() if v}

def is_lyndon(w): return bool(w) and all(w<w[i:] for i in range(1,len(w)))

def words_weight(alphabet,n,prefix=()):
 if not n:yield prefix;return
 for a in alphabet:
  if a[0]<=n:yield from words_weight(alphabet,n-a[0],prefix+(a,))

@lru_cache(None)
def lie(w):
 if len(w)==1:return eulerian(w[0])
 for i in range(1,len(w)):
  if is_lyndon(w[i:]):
   u=lie(w[:i]);v=lie(w[i:]);return plus(mul(u,v),mul(v,u),-1)
 raise AssertionError(w)

def coordinates(basis,vector,ambient):
 m=len(basis); rows=[[b.get(x,F(0)) for b in basis]+[vector.get(x,F(0))] for x in ambient]
 row=0;piv=[]
 for j in range(m):
  i=next((i for i in range(row,len(rows)) if rows[i][j]),None)
  check(i is not None,'basis independence')
  rows[row],rows[i]=rows[i],rows[row]; c=rows[row][j];rows[row]=[x/c for x in rows[row]]
  for i in range(len(rows)):
   if i!=row and rows[i][j]:
    c=rows[i][j];rows[i]=[x-c*y for x,y in zip(rows[i],rows[row])]
  piv.append(row);row+=1
 check(all(any(r[:-1]) or not r[-1] for r in rows),'primitive vector in Lie span')
 return [rows[r][-1] for r in piv]

def determinant(matrix):
 A=[row[:] for row in matrix];N=len(A);d=F(1)
 for k in range(N):
  i=next((i for i in range(k,N) if A[i][k]),None)
  check(i is not None,'primitive determinant nonzero')
  if i!=k:A[k],A[i]=A[i],A[k];d=-d
  pivot=A[k][k];d*=pivot
  for i in range(k+1,N):
   q=A[i][k]/pivot
   for j in range(k+1,N):A[i][j]-=q*A[k][j]
 return d

def matching(matrix):
 N=len(matrix)
 def feasible(rows,available):
  owner={}
  def aug(r,seen):
   for c in sorted(available):
    if matrix[r][c] and c not in seen:
     seen.add(c)
     if c not in owner or aug(owner[c],seen):owner[c]=r;return True
   return False
  return all(aug(r,set()) for r in rows)
 ans=[];available=set(range(N))
 for r in range(N):
  c=next((c for c in sorted(available) if matrix[r][c] and feasible(range(r+1,N),available-{c})),None)
  check(c is not None,'Hall support matching');ans.append(c);available.remove(c)
 return ans

def component_coproduct_walsh(n,g):
 out=defaultdict(F)
 for bits in range(1<<n):
  if any((g>>k&1) and ((bits>>i&1)!=(bits>>j&1)) for k,(i,j) in enumerate(pairs(n))):continue
  A=tuple(i for i in range(n) if bits>>i&1); B=tuple(i for i in range(n) if not bits>>i&1)
  # Use actual induced edge masks, not canonical graph masks: canonicalization
  # may change the orientation-line sign. The subsets retain increasing labels.
  def induced(vertices):
   oldix={p:k for k,p in enumerate(pairs(n))};ans=0
   for k,(i,j) in enumerate(pairs(len(vertices))):
    if g>>oldix[(vertices[i],vertices[j])]&1:ans|=1<<k
   return ans
  left=walsh(len(A),induced(A));right=walsh(len(B),induced(B))
  for u,x in left.items():
   for v,y in right.items():out[(u,v)]+=x*y
 return {k:v for k,v in out.items() if v}

def converse(poly):
 ans={}
 for (n,t),v in poly.items():ans[(n,catalog(n,True)[1][t^((1<<len(pairs(n)))-1)])]=v
 return ans

def run(max_n,graph_census_n):
 allrows=[];alphabet=[]
 for n in range(1,max(max_n,graph_census_n)+1):
  reps,_,_=catalog(n,True);s=[t for t in reps if strong(n,t)]
  cr=[g for g in catalog(n,False)[2] if connected(n,g)]
  ce=sum(g.bit_count()%2==0 for g in cr);co=len(cr)-ce
  print(f'CENSUS n={n} tournaments={len(reps)} strong={len(s)} connected_even={len(cr)} even_edges={ce} odd_edges={co}')
  if n>max_n:continue
  for g in catalog(n,False)[2]:
   check(coproduct(walsh(n,g))==component_coproduct_walsh(n,g),'full Walsh coalgebra identity')
  alphabet.extend((n,t) for t in s)
  ws=sorted(w for w in words_weight(alphabet,n) if is_lyndon(w))
  check(len(ws)==len(cr),'weighted Lyndon count')
  basis=[lie(w) for w in ws]
  for b in basis:check(not coproduct(b,True),'Lyndon primitive')
  for a in ((n,t) for t in s):
   check(not coproduct(eulerian(a),True),'Eulerian primitive')
   check(eulerian(a).get(a)==1,'Eulerian strong leading coefficient')
  ambient=[(n,t) for t in reps];cols=[]
  for g in cr:
   f=walsh(n,g);check(not coproduct(f,True),'connected Walsh primitive')
   check(converse(f)=={k:(-1)**g.bit_count()*v for k,v in f.items()},'converse edge parity')
   cols.append(coordinates(basis,f,ambient))
  mat=[list(row) for row in zip(*cols)] if cols else []
  det=determinant(mat)
  check(abs(det.numerator)&(abs(det.numerator)-1)==0 and det.denominator&(det.denominator-1)==0,'dyadic primitive determinant')
  print('DETERMINANT',n,det)
  match=matching(mat)
  for i,w in enumerate(ws):
   if n<=5:
    for j,g in enumerate(cr):
     if mat[i][j]:check((len(w)-1-g.bit_count())%2==0,'primitive parity support')
   print('MATCH',n,'word',w,'graph',cr[match[i]],'edges',cr[match[i]].bit_count(),'coefficient',mat[i][match[i]])
  allrows.append({'n':n,'strong_types':s,'connected_graph_masks':cr,'lyndon_words':ws,'matrix':[[str(x) for x in row] for row in mat],'matching':match,'determinant':str(det),'eulerian_generators':{str(t):{str(k):str(v) for k,v in eulerian((n,t)).items()} for t in s}})
 check(not strong(2,0),'no strong order two')
 check(bool(coproduct(walsh(2,0),True)),'disconnected Walsh image need not be primitive')
 check(walsh(3,0)!={(3,0):F(1)},'Walsh is not ordinal-join multiplicative')
 print(f'ALL CHECKS PASSED: {CHECKS}')
 return allrows

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--max-n',type=int,default=5);p.add_argument('--graph-census-n',type=int,default=6);p.add_argument('--json',default='primitive_walsh.json');a=p.parse_args()
 if not 1<=a.max_n<=6 or not a.max_n<=a.graph_census_n<=6:raise SystemExit('finite atlas bounds: 1 <= max_n <= graph_census_n <= 6')
 data=run(a.max_n,a.graph_census_n)
 with open(a.json,'w') as f:json.dump(data,f,indent=2);f.write('\n')
