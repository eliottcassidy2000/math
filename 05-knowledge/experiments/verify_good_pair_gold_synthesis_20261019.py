#!/usr/bin/env python3
"""Independent finite controls for the attributed good-pair Gold synthesis."""
from collections import defaultdict,Counter
from itertools import combinations,product
from functools import lru_cache
import hashlib,json

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def close(p):
 p=list(p)
 for k in range(len(p)):
  for v in range(len(p)):
   if p[v]>>k&1:p[v]|=p[k]
 return tuple(p)
def posets(n):
 edges=list(combinations(range(n),2));seen=set()
 for mask in range(1<<len(edges)):
  p=[0]*n
  for k,(u,v) in enumerate(edges):
   if mask>>k&1:p[v]|=1<<u
  q=close(p)
  if q not in seen:seen.add(q);yield q
def extensions(p):
 out=[];n=len(p)
 def rec(mask,w):
  if len(w)==n:out.append(tuple(w));return
  for v in range(n):
   if not mask>>v&1 and not p[v]&~mask:rec(mask|1<<v,w+[v])
 rec(0,[]);return out

def check(q):
 E=sum(q);require(E>0 and 2*q[0]<=E and all(q[i]>=q[i+1]>=0 for i in range(len(q)-1)),'Invalid law')
 if 2*q[0]==E:cuts=(0,0)
 else:
  prefix=0
  for ell,x in enumerate(q):
   if 2*(prefix+x)>=E:break
   prefix+=x
  require(1<=ell<len(q)-1 and x<=prefix,'Crossing boundary')
  cuts=(ell-1,ell) if 2*prefix+x>=E else (ell,ell-1)
 C=[[0,0],[0,0]]
 for i,w in enumerate(q):C[i<=cuts[0]][i<=cuts[1]]+=w
 for row in C:
  for cell in row:
   if cell:require(E>=sum(row)+cell,'Gold failed')
 return cuts
stats=Counter();digest=hashlib.sha256()
for m in range(2,8):
 for q in product(range(5),repeat=m):
  if not sum(q) or 2*q[0]>sum(q) or any(q[i]<q[i+1] for i in range(m-1)):continue
  check(q);stats['integer_monotone_laws']+=1
for n in range(2,7):
 for p in posets(n):
  stats['posets']+=1;words=extensions(p);wordset=set(words)
  U=tuple(sum(1<<v for v in range(n) if p[v]>>u&1) for u in range(n))
  for dual in (False,True):
   D,V=(U,p) if dual else (p,U)
   oriented=[tuple(reversed(w)) if dual else w for w in words];allowed=set(oriented)
   def ischain(T):return all(D[a]>>b&1 or D[b]>>a&1 for a,b in combinations(T,2))
   for a in range(n):
    for b in range(n):
     if a==b or (D[a]|V[a])>>b&1 or D[a]&~D[b]:continue
     T=[v for v in range(n) if (V[b]&~V[a])>>v&1]
     if not ischain(T):continue
     T.sort(key=lambda v:D[v].bit_count());chain=[b]+T;S=set(chain+[a]);fibres=defaultdict(lambda:[0]*(len(chain)+1))
     for w in oriented:
      r={v:i for i,v in enumerate(w)};j=sum(r[c]<r[a] for c in chain)
      fibres[tuple(-1 if v in S else v for v in w)][j]+=1
      if j:
       z=chain[j-1];out=list(w);out[r[a]],out[r[z]]=out[r[z]],out[r[a]]
       require(tuple(out) in allowed,'Invalid swap')
       rr={v:i for i,v in enumerate(out)}
       require(sum(rr[c]<rr[a] for c in chain)==j-1,'Wrong target bin')
       stats['swap_checks']+=1
     q=[0]*(len(chain)+1);weighted=q.copy()
     for sk,f in fibres.items():
      require(all(f[i]>=f[i+1] for i in range(len(f)-1)),'Skeleton monotonicity')
      stats['skeleton_laws']+=1
      weight=1+sum((i+1)*(v+2) for i,v in enumerate(sk))%7
      for i,x in enumerate(f):q[i]+=x;weighted[i]+=weight*x
     stats['structural_orientations']+=1
     if 2*q[0]<=sum(q):
      cuts=check(q);stats['good_orientations']+=1;digest.update(repr((p,dual,a,b,q,cuts)).encode())
     if 2*weighted[0]<=sum(weighted):check(weighted);stats['weighted_good_orientations']+=1
     other=[v for v in range(n) if (V[a]&~V[b])>>v&1]
     if a<b and D[a]==D[b] and ischain(other):
      # At least one orientation has its forward count at most one half.
      aa,bb=(a,b) if 2*q[0]<=sum(q) else (b,a)
      extras=[v for v in range(n) if (V[bb]&~V[aa])>>v&1];extras.sort(key=lambda v:D[v].bit_count())
      ch=[bb]+extras;qq=[0]*(len(ch)+1)
      for w in oriented:
       r={v:i for i,v in enumerate(w)};qq[sum(r[c]<r[aa] for c in ch)]+=1
      check(qq);stats['very_good_pairs']+=1
# Distinguishing full-upset hypothesis: labels u,v,s,t,w.
p=(0,0,1,5,3);w=extensions(p);N=sum(z.index(1)<z.index(0) for z in w)
require(2*N<len(w),'Example orientation condition')
print(json.dumps({'status':'EXACT_CONTROLS_PASSED','scope':'Symbolic arbitrary-size theorem; finite naturally labelled controls through six vertices, both dual orientations; nonnegative integer laws of length two through seven with entries zero through four.','stats':dict(stats),'digest':digest.hexdigest(),'full_upset_distinction':{'pred':p,'extensions':len(w),'v_before_u_count':N}},indent=2))
