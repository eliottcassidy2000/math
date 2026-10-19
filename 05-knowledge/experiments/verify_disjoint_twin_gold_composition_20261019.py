#!/usr/bin/env python3
"""Exact skeleton and correlated two-gadget controls; no independence assumption."""
from collections import defaultdict,Counter
from itertools import combinations,permutations
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
stats=Counter();digest=hashlib.sha256()
for n in range(3,7):
 for p in posets(n):
  stats['posets']+=1;U=tuple(sum(1<<v for v in range(n) if p[v]>>u&1) for u in range(n));words=extensions(p);ranks=[{v:i for i,v in enumerate(w)} for w in words];gadgets=[]
  for dual in [False,True]:
   D,V=(U,p) if dual else (p,U)
   for u in range(n):
    for v in range(n):
     if u==v or D[u]!=D[v] or V[v]&~V[u]:continue
     extra=V[u]&~V[v]
     if extra.bit_count()!=1:continue
     t=extra.bit_length()-1;triple=frozenset([u,v,t]);fibres=defaultdict(lambda:[0,0,0])
     for w,r in zip(words,ranks):
      rr={a:n-1-b for a,b in r.items()} if dual else r
      order=tuple(sorted(triple,key=rr.get));index={(u,v,t):0,(u,t,v):1,(v,u,t):2}[order]
      skeleton=tuple(-1 if a in triple else a for a in w);fibres[skeleton][index]+=1
     for skeleton,(A,B,C) in fibres.items():
      require(A==C and B<=A,'Skeleton invariance certificate failed')
      stats['individual_skeletons']+=1;digest.update(repr((p,dual,u,v,t,skeleton,A,B,C)).encode())
     stats['gadgets']+=1
     # In the dual, both orientation bits are reversed; counts are unchanged.
     gadgets.append((triple,(u,v),(v,t)))
  for g,h in permutations(gadgets,2):
   if g[0]&h[0]:continue
   stats['ordered_disjoint_gadget_pairs']+=1
   comparisons=[g[1],g[2],h[1],h[2]];C=Counter(tuple(r[a]<r[b] for a,b in comparisons) for r in ranks);E=len(words)
   for path,c4 in C.items():
    c1=sum(c for k,c in C.items() if k[:1]==path[:1]);c2=sum(c for k,c in C.items() if k[:2]==path[:2]);c3=sum(c for k,c in C.items() if k[:3]==path[:3])
    require(E>=c1+c2,'First Gold block failed');require(c2>=c3+c4,'Conditional second Gold block failed');require(4*c4<=E,'Quarter-survivor bound failed')
    stats['four_answer_paths']+=1
# Actual independent-gadget obstruction to taking just the first comparison from each.
p=(0,0,1,0,0,8);w=extensions(p);r=[{a:i for i,a in enumerate(v)} for v in w]
E=len(w);heavy=sum(q[0]<q[1] for q in r);both=sum(q[0]<q[1] and q[3]<q[4] for q in r)
require((E,heavy,both)==(180,120,80),'Naive-composition control mismatch');require(heavy+both>E,'Naive-composition failure absent')
# A correlated example: the two B events cannot coexist.
p=close((0,0,(1<<0)|(1<<3)|(1<<4),0,0,(1<<0)|(1<<1)|(1<<3)))
w=extensions(p);r=[{a:i for i,a in enumerate(v)} for v in w]
B1=sum(q[0]<q[2]<q[1] for q in r);B2=sum(q[3]<q[5]<q[4] for q in r);joint=sum(q[0]<q[2]<q[1] and q[3]<q[5]<q[4] for q in r)
require(len(w)==60 and B1>0 and B2>0 and joint==0,'Correlated control failed')
print(json.dumps({'status':'EXACT_CONTROLS_PASSED','scope':'Finite controls through six vertices; arbitrary finite r follows from skeleton-preserving injections and induction.','stats':dict(stats),'digest':digest.hexdigest(),'naive_crossed_first_comparisons':{'extensions':E,'first_survivors':heavy,'second_survivors':both,'sum':heavy+both},'correlated_example':{'pred':p,'extensions':len(w),'B1':B1,'B2':B2,'joint_B1_B2':joint}},indent=2))
