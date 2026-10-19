#!/usr/bin/env python3
"""Independent finite controls for a proof valid at every finite size."""
from itertools import combinations,permutations
from pathlib import Path
from collections import Counter
import json,hashlib

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def close(p):
 p=list(p)
 for k in range(len(p)):
  for j in range(len(p)):
   if p[j]>>k&1:p[j]|=p[k]
 return tuple(p)
def posets(n):
 edges=list(combinations(range(n),2));seen=set()
 for mask in range(1<<len(edges)):
  p=[0]*n
  for k,(a,b) in enumerate(edges):
   if mask>>k&1:p[b]|=1<<a
  q=close(p)
  if q not in seen:seen.add(q);yield q
def extensions(p):
 n=len(p);out=[]
 def rec(mask,w):
  if len(w)==n:out.append(tuple(w));return
  for v in range(n):
   if not mask>>v&1 and not p[v]&~mask:rec(mask|1<<v,w+[v])
 rec(0,[]);return out
stats=[];digest=hashlib.sha256()
for n in range(3,7):
 total=cases=dual_cases=0
 for p in posets(n):
  total+=1;U=tuple(sum(1<<j for j in range(n) if p[j]>>i&1) for i in range(n));words=extensions(p)
  for reverse in [False,True]:
   D,V=(U,p) if reverse else (p,U)
   for u in range(n):
    for v in range(n):
     if u==v or D[u]!=D[v] or V[v]&~V[u]:continue
     extra=V[u]&~V[v]
     if extra.bit_count()!=1:continue
     t=extra.bit_length()-1;C=Counter()
     for word in words:
      r={a:i for i,a in enumerate(reversed(word) if reverse else word)}
      triple=tuple(sorted([u,v,t],key=r.get));C[triple]+=1
     A=C[u,v,t];B=C[u,t,v];T=C[v,u,t];E=len(words)
     require(A==T,'Twin swap equality failed');require(0<B<=A,'Cycle injection inequality failed');require(E==2*A+B,'Three-order exhaustion failed')
     require(min(E-(A+B)-A,E-(A+B)-B,E-T-T)>=0,'Gold path inequality failed')
     # Verify the two explicit maps directly on every relevant word.
     wordset=set(tuple(reversed(w)) if reverse else w for w in words)
     images=set()
     for word in wordset:
      r={a:i for i,a in enumerate(word)}
      if r[u]<r[t]<r[v]:
       image=list(word);image[r[u]]=v;image[r[t]]=u;image[r[v]]=t;image=tuple(image)
       require(image in wordset,'Three-cycle image invalid');require(image not in images,'Three-cycle collision');images.add(image)
      if r[u]<r[v]<r[t]:
       image=list(word);image[r[u]],image[r[v]]=image[r[v]],image[r[u]]
       require(tuple(image) in wordset,'Twin-swap image invalid')
     cases+=not reverse;dual_cases+=reverse;digest.update(repr((p,reverse,u,v,t,A,B,T)).encode())
 stats.append({'vertices':n,'naturally_labeled_posets':total,'extra_successor_markings':cases,'dual_extra_predecessor_markings':dual_cases})
 require(total=={3:7,4:40,5:357,6:4824}[n],'Poset census mismatch')
result={'status':'EXACT_FINITE_CONTROLS_PASSED','scope':'All naturally labeled posets through six vertices; controls only. Arbitrary finite size is proved by the injections in the companion theorem.','stats':stats,'certificate_digest':digest.hexdigest()}
print(json.dumps(result,indent=2))
