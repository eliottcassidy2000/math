"""Finite reduced-core engine graded by excess plus number of leaves.
Uses candidate extension under independent audit of source local run decomposition.
No assertion of a uniform remainder for the infinite expansion.
"""
import itertools, functools
import numpy as np
from loop_engine import automorphisms
from loop_engine_fast import CoreB,basis

def decorated_cores(weight):
 found=[]
 for V in range(2 if weight==1 else 1,2*weight+1):
  pairs=[(a,b) for a in range(V) for b in range(a,V)]
  for ell in range(0,V+1):
   e=weight-ell;E=V+e
   if e< -1 or E<1:continue
   seen=set()
   for ids in itertools.combinations_with_replacement(range(len(pairs)),E):
    edges=[pairs[i] for i in ids];deg=[0]*V
    for a,b in edges:deg[a]+=1;deg[b]+=1
    if deg.count(1)!=ell or any(d not in (1,) and d<3 for d in deg):continue
    # No isolated vertices or degree two vertices: the latter are suppressed.
    reached={0}
    while True:
     more=reached|{b for a,b in edges if a in reached}|{a for a,b in edges if b in reached}
     if more==reached:break
     reached=more
    if len(reached)!=V:continue
    canon=min(tuple(sorted((min(p[a],p[b]),max(p[a],p[b])) for a,b in edges)) for p in itertools.permutations(range(V)))
    if canon in seen:continue
    seen.add(canon)
    found.append({'V':V,'edges':list(canon),'excess':e,'leaves':ell,'weight':weight,'automorphisms':automorphisms(V,list(canon))})
 # weight1 includes figureeight with one internal degree4 vertex.
 if weight==1:
  edges=[(0,0),(0,0)]
  found.insert(0,{'V':1,'edges':edges,'excess':1,'leaves':0,'weight':1,'automorphisms':8})
 return found

@functools.lru_cache(None)
def compiled(weight):
 return [(d,CoreB(d['V'],d['edges'])) for d in decorated_cores(weight)]

@functools.lru_cache(None)
def compiled_reduced(weight):
 import copy
 change=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-1,0]],dtype=np.int64)
 out=[]
 for d,old in compiled(weight):
  c=copy.copy(old);W=c.Wb
  for _ in d['edges']:W=np.tensordot(W,change,axes=([0],[0]))
  c.Wb=W;c.pathb=None;out.append((d,c))
 return out

def evaluate(M,weight,z=1.,detail=False,reduced=True):
 B=basis(M,z)
 cs=compiled_reduced(weight) if reduced else compiled(weight)
 if reduced:B=B[:4]
 vals=[(d,c.value_b(B)) for d,c in cs]
 return vals if detail else sum(x for _,x in vals)

if __name__=='__main__':
 import json,collections,time
 out={}
 for w in (1,2):
  t=time.time();cs=compiled(w)
  out[w]={'count':len(cs),'by_type':dict(collections.Counter(str((d['excess'],d['leaves'])) for d,c in cs)),'cores':[d for d,c in cs],'seconds':time.time()-t}
 print(json.dumps(out,indent=2))
