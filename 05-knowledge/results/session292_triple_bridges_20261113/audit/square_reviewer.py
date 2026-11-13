#!/usr/bin/env python3
"""Independent reviewer: no imports from producer scripts."""
from itertools import permutations,combinations
from math import factorial
import json
from pathlib import Path

def ck(x,s):
 if not x:raise RuntimeError(s)
def plus(a,b):
 z=a.copy()
 for k,v in b.items():z[k]=z.get(k,0)+v
 return {k:v for k,v in z.items() if v}
def times(a,b):
 z={}
 for u,x in a.items():
  for v,y in b.items():
   if u.isdisjoint(v):z[u|v]=z.get(u|v,0)+x*y
 return {k:v for k,v in z.items() if v}
def pw(a,k):
 z={frozenset():1}
 for _ in range(k):z=times(z,a)
 return z
models=orders=0
for n in range(1,5):
 pairs=list(combinations(range(n),2))
 for code in range(2**len(pairs)):
  arcs=set((j,i) if code>>b&1 else (i,j) for b,(i,j) in enumerate(pairs))
  direct=sum(all((p[i],p[i+1]) in arcs for i in range(n-1)) for p in permutations(range(n)))
  for perm in permutations(range(n)):
   q={(i,j):{frozenset([j]):1} for i,j in arcs}
   for j in range(n):q[(n,j)]={frozenset([j]):1};q[(j,n+1)]={frozenset():1}
   live=set(range(n+2))
   for k in perm:
    live.remove(k); old=q.copy()
    for i in live:
     ck(not times(old.get((i,k),{}),old.get((k,k),{})), 'destination annihilation')
     for j in live:q[(i,j)]=plus(old.get((i,j),{}),times(old.get((i,k),{}),old.get((k,j),{})))
   ck(q.get((n,n+1),{}).get(frozenset(range(n)),0)==direct,'Hamiltonian coefficient')
   orders+=1
  models+=1
sharp=[]
for h in range(1,8):
 D={frozenset(s):(-1)**len(s) for k in range(h+1) for s in combinations(range(h),k)}
 r={s:-c for s,c in D.items() if s}
 ck(pw(r,h)[frozenset(range(h))]==factorial(h),'triangle coefficient')
 inverse={frozenset():1}
 for j in range(1,h+1):inverse=plus(inverse,pw(r,j))
 ck(all(c==1 for c in inverse.values()) and len(inverse)==2**h,'full inverse')
 sharp.append(h)
B=[1,2,3,4,6,8,9,12];W=[-216,1188,-2464,2079,-924,594,-264,7]
ck([sum(w*b**j for w,b in zip(W,B)) for j in range(8)]==[0]*7+[3991680],'eight-mode witness')
R={'status':'PASS','all_tournaments_through_4':models,'all_pivot_orders':orders,'triangle_blocks':sharp,'eight_mode_witness':'PASS','scope':'Independent literal Hamiltonian permutation counts versus occupancy elimination, explicit determinant-denominator algebra, exact delayed-profile witness. Does not establish physical replacement or coverage.'}
Path(__file__).with_name('square_reviewer.json').write_text(json.dumps(R,indent=2)+'\n');print(json.dumps(R,indent=2))
