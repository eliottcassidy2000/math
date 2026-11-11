#!/usr/bin/env python3
"""Independent literal control on all four-vertex tournaments and root cuts."""
from itertools import combinations,permutations
from pathlib import Path
import json

def require(x,msg):
 if not x:raise RuntimeError(msg)
def paths(T,S):
 return [p for p in permutations(S) if all(p[i+1] in T[p[i]] for i in range(len(p)-1))]
def transitive(T,S):
 return not any(b in T[a] and c in T[b] and a in T[c] for a,b,c in permutations(S,3))
def review(T,z):
 V=set(range(len(T)))-{z};A=T[z];B=V-A
 adj={v:set() for v in V}
 for a in A:
  for b in B:
   if b in T[a]:adj[a].add(b);adj[b].add(a)
 H=len(paths(T,range(len(T))));D=0;defect=0
 for bits in range(1<<len(V)):
  L={v for i,v in enumerate(sorted(V)) if bits>>i&1};R=V-L;S=(A-L)|(B&L)
  eps=int(all(v in S or adj[v]&S for v in V));D+=eps
  left=sum(not p or p[-1] in B for p in paths(T,L))
  right=sum(not p or p[0] in A for p in paths(T,R))
  require(left*right>=eps,('fibre',T,z,L));defect+=left*right-eps
 require(H-D==defect,('sum',T,z))
 require((H==D)==transitive(T,V),('equality',T,z))
 require(H==D or H>=D+2,('strict gap',T,z))
 return H,D
pairs=list(combinations(range(4),2));checks=0
for mask in range(1<<len(pairs)):
 T=[set() for _ in range(4)]
 for k,(a,b) in enumerate(pairs):
  if mask>>k&1:T[a].add(b)
  else:T[b].add(a)
 for z in range(4):review(T,z);checks+=1
# Fixed 2K2 graph on shores {0,1},{2,3}; apex4.
h=[]
for x,y in [(0,0),(0,1),(1,0),(1,1)]:
 T=[set() for _ in range(5)];T[4]|={0,1};T[2].add(4);T[3].add(4)
 T[x].add(1-x);T[2+y].add(3-y)
 for a in (0,1):
  for b in (2,3):
   if b==a+2:T[a].add(b)
   else:T[b].add(a)
 H,D=review(T,4);require(D==9 and H>D,'2K2');h.append(H)
out={'status':'PASS','four_vertex_root_checks':checks,'four_vertex_cut_checks':checks*8,'twoK2_H_values':h,'twoK2_domination':9}
Path(__file__).with_name('bipartite_defect_tiny_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
