"""Recovered C2 witness; independent construction and meet-in-middle H counts.
Source eea15f5d968c8e493a809862bbfec3989d09c436 / moving_selector_atoms.
No producer imports or receipt-dependent calculations.
"""
from itertools import permutations,combinations
from functools import lru_cache
import json

def need(b,s):
 if not b:raise RuntimeError(s)
def relation(rows):return {(u,v) for u,row in enumerate(rows) for v in range(len(rows)) if row>>v&1}
def strong(E,n):
 for a in range(n):
  reach={a}
  while True:
   new=reach|{v for u,v in E if u in reach}
   if new==reach:break
   reach=new
  if len(reach)!=n:return False
 return True
def substitute(children):
 offsets=[0,len(children[0]),len(children[0])+len(children[1])];E=set()
 for i,rows in enumerate(children):
  E|={(u+offsets[i],v+offsets[i]) for u,v in relation(rows)}
  j=(i+1)%3
  E|={(u+offsets[i],v+offsets[j]) for u in range(len(rows)) for v in range(len(children[j]))}
 return E
def hp_middle(E,n=13):
 full=(1<<n)-1
 @lru_cache(None)
 def end(mask,v):
  if mask==1<<v:return 1
  rest=mask^(1<<v)
  return sum(end(rest,u) for u in range(n) if rest>>u&1 and (u,v) in E)
 @lru_cache(None)
 def start(mask,v):
  if mask==1<<v:return 1
  rest=mask^(1<<v)
  return sum(start(rest,w) for w in range(n) if rest>>w&1 and (v,w) in E)
 total=0
 for C in combinations(range(n),n//2):
  left=sum(1<<v for v in C);right=full^left
  for u in C:
   h=end(left,u)
   if h:
    total+=h*sum(start(right,v) for v in range(n) if right>>v&1 and (u,v) in E)
 return total
A=(10,12,9,16,7);D=(2,4,9,3);T=(14,12,8,0)
EA=relation(A);need(strong(EA,5),'A strong');need(strong(relation(D),4),'D strong')
G=[g for g in permutations(range(5)) if {(g[u],g[v]) for u,v in EA}==EA]
need(len(G)==3,'Aut A exactlyC3')
U=(3,0,1,2,4);V=(3,0,2,1,4)
for E in (EA,{(v,u) for u,v in EA}):
 maps=[g for g in permutations(range(5)) if {(g[u],g[v]) for u,v in E}==E and g[3]==3 and g[4]==4]
 need(all(any(g[u]==v for g in maps) for u,v in zip(U,V)),'all five rooted cursor probes')
 need(not any(all(g[u]==v for u,v in zip(U,V)) for g in maps),'no common transporter')
results=[]
for kind,child,expected in [('transitive',T,(2483993,2491877)),('all_strong',D,(2899675,2906091))]:
 E=substitute((A,child,child));need(strong(E,13),'source strong')
 for middle in permutations(range(3)):
  word=(3,)+middle+(4,);P=[]
  for i,v in enumerate(word):
   P.append(v)
   if i<4:P.extend((5+i,9+i))
  need(len(set(P))==13,'Hamilton vertex set')
  seams=set(zip(P,P[1:]));need(seams<=E,'legal native Hamilton guard')
  B={(u,v) if (u,v) in seams else (v,u) for u,v in E}
  need(strong(B,13),'partner strong')
  parity=sum(middle[i]>middle[j] for i in range(3) for j in range(i+1,3))%2
  h=hp_middle(B);need(h==expected[parity],'H mismatch')
  results.append({'kind':kind,'middle':middle,'parity':parity,'H':h})
print(json.dumps({'source_pin':'eea15f5d968c8e493a809862bbfec3989d09c436','cursor_sheets':2,'cursor_positions_each':5,'child_aut':3,'method':'6/7 meet-in-middle subset Hamilton paths','cases':results},sort_keys=True))
