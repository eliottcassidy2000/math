from itertools import product
from functools import lru_cache
from collections import defaultdict
import json
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
m=7;words=list(product(range(3),repeat=2));vertices=list(product(range(m),words));C={1,3,5}
def depth(i,j):return 2 if i in C and j in C else 1
def outer(i,j):return (j-i)%7 in (1,2,3)
def inner(u,v):
 for a,b in zip(u,v):
  if a!=b:return (b-a)%3==1
 return False
def arc(x,y):
 i,u=x;j,v=y
 if i==j:return inner(u,v)
 equal=u[:depth(i,j)]==v[:depth(i,j)]
 return equal if outer(i,j) else not equal
need(all(sum(arc(x,y) for y in vertices)==31 for x in vertices),'actual regular63 host')
path=[(i,w) for i in reversed(range(7)) for w in words]
need(len(set(path))==63 and all(arc(x,y) for x,y in zip(path,path[1:])),'actual reverse-cell lex path')
def component(i):return 'C135' if i in C else str(i)
def slot(v):return (component(v[0]),v[1][0])
slots=sorted({slot(v) for v in path});idx={s:i for i,s in enumerate(slots)};G=[set() for _ in slots];loops=set()
for x,y in zip(path,path[1:]):
 a,b=idx[slot(x)],idx[slot(y)]
 if a==b:loops.add(a)
 else:G[a].add(b);G[b].add(a)
need(len(slots)==15 and len(loops)==15,'slots and available singleton flags')
@lru_cache(None)
def matching(mask):
 if not mask:return ()
 a=(mask&-mask).bit_length()-1;rest=mask^(1<<a);best=matching(rest)
 for b in G[a]:
  if rest>>b&1:
   candidate=((a,b),)+matching(rest^(1<<b))
   if len(candidate)>len(best):best=candidate
 return best
match=matching((1<<len(slots))-1)
need(len(match)==6,'exact maximum matching')
colors={};side=[set(),set()]
for a in range(len(slots)):
 if a in colors:continue
 colors[a]=0;todo=[a]
 while todo:
  u=todo.pop();side[colors[u]].add(u)
  for v in G[u]:
   if v in colors:need(colors[v]!=colors[u],'slot bipartition')
   else:colors[v]=1-colors[u];todo.append(v)
need(sorted(map(len,side))==[6,9],'imbalanced bipartition')
# Turn matching plus singleton loops into actual indexed path-arc observations.
needed=[frozenset((a,b)) for a,b in match]
used={a for pair in match for a in pair}
needed += [frozenset((a,)) for a in range(len(slots)) if a not in used]
flags=[]
for requirement in needed:
 flags.append(next(i for i,(x,y) in enumerate(zip(path,path[1:])) if frozenset((idx[slot(x)],idx[slot(y)]))==requirement))
need(len(flags)==9 and {idx[slot(v)] for i in flags for v in path[i:i+2]}==set(range(15)),'actual9flags cover all slots')
print(json.dumps({'status':'PASS','checks':checks,'vertices':63,'outdegree':31,'full_resolution_component':[1,3,5],'declared_chart_order':3**16,'owner_base':15,'slot_bipartition_sizes':sorted(map(len,side)),'maximum_matching':[[slots[a],slots[b]] for a,b in match],'minimum_indexed_arc_flags':9,'actual_flag_positions_zero_based':flags,'scope':'specified named-cell group orbit of actual path; not all paths or full uncolored group'},indent=2))
