from itertools import combinations
from math import isqrt

def require(test, message):
 if not test: raise ValueError(message)
def square_edge(x,y):
 s=x+y
 return x!=y and isqrt(s)**2 == s
def edges(vertices):
 return {tuple(pair) for pair in combinations(sorted(vertices),2) if square_edge(*pair)}
def components(vertices, edge_set):
 todo=set(vertices); out=[]
 while todo:
  v=min(todo); todo.remove(v); reached={v}; queue=[v]
  for x in queue:
   for y in sorted(todo.copy()):
    if tuple(sorted((x,y))) in edge_set:
     todo.remove(y); reached.add(y); queue.append(y)
  out.append(sorted(reached))
 return sorted(out,key=lambda c:(len(c),c))

V=set(range(1,19))-{13,14,15}
source=edges(V)
comps=components(V,source)
S={len(c):c for c in comps}
A=[2,3,5,7,11,13]                  # Fi22 support
B=[17,19,23,31,47]                # Baby Monster support minus Fi22
C=[29,41,59,71]                   # Monster support minus Baby Monster
T={6:A,4:C,5:B}
require(comps==[[4,5,11,12],[2,7,9,16,18],[1,3,6,8,10,17]],'components changed')
f={x:y for size in S for x,y in zip(S[size],T[size])}
require(len(f)==len(set(f.values()))==15,'mapping is not bijective')
pairs=set(combinations(sorted(V),2))
back={pair for pair in pairs if square_edge(f[pair[0]],f[pair[1]])}
target=edges(set(f.values()))
require(components(set(f.values()),target)==[sorted(f.values())], 'prime target not connected')
preserved=source&back; lost=source-back; created=back-source
require((len(source),len(target),len(preserved),len(lost),len(created))==(12,14,5,7,9),'edge census changed')
require(not square_edge(f[1],f[3]) and square_edge(1,3),'hostile changed')
within={k:edges(T[k]) for k in T}
require({k:len(v) for k,v in within.items()}=={6:3,4:2,5:2},'within-shell census changed')
cross=target-set().union(*within.values())
require(cross=={(2,23),(2,47),(5,31),(5,59),(7,29),(13,23),(23,41)},'cross-shell census changed')
for size in (6,4,5):
 se=edges(S[size]); require(len(se)==size-1,'source component not a tree')
 require(len(se)!=len(within[size]),'unexpected shell isomorphism')

def mapped(E):return [(x,y,f[x],f[y]) for x,y in sorted(E)]
print('source components:',comps)
print('rank map:',sorted(f.items()))
print('source edges:',sorted(source))
print('prime target edges:', sorted(target))
print('preserved:',mapped(preserved))
print('lost:',mapped(lost))
print('created:',mapped(created))
print('within-shell edge totals:',{k:len(v) for k,v in within.items()})
print('cross-shell target edges:',sorted(cross))
print('source/target components:',comps,components(set(f.values()),target))
print('PASS: full pair census; no edge-count-preserving isomorphism exists.')
