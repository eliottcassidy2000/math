"""Independent edge-mask union-find native census; paths generate all legal sources."""
from itertools import combinations, permutations
from collections import Counter
import json

def need(x,msg):
 if not x:raise RuntimeError(msg)
out=[]
for n in range(3,7):
 pairs=list(combinations(range(n),2));index={e:i for i,e in enumerate(pairs)};N=len(pairs);parent=list(range(1<<N));degree=[0]*(1<<N);edges=0
 def root(x):
  while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
  return x
 for P in permutations(range(n)):
  fixed=0;keep=0
  for u,v in zip(P,P[1:]):
   k=index[tuple(sorted((u,v)))];keep|=1<<k
   if u>v:fixed|=1<<k
  free=[k for k in range(N) if not keep>>k&1];delta=((1<<N)-1)^keep
  # Force one free bit0 to enumerate each undirected native edge exactly once.
  states=[fixed]
  for k in free[:-1]:states += [x|(1<<k) for x in states]
  for x in states:
   y=x^delta;rx,ry=root(x),root(y)
   if rx!=ry:parent[rx]=ry
   degree[x]+=1;degree[y]+=1;edges+=1
 sizes=Counter(root(x) for x in range(1<<N));hist=Counter(sizes.values())
 if n==3:need(sorted(sizes.values())==[4,4] and Counter(degree)=={1:6,3:2},'two K1,3')
 if n==4:need(list(sizes.values())==[64],'native n4 connected')
 if n==5:need(sorted(sizes.values())==[512,512],'native n5 two sheets')
 out.append({'n':n,'vertices':1<<N,'edges':edges,'component_size_histogram':dict(sorted(hist.items()))})
print(json.dumps({'status':'FINITE-EXACT','method':'all ordered paths, all free nonpath orientations; direct union-find','census':out},indent=2,sort_keys=True))
