from functools import lru_cache
from itertools import combinations
from math import isqrt
from collections import deque

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def graph(n):
 g=[0]*n
 for a,b in combinations(range(n),2):
  q=isqrt(a+b+2)
  if q*q==a+b+2:g[a]|=1<<b;g[b]|=1<<a
 return g
def hamiltons(g):
 n=len(g);full=(1<<n)-1
 @lru_cache(None)
 def recurse(last,seen):
  if seen==full:return 1
  bits=g[last]&~seen;count=0
  while bits:
   bit=bits&-bits;bits-=bit;count+=recurse(bit.bit_length()-1,seen|bit)
  return count
 return sum(recurse(v,1<<v) for v in range(n))
def centers(g):
 n=len(g);e=[]
 for s in range(n):
  q=deque([s]);d={s:0}
  while q:
   v=q.popleft();bits=g[v]
   while bits:
    bit=bits&-bits;bits-=bit;u=bit.bit_length()-1
    if u not in d:d[u]=d[v]+1;q.append(u)
  require(len(d)==n,'disconnected')
  e.append(max(d.values()))
 m=min(e);return m,tuple(i+1 for i,x in enumerate(e) if x==m)
def components(g,subset):
 left=set(subset);result=[]
 while left:
  s=left.pop();p=[s];c={s}
  for u in p:
   for v in list(left):
    if g[u] & (1<<v):left.remove(v);c.add(v);p.append(v)
  result.append(c)
 return result
for n in range(1,19):
 g=graph(n); count=hamiltons(g);expected=1 if n==1 else 2 if n in {15,16,17} else 0
 require(count==expected,'Hamilton counts '+str(n));print(n,count)
for n in range(14,23):
 r,c=centers(graph(n));print('center',n,r,c)
 expected={14:(4,5),15:(4,5),16:(5,),17:(4,5),18:(4,5),19:(4,5),20:(4,),21:(4,5),22:(14,22)}
 require(c==expected[n],'center '+str(n))
g=graph(12); cs=components(g,set(range(12)));require(sorted(map(len,cs))==[3,4,5],'G12 component sizes')
# Component-induced eccentricities, by BFS through each component
out=[]
for comp in cs:
 ecc={}
 for s in comp:
  q=deque([s]);d={s:0}
  while q:
   v=q.popleft()
   for u in comp:
    if u not in d and g[v]&(1<<u):d[u]=d[v]+1;q.append(u)
  ecc[s]=max(d.values())
 z=min(ecc.values());out += [v+1 for v in comp if ecc[v]==z]
require(set(out)=={3,4,5,7},'G12 centers')
center_edges={(u+1,v+1) for u,v in combinations([x-1 for x in sorted(out)],2) if g[u]&(1<<v)}
require(center_edges=={(4,5)},'G12 center-induced edges')
for n,leaf in [(14,{8,9,10}),(15,{8,9}),(16,{8,16}),(17,{16,17}),(18,{16,17,18})]:
 g=graph(n); actual={i+1 for i,x in enumerate(g) if x.bit_count()==1}; edgecount=sum(x.bit_count() for x in g)//2
 require(actual==leaf,'leaf check '+str(n));require(len(components(g,set(range(n))))==1,'connectivity');require(edgecount==(n-1 if n==14 else n),'cycle rank '+str(n))
 print('graph',n,'leaves',actual,'edges',edgecount)
print('verified Hamilton counts, G12 center-induced graph, G14-22 center sequence, connectivity, cycle ranks and leaves')
p=(8,1,15,10,6,3,13,12,4,5,11,14,2,7,9)
for n,path in [(15,p),(16,p+(16,)),(17,(17,)+p+(16,))]:
 g=graph(n);require(set(path)==set(range(1,n+1)) and len(path)==n,'path covers G'+str(n))
 es={frozenset((a,b)) for a,b in zip(path,path[1:])}
 require(all(g[a-1]&(1<<(b-1)) for a,b in zip(path,path[1:])),'path edge G'+str(n))
 require(frozenset((4,5)) in es and frozenset((1,3)) not in es,'trichotomy')
 require(g[0]&(1<<2) and not g[1]&(1<<3),'ambient pair labels')
 # Since exactly two ordered paths were independently counted, this witness and its reversal exhaust all.
 print('G%d unique up to reversal: 4–5 forced, 1–3 excluded, 2–4 absent'%n)
