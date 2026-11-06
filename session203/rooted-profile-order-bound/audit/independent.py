import itertools as it, random, json
from pathlib import Path
from functools import lru_cache
P=Path(__file__).parent
checks=0;positive=0;strong=0

def need(x,msg):
 if not x:raise RuntimeError(msg)

def profile(a,r):
 n=len(a);h=s=L=0
 # A cover is uniquely a full permutation split immediately after root.
 for p in it.permutations(range(n)):
  valid=[bool(a[x]>>y&1) for x,y in zip(p,p[1:])]
  if all(valid):h+=1;s+=p[0]==r
  k=p.index(r)
  if all(valid[:k]) and all(valid[k+1:]):L+=1
 return h,s,L,h+L-s

def comps(a):
 n=len(a);reach=[]
 for i in range(n):
  seen={i};stack=[i]
  while stack:
   u=stack.pop()
   for v in range(n):
    if a[u]>>v&1 and v not in seen:seen.add(v);stack.append(v)
  reach.append(seen)
 left=set(range(n));cs=[]
 while left:
  v=min(left);c={u for u in left if u in reach[v] and v in reach[u]};cs.append(c);left-=c
 return sorted(cs,key=lambda c:-sum(a[next(iter(c))]>>v&1 for d in cs if d!=c for v in [next(iter(d))]))

def induce(a,vs):return tuple(sum(((a[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs)

@lru_cache(None)
def allprofiles(a):return tuple(profile(a,r) for r in range(len(a)))
def check(a):
 global checks,positive,strong
 cs=comps(a)
 for r,p in enumerate(allprofiles(a)):
  checks+=1;h,s,L,B=p
  if not s:continue
  positive+=1;need(r in cs[0],('initial',a,r))
  keep=sorted(set.union(cs[0],*(c for c in cs[1:] if len(c)>1)))
  pp=allprofiles(induce(a,keep))[keep.index(r)]
  need(p==pp,('preservation',a,r,p,pp))
  need(len(keep)<=s+B.bit_length(),('cap',a,r,p))
  c0=sorted(cs[0]);p0=allprofiles(induce(a,c0))[c0.index(r)];q=1
  for c in cs[1:]:q*=allprofiles(induce(a,sorted(c)))[0][0]
  need(p==tuple(q*x for x in p0),('product',a,r))
  if len(cs)==1:
   strong+=1;need(L*(1<<s)>=1<<(len(a)-1),('strong',a,r))

def graph(n,code):
 a=[0]*n
 for k,(i,j) in enumerate(it.combinations(range(n),2)):
  u,v=(i,j) if code>>k&1 else (j,i);a[u]|=1<<v
 return tuple(a)
counts={}
for n in range(1,6):
 counts[n]=1<<(n*(n-1)//2)
 for k in range(counts[n]):check(graph(n,k))
rng=random.Random(203137)
for n in [6,7]:
 codes=rng.sample(range(1<<(n*(n-1)//2)),128);counts[n]=len(codes)
 for code in codes:check(graph(n,code))
print(json.dumps({'status':'PASS','method':'Full permutation enumeration; split after root counts literal covers; DFS mutual reachability SCCs; all roots per graph','graphs_by_order':counts,'rooted_graph_checks':checks,'positive_roots':positive,'strong_roots':strong,'seed':203137,'scope':'All labeled tournaments orders 1–5, 128 distinct seeded samples each of orders 6 and 7. Not an all-order or occupancy proof.'},indent=2))
