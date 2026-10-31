from math import isqrt
from collections import deque, Counter

def graph(n,odd=False):
 a={x:set() for x in range(1,n+1)}
 for x in a:
  for k in range(isqrt(x)+1,isqrt(x+n)+1):
   if odd and k%2==0:continue
   y=k*k-x
   if x!=y:a[x].add(y)
 return a

def comps(a):
 seen=set();out=[]
 for v in a:
  if v in seen:continue
  c={v};seen.add(v);q=[v]
  while q:
   for w in a[q.pop()]:
    if w not in seen:seen.add(w);c.add(w);q.append(w)
  out.append(sorted(c))
 return out

def girth(a,odd=False):
 ans=None
 for s in a:
  d={s:0};parent={s:None};q=deque([s])
  while q:
   v=q.popleft()
   for w in a[v]:
    if w not in d:d[w]=d[v]+1;parent[w]=v;q.append(w)
    elif parent[v]!=w and parent[w]!=v:
     z=d[v]+d[w]+1
     if (not odd or z%2) and (ans is None or z<ans):ans=z
 return ans

def cuts(a):
 t=0;tin={};low={};aps=set();bridges=[]
 def dfs(v,par):
  nonlocal t
  tin[v]=low[v]=t;t+=1;children=0
  for w in sorted(a[v]):
   if w==par:continue
   if w in tin:low[v]=min(low[v],tin[w])
   else:
    dfs(w,v);low[v]=min(low[v],low[w]);children+=1
    if low[w]>tin[v]:bridges.append(sorted((v,w)))
    if par is not None and low[w]>=tin[v]:aps.add(v)
  if par is None and children>1:aps.add(v)
 for v in a:
  if v not in tin:dfs(v,None)
 return sorted(aps),sorted(bridges)

def invariants(n):
 a=graph(n);d=Counter(map(len,a.values()));c=comps(a);aps,bridges=cuts(a)
 return dict(n=n,edges=sum(k*v for k,v in d.items())//2,components=len(c),component_sizes=sorted(map(len,c),reverse=True),min_degree=min(d),max_degree=max(d),degree_histogram=dict(sorted(d.items())),leaves=d[1],articulation_points=aps,bridges=bridges,girth=girth(a),odd_girth=girth(a,True),odd_root_components=len(comps(graph(n,True))),arrival_roots=list(range(isqrt(n)+1,isqrt(2*n-1)+1)))
