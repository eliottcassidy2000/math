#!/usr/bin/env python3
"""Independent finite cross-checks of orientation moments and vertex-insertion recursion."""
from itertools import combinations,permutations
from functools import lru_cache
from math import comb
from orientation import moments, graph, paths

def require(p,msg):
 if not p: raise ValueError(msg)

def paths_permutation(n,E):
 E=set(E)
 return [p for p in permutations(range(1,n+1)) if p[0]<p[-1] and all(tuple(sorted(x)) in E for x in zip(p,p[1:]))]

def verify_graph(n,E,brute_orientation=True):
 P=paths_permutation(n,E)
 events,B,C=moments(E,P)
 # Independent direct orientation enumeration.
 if brute_orientation:
  idx={e:i for i,e in enumerate(E)}; hist=[0]*(2*len(P)+1)
  for code in range(1<<len(E)):
   count=0
   for p in P:
    for path in (p,p[::-1]):
     if all(bool(code>>idx[tuple(sorted((a,b)))]&1)==(a<b) for a,b in zip(path,path[1:])):count+=1
   hist[count]+=1
  require(hist==C,'orientation mismatch n=%s E=%s'%(n,E))
 # F(S,x): oriented spanning paths of S ending in x, via endpoint recurrence.
 adj={x:{y for e in E if x in e for y in e if y!=x} for x in range(1,n+1)}
 @lru_cache(None)
 def f(S,x):
  if not (S>>(x-1))&1:return 0
  T=S^(1<<(x-1))
  if not T:return 1
  return sum(f(T,y) for y in adj[x] if (T>>(y-1))&1)
 old=(1<<(n-1))-1; neigh=adj[n]; v1=2*sum(f(old,x) for x in neigh)
 v2=0
 # ordered nonempty partitions of the old vertices, no double-count correction
 S=(old-1)&old
 while S:
  T=old^S
  if T:
   v2+=sum(f(S,x)*f(T,y) for x in neigh for y in neigh if (S>>(x-1))&1 and (T>>(y-1))&1)
  S=(S-1)&old
 require(v1+v2==2*len(P),'insertion recurrence mismatch n=%s E=%s'%(n,E))
 return len(P)

if __name__=='__main__':
 count=0
 for n in range(2,6):
  potential=list(combinations(range(1,n+1),2))
  for bit in range(1<<len(potential)):
   E=[e for i,e in enumerate(potential) if (bit>>i)&1]
   verify_graph(n,E,brute_orientation=(n<=4 or bit==(1<<len(potential))-1)); count+=1
 print('Verified insertion for',count,'labelled graphs n=2..5, orientation histogram for all n<=4 and K5.')
 # Edge cases/path enumeration independent in target cases through direct DFS from each oriented start.
 for n in (15,16,17,18,23,24,25,26,27):
  E,P=paths(n); adj={v:set() for v in range(1,n+1)}
  for a,b in E:adj[a].add(b);adj[b].add(a)
  @lru_cache(None)
  def f(mask,x):
   if mask==1<<(x-1):return 1
   old=mask^(1<<(x-1))
   return sum(f(old,y) for y in adj[x] if old>>(y-1)&1)
  full=(1<<n)-1
  actual=sum(f(full,x) for x in adj)
  require(actual==2*len(P),'target independent DP mismatch at n=%d'%n)
  print('n=%d directed count=%d, independent endpoint DP states=%d'%(n,actual,f.cache_info().currsize))
 print('Independent checks passed.')
