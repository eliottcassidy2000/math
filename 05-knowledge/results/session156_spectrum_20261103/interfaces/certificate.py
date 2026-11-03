def require(x):
 if not x: raise AssertionError("interface check failed")
#!/usr/bin/env python3
"""Finite controls for the three exact clone families; Python standard library."""
from itertools import permutations, combinations, product
from collections import Counter
from functools import lru_cache
from math import isqrt
T=[62,60,120,16,96,8,43]

def edge(a,u,v): return (a[u]>>v)&1

def tournament(a):
 return all(not edge(a,u,u) for u in range(len(a))) and all(edge(a,u,v)+edge(a,v,u)==1 for u,v in combinations(range(len(a)),2))

def transitive(a,delete):
 S=[i for i in range(len(a)) if i not in delete]
 return not any(edge(a,x,y) and edge(a,y,z) and edge(a,z,x) for x,y,z in permutations(S,3))

def blowup(a,r,k):
 vs=[v for v in range(len(a)) if v!=r]+[r]*k
 b=[0]*len(vs)
 for i,j in combinations(range(len(vs)),2):
  if vs[i]==vs[j] or edge(a,vs[i],vs[j]): b[i]|=1<<j
  else: b[j]|=1<<i
 return b

def hp(a):
 n=len(a)
 @lru_cache(None)
 def d(mask,v):
  if mask==1<<v: return 1
  return sum(d(mask^(1<<v),u) for u in range(n) if mask>>u&1 and edge(a,u,v))
 return sum(d((1<<n)-1,v) for v in range(n))

def tally(a,r):
 U=[v for v in range(len(a)) if v!=r]; out=Counter()
 for p in permutations(U):
  valid=[edge(a,r,p[0])]+[edge(a,p[i],r) and edge(a,r,p[i+1]) for i in range(len(U)-1)]+[edge(a,p[-1],r)]
  defects=[i+1 for i in range(len(U)-1) if not edge(a,p[i],p[i+1])]
  if all(valid[s] for s in defects): out[(sum(valid),len(defects))]+=1
 return out

def apex_counts(bound):
 # Exact by S155's core reduction and H>=2^(m-2)+1.
 out={1}; m=2
 while 2**(m-2)+1<=bound:
  for mid in product((0,1),repeat=m-2):
   h,a,p=1,0,0
   for b in (1,*mid,0):
    if b: h,a,p=h,2*a+1,1
    else: h,a,p=h+a+p,2*a,0
   if h<=bound: out.add(h)
  m+=1
 return out

def prime(x): return x>=2 and all(x%d for d in range(2,isqrt(x)+1))

require(tournament(T) and hp(T)==35)
require(all(any(edge(T,u,v) for u in range(7) if mask>>u&1 for v in range(7) if not mask>>v&1) for mask in range(1,127)))
require(all(not transitive(T,{v}) for v in range(7)))
config={0:({(1,0):3,(2,0):14,(2,1):4},lambda k:18*2**k-1,{5,6}),3:({(1,1):5,(2,0):11,(2,1):8},lambda k:19*2**k-3,{5,6}),5:({(1,0):7,(1,1):4,(2,0):12},lambda k:12*2**k+11,{2,4})}
for r,(expected,f,deletion) in config.items():
 require(transitive(T,deletion) and r not in deletion)
 require(tally(T,r)==expected)
 for k in range(1,6):
  a=blowup(T,r,k); require(tournament(a) and hp(a)==f(k))
  # Clones are listed last; map deleted old vertices to the new labels.
  keep=[v for v in range(7) if v!=r]
  require(transitive(a,{keep.index(v) for v in deletion}))
 print('root',r,'tally',sorted(expected.items()),'H(1..5)',[f(k) for k in range(1,6)])

spec=apex_counts(4861)
for r,k in [(0,2),(3,3),(3,4),(3,6),(3,7),(3,8),(5,2)]:
 h=config[r][1](k); require(h not in spec)
 print('one-apex excluded:',r,k,h)
for p in [29,41,43,347,59,71,149,1213,4861]: require(prime(p))
for r,k,p in [(0,4,41),(3,4,43),(3,7,347),(5,4,29)]: require(config[r][1](k)==7*p)
print('all controls passed')
