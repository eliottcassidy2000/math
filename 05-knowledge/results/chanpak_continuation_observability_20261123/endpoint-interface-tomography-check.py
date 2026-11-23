from itertools import permutations, combinations
from functools import lru_cache
from random import Random

def H(T):
 n=len(T);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for A in range(1,1<<n):
  for i in range(n):
   if A>>i&1 and (m:=dp[A][i]):
    for j in range(n):
     if not A>>j&1 and T[i]>>j&1:dp[A|1<<j][j]+=m
 return sum(dp[-1])

def orient(T,u,v):
 T[u]|=1<<v;T[v]&=~(1<<u)

def chart(n,rng):
 T=[0]*n
 for u,v in combinations(range(n),2):orient(T,*( (u,v) if rng.randrange(2) else (v,u)))
 return T

def scaffold(n,p,i,mark=None,extra=0):
 O=[i]+[u for u in range(n) if u not in (p,i)]
 if mark:
  a,b=mark
  if i not in mark:O=[i,a,b]+[u for u in O if u not in (*mark,i)]
  else:O=[i,next(u for u in mark if u!=i)]+[u for u in O if u not in mark]
 path=[];new=n
 for k,u in enumerate(O):
  path.append(u)
  if k+1<len(O) and (not mark or {u,O[k+1]}!=set(mark)):
   path.append(new);new+=1
 for _ in range(1+bool(mark)+extra):
  path.append(new);new+=1
 assert new==2*n-1+extra
 return path

def test(T,p,i,rng,mark):
 n=len(T);path=scaffold(n,p,i,mark)
 fixed=chart(2*n-1,rng)
 for u in range(n):fixed[u]=(fixed[u]&~((1<<n)-1))|T[u]
 actual=0;edges=list(zip(path,path[1:]))
 eps=1
 if mark:
  seam=next((u,v) for u,v in edges if {u,v}==set(mark))
  edges.remove(seam)
  eps=1 if T[seam[0]]>>seam[1]&1 else -1
 for mask in range(1<<len(edges)):
  cur=fixed[:]
  for j,(u,v) in enumerate(edges):orient(cur, *((u,v) if mask>>j&1 else (v,u)))
  x=H(cur)
  if mark:
   a,b=mark
   cur[a]^=1<<b;cur[b]^=1<<a
   x-=H(cur)
  actual+=x*(-1)**(len(edges)-mask.bit_count())
 assert actual==eps*(2*(fixed[p]>>i&1)+2*(fixed[path[-1]]>>p&1)-2)

def strong_cycle(n,path):
 forbidden={frozenset((u,v)) for u,v in zip(path,path[1:])}
 allowed=lambda u,v:frozenset((u,v)) not in forbidden
 gaps={u:[v for v in range(n,2*n) if allowed(u,v) and allowed(v,(u+1)%n)] for u in range(n)}
 mt={}
 def aug(u,seen):
  for v in gaps[u]:
   if v in seen:continue
   seen.add(v)
   if v not in mt or aug(mt[v],seen):mt[v]=u;return True
  return False
 for u in range(n):assert aug(u,set())
 match={u:v for v,u in mt.items()}
 l=[x for u in range(n) for x in (u,match[u])]
 for k in range(1,len(l),2):
  if allowed(l[k],2*n) and allowed(2*n,l[(k+1)%(2*n)]):
   l.insert(k+1,2*n);break
 assert len(l)==2*n+1
 assert all(allowed(u,v) and (u>=n or v>=n) for u,v in zip(l,l[1:]+l[:1]))
 return l

r=Random(157)
controls=0;strong_count=0
for n in range(2,5):
 for _ in range(3):
  T=chart(n,r)
  for p,i in permutations(range(n),2):
   test(T,p,i,r,None);controls+=1
   if n>=3:
    mark=tuple(v for v in range(n) if v!=p)[:2]
    test(T,p,i,r,mark);controls+=1
for n in range(8,24):
 for _ in range(20):
  p,i=r.sample(range(n),2)
  for mark in (None,tuple(r.sample([v for v in range(n) if v!=p],2))):
   l=strong_cycle(n,scaffold(n,p,i,mark,2))
   strong_count+=1
print(controls,"literal count/gap calibrations;",strong_count,"strong probe scaffold cycles")

from collections import Counter
def covers(T,colors):
 n=len(T);out=[Counter() for _ in range(n+1)]
 for pp in permutations(range(n)):
  for sig in range(1<<(n-1)):
   if any(not sig>>i&1 and not T[u]>>v&1
          for i,(u,v) in enumerate(zip(pp,pp[1:]))):continue
   cuts=[0]+[i+1 for i in range(n-1) if sig>>i&1]+[n]
   key=tuple((colors[pp[u]],colors[pp[v-1]]) for u,v in zip(cuts,cuts[1:]))
   out[len(key)][key]+=1
 return out
def contract(C,U,ph):
 m=len(U);value=0
 for pv in permutations(range(m)):
  for R in range(1<<(m+1)):
   if any(not R>>(i+1)&1 and not U[u]>>v&1 for i,(u,v) in enumerate(zip(pv,pv[1:]))):continue
   g=[i for i in range(m+1) if R>>i&1]
   if len(g)>=len(C):continue
   for be,v in C[len(g)].items():
    if all((i==0 or ph[pv[i-1]][a]) and (i==m or not ph[pv[i]][b]) for i,(a,b) in zip(g,be)):
     value+=v
 return value
for _ in range(60):
 n=r.randrange(2,6);m=r.randrange(5);T=chart(n,r);U=chart(m,r)
 colors=[r.randrange(3) for i in range(n)]
 ph=[[bool(r.randrange(2)) for j in range(3)] for i in range(m)]
 Z=T+[v<<n for v in U]
 for i,a in enumerate(colors):
  for j in range(m):
   orient(Z, *((n+j,i) if ph[j][a] else (i,n+j)))
 assert H(Z)==contract(covers(T,colors),U,ph)
print("60 literal colored contractions")
