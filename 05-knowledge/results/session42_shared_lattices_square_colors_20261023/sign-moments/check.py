import itertools, math, json
from collections import Counter
import numpy as np

def need(b,msg):
 if not b: raise RuntimeError(msg)
def fwht(a):
 h=1
 while h<len(a):
  r=a.reshape(-1,2*h);x=r[:,:h].copy();y=r[:,h:].copy();r[:,:h]=x+y;r[:,h:]=x-y;h*=2
 return a

def extremal(n,mask,E):
 deg=[0]*n;adj=[[] for _ in range(n)]
 for k,(a,b) in enumerate(E):
  if mask>>k&1:deg[a]+=1;deg[b]+=1;adj[a].append(b);adj[b].append(a)
 active=[v for v in range(n) if deg[v]]
 if sorted(deg[v] for v in active)==[1,1,2]:return True
 if sorted(deg[v] for v in active)!=[2,2,2,2]:return False
 seen={active[0]};todo=list(seen)
 while todo:
  for w in adj[todo.pop()]:
   if w not in seen:seen.add(w);todo.append(w)
 return len(seen)==4
out=[]; arrays={}
for n in range(2,8):
 E=list(itertools.combinations(range(n),2));a=np.zeros(1<<len(E),dtype=np.int64)
 for p in itertools.permutations(range(n)):
  rank=[0]*n
  for k,v in enumerate(p):rank[v]=k
  mask=sum(1<<k for k,(i,j) in enumerate(E) if rank[i]>rank[j]);a[mask]+=1
 fwht(a);arrays[n]=a
 fac=math.factorial(n);need(int(a[0])==fac,'empty moment')
 maxabs=int(np.max(np.abs(a[1:])));need(3*maxabs<=fac,'one third bound')
 equality=[int(x) for x in np.flatnonzero(3*np.abs(a)==fac)]
 need(all(extremal(n,x,E) for x in equality),'unexpected equality shape')
 expected=3*math.comb(n,3)+(3*math.comb(n,4) if n>=4 else 0)
 need(len(equality)==expected,'equality count')
 # Reversal proves every odd-edge coefficient zero, checked exhaustively.
 need(all(a[x]==0 for x in range(1,len(a)) if x.bit_count()%2),'odd moment')
 recurrence=0
 if n>=3:
  # Exhaustive to6; deterministic spread to7. Integer signed sums satisfy S_n=sum eps S_(n-1).
  tests=range(len(a)) if n<=6 else sorted({0,len(a)-1,*((i*104729)%len(a) for i in range(1000))})
  low=list(itertools.combinations(range(n-1),2));idx={e:k for k,e in enumerate(low)}
  for mask in tests:
   total=0
   for v in range(n):
    rem=[x for x in range(n) if x!=v];rename={x:k for k,x in enumerate(rem)};small=0;sign=1
    for k,(i,j) in enumerate(E):
     if not(mask>>k&1):continue
     if i==v or j==v:
      if j==v:sign=-sign
     else:small |= 1<<idx[rename[i],rename[j]]
    total+=sign*int(arrays[n-1][small])
   need(total==a[mask],'first vertex recurrence');recurrence+=1
 out.append(dict(n=n,graphs=len(a),permutations=fac,max_abs_numerator=maxabs,denominator=fac,equality_count=len(equality),recurrence_checks=recurrence))
print(json.dumps({'status':'PASS','scope':'all graph sign moments through n7 via exact int64 FWHT; theorem proved separately','orders':out},indent=2))
