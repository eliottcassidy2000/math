#!/usr/bin/env python3
import json,math
from functools import lru_cache

def chain(parts):
 n=sum(parts)+1;a=[sum(1<<j for j in range(i+1,n)) for i in range(n)];s=0
 for p in parts:
  a[s]^=1<<(s+p);a[s+p]|=1<<s;s+=p
 return a

def hp(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   c=dp[s][v]
   if not c:continue
   z=a[v]&~s
   while z:
    b=z&-z;z-=b;dp[s|b][b.bit_length()-1]+=c
 return sum(dp[-1])

def comps(t):
 if t==0:yield ();return
 for p in range(2,t+1):
  for rest in comps(t-p):yield (p,)+rest

physical=0
for t in range(11):
 for ps in comps(t):
  assert hp(chain(ps))==math.prod(2**(p-1)+1 for p in ps)
  physical+=1

def primes(n):
 d=2
 while d*d<=n:
  if n%d==0:
   yield d
   while n%d==0:n//=d
  d+=1
 if n>1:yield n
primitive=[]
for m in range(1,21):
 q=2**m+1
 good=[p for p in primes(q) if pow(2,2*m,p)==1 and all(pow(2,j,p)!=1 for j in range(1,2*m))]
 assert bool(good)==(m!=3)
 primitive.append(dict(m=m,q=q,primitive_primes=good))
assert 2**3+1==(2**1+1)**2

def recognize(h):
 if h<1:return None
 if h==1:return (1,())
 rem=h;cost=1;f=[]
 for m in range((h-1).bit_length()-1,0,-1):
  if m==3:continue
  q=(1<<m)+1;e=0
  while rem%q==0:rem//=q;cost+=m+1;e+=1
  if e:f.append((m,e))
 return (cost,tuple(f)) if rem==1 else None

# Independent ascending-product closure, tracking every possible weighted cost.
cap=20000;generators=[(2**m+1,m+1) for m in range(1,cap.bit_length()) if m!=3 and 2**m+1<=cap]
costs={1:{0}}
for h in range(1,cap+1):
 if h not in costs:continue
 for q,w in generators:
  if h*q<=cap:costs.setdefault(h*q,set()).update(c+w for c in costs[h])
for h in range(1,cap+1):
 r=recognize(h)
 assert (r is not None)==(h in costs),h
 if r:
  assert costs[h]=={r[0]-1},h
  assert math.prod((2**m+1)**e for m,e in r[1])==h

# Canonical partitions and physical value injectivity through weight24.
products={}
def partition_visit(rem,start,ps=()):
 if rem==0:
  v=math.prod(2**(p-1)+1 for p in ps)
  assert v not in products or products[v]==ps,(v,ps,products.get(v))
  products[v]=ps;return 1
 count=0
 for p in range(max(2,start),rem+1):
  if p!=4:count+=partition_visit(rem-p,p,ps+(p,))
 return count
exact=[partition_visit(t,2) for t in range(25)]
for h,ps in products.items():assert recognize(h)[0]==sum(ps)+1

N=100;p=[0]*(N+1);p[0]=1
for w in range(1,N+1):
 for t in range(w,N+1):p[t]+=p[t-w]
q=[0]*(N+1);q[0]=1
for w in range(2,N+1):
 if w==4:continue
 for t in range(w,N+1):q[t]+=q[t-w]
def at(t):return p[t] if t>=0 else 0
for t in range(N+1):
 assert q[t]==at(t)-at(t-1)-at(t-4)+at(t-5)
 assert sum(q[:t+1])==at(t)-at(t-4)
 if t<=24:assert exact[t]==q[t]
for M in range(3,100):assert sum(m+1 for m in range(1,M+1) if m!=3)==M*(M+3)//2-4
print(json.dumps({'status':'PASS','physical_chain_compositions_through_order11':physical,'primitive_controls':primitive,'distinct_canonical_products_through_weight24':len(products),'greedy_inputs_checked':cap,'greedy_members_found':len(costs),'partition_rows':[{'n':n,'chain_count':q[n-1],'padded_count':sum(q[:n])} for n in range(1,31)]},indent=2))
