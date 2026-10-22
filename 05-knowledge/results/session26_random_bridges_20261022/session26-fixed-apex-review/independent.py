import itertools as it, math, json, hashlib
from fractions import Fraction as F
from collections import Counter
from functools import lru_cache

def need(test, msg):
 if not test: raise RuntimeError(msg)
def adj(n,bits,pairs):
 a=[0]*n
 for t,(i,j) in enumerate(pairs):
  if bits>>t&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return tuple(a)
def data(a):
 n=len(a);full=(1<<n)-1
 trans=[False]*(1<<n);trans[0]=True
 for s in range(1,1<<n):
  trans[s]=any(a[i]&(s^(1<<i))==(s^(1<<i)) and trans[s^(1<<i)] for i in range(n) if s>>i&1)
 tau=n-max(s.bit_count() for s in range(1<<n) if trans[s])
 dp=[[0]*n for _ in trans]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   if not dp[s][i]:continue
   targets=a[i]&~s&full
   while targets:
    q=targets&-targets;j=q.bit_length()-1;targets-=q;dp[s|q][j]+=dp[s][i]
 return trans,tau,sum(dp[-1])
def canon(a):
 n=len(a);ps=list(it.permutations(range(n)))
 codes=[tuple((a[p[i]]>>p[j])&1 for i in range(n) for j in range(i+1,n)) for p in ps]
 c=min(codes)
 return c,codes.count(c)
rows=[]
for n in range(2,6):
 pairs=list(it.combinations(range(n),2)); classes={}
 for bits in range(1<<len(pairs)):
  a=adj(n,bits,pairs);c,aut=canon(a)
  if c not in classes:classes[c]=(a,aut,data(a))
 for k in range(1,n):
  m=n-k;free=[p for p in pairs if p[1]>=m];cnt=Counter();bad=0;coll=0;sumh=0
  for bits in range(1<<len(free)):
   a=[sum(1<<j for j in range(i+1,m)) if i<m else 0 for i in range(n)]
   for t,(i,j) in enumerate(free):
    if bits>>t&1:a[i]|=1<<j
    else:a[j]|=1<<i
   a=tuple(a);tr,t,h=data(a);c,aut=canon(a);cnt[c]+=1
   f=sum(tr[s] for s in range(1<<n) if s.bit_count()==m);bad+=f>1
   coll+=len({a[i]&((1<<m)-1) for i in range(m,n)})<k
   need(h<=math.factorial(k)*(k+1)**m,'H upper');sumh+=h
  for c,(a,aut,(tr,t,h)) in classes.items():
   f=sum(tr[s] for s in range(1<<n) if s.bit_count()==m)
   need(cnt[c]==F(math.factorial(k)*f,aut),'per-class multiplicity')
  bound=sum(math.comb(k,r)*math.comb(m,r)*F(m-r+1,2**(m-r))**r for r in range(1,min(k,m)+1))
  need(F(bad,1<<len(free))<=bound,'uniqueness bound')
  need(F(coll,1<<len(free))<=F(math.comb(k,2),2**m),'collision bound')
  v=sum((F(h,aut) for a,aut,(tr,t,h) in classes.values() if t<=k),F(0))
  gap=F(sumh,math.factorial(k))-v
  need(0<=gap<=bad*(k+1)**m,'mass correction')
  rows.append(dict(n=n,k=k,classes=len(cnt),presentations=sum(cnt.values()),nonunique=bad,union_bound=str(bound),weighted_mass=str(v),gap=str(gap)))
# Every Hamilton-labelled tournament with distinguished path 0->...->n-1.
# This exhausts all tournament/path pairs up to simultaneous relabelling.
transfer=[]
for n in range(1,8):
 free=[(i,j) for i in range(n) for j in range(i+2,n)];hist=Counter()
 for bits in range(1<<len(free)):
  a=[(1<<(i+1)) if i+1<n else 0 for i in range(n)]
  for t,(i,j) in enumerate(free):
   if bits>>t&1:a[i]|=1<<j
   else:a[j]|=1<<i
  b=[(1<<(i+1)) if i+1<n else 0 for i in range(n)]
  for i,j in free:
   if a[i]>>j&1:b[j]|=1<<i
   else:b[i]|=1<<j
  ta=data(a)[1];tb=data(b)[1]
  need(ta+tb>=n//3,'lower transfer');need(tb<=(n+ta)//2,'upper transfer')
  hist[ta,tb]+=1
 transfer.append(dict(n=n,path_presentations=1<<len(free),feedback_pairs={str(k):v for k,v in sorted(hist.items())}))
print(json.dumps(dict(status='PASS',class_checks=rows,transfer_checks=transfer),indent=2))
