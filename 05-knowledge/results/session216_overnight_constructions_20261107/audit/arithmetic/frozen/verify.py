"""Replay exact coefficient certificates and independent literal subset-DP checks."""
import json,itertools,pathlib,hashlib
BASE=pathlib.Path(__file__).resolve().parent
W=json.loads((BASE/'interval_witnesses.json').read_text())
def require(b,msg):
 if not b:raise RuntimeError(msg)
def valid(a):
 n=len(a)
 return all(not(a[i]>>i&1) for i in range(n)) and all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i in range(n) for j in range(i+1,n))
def coefficients(a,v):
 n=len(a);c=[0]*(n+1)
 for p in itertools.permutations(i for i in range(n) if i!=v):
  forced=optional=0
  for g in range(n):
   l=p[g-1] if g else None;r=p[g] if g<n-1 else None
   e=l is None or r is None or bool(a[l]>>r&1)
   f=(l is None or bool(a[l]>>v&1)) and (r is None or bool(a[v]>>r&1))
   if not(e or f):break
   if f:
    if e:optional+=1
    else:forced+=1
  else:
   from math import comb
   for j in range(optional+1):c[forced+j]+=comb(optional,j)
 return c
def H_all(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)];hs=[0]*(1<<n);hs[0]=1
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for v,x in enumerate(dp[s]):
   if not x:continue
   z=a[v]&~s
   while z:
    t=z&-z;w=t.bit_length()-1;dp[s|t][w]+=x;z-=t
  hs[s]=sum(dp[s])
 return hs
def substitute(a,v,u):
 old=[x for x in range(len(a)) if x!=v];m=len(old);b=[0]*(m+len(u))
 for i,x in enumerate(old):
  for j,y in enumerate(old):
   if a[x]>>y&1:b[i]|=1<<j
  for j in range(len(u)):
   if a[x]>>v&1:b[i]|=1<<(m+j)
   else:b[m+j]|=1<<i
 for i,x in enumerate(u):b[m+i]|=x<<m
 return b
def enum(n):
 for bits in range(1<<(n*(n-1)//2)):
  a=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if bits>>k&1:a[i]|=1<<j
    else:a[j]|=1<<i
    k+=1
  yield a
receipts=[];checks=0
seeds=[u for n in range(1,5) for u in enum(n)]
for w in W:
 a,v=w['rows'],w['root'];require(valid(a),'bad core');cy=w['hamilton_cycle'];require(sorted(cy)==list(range(len(a))) and all(a[x]>>y&1 for x,y in zip(cy,cy[1:]+cy[:1])),'bad Hamilton cycle');c=coefficients(a,v)
 require(c[1]==w['h'] and c[2]==36 and all(z==0 for z in c[3:]),'coefficient mismatch')
 for u in seeds:
  hs=H_all(u);h=hs[-1];full=len(hs)-1;q=sum(hs[s]*hs[full^s] for s in range(full+1))
  out=H_all(substitute(a,v,u))[-1]
  require(out==36*q+(w['h']-72)*h,'arbitrary seed identity');checks+=1
 for k in range(1,7):
  u=[sum(1<<j for j in range(i+1,k)) for i in range(k)]
  out=H_all(substitute(a,v,u))[-1]
  require(out==36*(1<<k)+w['h']-72,'clone identity');checks+=1
 receipts.append({'h':w['h'],'coefficients':c})
require([w['h'] for w in W]==list(range(65,86,2)),'interval')
print(json.dumps({'status':'PASS','literal_checks':checks,'seed_count':len(seeds),'coefficient_certificates':receipts},indent=2))
