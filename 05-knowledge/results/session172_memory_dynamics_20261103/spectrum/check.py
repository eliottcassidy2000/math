#!/usr/bin/env python3
"""Independent literal subset dynamic programming; no repository imports."""
import itertools as it,json,math
B=[(16,1,3,7,14),(20,1,2,7,14)]
def require(ok,label):
 if not ok:raise RuntimeError(label)
def paths(a):
 n=len(a);h=[0]*(1<<n);h[0]=1;dp=[[0]*n for _ in h]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for j,v in enumerate(dp[s]):
   h[s]+=v
   if not v:continue
   cand=a[j]&~s
   while cand:
    b=cand&-cand;cand-=b;k=b.bit_length()-1;dp[s|b][k]+=v
 return h
def cover(a):
 h=paths(a);n=len(a);p=[[0]*(n+1) for _ in h];p[0][0]=1
 for s in range(1,len(h)):
  least=s&-s;rest=s^least;t=rest
  while True:
   A=t|least;rem=s^A
   if h[A]:
    for k in range(1,rem.bit_count()+2):p[s][k]+=h[A]*p[rem][k-1]
   if not t:break
   t=(t-1)&rest
 return p[-1]
def evalF(p,u):
 ans=0;fall=1
 for k,c in enumerate(p):
  if k:fall*=u-k+1
  ans+=c*fall
 return ans
def F(u):return u**5+6*u**3+2*u
def ext(a):
 n=len(a);return tuple(v|((1<<n) if i else 0) for i,v in enumerate(a))+(1,)
def ordinal(parts):
 n=sum(map(len,parts));a=[];start=0
 for p in parts:
  end=start+len(p);later=((1<<n)-1)^((1<<end)-1)
  a.extend((v<<start)|later for v in p);start=end
 return a
def wrap(parts):
 a=ordinal(parts);n=len(a);m=len(parts[0])
 return tuple(v|((1<<n) if i>=m else 0) for i,v in enumerate(a))+((1<<m)-1,)
def strong(a,deleted=None):
 mask=((1<<len(a))-1)^(0 if deleted is None else 1<<deleted)
 for start in range(len(a)):
  if not(mask>>start&1):continue
  seen=1<<start
  for _ in a:
   old=seen
   for i in range(len(a)):
    if seen>>i&1:seen|=a[i]&mask
   if seen==old:break
  if seen!=mask:return False
 return True
def formula(ps,u):
 f=[evalF(p,u) for p in ps];lo=[evalF(p,u-1) for p in ps];hi=[evalF(p,u+1) for p in ps]
 return u*((lo[0]-f[0])*math.prod(lo[1:])+(3*f[0]-hi[0]-lo[0])*math.prod(f[1:])+(hi[0]-f[0])*math.prod(hi[1:]))
ps=[cover(a) for a in B];qs=[cover(ext(a)) for a in B]
require(ps[0]==ps[1]==[0,9,33,31,10,1],'seed profile')
require([[evalF(p,u) for u in [1,2]] for p in qs]==[[17,224],[19,256]],'query profiles')
records=[]
for bits in it.product(range(2),repeat=2):
 parts=[B[b] for b in bits]
 for address in [None,0,1]:
  pp=parts.copy();profiles=[ps[b] for b in bits]
  if address is not None:pp[address]=ext(pp[address]);profiles[address]=qs[bits[address]]
  a=wrap(pp);p=cover(a)
  require(strong(a),'strong wrapper')
  require([i for i in range(len(a)) if not strong(a,i)]==[len(a)-1],'unique articulation')
  for u in range(len(a)+1):require(evalF(p,u)==formula(profiles,u),'full polynomial formula')
  M=2
  expect=5787 if address is None else (207+30*bits[0])*84-(173+26*bits[0])*9 if address==0 else (16800+2400*bits[1])-(969+114*bits[1])
  require(p[1]==expect,'local read formula')
  records.append({'bits':bits,'address':address,'order':len(a),'H':p[1],'profile':p})
# Labelled isometry is inherited literally: all wrapper arcs are fixed.
for bits in it.product(range(2),repeat=2):
 for other in it.product(range(2),repeat=2):
  a=wrap([B[b] for b in bits]);b=wrap([B[b] for b in other])
  require(sum((x^y).bit_count() for x,y in zip(a,b))//2==sum(x!=y for x,y in zip(bits,other)),'arc isometry')
print(json.dumps({'status':'verified','literal_records':records,'query_profiles':qs,'examples':[{'M':M,'order':5*M+1,'H':75*84**(M-1)-57*9**(M-1)} for M in range(2,9)]},indent=2,sort_keys=True))
