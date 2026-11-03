"""Literal tournament counting recurrences; no dependencies."""
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
def ordinal(parts):
 n=sum(map(len,parts));a=[];start=0
 for p in parts:
  end=start+len(p);later=((1<<n)-1)^((1<<end)-1)
  a.extend((v<<start)|later for v in p);start=end
 return a
def strong(a):
 mask=(1<<len(a))-1
 for start in range(len(a)):
  seen=1<<start
  for _ in a:
   old=seen
   for i in range(len(a)):
    if seen>>i&1:seen|=a[i]
   if seen==old:break
  if seen!=mask:return False
 return True
def dual(a):return tuple(((1<<len(a))-1)^v^(1<<i) for i,v in enumerate(a))
def wrap(a,r):
 n=len(a);old=ordinal([a,dual(a)]);z=((1<<r)-1)<<(2*n)
 return tuple(v|(z if i>=n else 0) for i,v in enumerate(old))+tuple(((1<<n)-1)|sum(1<<(2*n+j) for j in range(i+1,r)) for i in range(r))
def induced(a,S):return tuple(sum(1<<j for j,y in enumerate(S) if a[x]>>y&1) for x in S)
def survive(a,D):return strong(induced(a,[i for i in range(len(a)) if i not in D]))
def f(j):return j**5+6*j**3+2*j
