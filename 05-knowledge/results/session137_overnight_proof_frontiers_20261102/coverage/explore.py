from functools import lru_cache
from itertools import permutations

def tournament(n,code):
 a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if code&1:a[i]|=1<<j
   else:a[j]|=1<<i
   code>>=1
 return tuple(a)
def extend(a,s):
 n=len(a);return tuple(a[i]|((1<<n) if not(s>>i&1) else 0) for i in range(n))+(s,)
def H(a):
 n=len(a)
 if not n:return 1
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for m in range(1,1<<n):
  for i,v in enumerate(dp[m]):
   if not v:continue
   z=a[i]&~m
   while z:
    b=z&-z;z-=b;j=b.bit_length()-1;dp[m|b][j]+=v
 return sum(dp[-1])
def strong(a):
 n=len(a)
 if n<3:return False
 def reach(rev):
  seen=1
  while True:
   new=seen
   for i in range(n):
    if seen>>i&1:
     if not rev:new|=a[i]
     else:new|=sum(1<<j for j in range(n) if a[j]>>i&1)
   if new==seen:return seen==(1<<n)-1
   seen=new
 return reach(False) and reach(True)
def literal(a):
 return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(a)-1)) for p in permutations(range(len(a))))
def endpoint_Q(a):
 n=len(a);s=[0]*n;e=[0]*n;Q=[[0]*n for _ in a]
 for p in permutations(range(n)):
  defects=[i for i in range(n-1) if not(a[p[i]]>>p[i+1]&1)]
  if not defects:
   s[p[0]]+=1;e[p[-1]]+=1
   for i in range(n-1):Q[p[i]][p[i+1]]+=1
  elif len(defects)==1:
   i=defects[0];Q[p[i]][p[i+1]]+=1
 return s,e,Q

if __name__=='__main__':
 for n in range(3,6):
  images=[];hh=set()
  for code in range(1<<(n*(n-1)//2)):
   a=tournament(n,code)
   if not strong(a):continue
   im={H(extend(a,s)) for s in range(1,(1<<n)-1)}
   images.append((code,H(a),im));hh|=im
  print('n',n,'parents',len(images),'union',sorted(hh),'intersection',sorted(set.intersection(*(t[2] for t in images))))
  print('pair_disjoint',next(((x[:2],y[:2]) for x in images for y in images if not x[2]&y[2]),None))

def endpoint_Q_dp(a):
 n=len(a);N=1<<n;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)]
 for i in range(n):E[1<<i][i]=S[1<<i][i]=1
 for m in range(1,N):
  for i in range(n):
   if not(m>>i&1):continue
   z=m^(1<<i)
   if not z:continue
   E[m][i]=sum(E[z][j] for j in range(n) if z>>j&1 and a[j]>>i&1)
   S[m][i]=sum(S[z][j] for j in range(n) if z>>j&1 and a[i]>>j&1)
 q=[[0]*n for _ in a]
 for m in range(1,N-1):
  mm=N-1-m
  for i in range(n):
   if not E[m][i]:continue
   for j in range(n):
    if S[mm][j]:q[i][j]+=E[m][i]*S[mm][j]
 return S[-1],E[-1],q

def c3_blowup(q):
 k=len(q);a=[0]*(3*k)
 for i in range(k):
  for r in range(3):
   v=3*i+r;a[v]|=1<<(3*i+(r+1)%3)
   for j in range(k):
    if q[i]>>j&1:a[v]|=7<<(3*j)
 return tuple(a)
