from itertools import combinations,permutations
from collections import defaultdict

def check(cond,msg):
 if not cond: raise AssertionError(msg)

def tournament(n,b):
 a=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if b>>k&1:a[j]|=1<<i
  else:a[i]|=1<<j
 return a

def hp_counts(a):
 n=len(a); m=1<<n; dp=[[0]*n for _ in range(m)]; h=[0]*m;h[0]=1
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,m):
  for v in range(n):
   if s>>v&1:
    rem=s^(1<<v)
    if rem:dp[s][v]=sum(dp[rem][u] for u in range(n) if rem>>u&1 and a[u]>>v&1)
    else:dp[s][v]=1
  h[s]=sum(dp[s])
 return h

def cycles(a):
 n=len(a);res=[]
 for k in range(3,n+1,2):
  for s in combinations(range(n),k):
   for p in permutations(s[1:]):
    q=(s[0],)+p
    if all(a[q[i]]>>q[(i+1)%k]&1 for i in range(k)):
     res.append((sum(1<<v for v in s),k))
 return res

def poly_cycle(a):
 n=len(a); cs=cycles(a); out=defaultdict(int)
 def dfs(i,used,c):
  if i==len(cs):out[n-used.bit_count()+c]+=2**c;return
  dfs(i+1,used,c)
  mask,_=cs[i]
  if not mask&used:dfs(i+1,used|mask,c+1)
 dfs(0,0,0)
 return dict(out),cs

def eval_poly(poly,r):return sum(coef*r**k for k,coef in poly.items())
def f_from_subsets(h,n,r):
 if r==0:return int(n==0)
 allm=(1<<n)-1
 def walk(s,k):
  if k==1:return h[s]
  total=0;t=s
  while True:
   total+=h[t]*walk(s^t,k-1)
   if not t:break
   t=(t-1)&s
  return total
 return walk(allm,r)

for n in range(1,6):
 for b in range(1<<(n*(n-1)//2)):
  a=tournament(n,b);h=hp_counts(a);pol,cs=poly_cycle(a)
  for r in range(1,5):check(eval_poly(pol,r)==f_from_subsets(h,n,r),(n,b,r,pol))
  check(all((k-n)%2==0 for k in pol),(n,b,pol))
  f2=eval_poly(pol,2)
  if n>=2:check(f2%4==0,(n,b,f2))
  if n>=4 and n%2==0:
   c=sum(k==n-1 for _,k in cs)
   check(f2%16==(8*c)%16,(n,b,f2,c))
  if n>=5 and n%2==1:
   c=sum(k==n for _,k in cs)
   check(f2%16==(4*c)%16,(n,b,f2,c))
 print('n=',n,'checked=',1<<(n*(n-1)//2),flush=True)
print('PASS')
