from math import comb, factorial

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def qint(n,negative=False):return [(-1)**i if negative else 1 for i in range(n)]
def dp_poly(n,edges):
 adj=[sum(1<<u for u in range(v+1,n) if (v,u) in edges) for v in range(n)]
 dp={0:[1]}
 for mask in range(1,1<<n):
  deg=mask.bit_count()*(mask.bit_count()-1)//2
  c=[0]*(deg+1)
  for v in range(n):
   if mask>>v&1:
    prev=mask^(1<<v); shift=(prev>>(v+1)).bit_count()
    sign=(-1)**((prev&adj[v]).bit_count())
    for i,a in enumerate(dp[prev]):c[i+shift]+=sign*a
  dp[mask]=c
 return dp[(1<<n)-1]
def threshold(r,j):
 n=r; edges=set()
 for h in range(j):
  for u in range(n):edges.add((u,n))
  n+=1
  if h<j-1:n+=1
 return n,edges
for r,j in [(3,2),(3,4),(3,6),(2,2),(4,3)]:
 n,E=threshold(r,j); actual=dp_poly(n,E); expected=[1]
 for i in range(1,r+1):expected=mul(expected,qint(i))
 for t in range(1,j+1):expected=mul(expected,qint(r+2*t-1,True))
 for t in range(1,j):expected=mul(expected,qint(r+2*t))
 if actual!=expected:raise ValueError(('factor',r,j))
 jets=[sum(a*comb(i,h) for i,a in enumerate(actual) if i>=h) for h in range(j+1)]
 if r%2 and (any(jets[:j]) or jets[j]*2**j!=(-1)**j*factorial(n)):raise ValueError(('jet',r,j,jets))
 print('r,j,n',r,j,n,'coefficients',len(actual),'q1jets',jets,'PASS')
for mask in (255,2015):
 pairs=[(i,j) for i in range(6) for j in range(i+1,6)]
 D=dp_poly(6,{p for h,p in enumerate(pairs) if mask>>h&1})
 jets=[sum(a*comb(i,h) for i,a in enumerate(D) if i>=h) for h in range(3)]
 print('mask',mask,'jets',jets,'D',D)
