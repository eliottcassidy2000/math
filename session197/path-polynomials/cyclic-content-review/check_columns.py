import json, math, random
from functools import reduce

def must(c,s):
 if not c:raise RuntimeError(s)
def multiply(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def base(k,odd):
 p=[0,1] if odd else [0,0,1]
 for j in range(1,k+1):p=multiply(p,[-j*j,0,1])
 return p
def poly(u,odd):
 out=[0]*(2*len(u)+(0 if odd else 1))
 for k,x in enumerate(u):
  for i,v in enumerate(base(k,odd)):out[i]+=x*v
 return out
def val(p,x):
 y=0
 for z in reversed(p):y=y*x+z
 return y
def newton(p):
 values=[val(p,x)for x in range(len(p))];f=[]
 while values:
  f.append(values[0]);values=[values[i+1]-values[i]for i in range(len(values)-1)]
 return f
def kernel(f):
 # Independently count the three possible starts of every cyclic run-color word.
 n=len(f)-1;K=[[0]*(n+1)for _ in range(n+1)]
 for r in range(1,n+2):
  for s in range(1,n+2):
   for t in range(1,n+1):
    if max(r,s,t)-min(r,s,t)>1:continue
    count=0
    for start in range(3):
     counts=[0,0,0]
     for j in range(r+s+t):counts[(start+j)%3]+=1
     count+=counts==[r,s,t]
    K[r-1][s-1]+=count*math.factorial(r)*math.factorial(s)*f[t]
 return K
def gcd(xs):return reduce(math.gcd,xs,0)
def bank(K,L):
 S=[[0]*(L+1)for _ in range(L+1)];S[0][0]=1
 for i in range(1,L+1):
  for j in range(1,i+1):S[i][j]=S[i-1][j-1]+j*S[i-1][j]
 return [sum(S[i][r]*K[r-1][s-1]*S[j][s]for r in range(1,L+1)for s in range(1,L+1))for i in range(1,L+1)for j in range(1,L+1)]
rng=random.Random(1971445);controls=columns=0
for odd in [False,True]:
 for d in range(1,11):
  for trial in range(60):
   u=[rng.randrange(-20,21)for _ in range(d)];u[0]=2*rng.randrange(-20,21)+1
   P=poly(u,odd);f=newton(P);K=kernel(f)
   gamma=[math.factorial(2*k+1)*(1 if odd else k+1)for k in range(d)]
   terms=[gamma[k]*u[k]for k in range(d)]
   g=gcd(f);G=gcd(math.factorial(r-1)*f[r]for r in range(1,len(f)));C=gcd(x for row in K for x in row)
   a=[math.factorial(k)**2*terms[k]for k in range(d)]
   must(g==gcd(terms)and G==gcd(math.factorial(k)*terms[k]for k in range(d))and C==gcd(a),'three exact central contents')
   must(C==gcd(bank(K,max(2,d-1))),'short Stirling bank content')
   must(G*G%(g*C)==0 and g**3%(G*G)==0 and g*g%C==0 and G**4%(C**3)==0,'valuation hierarchy')
   y=a[1]if d>1 else 0
   must((K[0][0],K[0][1])==((3*u[0]+y,2*u[0]+2*y)if odd else(5*u[0]+y,6*u[0]+2*y)),'exact boundary matrix')
   controls+=1
 for k in range(12):
  f=newton(base(k,odd));K=kernel(f);gamma=math.factorial(2*k+1)*(1 if odd else k+1);div=math.factorial(k)**2*gamma
  must(all(x%div==0 for row in K for x in row),'integral normalized kernel column')
  if k>=2:must(K[k-1][k-1]==div,'unit diagonal beyond boundary')
  columns+=1
host=kernel(newton([0,1,0,1]));must(gcd(x for row in host for x in row)==4,'even-value counterexample')
print(json.dumps({'status':'PASS','signed_odd_value_profiles':controls,'normalized_kernel_columns':columns,'even_value_counterexample':4,'short_bank_max_size':9},sort_keys=True,indent=2))
