"""Independent residue-based arity recurrence and named-multiplier audit."""
from fractions import Fraction as Q
from math import gcd,comb,prod
from itertools import combinations
import json
checks=0

def check(a,b,tag):
 global checks
 checks+=1
 if a!=b:raise RuntimeError((tag,a,b))
def factor(n):
 out=[];p=2
 while p*p<=n:
  a=0
  while n%p==0:n//=p;a+=1
  if a:out.append((p,a))
  p+=1
 if n>1:out.append((n,1))
 return out

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def formula(n,r):
 fs=factor(n);z=Q((n-1)**r)
 for mask in range(1,2**len(fs)):
  lam=Q(n);pol=Q(1)
  for i,(p,a) in enumerate(fs):
   if mask>>i&1:
    lam*=Q(p-1,p)
    pol*=sum((Q(prod(range(r,r+t)),prod(range(1,t+1))*p**t) for t in range(a)),Q(0))
  z+=(-1)**mask.bit_count()*lam**r*pol
 return z
for n in range(2,61):
 fs=factor(n);dist=[0]*n;dist[1]=1;values=[]
 for r in range(15):
  values.append(dist[0]);check(dist[0],formula(n,r),('formula',n,r))
  new=[0]*n
  for a,v in enumerate(dist):
   for x in range(1,n):new[a*x%n]+=v
  dist=new
 if len(fs)==1 and fs[0][1]==1:continue
 annih=[-(n-1),1]
 for mask in range(1,2**len(fs)):
  lam=Q(n);depth=0
  for i,(p,a) in enumerate(fs):
   if mask>>i&1:lam*=Q(p-1,p);depth+=a-1
  for _ in range(depth+1):annih=mul(annih,[-lam,1])
 for r in range(len(values)-len(annih)+1):check(sum(c*values[r+i] for i,c in enumerate(annih)),0,('recurrence',n,r))
 divisors=[d for d in range(1,n+1) if n%d==0]
 for a,b in combinations(divisors,2):
  check(any((a*x%n==0)!=(b*x%n==0) for x in range(1,n)),True,('named distinguish',n,a,b))
 order=2**len(fs)+2**(len(fs)-1)*sum(a-1 for p,a in fs)
 check(len(divisors)==order,sum(a>1 for p,a in fs)<=1,('dimension gap',n))
# Arbitrary correlated CRT source, with signed weights and the zero residue retained.
for n in [4,6,8,12,18,36]:
 fs=factor(n)
 for weights in [[(x%3)-1 for x in range(n)],[int(x==0)+int(x==1)*2-int(x==2) for x in range(n)]]:
  annih=[1]
  for mask in range(2**len(fs)):
   selected=[(p,a) for i,(p,a) in enumerate(fs) if mask>>i&1]
   lam=sum(weights[x] for x in range(n) if all(x%p for p,a in selected))
   exponent=1+sum(a-1 for p,a in selected)
   for _ in range(exponent):annih=mul(annih,[-lam,1])
  state=[0]*n;state[1]=1;values=[]
  for _ in range(len(annih)+5):
   values.append(state[0]);new=[0]*n
   for x,v in enumerate(state):
    for y,w in enumerate(weights):new[x*y%n]+=v*w
   state=new
  for r in range(6):check(sum(a*values[r+j] for j,a in enumerate(annih)),0,('weighted ring',n,weights,r))
# Divisor probes invert all source masses, with total mass supplied.
def mobius(k):
 f=factor(k)
 return 0 if any(a>1 for p,a in f) else (-1)**len(f)
for n in range(4,61):
 fs=factor(n)
 if len(fs)==1 and fs[0][1]==1:continue
 ds0=[d for d in range(1,n+1) if n%d==0]
 masses={d:Q((3*d)%7-3,5) for d in ds0}
 probes={d:sum((m for g,m in masses.items() if g%d==0),Q(0)) for d in ds0}
 for g in ds0:
  recovered=sum((mobius(d//g)*probes[d] for d in ds0 if d%g==0),Q(0))
  check(recovered,masses[g],('Mobius restoration',n,g))
# All-order hidden mode certificate, with independently built residue multiplication.
ds=[1,2,3,4,6,9,12,18,36];h=[0,-3,2,3,1,-2,-3,2,0]
M=[[sum(gcd(d*x,36)==e for x in range(1,36)) for e in ds] for d in ds]
check([sum(h[i]*M[i][j] for i in range(9)) for j in range(9)],[12*v for v in h],'eigenmode')
check(sum(h[i] for i,d in enumerate(ds) if d*2%36==0),2,'named separator')
check(h[-1],0,'unobserved zero')
print(json.dumps({'status':'PASS','checks':checks,'moduli_through':60,'producer_imports':False},sort_keys=True))
