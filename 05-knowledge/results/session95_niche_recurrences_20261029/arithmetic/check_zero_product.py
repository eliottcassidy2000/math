"""Exact independent controls for zero-product arity recurrences (stdlib only)."""
from fractions import Fraction as F
from itertools import product, combinations
from math import gcd, prod, comb
import json

checks=0

def require(v,msg):
 global checks
 checks+=1
 if not v: raise RuntimeError(msg)

def factors(n):
 out=[];p=2
 while p*p<=n:
  if n%p==0:
   a=0
   while n%p==0:n//=p;a+=1
   out.append((p,a))
  p+=1
 if n>1:out.append((n,1))
 return out

def divisors(n):return [d for d in range(1,n+1) if n%d==0]

def ps(r,p,a):
 return sum((F(prod(range(r,r+t)),prod(range(1,t+1))*p**t) for t in range(a)),F(0))

def formula(n,r):
 fs=factors(n);out=F((n-1)**r)
 for mask in range(1,1<<len(fs)):
  lam=n;pol=F(1)
  for i,(p,a) in enumerate(fs):
   if mask>>i&1:lam=lam//p*(p-1);pol*=ps(r,p,a)
  out+=(-1)**mask.bit_count()*lam**r*pol
 require(out.denominator==1,('integrality',n,r))
 return int(out)

def multiply_poly(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def recurrence(n,weights=None,floor=False):
 fs=factors(n);q=[F(1)]
 for mask in range(1<<len(fs)):
  ds=sum(a-1 for i,(p,a) in enumerate(fs) if mask>>i&1)
  if weights is None:
   lam=n-1 if mask==0 else n
   if mask:
    for i,(p,a) in enumerate(fs):
     if mask>>i&1:lam=lam//p*(p-1)
   if mask==0 and floor:lam=n*(n-1)//2
  else:
   lam=sum(w for x,w in enumerate(weights) if all(x%p for i,(p,a) in enumerate(fs) if mask>>i&1))
  for _ in range(ds+1):q=multiply_poly(q,[-lam,1])
 return q

def apply_recurrence(q,seq):
 for r in range(len(seq)-len(q)+1):
  require(sum(q[j]*seq[r+j] for j in range(len(q)))==0,('recurrence',r,q,seq))

def dp(n,rmax,weights=None,initial=None):
 ds=divisors(n)
 if weights is None:weights=[0]+[1]*(n-1)
 v={d:F(0) for d in ds};v[1]=1
 if initial is not None:v={d:F(initial.get(d,0)) for d in ds}
 out=[]
 # Aggregate by gcd first, preserving total source mass.
 gw={d:sum(weights[x] for x in range(n) if gcd(x,n)==d) for d in ds}
 for r in range(rmax+1):
  out.append(v[n]);nv={d:F(0) for d in ds}
  for d in ds:
   for e in ds:nv[gcd(d*e,n)]+=v[d]*gw[e]
  v=nv
 return out

def rank(rows):
 a=[list(map(F,row)) for row in rows];r=0
 for col in range(len(a[0])):
  j=next((j for j in range(r,len(a)) if a[j][col]),None)
  if j is None:continue
  a[r],a[j]=a[j],a[r];z=a[r][col];a[r]=[x/z for x in a[r]]
  for j in range(r+1,len(a)):
   z=a[j][col]
   if z:a[j]=[x-z*y for x,y in zip(a[j],a[r])]
  r+=1
  if r==len(a):break
 return r

# Literal tuples, not valuation-state calculations.
entries=0
for n in range(2,14):
 for r in range(1,5):
  z=s=0
  for xs in product(range(1,n),repeat=r):
   x=prod(xs);z+=x%n==0;s+=x//n;entries+=1
  require(z==formula(n,r),('literal-zero',n,r))
  rhs=F((n*(n-1)//2)**r,n)-F((n-1)**r,2)+F(z,2)
  require(s==rhs,('literal-floor',n,r))

# Full divisor automata, all composite and prime moduli through 100.
for n in range(2,101):
 seq=dp(n,12)
 for r,z in enumerate(seq):require(z==formula(n,r),('dp',n,r))
 fs=factors(n)
 if len(fs)==1 and fs[0][1]==1:
  require(not any(seq),('prime-zero',n));continue
 q=recurrence(n)
 d=(1<<len(fs))+(1<<(len(fs)-1))*sum(a-1 for p,a in fs)
 require(len(q)-1==d,('order',n))
 gap=prod(a+1 for p,a in fs)-d
 expected=sum((1<<(len(fs)-mask.bit_count()))*prod(a-1 for i,(p,a) in enumerate(fs) if mask>>i&1) for mask in range(1<<len(fs)) if mask.bit_count()>=2)
 require(gap==expected and (gap==0)==(sum(a>1 for p,a in fs)<=1),('dimension-gap',n))
 seq=[formula(n,r) for r in range(2*d+4)]
 apply_recurrence(q,seq)
 floors=[F((n*(n-1)//2)**r,n)-F((n-1)**r,2)+F(formula(n,r),2) for r in range(1,2*d+5)]
 apply_recurrence(recurrence(n,floor=True),floors)
 ds=divisors(n)
 for u,v in combinations(ds,2):
  require(any((gcd(u*x,n)==n)!=(gcd(v*x,n)==n) for x in ds if x<n),('distinguish',n,u,v))

# Hankel rank independently detects exact scalar minimal order.
hankel=[]
for n in [4,6,8,9,10,12,18,24,30,36,60,72]:
 q=recurrence(n);d=len(q)-1;z=dp(n,2*d+2)
 rk=rank([[z[i+j] for j in range(d+1)] for i in range(d+1)])
 require(rk==d,('Hankel',n,rk,d));hankel.append([n,d])

# Arbitrary positive, signed, and rational source weights, including residue zero.
for n in range(2,31):
 banks=[[F((x*x+3*x+1)%7) for x in range(n)],
        [F((x*5+2)%9-4,3) for x in range(n)],
        [F(int(x==0)) for x in range(n)],
        [F(int(x==1)) for x in range(n)]]
 for w in banks:
  q=recurrence(n,w);z=dp(n,2*len(q)+2,w)
  apply_recurrence(q,z)
  # Nontrivial initial source supported on divisor states.
  initial={d:F((d*3)%5-2,2) for d in divisors(n)}
  apply_recurrence(q,dp(n,2*len(q)+2,w,initial))

# Classical Moebius inversion of the full one-step probe bank.
def moebius(n):
 fs=factors(n)
 return 0 if any(a>1 for p,a in fs) else (-1)**len(fs)
for n in range(2,101):
 ds=divisors(n)
 for mode in range(2):
  masses={g:F((3*g+2)%11-(5 if mode else 0),7) for g in ds}
  probes={1:sum(masses.values())}
  for d in ds:
   if d>1:
    probes[d]=sum(m for g,m in masses.items() if gcd(g*(n//d),n)==n)
    require(probes[d]==sum(m for g,m in masses.items() if g%d==0),('probe-zeta',n,d))
  for g in ds:
   recovered=sum(moebius(d//g)*probes[d] for d in ds if d%g==0)
   require(recovered==masses[g],('probe-inverse',n,g,mode))

# All-order hostile verified by an exact left-eigenvector equation.
n=36;ds=divisors(n);h=[0,-3,2,3,1,-2,-3,2,0]
require(ds==[1,2,3,4,6,9,12,18,36],'state-order')
M=[[sum(gcd(d*x,n)==e for x in range(1,n)) for e in ds] for d in ds]
require([sum(h[i]*M[i][j] for i in range(len(ds))) for j in range(len(ds))]==[12*x for x in h],'hostile-eigenvector')
require(sum(h)==0 and h[-1]==0,'hostile-mass-and-output')
require(sum(h[i] for i,d in enumerate(ds) if gcd(2*d,n)==n)==2,'named-multiplier')
require(sum(max(x,0) for x in h)==8,'probability-normalization')

print(json.dumps({'status':'PASS','checks':checks,'literal_tuple_entries':entries,'hankel_orders':hankel,'mod36_hostile_eigenvalue':12,'mod36_named_probability_gap':'1/4'},sort_keys=True))
