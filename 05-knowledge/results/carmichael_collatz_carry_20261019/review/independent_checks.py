#!/usr/bin/env python3
"""Independent bounded review of exact inheritance propositions."""
import itertools,json,math
from fractions import Fraction
from collections import Counter
checks=0
def require(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)

def factor(n):
 fs=[]
 for p in range(2,math.isqrt(n)+1):
  if n%p==0:
   e=0
   while n%p==0:n//=p;e+=1
   fs.append((p,e))
 if n>1:fs.append((n,1))
 return fs

def carmichael(n):
 fs=factor(n)
 return len(fs)>=2 and all(e==1 and (n-1)%(p-1)==0 for p,e in fs)

def mu(n):
 fs=factor(n)
 return 0 if any(e>1 for p,e in fs) else (-1)**len(fs)
def phi(n):
 return sum(math.gcd(i,n)==1 for i in range(n))

def apply(word,x):
 for b in word:x=(3*x+1)/2 if b else x/2
 return x

def carry(word):
 p=sum(word)
 return sum((2**j)*3**sum(word[j+1:]) for j,b in enumerate(word) if b)

def compose(f,g,m):return(f[0]*g[0]%m,(f[0]*g[1]+f[1])%m)
def inv(f,m):
 a=pow(f[0],-1,m);return(a,-a*f[1]%m)

kernel_n=0;core_pairs=0
for n in range(2,801):
 fs=factor(n);units=[a for a in range(1,n) if math.gcd(a,n)==1]
 image={pow(a,n-1,n) for a in units}
 expect=math.prod(p**(e-1)*(p-1)//math.gcd(n-1,p-1) for p,e in fs)
 require(len(image)==expect,('Fermat image index',n))
 require(len(units)==expect*sum(pow(a,n-1,n)==1 for a in units),('first isomorphism',n));kernel_n+=1
primes=[n for n in range(2,251) if factor(n)==[(n,1)]]
for n in range(2,251):
 fs=factor(n)
 if any(e!=1 for p,e in fs):continue
 lam=math.lcm(*(p-1 for p,e in fs))
 for p in primes:
  if p>n or n%p==0:continue
  criterion=(n-1)%(p-1)==0 and (n*p-1)%lam==0
  require(carmichael(n*p)==criterion,('core completion',n,p));core_pairs+=1
  if math.gcd(n,lam)>1:require(not criterion,('hereditary gcd',n,p))

moduli=[]
for m in range(5,61):
 if math.gcd(m,6)!=1:continue
 h=pow(2,-1,m);f0=(h,0);f1=(3*h%m,h)
 comm=compose(compose(f0,f1,m),compose(inv(f0,m),inv(f1,m),m),m)
 require(comm==(1,-pow(4,-1,m)%m),('commutator',m))
 require({comm[1]*j%m for j in range(m)}==set(range(m)),('translation generator',m))
 moduli.append(m)

Z={};N={};P={};words=0;native=0
for q in range(2,13):
 zero=Counter();orb={};pr={}
 for w in itertools.product((0,1),repeat=q):
  p=sum(w)
  if p in(0,q):continue
  words+=1
  b=apply(w,Fraction(0));a=apply(w,Fraction(1))-b;x=b/(1-a)
  B=carry(w);D=2**q-3**p
  require(b*2**q==B,('direct affine carry',w))
  require(x==Fraction(B,D),('fixed point',w))
  rot=w[1:]+w[:1]
  require(2*carry(rot)==3**w[0]*B+w[0]*D,('rotation',w))
  if x.denominator==1:
   zero[p]+=1;key=min(w[j:]+w[:j] for j in range(q));orb.setdefault(p,set()).add(key)
   period=next(k for k in range(1,q+1) if q%k==0 and w==w[:k]*(q//k))
   if period==q:pr.setdefault(p,set()).add(key)
   y=x
   for s in w:
    require(y.denominator==1 and y.numerator%2==s,('native branch',w,y))
    y=(3*y+1)/2 if s else y/2
   require(y==x,('native cycle',w));native+=1
 for p in range(1,q):
  Z[q,p]=zero[p];N[q,p]=len(orb.get(p,set()));P[q,p]=len(pr.get(p,set()))
  ds=[d for d in range(1,math.gcd(q,p)+1) if q%d==0 and p%d==0]
  require(q*N[q,p]==sum(phi(d)*Z[q//d,p//d] for d in ds),('Burnside',q,p))
  require(q*P[q,p]==sum(mu(d)*Z[q//d,p//d] for d in ds),('primitive Mobius',q,p))

# Verify the additive character recursion exactly in a finite-field root-of-unity model.
fourier_cases=0
for m,ell in [(5,11),(7,29),(11,23),(13,53)]:
 z=next(z for z in range(2,ell) if pow(z,m,ell)==1 and all(pow(z,d,ell)!=1 for d in range(1,m)))
 for a in range(m):
  for q in range(2,9):
   C={(0,j):1 for j in range(m)}
   for j in range(q):
    out={}
    for k in range(j+2):
     for freq in range(m):out[k,freq]=(C.get((k,freq),0)+pow(z,freq*2**j,ell)*C.get((k-1,3*freq%m),0))%ell
    C=out
   for p in range(1,q):
    direct=sum(pow(z,a*carry(w),ell) for w in itertools.product((0,1),repeat=q) if sum(w)==p)%ell
    require(C[p,a]==direct,('Fourier carry recurrence',m,ell,q,p,a));fourier_cases+=1

occupancies=0
for m in range(2,7):
 for cs in itertools.product(range(3),repeat=m):
  W=sum(cs)
  if not W:continue
  a,r=divmod(W,m);E=m*sum(x*x for x in cs)-W*W
  require(E>=r*(m-r),('occupancy lower bound',cs))
  require((E==r*(m-r))==all(x in(a,a+1) for x in cs),('occupancy equality',cs))
  if W<m:
   tv=sum(abs(Fraction(c,W)-Fraction(1,m)) for c in cs)/2
   require(tv>=1-Fraction(W,m),('sparse TV',cs))
  occupancies+=1
print(json.dumps({'status':'PASS','explicit_checks':checks,'kernel_integers':kernel_n,'core_prime_pairs':core_pairs,'commutator_moduli':moduli,'mixed_words':words,'native_integral_words':native,'Fourier_recurrence_cases':fourier_cases,'occupancy_vectors':occupancies,'primitive_mixed_shapes':[list(k) for k,v in P.items() if v],'repeated_necklace_hostile':{'q':4,'p':2,'zero_words':Z[4,2],'orbits':N[4,2],'primitive_orbits':P[4,2]},'scope':'finite exact controls; fixed-m asymptotic reviewed as conditional proof, not inferred computationally'},indent=2))
