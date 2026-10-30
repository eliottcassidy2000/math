"""Exact weighted-cycle determinants, modular images, and Collatz word tests."""
from itertools import product
from math import gcd
from fractions import Fraction
import json

def require(ok, msg):
    if not ok: raise RuntimeError(msg)
def det(A):
    n=len(A)
    if n==0:return 1
    a=[list(map(Fraction,r)) for r in A];d=Fraction(1)
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[j],a[pivot]=a[pivot],a[j];d=-d
        v=a[j][j];d*=v
        for k in range(j,n):a[j][k]/=v
        for i in range(j+1,n):
            v=a[i][j]
            for k in range(j,n):a[i][k]-=v*a[j][k]
    require(d.denominator==1,'integer determinant')
    return d.numerator

def matrix(p,u,a):
    k=len(a);M=[[0]*k for _ in a]
    for i,v in enumerate(a):M[i][i]-=u;M[i][(i+1)%k]+=p**v
    return M

def numerators(p,u,a,b):
    k=len(a);out=[]
    for i in range(k):
        s=N=0
        for j in range(k):
            N+=u**(k-1-j)*p**s*b[(i+j)%k];s+=a[(i+j)%k]
        out.append(N)
    return out

def period(v):
    return next(d for d in range(1,len(v)+1) if len(v)%d==0 and all(v[i]==v[i%d] for i in range(len(v))))

matrices=minor_checks=0
for p in range(1,6):
 for u in range(-5,6):
  if gcd(p,u)!=1:continue
  for k in range(1,6):
   for a in product((1,2),repeat=k):
    M=matrix(p,u,a);D=p**sum(a)-u**k
    if D==0:require(p==1 and (u==1 or (u==-1 and k%2==0)),'complete singular classification')
    require(det(M)==(-1)**(k-1)*D,('det',p,u,a))
    if k>1:
     first=det([r[:k-1] for r in M[:k-1]])
     shifted=det([r[1:] for r in M[:k-1]])
     require(first==(-u)**(k-1) and shifted==p**sum(a[:-1]),'displayed minors')
     require(gcd(first,shifted)==1,'Smith invariant factors before last')
     minor_checks+=1
    matrices+=1

identities=integral=0
for p,u in [(2,3),(3,2),(2,-3),(1,1),(1,-1)]:
 for k in range(1,5):
  for a in product((1,2),repeat=k):
   D=p**sum(a)-u**k
   for b in product((-1,0,1),repeat=k):
    Ns=numerators(p,u,a,b)
    require(all(p**a[i]*Ns[(i+1)%k]-u*Ns[i]==D*b[i] for i in range(k)),'adjugate recurrence')
    if D:
     require((Ns[0]%D==0)==all(N%D==0 for N in Ns),'one congruence sufficient')
     integral+=int(Ns[0]%D==0)
    else:
     require((Ns[0]==0)==all(N==0 for N in Ns),'singular compatibility')
     for initial in (-2,0,3):
      value=initial
      for rhs in b:value=u*value+rhs
      require((value==initial)==(Ns[0]==0),'singular free parameter')
     integral+=int(Ns[0]==0)
    identities+=1

modular_b=modular_systems=0
for p,u in ((2,3),(1,1),(1,-1)):
 for k in range(1,4):
  for a in product((1,2,3),repeat=k):
   D=p**sum(a)-u**k;M=matrix(p,u,a)
   for m in range(2,7):
    image={}
    for x in product(range(m),repeat=k):
     b=tuple(sum(M[i][j]*x[j] for j in range(k))%m for i in range(k))
     image[b]=image.get(b,0)+1
    g=gcd(D,m)
    require(image.get((0,)*k)==g,'modular kernel count')
    for b in product(range(m),repeat=k):
     N=numerators(p,u,a,b)[0]
     require((N%g==0)==(b in image),'modular obstruction')
     if b in image:require(image[b]==g,'modular fiber count')
     modular_b+=1
    modular_systems+=1

word_count=integer_words=positive_words=negative_words=0;cycles=set()
for k in range(1,8):
 for a in product(range(1,5),repeat=k):
  D=2**sum(a)-3**k;N=numerators(2,3,a,(1,)*k)[0]
  require(D%2 and N%2 and gcd(D,3)==1 and gcd(N,3)==1,'odd/unit arithmetic')
  word_count+=1
  if N%D:continue
  xs=tuple(x//D for x in numerators(2,3,a,(1,)*k));integer_words+=1
  require(all(x%2 for x in xs),'odd quotient')
  require(all((3*xs[i]+1)==2**a[i]*xs[(i+1)%k] for i in range(k)),'exact word')
  require(all((x>0)==(D>0) for x in xs),'sign criterion')
  require(period(xs)==period(a),'primitive word iff primitive orbit')
  positive_words+=D>0;negative_words+=D<0
  d=period(xs);v=xs[:d];cycles.add(min(v[i:]+v[:i] for i in range(d)))

w=(1,1,1,2,1,1,4);v=(1,1,1,1,1,2,4)
require(sorted(w)==sorted(v),'same multiset hostile')
D=2**sum(w)-3**len(w)
Nw=numerators(2,3,w,(1,)*len(w))[0];Nv=numerators(2,3,v,(1,)*len(v))[0]
require((D,Nw,Nv)==(-139,2363,2123),'ordered carry hostile')
require(Nw%D==0 and Nv%D!=0,'integrality changes under noncyclic permutation')
require((-7,-5) in cycles and (1,) in cycles and (-1,) in cycles,'cycle controls')
# Coprimality is essential to the cyclic Smith conclusion.
require(matrix(2,2,(1,1))==[[-2,2],[2,-2]],'exceptional singular control')
require(gcd(*[abs(z) for r in matrix(2,2,(1,1)) for z in r])==2,'nontrivial first invariant')
print(json.dumps({'status':'PASS','matrix_cases':matrices,'coprime_minor_cases':minor_checks,'numerator_identity_cases':identities,'integral_rhs_cases':integral,'modular_systems':modular_systems,'modular_rhs_cases':modular_b,'collatz_words':word_count,'integral_collatz_words':integer_words,'positive_words':positive_words,'negative_words':negative_words,'primitive_integer_cycles_in_bounded_scan':[list(v) for v in sorted(cycles)],'ordered_multiset_hostile':{'words':[w,v],'D':D,'numerators':[Nw,Nv]},'scope':'exact bounded checks; no positive Collatz cycle exclusion beyond this word box'},sort_keys=True,indent=2))
