#!/usr/bin/env python3
"""Independent finite-conditioning audit; checks survive python -O."""
from fractions import Fraction as F
from itertools import product
from math import prod, gcd, isqrt

def require(ok, message):
    if not ok: raise RuntimeError(message)

def kernel(p,k,v):
    alpha=1-F(1,p)
    A=B=F(0)
    for t in range(k+1):
        weight=alpha/F(p**t) if t<k else F(1,p**k)
        a=1+t*alpha+F(t*(t+1),2)*alpha**2 if t<=v else F((v+1)*(v+2),2)*alpha**2
        b=1+t*alpha**2 if t<=v else (v+1)*alpha**2+F(1,p)
        A+=weight*a; B+=weight*b
    return A,B,1+F(1,p**3)-(F(1,p*p) if v>0 else 0)

def mass(factors,v):
    branches=[prod(kernel(p,k,t)[j]*F((p-1)*p**(k-t-1) if t<k else 1,p**k) for (p,k),t in zip(factors,v)) for j in range(3)]
    return sum(w*z for w,z in zip([F(9,16),F(6,16),F(1,16)],branches))

steps=0
for p in [3,5,7,11,17,31,101]:
    for k in range(1,11):
        vals=[kernel(p,k,v) for v in range(k+1)]
        require(all(x>0 for row in vals for x in row),'positivity')
        for v in range(k):
            A,B,C=vals[v]; a,b,c=vals[v+1]
            require(a*B>A*b and b*C>B*c,'strict likelihood ratio')
            alpha=1-F(1,p); d=1-F(1,p*p)
            require((a-A)/(b-B)==(v+2+F(p,p*p-1) if v<k-1 else k+1/alpha),'increment quotient')
            steps+=1
        for j in range(3):
            require(sum(F((p-1)*p**(k-v-1) if v<k else 1,p**k)*row[j] for v,row in enumerate(vals))==1,'normalization')
print('Independent finite-U sums: positivity, normalization, MLR and endpoint algebra:',steps,'steps')

# Independent direct integer counts, retaining four shared labels.
q=15
counts=[[0]*q for _ in range(4)]
for d,u,D in product(range(q),repeat=3):
    g=gcd(gcd(d,u),q); h=gcd(gcd(D,u),q)
    for a,b,r in product(range(2),range(2),range(q)):
        if (D+a)%g==0 and (r+b*d)%h==0: counts[2*a+b][r]+=g*h
for r in range(q):
    expected=[prod(kernel(p,1,int(r%p==0))[j] for p in (3,5)) for j in (0,1,1,2)]
    require([F(c[r],q**3) for c in counts]==expected,'branchwise CRT')
print('Independent integer enumeration matches branchwise CRT at every residue modulo 15')

factors=[(3,2),(5,2),(7,1)]
U=list(product(*(range(k+1) for p,k in factors)))
m={v:mass(factors,v) for v in U}
require(sum(m.values())==1,'joint normalization')
for x,y in product(U,repeat=2):
    low=tuple(min(a,b) for a,b in zip(x,y)); high=tuple(max(a,b) for a,b in zip(x,y))
    delta=m[low]*m[high]-m[x]*m[y]
    require(delta>0 if low not in (x,y) else delta==0,'strict MTP2')
    if low==x and x!=y:
        def density(v):
            multiplicity=prod((p-1)*p**(k-t-1) if t<k else 1 for (p,k),t in zip(factors,v))
            return m[v]/multiplicity
        require(density(x)<density(y),'density monotonicity')
# Every single-coordinate conditioning and each threshold pair.
cases=0
for fixed in range(3):
    rem=[i for i in range(3) if i!=fixed]
    for value in range(factors[fixed][1]+1):
        selected={v:z for v,z in m.items() if v[fixed]==value}; total=sum(selected.values())
        for t in range(1,factors[rem[0]][1]+1):
            for s in range(1,factors[rem[1]][1]+1):
                a=sum(z for v,z in selected.items() if v[rem[0]]>=t)
                b=sum(z for v,z in selected.items() if v[rem[1]]>=s)
                both=sum(z for v,z in selected.items() if v[rem[0]]>=t and v[rem[1]]>=s)
                require(both*total>a*b,'strict conditional association')
                cases+=1
print('Independent strict MTP2:',len(U)**2,'pairs; strict conditional threshold covariances:',cases)

# Independently rebuild rigorous product interval using a set sieve and decimals
# represented only as integers; no producer imports and no floating arithmetic.
N=1_000_000
composite=bytearray(N+1)
for p in range(2,isqrt(N)+1):
    if not composite[p]:
        for n in range(p*p,N+1,p): composite[n]=1
S=10**25
lower=[S]*3; upper=[S]*3
for p in range(3,N+1,2):
    if composite[p]: continue
    for j,num in enumerate((p**3-2*p+1,p**3-p+1,p**3+1)):
        den=p**3
        lower[j]=lower[j]*num//den
        upper[j]=-((-upper[j]*num)//den)
lo=F(9*lower[0]+6*lower[1]+lower[2],16*S)-F(3,2*N)
hi=F(9*upper[0]+6*upper[1]+upper[2],16*S)
require(F(771866451842,10**12)<lo<hi<F(771867951843,10**12),'published rigorous enclosure')
print('Published interval independently certified: 0.771866451842 < Csharp < 0.771867951843')
print('PASS (all checks explicit, including under -O)')
