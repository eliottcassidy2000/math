#!/usr/bin/env python3
"""Portable exact rational checks; Python standard library only."""
from fractions import Fraction as Q
from itertools import product
from math import gcd, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def kernels(p,k,v):
    a=1-Q(1,p); d=1-Q(1,p*p); b=1+Q(1,p)-Q(1,p*p); c=1+Q(2,p)
    if v<k:
        z=Q(1,p**(v+1)); A=c-z*(c+(v+1)*d); B=b-z*d
    else:
        A=c-Q(2+k*a,p**(k+1)); B=b-Q(a,p**(k+1))
    C=1+Q(1,p**3)-(Q(1,p*p) if v else 0)
    return A,B,C


def valuation(r,p,k):
    v=0
    while v<k and r%p==0:
        v+=1; r//=p
    return v


def direct(q):
    out=[[Q(0) for _ in range(q)] for _ in range(4)]
    for d,u,D in product(range(q),repeat=3):
        g=gcd(gcd(d,u),q); h=gcd(gcd(D,u),q)
        for a,b in product(range(2),repeat=2):
            if (D+a)%g: continue
            for r in range((-b*d)%h,q,h):
                out[2*a+b][r]+=Q(g*h,q**3)
    return out


def main():
    primes=[3,5,7,11,13,17,19,23,29,31,37,41,43]
    n=0
    for p in primes:
        for k in range(1,7):
            vals=[kernels(p,k,v) for v in range(k+1)]
            for x,y in zip(vals,vals[1:]):
                A,B,C=x; Ap,Bp,Cp=y
                require(Ap*B>A*Bp and Bp*C>B*Cp, f"Likelihood ratio failure: {p=}, {k=}, {x=}, {y=}")
                require(Ap>A and Bp>B, f"Kernel monotonicity failure: {p=}, {k=}, {x=}, {y=}")
                n+=1
            for j in range(3):
                require(sum(Q((p-1)*p**(k-v-1) if v<k else 1,p**k)*vals[v][j]
                            for v in range(k+1))==1, f"Normalization failure: {p=}, {k=}, {j=}")
    print('Exact prime-power likelihood-ratio adjacent checks:',n)
    for p,k in [(3,1),(3,2),(3,3),(5,1),(5,2),(7,1)]:
        q=p**k; out=direct(q)
        for r in range(q):
            A,B,C=kernels(p,k,valuation(r,p,k))
            require([out[j][r] for j in range(4)]==[A,B,B,C], f"Direct kernel mismatch: {q=}, {r=}")
        print('Independent direct branch enumeration passed:',q)
    for factors in [[(3,1),(5,1)],[(3,2),(5,1)],[(3,2),(5,2),(7,1)]]:
        universe=list(product(*(range(k+1) for p,k in factors)))
        densities={}
        masses={}
        for v in universe:
            branch=[prod(kernels(p,k,t)[j] for (p,k),t in zip(factors,v)) for j in range(3)]
            L=(9*branch[0]+6*branch[1]+branch[2])/16
            density=prod(Q((p-1)*p**(k-t-1) if t<k else 1,p**k) for (p,k),t in zip(factors,v))
            densities[v]=L; masses[v]=L*density
        require(sum(masses.values())==1, f"Joint normalization failure: {factors=}")
        for x,y in product(universe,repeat=2):
            meet=tuple(map(min,x,y)); join=tuple(map(max,x,y))
            delta=masses[meet]*masses[join]-masses[x]*masses[y]
            require(delta>=0, f"MTP2 failure: {factors=}, {x=}, {y=}")
            if meet!=x and meet!=y:
                require(delta>0, f"Strict MTP2 failure: {factors=}, {x=}, {y=}")
            if x!=y and meet==x:
                require(densities[x]<densities[y], f"Density monotonicity failure: {factors=}, {x=}, {y=}")
        print('Exact strict MTP2 and valuation monotonicity:',factors,len(universe)**2,'pairs')
    out=direct(15)
    joint=(9*out[0][0]+3*out[1][0]+3*out[2][0]+out[3][0])/16
    marg=[]
    for p in [3,5]:
        A,B,C=kernels(p,1,1); marg.append((9*A+6*B+C)/16)
    require(joint-prod(marg)==Q(1,75), "Hostile CRT witness mismatch")
    print('Hostile product-of-mixtures gap at q=15:',joint-prod(marg))
    # Descriptive approximate value with a rigorous analytic truncation bound.
    N=1_000_000
    sieve=bytearray(b'\x01')*(N+1); sieve[:2]=b'\x00\x00'
    for p in range(2,int(N**.5)+1):
        if sieve[p]: sieve[p*p:N+1:p]=b'\x00'*((N-p*p)//p+1)
    A=B=C=1.0
    scale=10**30
    lo=[scale]*3; hi=[scale]*3
    for p in range(3,N+1,2):
        if sieve[p]:
            A*=1-2/p**2+1/p**3; B*=1-1/p**2+1/p**3; C*=1+1/p**3
            den=p**3
            for j,num in enumerate([den-2*p+1,den-p+1,den+1]):
                lo[j]=lo[j]*num//den
                hi[j]=(hi[j]*num+den-1)//den
    estimate=(9*A+6*B+C)/16
    print('Descriptive sharp-floor partial product:',format(estimate,'.12f'))
    print('Analytic truncation error bound:',Q(3,2*N),'(floating rounding separate)')
    lower=Q(9*lo[0]+6*lo[1]+lo[2],16*scale)-Q(3,2*N)
    upper=Q(9*hi[0]+6*hi[1]+hi[2],16*scale)
    digits=10**12
    lower_dec=(lower.numerator*digits)//lower.denominator
    upper_dec=-((-upper.numerator*digits)//upper.denominator)
    print('Certified outward interval:',f'{lower_dec//digits}.{lower_dec%digits:012d}',
          f'{upper_dec//digits}.{upper_dec%digits:012d}')
    require(Q(771866,10**6)<lower<upper<Q(771869,10**6), "Certified enclosure failure")
    print('PASS')

if __name__=='__main__': main()
