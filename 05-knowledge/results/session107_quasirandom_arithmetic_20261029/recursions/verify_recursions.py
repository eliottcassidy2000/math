#!/usr/bin/env python3
"""Standalone exact arithmetic audit; explicit checks survive python -O."""
import json
from math import gcd, isqrt, prod
from fractions import Fraction

def require(ok, label):
    if not ok:
        raise RuntimeError('FAILED: '+label)

def prime(n):
    return n >= 2 and all(n % d for d in range(2,isqrt(n)+1))

def f(x): return x*x-x+1

def orbit(m):
    x=2 % m; seen={}; a=[]
    while x not in seen:
        seen[x]=len(a); a.append(x); x=f(x)%m
    return a,seen[x]

def iterate(x,n,m):
    for _ in range(n): x=f(x)%m
    return x

def valuation(x,p):
    if not x: raise ValueError('zero valuation unspecified')
    v=0
    while x%p==0: x//=p;v+=1
    return v

s=[2]
for _ in range(6):s.append(f(s[-1]))
require(s[:6]==[2,3,7,43,1807,3263443],'Sylvester terms')
require(1807==13*139 and all(prime(p) for p in [13,139,3263443]),'factorization and next prime')
for n in range(1,len(s)):
    require(s[n]-1==prod(s[:n]),'product recursion')
    require(sum((Fraction(1,x) for x in s[:n]),Fraction())==1-Fraction(1,s[n]-1),'Egyptian identity')
    for k in range(n): require(gcd(s[n],s[k])==1,'pairwise coprime')
require(orbit(19)==([2,3,7,5],0),'requested cycle identified modulo19')
require(orbit(13)==([2,3,7,4,0,1],5),'13 absorption')
require(orbit(139)==([2,3,7,43,0,1],5),'139 absorption')
require(1805==5*19**2,'prime-square return defect')
periods=[]
for k in range(1,7):
    m=19**k; a,j=orbit(m); actual=len(a)-j
    expected=4 if k<=2 else 12*19**(k-3)
    require(j==0 and actual==expected,'period law level '+str(k))
    periods.append(actual)
base_cycle=set(orbit(19**3)[0])
for k in range(3,7):
    cycle=orbit(19**k)[0]
    counts={c:0 for c in base_cycle}
    for x in cycle:
        require(x%(19**3) in counts,'allowed low-digit support')
        counts[x%(19**3)]+=1
    require(set(counts.values())=={19**(k-3)},'exact full-fiber uniformity')
# Exact Fourier nonvanishing at every low-digit frequency, no floating point.
# Reduce sum X^e modulo Phi_(19^t)=sum_{j=0}^{18} X^(j*19^(t-1)).
for b in range(19**3):
    common=gcd(b,19**3); q=19**3//common
    if q==1:
        require(len(base_cycle)==12,'trivial Fourier frequency')
        continue
    unit=b//common; block=q//19; remainder={}
    for c in base_cycle:
        e=(unit*c)%q
        if e<18*block:
            remainder[e]=remainder.get(e,0)+1
        else:
            r=e-18*block
            for j in range(18):
                ee=r+j*block; remainder[ee]=remainder.get(ee,0)-1
    require(any(remainder.values()),'nonzero cyclotomic remainder frequency '+str(b))
x=2; derivative=1
for _ in range(4): derivative*=2*x-1; x=f(x)
require(derivative==16575 and derivative%19==7,'multiplier')
require(pow(7,3,19)==1 and 7!=1 and pow(7,2,19)!=1,'multiplier order3')
require(iterate(2,12,19**4)-2==17*19**3,'nonzero h defect')
y=2; dh=1
for _ in range(12):
    dh=dh*(2*y-1)%19; y=f(y)%19
require(dh==1,'h tangent to identity')
require([x for x in range(19) if f(x)%19==0]==[8,12],'unreachable roots mod19')
require(valuation(iterate(2,12,19**4)-2,19)==3,'exact valuation')
require(valuation(pow(2,18,19**3)-1,19)==1,'19 not base2 Wieferich')
D=[d for d in range(1,1807) if 1806%d==0]
selected=[d+1 for d in D if prime(d+1)]
require(selected==[2,3,7,43],'ALL divisor selector audit')
require(1806%12!=0 and 1806%138!=0,'selector excludes 13 and139')
require(prime(19) and 19%6==1 and 0 not in orbit(19)[0], 'splitting insufficient hostile')
fermat=[2**(2**n)+1 for n in range(6)]
for n in range(1,len(fermat)):
    require(fermat[n]==(fermat[n-1]-1)**2+1,'Fermat quadratic recursion')
    require(fermat[n]-2==prod(fermat[:n]),'Fermat product recursion')
for n,p in enumerate(fermat[:5]):
    require(prime(p),'Fermat primality')
    d=2**(n+1)
    require(2**d-1==p*(p-2),'Fermat non-Wieferich factorization')
require(fermat[5]==641*6700417,'Fermat composite hostile')
# Local necessary constraint, deliberately not an iff.
for x in range(2,100):
    z=f(x)
    for p in range(2,isqrt(z)+1):
        if z%p==0 and prime(p):
            require(p==3 or p%6==1,'Phi6 prime divisor constraint')
report={'status':'PASS','Sylvester':s,'selector':selected,'cycle19':orbit(19),'periods19_k1_to6':periods,'low_digit_Fourier_frequencies_exactly_nonzero':19**3,'g_derivative':derivative,'h_defect_mod19power4':17*19**3,'Fermat':fermat,'hostiles':['19 admissible mod6 but never divisor','19 dynamical square delay but not base2 Wieferich','1807 composite followed by3263443 prime','all divisors needed for selector fixed point']}
print(json.dumps(report,indent=2))
