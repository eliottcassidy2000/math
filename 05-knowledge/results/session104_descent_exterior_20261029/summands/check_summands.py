#!/usr/bin/env python3
"""Exact finite controls; all guards stay active under python -O."""
from fractions import Fraction
from math import gcd, isqrt
import json
from pathlib import Path

require_calls = 0
def require(test, message):
    global require_calls
    require_calls += 1
    if not test:
        raise RuntimeError(message)

def positive_int(x, name):
    if type(x) is not int or x < 1:
        raise ValueError(name + ' must be a positive integer')

def defect(a, k, p, d=1):
    positive_int(k, 'k')
    positive_int(d, 'd')
    if type(p) is not int or p not in (0, 1, 2):
        raise ValueError('unsupported power')
    if type(a) not in (int, Fraction):
        raise ValueError('exact base required')
    return sum((a+d*i)**p for i in range(k+1))-sum((a+d*i)**p for i in range(k+1,2*k+1))

def rad(n):
    positive_int(n, 'n')
    r, p = 1, 2
    while p*p <= n:
        if n % p == 0:
            r *= p
            while n % p == 0:
                n //= p
        p += 1
    return r*n

checks = 0
for k in range(1,101):
    for d in (1,2,3,7):
        for a in (-k*d, 0, 1, d*k*k, d*k*(2*k+1), d*k*(2*k+1)-1, d*k*(2*k+1)+1):
            require(defect(a,k,1,d)==a-d*k*k, 'linear polynomial')
            require(defect(a,k,2,d)==(a+d*k)*(a-d*k*(2*k+1)), 'square polynomial')
            checks += 2
    for a in (0, k*k, k*(2*k+1)):
        require(defect(a+1,k,2)-defect(a,k,2)==2*defect(a,k,1)+1, 'forward difference')
        require(defect(Fraction(2*a+1,2),k,2)-defect(Fraction(2*a-1,2),k,2)==2*defect(a,k,1), 'centered difference')
        checks += 2
    require(defect(2*k*k,k,1,2)==0, 'even identity')
    require(defect(k*(2*k+1),k,2)==0, 'square identity')
    require(defect(k*(2*k+1)-1,k,2)!=0, 'wrong start hostile')
    checks += 3

pairs_total=0
for n in range(1,401):
    pairs={(a,b) for a in range(1,n) for b in (n-a,) if a<b}
    require(len(pairs)==(n-1)//2,'pair fibre count')
    arcs={a for pair in pairs for a in pair}
    require(arcs=={x for x in range(1,n) if n!=2*x},'strict shadow')
    require(len(arcs)==2*len(pairs),'arc versus pair count')
    weak={(a,n-a) for a in range(1,n) if a<=n-a}
    require(len(weak)-len(pairs)==int(n%2==0),'additive diagonal')
    pairs_total+=len(pairs)
    require(pairs_total==(n-1)**2//4,'cumulative pair count')
    if n>=2:
        mw={(a,n//a) for a in range(1,n+1) if n%a==0 and a<=n//a}
        ms={(a,b) for a,b in mw if a<b}
        require(len(mw)-len(ms)==int(isqrt(n)**2==n),'square diagonal')
    checks += 6

live={2,3}
for n in range(1,1001):
    if any(a in live and n-a in live for a in range(1,(n+1)//2)):
        live.add(n)
require(live==set(range(1,1001))-{1,4,6},'strict closure')
require(3**2+4**2==5**2 and defect(3,1,2)==0,'Pythagorean startup')
require(rad(9*16*25)==30 and rad(1*2*3)==6,'abc normalized examples')
for n in range(1,101):
    require(rad(n*n)==rad(n),'square prime support')
    require(rad(2*n)==rad(n)*(1 if n%2==0 else 2),'double prime support')
    checks += 2
A,B,C=3,4,5
legs=[]
for j in range(7):
    require(A*A+B*B==C*C,'Gaussian Pythagorean identity')
    require(gcd(abs(A),abs(B))==1 and A%2==1 and B%2==0,'Gaussian primitivity')
    require(C==5**(2**j),'hypotenuse power')
    require(all(gcd(abs(A),abs(old))==1 for old in legs),'fresh odd-leg support')
    require(gcd(abs(A),5)==1 and abs(A)>1,'leg prime absent from hypotenuse')
    legs.append(A)
    A,B,C=A*A-B*B,2*A*B,C*C
    checks += 5
hostiles=0
for args in [(1,0,1),(1,True,1),(1,1,3),(1.0,1,1),(1,1,1,0),(1,1,1,True)]:
    try:
        defect(*args)
    except ValueError:
        hostiles+=1
    else:
        raise RuntimeError('illegal input accepted')
require(hostiles==6,'hostile guards')
# Deliberate false predicate demonstrates the check mechanism itself remains live.
try:
    require(False,'guard self-test')
except RuntimeError:
    pass
else:
    raise RuntimeError('guards inactive')
result={'status':'PASS','explicit_require_calls':require_calls,'invalid_input_controls':hostiles,'square_and_even_family_k_max':100,'pair_fibres_n_max':400,'closure_n_max':1000,'Gaussian_iterates':7,'historical_three_piece_partition':'recovered separately in THREE_CHAINS_RECOVERED.md','all_k_status':'proved algebraically in SUMMAND_SQUARES.md'}
print(json.dumps(result,indent=2,sort_keys=True))
