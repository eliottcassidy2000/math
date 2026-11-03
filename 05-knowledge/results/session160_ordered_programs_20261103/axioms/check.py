#!/usr/bin/env python3
from math import isqrt

def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))

def in_m(n):
    return n > 0 and n % 2 == 1 and n not in (7, 21)

def proposed_atom(n):
    return (n % 2 == 1 and prime(n) and n != 7 or
            n % 7 == 0 and n//7 != 3 and n//7 % 2 == 1 and prime(n//7) or
            n in (63,343))

def prime_factors(n):
    out=[]
    while n > 1:
        p=next((d for d in range(2,isqrt(n)+1) if n%d==0),n)
        out.append(p); n//=p
    return out

def factors(n):
    a=b=0
    while n%7==0: n//=7; a+=1
    while n%3==0: n//=3; b+=1
    r=prime_factors(n)
    if a==0: return [3]*b+r
    if a%2==0: return [49]*(a//2)+[3]*b+r
    if a>=3: return [343]+[49]*((a-3)//2)+[3]*b+r
    if b>=2: return [63]+[3]*(b-2)+r
    if not r: raise ValueError('forbidden 7 or 21')
    return [7*r[0]]+[3]*b+r[1:]

def weight2(n):
    f=prime_factors(n)
    return sum(1 if p==7 else 2 for p in f)

for n in range(1,100002,2):
    if not in_m(n): continue
    fs=factors(n)
    prod=1
    for x in fs:
        if not proposed_atom(x): raise SystemExit(('nonatom factor', n,x))
        prod*=x
    if prod!=n: raise SystemExit(('bad product',n,fs))
    actual=n>1 and not any(n%d==0 and in_m(d) and in_m(n//d) and d>1 and n//d>1 for d in range(3,isqrt(n)+1,2))
    if actual!=proposed_atom(n): raise SystemExit(('bad atom classification',n,actual))
    if actual and not (2<=weight2(n)<=5): raise SystemExit(('bad weight',n))
if 63*63 != 49*3**4 or weight2(63)!=5 or weight2(49)!=2:
    raise SystemExit('bad equality')
print('OK: odd integers through 100001; atom classification, factorization, weight and equality')
