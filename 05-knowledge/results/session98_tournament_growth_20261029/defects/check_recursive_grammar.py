#!/usr/bin/env python3
"""Cycle insertion hierarchy, independently compared with partition census."""
from fractions import Fraction
from functools import lru_cache
from math import gcd,lcm
import json
from check_defects import mul,factor,census,require

def recurrence(N,D):
    unit=(1,)+(0,)*D
    def power(a,k):
        out=unit
        for _ in range(k):out=tuple(mul(out,a,D))
        return out
    I={l:power(factor(l,D),(l-1)//2) for l in range(1,N+1)}
    for l in I:
        if l%2==0:I[l]=(0,)*(l//2)+I[l]
        I[l]=I[l][:D+1]
    K={(l,j):power(factor(lcm(l,j),D),gcd(l,j)) for l in range(1,N+1) for j in range(1,N+1)}
    @lru_cache(None)
    def A(n,w):
        if n==0:return unit
        out=[Fraction(0)]*(D+1)
        for l in range(1,n+1):
            front=mul(w[l-1],I[l],D)
            if not any(front):continue
            new=tuple(tuple(mul(w[j-1],K[l,j],D)) for j in range(1,n-l+1))
            term=mul(front,A(n-l,new),D)
            out=[x+Fraction(y,n) for x,y in zip(out,term)]
        return tuple(out)
    result=A(N,tuple(unit for _ in range(N)))
    require(all(Fraction(x).denominator==1 for x in result),'recursive integrality')
    return list(map(int,result)),A.cache_info().currsize

def main():
    out=[]
    for n in range(9):
        d=min(4,n*(n-1)//2);got,states=recurrence(n,d)
        require(got==census(n,d),'recursive grammar vs partition census')
        out.append({'n':n,'D':d,'counts':got,'states':states})
    print(json.dumps({'status':'PASS','comparisons':out},sort_keys=True))
if __name__=='__main__':main()
