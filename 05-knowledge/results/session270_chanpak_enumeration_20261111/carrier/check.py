#!/usr/bin/env python3
"""Exact small-universe controls for the common marked carrier."""
from fractions import Fraction
from math import gcd, factorial, comb
from collections import Counter
from itertools import permutations
import json
from pathlib import Path

def partitions(n,lo=1):
    if not n:
        yield ()
    for p in range(lo,n+1,2):
        for tail in partitions(n-p,p): yield (p,)+tail

def count(k,m):
    out=Fraction(0)
    for p in partitions(k):
        e=sum((a-1)//2 for a in p)+sum(gcd(a,b) for i,a in enumerate(p) for b in p[i+1:])
        z=1
        for a,r in Counter(p).items(): z*=a**r*factorial(r)
        out+=Fraction(1<<(e+m*len(p)),z)
    if out.denominator!=1:raise ArithmeticError('nonintegral class count')
    return int(out)

def path_mass(k,m):
    power=k*(k-1)//2-2*k
    x=Fraction(factorial(k)*sum(comb(k+1,j)*(2**k*j)**m for j in range(k+2)))
    x*=2**power if power>=0 else Fraction(1,2**(-power))
    if x.denominator!=1:raise ArithmeticError('nonintegral path mass')
    return int(x)

def universe(k,m):
    n=k+m;pairs=[(a,b) for a in range(k) for b in range(a+1,k)]
    pairs += [(a,k+b) for a in range(k) for b in range(m)]
    for mask in range(1<<len(pairs)):
        out=[0]*n
        for i,(a,b) in enumerate(pairs):
            if mask>>i&1:out[a]|=1<<b
            else:out[b]|=1<<a
        for a in range(k,n):
            for b in range(a+1,n):out[a]|=1<<b
        yield out

def canon(out,k):
    n=len(out);best=None
    for p in permutations(range(k)):
        order=p+tuple(range(k,n));code=0
        for i,a in enumerate(order):
            for b in order[i+1:]:code=(code<<1)|((out[a]>>b)&1)
        best=code if best is None else min(best,code)
    return best

def hamilton(out):
    n=len(out);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    full=(1<<n)-1
    for mask in range(1,1<<n):
        for v in range(n):
            ways=dp[mask][v]
            if not ways:continue
            nxt=out[v]&(full^mask)
            while nxt:
                bit=nxt&-nxt;nxt-=bit;w=bit.bit_length()-1
                dp[mask|bit][w]+=ways
    return sum(dp[full])

def run():
    rows=[]
    for k in range(1,5):
        for m in range(4 if k<=3 else 2):
            classes=set();mass=0;total=0
            for out in universe(k,m):
                classes.add(canon(out,k));mass+=hamilton(out);total+=1
            expected=count(k,m);h_expected=path_mass(k,m)
            if (len(classes),mass)!=(expected,h_expected):raise AssertionError((k,m,len(classes),mass,expected,h_expected))
            rows.append({'k':k,'m':m,'labelled_structures':total,'marked_classes':expected,'path_mass':mass})
    large={'B3_1000_bits':count(3,1000).bit_length(),'B4_1000_bits':count(4,1000).bit_length(),'Z2_1000_bits':path_mass(2,1000).bit_length()}
    result={'status':'FINITE-EXACT','controls':rows,'large_queries':large}
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run()
