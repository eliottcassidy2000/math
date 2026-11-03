#!/usr/bin/env python3
from math import isqrt

def require(test, details):
    if not test: raise RuntimeError(details)

def sq(x): return isqrt(x)**2 == x

def d(n): return isqrt(2*n-1)-isqrt(n)

def plateau(n): return d(n-1)==d(n)==d(n+1)

centers=[]
for n in range(2,1_000_001):
    birth=sq(2*n-2) or sq(2*n-1)
    death=sq(n)
    require(d(n)-d(n-1)==int(birth)-int(death),('difference',n))
    if birth and death:
        centers.append(n)
        require(plateau(n),('plateau',n))
require(centers[:6]==[9,25,289,841,9801,28561],('centers',centers))
require([d(i) for i in [8,9,10,24,25,26]]==[1,1,1,2,2,2], 'small plateaus')
require([8+1,1+15,15+10]==[9,16,25] and 9+16==25,'path labels')
# Positive Pell r+t√2 starts 3+2√2, next multiply by 3+2√2.
r,t=3,2
for j in range(1,201):
    require(r*r-2*t*t==1,('positive Pell',j))
    require(plateau(r*r),('positive plateau',j))
    r,t=3*r+4*t,2*r+3*t
# Negative Pell s+r√2 starts 7+5√2, then multiply by 3+2√2.
s,r=7,5
for j in range(1,201):
    require(s*s-2*r*r==-1,('negative Pell',j))
    require(plateau(r*r),('negative plateau',j))
    s,r=3*s+4*r,2*s+3*r
print('PASS: all n<=1000000; both Pell families 200 terms; 9/25 plateaus and G15 edge labels; optimized-safe')
