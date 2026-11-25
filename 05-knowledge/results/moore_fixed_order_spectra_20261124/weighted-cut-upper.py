#!/usr/bin/env python3
"""Exact diagonal-scaling / weighted Bregman bound, no graph is needed."""
from fractions import Fraction as F
from math import isqrt, factorial

n=3250
m=130
w=[F(z,1000) for z in (1000,1000,1000,1048,1208,1397,1619,1876,2176,2527,2934,3408,3960,4600,5343,6207,7208,8367,9707,11252,13022,15038,17288,19661,20431)]
ap=4
M=len(w)
assert M*m==n and w==sorted(w)
e0=F(27182818,10**7)
assert sum(F(1,factorial(i)) for i in range(13))>e0

def f(i,j):
    return w[j]/w[i]-w[i]/w[j] if i<j else F(0)
def sqrt_up(x,S=1000):
    return F(isqrt((x.numerator*S*S)//x.denominator)+1,S)
def rectangle(i,j):
    x,y=(i+1)*m,(M-j)*m
    return F(115*x*y,2*n)+F(15,2)*sqrt_up(F(x*y*(n-x)*(n-y),n*n))
s=M*m*(m-1)//2 + sum(m*(a/ap+ap/a) for a in w)
for i in range(M):
    for j in range(i+1,M):
        mix=f(i,j)-f(i+1,j)-f(i,j-1)+f(i+1,j-1)
        assert mix>=0
        s+=m*m*w[i]/w[j]+mix*rectangle(i,j)
row_correction=6*(ap+m*sum(w))
extra=sum(f(i-1,i) for i in range(1,M))
for K in (0,extra):
    base=(s+row_correction+K)/((n+1)*e0)
    # Round upward to one thousandth before powering to keep integers small.
    base=F(-(-base.numerator*1000//base.denominator),1000)
    print(base)
    for b in (27049,27050,27051,27052):
        print(b,base.numerator**(n+1)<(base.denominator**(n+1)<<b))


# Nonuniform AM-GM can exploit the remaining variation between bins.
w=[F(x,1000) for x in (1000,1000,1001,1008,1044,1146,1303,1517,1760,2065,2412,2820,3299,3858,4518,5288,6189,7233,8477,9915,11584,13508,15907,18467,20748)]
z=[F(x,1000) for x in (788,862,953,1053,1152,1234,1287,1291,1309,1300,1302,1304,1303,1306,1305,1304,1301,1301,1296,1289,1280,1271,1213,1134,1032)]
def g(i,j):
    return z[j]*w[j]/w[i]-z[i]*w[i]/w[j] if i<j else F(0)
def rect_bounds(i,j):
    x,y=(i+1)*m,(M-j)*m
    mean=F(115*x*y,2*n)
    d=F(15,2)*sqrt_up(F(x*y*(n-x)*(n-y),n*n))
    return max(F(0),mean-d),mean+d
s=6*ap+sum(m*ap/w[i]+z[i]*(m*(m-1)//2+m*w[i]/ap+6*m*w[i]) for i in range(M))
for i in range(M):
    for j in range(i+1,M):
        mix=g(i,j)-g(i+1,j)-g(i,j-1)+g(i+1,j-1)
        lower,upper=rect_bounds(i,j)
        s+=m*m*z[i]*w[i]/w[j]+mix*(upper if mix>=0 else lower)
from math import prod
Z=prod(z)**m
extra=sum(g(i-1,i) for i in range(1,M))
for K in (0,extra):
    base=(s+K)/((n+1)*e0)
    base=F(-(-base.numerator*1000//base.denominator),1000)
    print("refined",base)
    for b in (27020,27021,27022):
        ok=base.numerator**(n+1)*Z.denominator<(base.denominator**(n+1)*Z.numerator<<b)
        print(b,ok)
        if b==27021:
            assert ok
    assert 5*base.numerator**(n+1)*Z.denominator < 4*(base.denominator**(n+1)*Z.numerator<<27021)
p=13511
assert 5*((1<<(p-1))*((1<<p)-1)-1)>4*(1<<27021)
