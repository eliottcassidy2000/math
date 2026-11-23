#!/usr/bin/env python3
"""Exact prefix splitting and dyadic/beta tape replay."""
from fractions import Fraction as F
from itertools import product

def g(w,x):
    return sum(x**(j-i-1) for i,a in enumerate(w) if a=='O' for j in range(i+1,len(w)) if w[j]=='I')
def b(w,x):
    return sum(x**i for i,a in enumerate(w) if a=='I')
def H(w):
    return 1+g(w,2)+g(w,0)
def Q(w):
    return F(2)**(len(w)+1)+F(2)**(len(w)-1)*(g(w,F(3,2))+g(w,F(1,2)))

k=0
for n in range(2,9):
    for mid in product('OI',repeat=n-2):
        w='O'+''.join(mid)+'I'
        for a in (1,2,5):
            L='O'*a+w;R='O'+'I'*a+w
            assert H(L)==H(w)+(2**a-1)*b(w,2)
            assert H(R)==H(w)+2**a*(1+b(w,2))
            x=F(3,2)
            assert g(L,x)==g(w,x)+2*(x**a-1)*b(w,x)
            assert g(R,x)==g(w,x)+2*(x**a-1)+x**a*b(w,x)
            # Literal count formulas produce integer heights.
            assert Q(L).denominator==Q(R).denominator==1
            k+=1
print(k,'binary and base-3/2 exact prefix tests')
