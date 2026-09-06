#!/usr/bin/env python3
"""Exact hostile probe of the two-minimal-rung trinomial seed claim.

Universe: normalized f=z^-a + t*z^b + z^c, 1<=a,c<=height,
-a<b<c, b!=0; reflection reduces b>0. Extremes nonzero; t!=0.
The coefficient gauge normalization preserves every moment zero.
"""
import argparse
from math import factorial, gcd
import sympy as s

t=s.symbols('t')

def require(condition,detail):
    if not condition:
        raise RuntimeError(detail)

def moment(a,b,c,m):
    terms={}
    for z in range(m+1):
        numerator=a*m-(a+c)*z
        if numerator<0 or numerator%(a+b):
            continue
        y=numerator//(a+b)
        x=m-y-z
        if x<0:
            continue
        terms[y]=factorial(m)//(factorial(x)*factorial(y)*factorial(z))
    return s.Poly.from_dict({(k,):v for k,v in terms.items()},t)

def probe(height):
    patterns=ambiguous=0
    for a in range(1,height+1):
        for c in range(2,height+1):
            for b in range(1,c):
                if gcd(gcd(a,b),c)!=1:
                    continue
                patterns+=1
                m0=next(m for m in range(1,a+c+1) if not moment(a,b,c,m).is_zero)
                p=moment(a,b,c,m0)
                if len(p.terms())<2:
                    continue
                ambiguous+=1
                q=moment(a,b,c,2*m0)
                common=s.gcd(p,q)
                if len(common.terms())>1:
                    print('HOSTILE',a,b,c,m0,p.as_expr(),q.as_expr(),common.as_expr(),flush=True)
                    return
    print('PASS',{'height':height,'primitive_patterns':patterns,'ambiguous_minimal_rungs':ambiguous},flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--height',type=int,default=20)
    args=parser.parse_args()
    require(moment(2,1,4,3).as_expr()==3*t**2+3,'first-return hostile control')
    require(moment(2,1,4,6).as_expr()==15*t**4+60*t**2+15,'doubled-return positive control')
    print('controls: first-rung cancellation t^2=-1, next moment -30; exact positive control passed')
    probe(args.height)
