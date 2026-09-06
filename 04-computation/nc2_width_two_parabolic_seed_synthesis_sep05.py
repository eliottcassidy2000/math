#!/usr/bin/env python3
"""Exact checks of the width-two moment / rational parabolic-order identity.

No dynamic basin numerics are used. The universal bound uses the cited Fatou
critical-point theorem, independently of this finite exact coefficient audit.
"""
from fractions import Fraction
from itertools import product
import sympy as s

z=s.symbols('z')

def require(condition,detail):
    if not condition:
        raise RuntimeError(detail)

def multiply(left,right):
    result={}
    for i,a in left.items():
        for j,b in right.items():
            result[i+j]=result.get(i+j,0)+a*b
    return {k:v for k,v in result.items() if v}

def check(coefficients,label=None):
    d=len(coefficients)-1
    r=s.Integer(coefficients[0])
    R=s.Poly.from_list(list(reversed(coefficients)),z,domain=s.QQ)
    require(r and R.degree()==d,('exact endpoint hypotheses',coefficients))
    laurent={i-2:Fraction(a) for i,a in enumerate(coefficients) if a}
    current={0:Fraction(1)}
    m=0
    for index in range(1,d+1):
        current=multiply(current,laurent)
        if current.get(0,0):
            m=index
            ct=s.Rational(current[0].numerator,current[0].denominator)
            break
    require(m,('bound failed',coefficients))
    # Independent rational composition, no root product / moments used here.
    powers=[s.Poly(1,z,domain=s.QQ)]
    for _ in range(d):
        powers.append(powers[-1]*R)
    H=sum((s.Poly(coefficients[k]*(-r*z)**k,z)*powers[d-k]
           for k in range(d+1)),s.Poly(0,z,domain=s.QQ))
    defect=r*r*powers[d-1]-H
    order=min(exponent[0] for exponent,coefficient in defect.terms() if coefficient)
    require(order==2*m,(coefficients,m,order))
    # T^2-z = z*defect/H and H(0)=r^(d+1).
    leading=defect.nth(order)/r**(d+1)
    require(leading==-2*ct/(m*r**m),(coefficients,m,leading,ct))
    # Remove forced critical poles from the derivative numerator.
    derivative_numerator=R-s.Poly(z,z)*R.diff()
    pole_factor=s.gcd(R,R.diff())
    critical=derivative_numerator.exquo(pole_factor)
    distinct_active=critical.sqf_part().degree()
    distinct_roots=R.sqf_part().degree()
    require(m<=distinct_active<=distinct_roots<=d,
            ('critical count',coefficients,m,distinct_active,distinct_roots,d))
    if label:
        print(label,{'degree':d,'first_moment':m,'constant_term':str(ct),
                     'iterate_defect_order':order+1,'iterate_leading':str(leading),
                     'distinct_nonpole_critical_points':distinct_active,
                     'distinct_roots_R':distinct_roots})
    return m

def main():
    print('FINITE-EXACT: direct Laurent convolution versus exact rational composition')
    check([1,1,0,1,-1],'canonical parity/cancellation hostile')
    check([2,-1,3,0,4],'nonunit endpoint / nonzero first moment')
    for d in (3,5,7,9,11,13):
        check([1]+[0]*(d-1)+[1],f'sharp binomial degree {d}')
    for d in (4,6,8,10):
        R=s.Poly((1+z)**d,z)
        check(list(reversed(R.all_coeffs())),f'repeated-root degree {d}')
    cases=0
    deepest={}
    for d in range(3,7):
        positions=[i for i in range(1,d) if i!=2]
        maximum=0
        for endpoints in product((-1,1),repeat=2):
            for interior in product((-1,0,1),repeat=len(positions)):
                coefficients=[0]*(d+1)
                coefficients[0],coefficients[d]=endpoints
                for i,a in zip(positions,interior):
                    coefficients[i]=a
                maximum=max(maximum,check(coefficients))
                cases+=1
        deepest[d]=maximum
    print('signed coefficient-box audit',{'cases':cases,'degrees':'3..6',
          'extremes':[-1,1],'nonconstant_interiors':[-1,0,1],
          'constant_coefficient':0,'max_first_return_by_degree':deepest})
    # One-sided boundary: R=1+z has every positive constant term zero and T^2=id.
    T=-z/(1+z)
    require(s.cancel(T.subs(z,T)-z)==0,'one-sided exact involution')
    print('one-sided excluded control: R=1+z, rational involution exact')
    print('PASS: every exact identity and strengthened critical-count bound verified')

if __name__=='__main__':
    main()
