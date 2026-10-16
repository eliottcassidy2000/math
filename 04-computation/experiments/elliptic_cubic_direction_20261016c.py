#!/usr/bin/env python3
"""Independent exact certificate for the odd-cubic width-three direction.

Uses rational arithmetic in Q[theta,j,z]/(f3,j^2+11,f5), dimension 30.
The companion proof specifies the real embedding theta in (4/5,1),
j=sqrt(-11), z=zeta11. No floating-point root is used.
"""
import os
os.environ['SYMPY_GROUND_TYPES'] = 'python'
import json
import sympy as s


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    z,j,a=s.symbols('z j theta')
    f3=2*a**3+2*a*a-2*a-1
    f5=z**5+(1-j)*z**4/2-z**3+z**2-(1+j)*z/2-1
    gb=s.groebner([f3,j*j+11,f5],z,j,a,domain=s.QQ)
    red=lambda f:gb.reduce(s.expand(f))[1]
    basis=[a**i*j**k*z**h for h in range(5) for k in range(2) for i in range(3)]
    monoms=[(h,k,i) for h in range(5) for k in range(2) for i in range(3)]
    def vector(f):
        p=s.Poly(red(f),z,j,a)
        return s.Matrix([p.coeff_monomial(m) for m in monoms])
    den=a*z*z-a*z+1
    num=z*z*(a-a*z+z*z)
    matrix=s.Matrix.hstack(*(vector(den*b) for b in basis))
    # Independent linear-algebra construction, not the Cramer numerator route.
    values=matrix.inv()*vector(num)
    unit=s.expand(sum(c*b for c,b in zip(values,basis)))
    coeff=[s.factor(unit.coeff(z,h)) for h in range(5)]
    check(red(den*unit-num)==0,'exact division in the 30-dimensional algebra')
    check(red(z**11-1)==0,'eleventh-root relation')
    conjugate=unit.subs({z:z**10,j:-j},simultaneous=True)
    check(red(unit*conjugate-1)==0,'unit norm exactly one')
    check(coeff[4]==0 and coeff[0]!=0 and coeff[3]!=0,'exact width three')
    check(all(f3.subs(a,r)!=0 for r in [1,-1,s.Rational(1,2),s.Rational(-1,2)]),
          'cubic rational-root certificate')
    check(f3.subs(a,s.Rational(4,5))<0<f3.subs(a,1),'selected real interval')
    # f3'=6a^2+4a-2 is increasing and positive throughout this interval.
    check(s.diff(f3,a).subs(a,s.Rational(4,5))>0,'interval monotonicity')
    print(json.dumps({'status':'FINITE-EXACT independent algebraic certificate',
                      'witness_algebra_dimension':30,
                      'real_root_interval':['4/5','1'],
                      'cubic':str(f3),'quintic':str(f5),
                      'basis_order':'z exponent 0..4, j exponent 0..1, theta exponent 0..2',
                      'unit_coefficients':[str(v) for v in coeff],
                      'common_coefficient_denominator':int(s.ilcm(*(v.q for v in values))),
                      'norm_one':True,'width':3,
                      'carrier':'Q(zeta55,theta), not Q(zeta55)'},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
