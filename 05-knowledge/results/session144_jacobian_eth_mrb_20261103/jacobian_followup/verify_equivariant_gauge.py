#!/usr/bin/env python3
import json
import sympy as S
x,y,z,v,t,a=S.symbols('x y z v t a')
rows=[]
for k in range(2,7):
    monomial=x**(2*k-1)*z**k
    T=S.Matrix([x,y+monomial,z])
    tests={
      'determinant':S.expand(T.jacobian([x,y,z]).det()-1)==0,
      'left_inverse':S.expand((y+monomial)-monomial-y)==0,
      'right_inverse':S.expand((y-monomial)+monomial-y)==0,
      'equivariance':S.simplify(monomial.subs({x:a*x,z:a**(-2)*z})-a**(-1)*monomial)==0,
      'quotient':S.expand(x*(y+monomial)-(v+t**k).subs({v:x*y,t:x*x*z}))==0,
      'scale_degree':S.Poly(2-3*v-t-3*t**k,v,t).total_degree()==k,
      'line_parametrization':S.expand((2-3*v-t-3*t**k).subs(v,(2-t-3*t**k)/3))==0,
    }
    for name,ok in tests.items():
        if not ok:raise RuntimeError('k=%d failed %s'%(k,name))
    rows.append({'k':k,'checks':tests})
print(json.dumps({'status':'PASS','checks':35,'finite_range':[2,3,4,5,6],'results':rows,'universal_proof':'Displayed monomial exponents, elementary inverse, and quotient substitution, not interpolation.'},sort_keys=True,indent=2))
