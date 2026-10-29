"""Symbolic discriminant audit: fibre geometry, not generic MW rank17."""
import json,sys
from pathlib import Path
import sympy as S
from named_family import family
T=S.symbols('T');out=[]
for i in [273,302,4239]:
 a1,a2,a3,a4,a6=family(i,T)[0];b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
 c4=S.Poly(b2*b2-24*b4,T);c6=S.Poly(-b2**3+36*b2*b4-216*b6,T);D=S.Poly(-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6,T)
 if c4**3-c6**2!=1728*D:raise ValueError('identity')
 out.append({'curve':i,'c4_degree':c4.degree(),'c6_degree':c6.degree(),'discriminant_degree':D.degree(),'gcd_D_derivative_degree':S.gcd(D,D.diff()).degree(),'gcd_c4_D_degree':S.gcd(c4,D).degree(),'D_leading_coefficient':str(D.LC()),'D_primitive_coefficients':[str(x)for x in D.primitive()[1].all_coeffs()]})
print(json.dumps(out,indent=2))
