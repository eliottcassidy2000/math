#!/usr/bin/env python3
"""Exact new ell=4,5 complementary collision certificates; no inherited scan."""
import sympy as s,json,time
from math import factorial
from pathlib import Path
r,t=s.symbols('r t'); z=2*r+1
out=[]
for ell in [4,5,6,7]:
 start=time.time()
 P=sum(s.Rational(factorial(ell),factorial(i)*factorial(3*(ell-i)))*s.prod(z+j for j in range(2*i+1,2*ell+1))*t**i for i in range(ell+1))
 Q=sum(s.Rational(factorial(2*ell),factorial(i)*factorial(3*(2*ell-i)))*s.prod(2*z+j for j in range(2*i+1,4*ell+1))*t**i for i in range(2*ell+1))
 result=s.Poly(s.resultant(P,Q,t),r)
 forced=s.prod((r+j)**j*(2*r+2*j+1)**j for j in range(1,ell+1))
 residual=s.Poly(s.cancel((-1)**ell*result.as_expr()/forced),r)
 assert all(c>0 for c in residual.all_coeffs())
 item={'ell':ell,'resultant_degree':result.degree(),'residual_degree':residual.degree(),'residual_coefficients_desc':[str(x) for x in residual.all_coeffs()],'positive_coefficients':True}
 # Independent literal Sylvester determinant controls. Polynomial identity
 # additionally exact through the symbolic resultant division itself.
 for rr in [0,1,2,3,7,19]:
  pc=s.Poly(P.subs(r,rr),t).all_coeffs(); qc=s.Poly(Q.subs(r,rr),t).all_coeffs(); M=s.zeros(3*ell)
  for i in range(2*ell):
   for j,c in enumerate(pc):M[i,i+j]=c
  for i in range(ell):
   for j,c in enumerate(qc):M[2*ell+i,i+j]=c
  assert M.det()==result.eval(rr)
 item['sylvester_controls']=[0,1,2,3,7,19]
 out.append(item); Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASS',ell,'degrees',result.degree(),residual.degree(),'positive',len(residual.all_coeffs()),'seconds',round(time.time()-start,2),flush=True)
print('PASS: four new all-r symbolic strips; general ell and physical positive-x noncancellation remain OPEN')
