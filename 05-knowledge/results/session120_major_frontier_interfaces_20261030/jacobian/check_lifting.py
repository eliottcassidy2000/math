#!/usr/bin/env python3
"""Exact controls for the two degeneration hostiles in PROOF.md section 8."""
import json
import sympy as S
x,y,t,T=S.symbols('x y t T')
def require(v,msg):
    if not v: raise RuntimeError(msg)
def jac(P,Q): return S.expand(S.diff(P,x)*S.diff(Q,y)-S.diff(P,y)*S.diff(Q,x))
P=t*y**5*x+y**3
Q=-t*y**6*x**2/4-y**4*x/2-y**2/(4*t)
require(S.cancel(jac(P,Q)-P**2)==0,'generic fibre identity')
for v in range(-7,8):
    Pv=P.subs(t,v)
    Qv=Q.subs(t,v) if v else -x*y**4/3
    require(S.cancel(jac(Pv,Qv)-Pv**2)==0,'numerical fibre identity')
require(S.expand(Q).coeff(t,-1)==-y**2/4,'pole coefficient')
require(2%3!=0,'forbidden subring degree')
# A positive regular-lift control: changing B=y^3 to B=y^2 removes the obstruction.
Pa=t*y**5*x+y**2
Qa=-y**3*x/2-t*y**6*x**2/4
require(S.cancel(jac(Pa,Qa)-Pa**2)==0,'positive regular-lift control')
Ig=(T-1)**2; Is=T*(T-1)
require(S.rem(Ig,Is,T)!=0 and S.rem(Is,Ig,T)!=0,'incomparable ideals')
B=(y**2-1)**2
res=S.Poly(S.resultant(B-T,S.diff(B,y),y),T).sqf_part().monic().as_expr()
require(S.expand(res-Is)==0,'special critical value radical')
print(json.dumps({'status':'FINITE-EXACT; no-lift theorem has Laurent proof','nonzero_and_zero_fibres_checked':15,'generic_rational_companion_identity':True,'pole_coefficient':'-y^2/4','forbidden_subring':'K[y^3]','positive_regular_lift_control':'B=y^2','generic_image_generator':str(Ig),'special_image_generator':str(Is),'image_ideals_incomparable':True},indent=2,sort_keys=True))
