#!/usr/bin/env python3
"""Independent ordinary-derivative covariance and symmetric-square audit."""
import sympy as s
import json
x=s.symbols('x');y=s.Function('y')(x)
sig=(1-8*x)/(8*(1+x));ds=s.diff(sig,x)
P=(24*x*x+14*x-1)/(x*(x+1)*(8*x-1));Q=(8*x+2)/(x*(x+1)*(8*x-1))
Pb=s.cancel(P.subs(x,sig)*ds-s.diff(sig,x,2)/ds);Qb=s.cancel(Q.subs(x,sig)*ds**2)
cov=[s.factor(Pb+2/(1+x)-P),s.factor(Qb+Pb/(1+x)-Q)]
assert cov==[0,0]
t=x*(1-8*x)/(1+x);T=s.cancel(t/s.diff(t,x))
def th(z):return T*s.diff(z,x)
U=(1+x)*y;T1=th(U);T2=th(T1);T3=th(T2)
A=s.expand((1-34*t+t*t)*T3+(-51*t+3*t*t)*T2+(-27*t+3*t*t)*T1+(-5*t+t*t)*U)
SS=[1,3*P,s.diff(P,x)+4*Q+2*P**2,2*s.diff(Q,x)+4*P*Q]
vs=[s.diff(y,x,3),s.diff(y,x,2),s.diff(y,x),y]
cs=[s.cancel(A.coeff(v)) for v in vs];factor=s.factor(cs[0])
res=[s.factor(s.cancel(a-factor*b)) for a,b in zip(cs,SS)]
assert res==[0,0,0,0]
assert s.factor(sig.subs(x,sig)-x)==0
assert s.factor(t.subs(x,sig)-t)==0
assert s.factor((1+x)*(1+sig))==s.Rational(9,8)
print(json.dumps({'Franel_covariance_residuals':list(map(str,cov)),'Apery_sym_square_pullback_residuals':list(map(str,res)),'pullback_multiplier':str(factor),'sigma_involution':True,'t_sigma_invariant':True,'gauge_product':'9/8'},indent=2))
