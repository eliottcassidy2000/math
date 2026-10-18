#!/usr/bin/env python3
"""Exact characteristic-zero differential identities used by Apéry descent.
Coefficients are rational functions; no numerical evaluation or sampled points.
"""
import json
import sympy as s
x=s.symbols('x');y=s.Function('y')(x)
h=1-7*x-8*x*x; A2=x*h;A1=s.diff(A2,x);A0=-8*x-2
sig=(1-8*x)/(8*(1+x));g=1/(1+x)
C2=A2*g*s.diff(sig,x)**2
C1=A2*(2*s.diff(g,x)*s.diff(sig,x)+g*s.diff(sig,x,2))+A1*g*s.diff(sig,x)
C0=A2*s.diff(g,x,2)+A1*s.diff(g,x)+A0*g
factor=s.factor(C2/A2.subs(x,sig))
rr=[s.factor(C1-factor*A1.subs(x,sig)),s.factor(C0-factor*A0.subs(x,sig))]
if rr!=[0,0]:raise RuntimeError(('Franel covariance',rr))
t=x*(1-8*x)/(1+x)
theta=lambda f:s.factor(t/s.diff(t,x)*s.diff(f,x))
f=(1+x)*y;ts=[f]
for k in range(3):ts.append(theta(ts[-1]))
LA=(1-34*t+t*t)*ts[3]+(-51*t+3*t*t)*ts[2]+(-27*t+3*t*t)*ts[1]+(-5*t+t*t)*ts[0]
coeff=[s.factor(s.expand(LA).coeff(s.diff(y,x,k))) for k in range(4)]
P=s.factor(A1/A2);Q=s.factor(A0/A2)
S=[s.factor(2*s.diff(Q,x)+4*P*Q),s.factor(s.diff(P,x)+4*Q+2*P*P),3*P,1]
res=[s.factor(coeff[k]-coeff[3]*S[k]) for k in range(4)]
if res!=[0]*4:raise RuntimeError(('symmetric square pullback',res))
print(json.dumps({'status':'EXACT rational-function identities','Franel_covariance_multiplier':str(factor),'Franel_residuals':list(map(str,rr)),'Apery_pullback_leading_factor':str(coeff[3]),'symmetric_square_residuals':list(map(str,res))},indent=2))
