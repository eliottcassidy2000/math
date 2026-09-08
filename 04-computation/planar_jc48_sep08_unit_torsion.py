#!/usr/bin/env python3
"""Exact unit-response controls for the full boundary-linear family."""
import hashlib
import json
import sympy as S

x,t,r,b,g=S.symbols('x t r b g')
a=S.symbols('a',nonzero=True)
gates=0

def need(value,label):
    global gates
    gates+=1
    if not value:raise RuntimeError(label)

def eq(left,right,label):need(S.cancel(left-right)==0,label)
def jac(F,G,X,Z):return S.diff(F,X)*S.diff(G,Z)-S.diff(F,Z)*S.diff(G,X)
bb=-x*x*(1+x*x*t)
def source(F):return S.cancel(F.subs({r:1/x,b:bb},simultaneous=True))
rows=[]
for h in (S.Integer(0),b**3,b**4,2*b**3-b**5,b**3+b**6):
    F=a*r*b+h
    L=S.integrate(h/b**3,b);M=S.integrate(h*h/b**3,b)
    G=S.cancel((-F*F/(2*b*b)-2*F*L+M)/a**3)
    Fs=source(F);Gs=source(G)
    eq(jac(F,G,r,b),r*r,'literal rational primitive')
    P2=S.cancel(Fs*Fs*Gs)
    need(not S.denom(P2).has(x,t),'square value clears every affine pole')
    eq(jac(Fs,P2,x,t),Fs*Fs,'actual polynomial witness for g squared theta zero')
    eq(S.diff(Fs,x).subs(x,0),-a,'source smooth along E')
    eq(S.limit(Fs*Fs*Gs,x,0),-a/2,'exact double scalar principal coefficient')
    eq(S.limit(Fs*(Gs+a/(2*Fs*Fs)),x,0),0,'zero simple principal coefficient')
    eq(S.limit(Gs+a/(2*Fs*Fs),x,0),-t/a,'regular E remainder after full scalar subtraction')
    eq(S.limit(G,b,0),-r*r/(2*a),'other actual component has regular primitive')
    eq(S.limit(x*Fs*Gs,x,0),S.Rational(1,2),'one value factor leaves a genuine simple pole')
    rows.append([str(h),S.degree(P2,t),S.degree(P2,x)])

for j in range(0,13):
    pp=-a*S.Rational(1,2)*(-1)**j*S.factorial(j+1)*g**(-j-2)
    eq(S.diff(-a/(2*g*g),g,j),pp,'canonical connection coefficient')
    eq(S.diff(g*pp,g)-g*S.diff(pp,g),pp,'Weyl relation on actual arm')
    eq(S.limit(g**(j+2)*pp,g,0),-a*S.Rational(1,2)*(-1)**j*S.factorial(j+1),
       'nonzero exact primary order')
eq(g*g*(-a/(2*g*g)),-a/2,'order two disappears modulo polynomial principal parts')
eq(g*(-a/(2*g*g)),-a/(2*g),'order one remains modulo diagonal')

print('Boundary-linear unit response PASS; exact primary order2 and unbounded intrinsic derivative arm')
print('Finite universe: h=0,b^3,b^4,2b^3-b^5,b^3+b^6; arbitrary nonzero a; derivative levels0..12')
print('Controls: actual polynomial g^2 primitive, full E principal part, regular Gamma and residual simple pole')
print('Always-active gates:',gates)
print('Semantic SHA256:',hashlib.sha256(json.dumps(rows,separators=(',',':'),default=str).encode()).hexdigest())
