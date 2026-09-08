#!/usr/bin/env python3
"""Exact controls for the unbounded DG boundary-linear exclusion.

The rational-map degree and incompatible special-fibre valuations are
proved in the companion note, not extrapolated from these finite controls.
"""
import hashlib
import json
import sympy as S

x,t,r,b,g=S.symbols('x t r b g')
a,lam=S.symbols('a lam', nonzero=True)
c0,c1,c2,c3,c4,c5=S.symbols('c0:6')
gates=0
records=[]

def check(value,label):
    global gates
    gates+=1
    if not value:
        raise RuntimeError(label)

def eq(left,right,label):
    check(S.cancel(left-right)==0,label)

def jac(F,G,X,Z):
    return S.diff(F,X)*S.diff(G,Z)-S.diff(F,Z)*S.diff(G,X)

bb=-x*x*(1+x*x*t)
vv=x*(1+x*x*t)
pull={r:1/x,b:bb}

def source(F):
    return S.cancel(F.subs(pull,simultaneous=True))

eq(source(r*b),-vv,'literal global rb')
eq(jac(1/r,-r*r-r**4*b,r,b),r*r,'boundary Jacobian multiplier')
eq(source(c0*r),c0/x,'only constant coefficient of A produces a pole')
for m in range(1,10):
    FF=source(r*b**m)
    check(S.denom(FF)==1,'all positive powers give global boundary-linear monomials')
    if m>=2:
        eq(S.diff(FF,x).subs(x,0),0,'higher multiplicity source critical x derivative')
        eq(S.diff(FF,t).subs(x,0),0,'higher multiplicity source critical t derivative')
    eq(S.limit(FF/x**(2*m-1),x,0),(-1)**m,'exact x valuation of rb^m')
for m in range(0,10):
    FF=bb**m
    eq(S.diff(FF,x).subs(x,0),0,'B(b) source critical x derivative')
    eq(S.diff(FF,t).subs(x,0),0,'B(b) source critical t derivative')

# Relative derivative and its twice-F necessary primitive.
A=S.Function('A')(b);B=S.Function('B')(b)
rhs=lam*(g-B)**2/A**3
eq(S.diff(rhs,g,2),2*lam/A**3,'twice independent-coordinate derivative')
for m in range(1,9):
    primitive=b**(1-3*m)/(1-3*m)
    eq(S.diff(primitive,b),1/b**(3*m),'one-pole rational primitive positive control')
for degree in range(2,10):
    for roots in range(2,degree+1):
        check(3*degree-1>3*degree-roots,'multiple poles violate rational-map zero bound')
eq(S.residue(1/(b*(b-1))**3,b,0),-6,'two-root nonzero residue hostile')

# The residue classification is a symbolic identity, not a search over B.
h=c1*b+c2*b*b+c3*b**3+c4*b**4+c5*b**5
eq(S.residue(lam*(g-h)**2/(a**3*b**3),b,0),
   lam*(c1*c1-2*g*c2)/a**3,'full residue classification')
h=c3*b**3+c4*b**4+c5*b**5
L=S.integrate(h/b**3,b);M=S.integrate(h*h/b**3,b)
Grel=lam/a**3*(-g*g/(2*b*b)-2*g*L+M)
eq(S.diff(Grel,b),lam*(g-h)**2/(a**3*b**3),'full symbolic rational primitive')
F=a*r*b+h+c0
G0=S.cancel(Grel.subs(g,F-c0))
Q=S.cancel(G0+lam*r*r/(2*a))
check(S.Poly(Q,r,b).degree(r)<=1,'remaining numerator global rb form')
eq(S.Poly(Q,r,b).coeff_monomial(r),0,'no bare r term in global correction')
check(S.denom(source(Q)).has(x)==False,'complete correction has no source pole')
eq(jac(F,G0,r,b),lam*r*r,'symbolic boundary-chart rational Jacobian')
eq(jac(source(F),source(G0),x,t),lam,'independent source-chart rational Jacobian')
eq(S.limit(x*x*source(G0),x,0),-lam/(2*a),'exact double pole at E')
eq(S.limit((source(F)-c0)/x,x,0),-a,'special fibre reduced at E')
eq(S.limit((F-c0)/b,b,0),a*r,'special fibre reduced at Gamma')
eq(S.limit(G0,b,0),-lam*r*r/(2*a),'primitive regular along Gamma')
for k in range(1,6):
    Hrep=1/(F-c0)**k
    eq(S.limit(b**k*Hrep,b,0),1/(a*r)**k,'every value pole reappears at Gamma')

# Arbitrarily high degrees are represented by an exact closed formula in k;
# this finite bank checks its literal algebra and sign at k=3,...,9.
for k in range(3,10):
    FF=r*b+b**k
    TT=r*r+S.Rational(2*k,k-2)*r*b**(k-1)+S.Rational(k*k,(k-2)*(k-1))*b**(2*k-2)
    GG=-TT/2
    eq(jac(FF,GG,r,b),r*r,'monomial-family boundary Jacobian')
    eq(S.limit(x*x*source(GG),x,0),S.Rational(-1,2),'monomial-family pole')
    check(S.denom(source(FF))==1,'monomial-family global first coordinate')
    records.append([k,str(S.expand(TT))])

FF=r*b+b**4;HH=b*b;GG=-r*r/2-2*r*b**3-S.Rational(4,3)*b**6
Fs=source(FF);Gs=source(GG)
eq(FF-HH*HH,r*b,'quartic global square prefix and linear remainder')
check(S.degree(Fs,t)==4,'actual quartic source degree')
eq(jac(Fs,Gs,x,t),1,'quartic rational mate source Jacobian')
eq(S.diff(Fs,x).subs(x,0),-1,'quartic source noncritical at E')
eq(S.diff(FF,r),b,'quartic critical-point first equation')
eq(S.diff(FF,b),r+4*b**3,'quartic critical-point second equation')
eq(FF.subs(r,0),b**4,'nonconstant boundary restriction')
eq(S.limit(x*x*Gs,x,0),S.Rational(-1,2),'quartic E principal part')
eq(GG.subs(b,0),-r*r/2,'quartic Gamma regularity')

print('DG boundary-linear carrier controls PASS; unbounded analytic proof required')
print('Universe: global monomials rb^m (m=1..9), B=b^m (m=0..9), arbitrary residue coefficients through degree5')
print('Positive controls: explicit rational primitives, both literal chart Jacobians, k=3..9 monomial family')
print('Hostiles: two-root residue, source-critical multiplicities, and source-critical-free global-root quartic with incompatible fibre poles')
print('Always-active gates:',gates)
print('Semantic SHA256:',hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
