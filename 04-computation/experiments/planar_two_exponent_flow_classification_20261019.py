#!/usr/bin/env python3
"""Exact interface controls for two-exponent flow classification."""
import json
import sympy as S
x,t,u,z,e,l,a,v=S.symbols('x t u z e l a v')
checks=0
def check(b,label):
 global checks
 if not bool(b):raise RuntimeError(label)
 checks+=1
def zero(b,label):check(S.cancel(b)==0,label)
def sub(b,X,T):return b.subs({x:X,t:T},simultaneous=True)
def deriv(b,H):return S.diff(b,x)*S.diff(H,t)-S.diff(b,t)*S.diff(H,x)
us=x*x*t
for m in range(2,7):
 for k in range(0,2*m+4):
  aa=k-2*m;H=x**k*t**m*(1+us);dx=S.diff(H,t);dt=-S.diff(H,x)
  zero(deriv(us,H)+aa*x*H,'general u derivative')
  check(S.cancel(dx/t**(m-1)).is_polynomial(x,t),'source x order')
  check(S.cancel(dt/t**m).is_polynomial(x,t),'source t order')
  zero(sub(H,z*x,z**-2*t)-z**aa*H,'general scaling')
 for d in range(5):
  R=u**d
  zero(a*u*(1+u)*S.diff(R,u)+(m+(m+1)*u)*R-((a*d+m)*u**d+(a*d+m+1)*u**(d+1)),'exactness monomial')
 B=sum(S.Rational((-1)**j,j-m+1)*u**(j-m+1) for j in range(m-1))
 zero(S.diff(B,u)+(-1)**m*(1/(1+u)-1/u)-1/(u**m*(1+u)),'resonant partial fractions')
 Up=v**m*(1+v)/(u**m*(1+u));r=(1+v)*u/(v*(1+u));h=-(B.subs(u,v)-B)/((-1)**m)
 D=lambda b:S.diff(b,u)+S.diff(b,v)*Up
 zero(D(h)-D(r)/r,'resonant logarithmic identity')
positive=[(2,1),(2,2),(3,1),(3,3),(4,2),(6,3)]
negative=[(3,2),(5,2),(5,3),(7,2),(7,4),(8,3)]
for m,aa in positive:
 r=m//aa;k=2*m+aa;H=x**k*t**m*(1+us)
 xi=z**(m+1)/(e-z**aa)**r;ti=(e-z**aa)**(2*r+1)/z**(2*m+aa+2);zs=x*us**r
 zero(sub(zs,xi,ti)-z,'positive inverse parameter')
 zero(sub(H,xi,ti)-e,'positive inverse E')
 zero(xi.subs({z:zs,e:H},simultaneous=True)-x,'positive inverse x')
 zero(ti.subs({z:zs,e:H},simultaneous=True)-t,'positive inverse t')
 zero(sub(deriv(zs,H),xi,ti)-z**(m+aa+2)/(e-z**aa)**r,'positive parameter velocity')
 P=-S.integrate(z**aa*(e*z**aa-1)**r,z)
 check(S.degree(P,z)==m+aa+1,'positive primitive degree')
 zero(S.diff(P,z)*(-1/(z**aa*(e*z**aa-1)**r))-1,'positive primitive derivative')
for m,b in negative:
 r=(m+1)//b;k=2*m-b;H=x**k*t**m*(1+us)
 xi=z/(e*z**b-1)**r;ti=(e*z**b-1)**(2*r-1)/z**2;zs=x/us**r
 zero(sub(zs,xi,ti)-z,'negative inverse parameter')
 zero(sub(H,xi,ti)-e,'negative inverse E')
 zero(xi.subs({z:zs,e:H},simultaneous=True)-x,'negative inverse x')
 zero(ti.subs({z:zs,e:H},simultaneous=True)-t,'negative inverse t')
 zero(sub(deriv(zs,H),xi,ti)+1/(z**(b-2)*(e*z**b-1)**(r-2)),'negative parameter velocity')
 P=-S.integrate(z**(b-2)*(e*z**b-1)**(r-2),z)
 check(S.degree(P,z)==m-b,'negative primitive degree')
 zero(S.diff(P,z)*(-1/(z**(b-2)*(e*z**b-1)**(r-2)))-1,'negative primitive derivative')
# The exceptional rational branch, without treating a primitive degree as a field degree.
H=x**4*t**3*(1+us);D=1-l*x**3*t*t;B=(1+us)*D*D-us
X=x*D/B**2;T=t*B**3/D**2
zero(sub(H,X,T)-H,'exception energy')
zero(sub(1/(x**3*t*t),X,T)-(1/(x**3*t*t)-l),'exception translation')
zero(deriv(1/(x**3*t*t),H)+1,'exception infinitesimal translation')
check(S.gcd(D,t*B)==1,'exception nonpolynomial denominator')
zero(D.subs(t,0)-1,'exception source unit D')
zero(B.subs(t,0)-1,'exception source unit B')
# Generic affine primitive obstruction: d-1 coefficient under translation.
beta=S.symbols('beta')
for d in range(2,9):
 P=z**d+3*z**(d-1)+z+1
 zero(S.Poly(S.expand(P.subs(z,z+beta)-P),z).coeff_monomial(z**(d-1))-d*beta,'primitive translation coefficient')
print(json.dumps({'status':'PASS','checks':checks,'positive_parameter_controls':positive,'negative_parameter_controls':negative,
 'scope':'Finite exact interface controls; universal loci follow from scaling, trace, pole, primitive and rational-translation proofs.'},indent=2,sort_keys=True))
