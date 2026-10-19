#!/usr/bin/env python3
"""Exact interface controls for the universal monomial-family proof."""
import json
import sympy as S
x,t,u,v,z,e,l,a,c,V,T=S.symbols('x t u v z e l a c V T')
checks=0
def check(a,label):
 global checks
 if not bool(a):raise RuntimeError(label)
 checks+=1
def zero(a,label):check(S.cancel(a)==0,label)
for k in range(13):
 H=x**k*t*t*(1+x*x*t);dx=S.diff(H,t);dt=-S.diff(H,x)
 zero(2*x*t*dx+x*x*dt+(k-4)*x*H,'literal u derivation')
 zero(dx.subs(t,0),'source x order')
 zero(dt.subs(t,0),'source t order zero')
 zero(S.diff(dt,t).subs(t,0),'source t order one')
 zero(H.subs({x:z*x,t:z**-2*t},simultaneous=True)-z**(k-4)*H,'Hamiltonian scaling')
 zero(dx.subs({x:z*x,t:z**-2*t},simultaneous=True)-z**(k-2)*dx,'vector x scaling')
 zero(dt.subs({x:z*x,t:z**-2*t},simultaneous=True)-z**(k-5)*dt,'vector t scaling')
for d in range(7):
 R=u**d
 op=a*u*(1+u)*S.diff(R,u)+(2+3*u)*R
 zero(op-((a*d+2)*u**d+(a*d+3)*u**(d+1)),'exactness operator monomial')
A,B=S.symbols('A B');R=A*u+B
op=S.Poly(a*u*(1+u)*S.diff(R,u)+(2+3*u)*R-a*u*(1+u),u)
zero(op.coeff_monomial(1)-2*B,'linear constant coefficient')
zero(op.coeff_monomial(u)-((a+2)*A+3*B-a),'linear u coefficient')
zero(op.coeff_monomial(u*u)-((a+3)*A-a),'linear square coefficient')
# Cubic resonance leaves the linear and quadratic coefficients incompatible.
r0,r1,r2,r3=S.symbols('r0 r1 r2 r3');R=r0+r1*u+r2*u*u+r3*u**3
op=S.Poly(-u*(1+u)*S.diff(R,u)+(2+3*u)*R+u*(1+u),u)
zero(op.coeff_monomial(1)-2*r0,'resonant constant')
zero(op.coeff_monomial(u)-(3*r0+r1+1),'resonant linear')
zero(op.coeff_monomial(u*u)-(2*r1+1),'resonant square')
# k3: use independent U symbol and the vector-field preservation derivative.
Up=V*V*(1+V)/(u*u*(1+u))
def D(g):return S.diff(g,u)+S.diff(g,V)*Up
h=1/V-1/u;r=(1+V)*u/(V*(1+u))
zero(D(h)-D(r)/r,'k3 logarithmic differential')
zero(r.subs(V,u/(1+c*u))-(1+(c+1)*u)/(1+u),'k3 constant h consequence')
# k6: rational field coordinates and primitive.
H=x**6*t*t*(1+x*x*t);zs=x**3*t
xi=z**3/(e-z*z);ti=(e-z*z)**3/z**8
zero(zs.subs({x:xi,t:ti},simultaneous=True)-z,'k6 inverse z')
zero(H.subs({x:xi,t:ti},simultaneous=True)-e,'k6 inverse E')
zero(xi.subs({z:zs,e:H},simultaneous=True)-x,'k6 inverse x')
zero(ti.subs({z:zs,e:H},simultaneous=True)-t,'k6 inverse t')
dz=S.diff(zs,x)*S.diff(H,t)-S.diff(zs,t)*S.diff(H,x)
zero(dz-x**10*t**3,'k6 source z derivative')
zero(dz.subs({x:xi,t:ti},simultaneous=True)-z**6/(e-z*z),'k6 rational z derivative')
primitive=1/(3*z**3)-e/(5*z**5)
zero(S.diff(primitive,z)*z**6/(e-z*z)-1,'k6 primitive derivative')
Q=(3*e-5*z*z-15*l*z**5)*V**5+5*z**5*V*V-3*e*z**5
zero(15*V**5*z**5*(primitive.subs(z,V)-primitive-l)-Q,'k6 cleared time equation')
remainder=S.rem((e-V*V)**3-T*V**8,Q,V)
zero(S.Poly(remainder,V).coeff_monomial(V**4)-3*e,'k6 T remainder obstruction')
zero((V**3/(e-V*V))**3*((e-V*V)**3/V**8)-V,'k6 joint field recovery')
zero((V**3/(e-V*V))*V*V+V**3-(V**3/(e-V*V))*e,'k6 X cubic')
# Coprime coefficient test symbolically modulo zeta^5=1.
y=S.symbols('y')
const=5*z**5*V**2-5*z*z*V**5-15*l*z**5*V**5
zero(S.rem(S.expand(const.subs(V,y*z)-5*z**7*(y*y-1-3*l*z**3)),y**5-1,y),'k6 root-of-unity coprimality interface')
print(json.dumps({'status':'PASS','checks':checks,'tested_k_interfaces':list(range(13)),
 'scope':'Finite exact interface controls; the all-k classification follows from the scaling, trace, exactness and exceptional-case proofs.'},indent=2,sort_keys=True))
