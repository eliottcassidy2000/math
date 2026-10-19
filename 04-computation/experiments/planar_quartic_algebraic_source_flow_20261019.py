#!/usr/bin/env python3
"""Exact algebraic and finite-jet controls for the quartic flow proof."""
import json
import sympy as S
x,t,u,e,v,C,L=S.symbols('x t u e v C L')
checks=0
def check(a,label):
 global checks
 if not bool(a):raise RuntimeError(label)
 checks+=1
def zero(a,label):check(S.cancel(a)==0,label)
def P(a):return a**3/S.Integer(3)+a**4/S.Integer(4)
H=x**5*t**2*(1+x*x*t);us=x*x*t
DX=S.diff(H,t);DT=-S.diff(H,x)
def delta(a):return S.expand(S.diff(a,x)*DX+S.diff(a,t)*DT)
zero(DX-(2*x**5*t+3*x**7*t**2),'source x derivation')
zero(DT-(-5*x**4*t**2-7*x**6*t**3),'source t derivation')
zero(delta(H),'energy invariant')
zero(delta(us)+x*H,'u derivative')
zero(delta(P(us))+H**2,'time coordinate derivative')
zero(delta(delta(P(us))),'time coordinate second derivative')
xi=e/(u*u*(1+u));ti=u**5*(1+u)**2/e**2
zero((x*x*t).subs({x:xi,t:ti},simultaneous=True)-u,'inverse u')
zero(H.subs({x:xi,t:ti},simultaneous=True)-e,'inverse E')
zero(xi.subs({u:us,e:H},simultaneous=True)-x,'inverse x')
zero(ti.subs({u:us,e:H},simultaneous=True)-t,'inverse t')
zero((P(v)-P(u)).subs(v,u),'simple root value')
zero(S.diff(P(v)-P(u),v).subs(v,u)-u*u*(1+u),'simple root derivative')
check(u*u*(1+u)!=0,'nonzero generic derivative')
quartic=v**4+S.Rational(4,3)*v**3-4*C
remainder=(4*C+S.Rational(16,81))*v**3+S.Rational(8,3)*C*v*v+S.Rational(4,9)*C*v-S.Rational(16,27)*C
zero(S.rem(v**5*(1+v)**2,quartic,v)-remainder,'coordinate T remainder')
zero(4*(P(v)-C)-quartic,'quartic normalization')
zero((e/(v*v*(1+v)))**2*(v**5*(1+v)**2/e**2)-v,'joint field recovery')
for time in [1,-1,2,S.Rational(1,3)]:
 c=P(u)-time*e*e
 check(S.Poly(c,e).degree()==2,'nonzero C')
 check(S.Poly(remainder.subs(C,c),v).degree()==3,'nonconstant remainder')
 zero(S.rem(v**5*(1+v)**2,quartic.subs(C,c),v)-remainder.subs(C,c),'scalar remainder')
# Finite source jets: derive first, truncate in t, then compose polynomial identities.
N=8
def trunc(a):
 return S.Add(*(c*x**i*t**j for (i,j),c in S.Poly(S.expand(a),x,t).terms() if j<N))
def jet(a,time):
 term=a;out=a
 for n in range(1,N+1):
  term=trunc(delta(term));out+=time**n*term/S.factorial(n)
 return trunc(out)
for time in [1,-1,2]:
 X=jet(x,time);T=jet(t,time);V=jet(us,time)
 zero(trunc(X*X*T-V),'jet joint u')
 zero(trunc(P(V)-P(us)+time*H*H),'jet exact time relation')
 zero(trunc(X**5*T*T*(1+X*X*T)-H),'jet energy')
 zero(X.subs(t,0)-x,'jet fixed line x')
 zero(T.subs(t,0),'jet fixed line t')
 zero(trunc(X*V*V*(1+V)-H),'jet coordinate x identity')
 zero(trunc(T*H*H-V**5*(1+V)**2),'jet coordinate t identity')
print(json.dumps({'status':'PASS','checks':checks,'source_t_cutoff':N,'joint_extension_degree_proved':4,
 'scope':'Finite exact controls; universal degree and nonrationality follow from the written prime-valuation and degree proof.'},indent=2,sort_keys=True))
