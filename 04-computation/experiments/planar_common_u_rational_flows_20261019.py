#!/usr/bin/env python3
"""Exact controls for common-u rational Hamiltonian flows; no numerical tests."""
import json
import sympy as S
x,t,z,s,L,M=S.symbols('x t z s L M')
u=x*x*t;q=1+u
checks=0
def check(a,label):
 global checks
 if not bool(a):raise RuntimeError(label)
 checks+=1
def zero(a,label):check(S.cancel(a)==0,label)
def flow(f,time):
 a=S.diff(f,z).subs(z,u);D=1-time*x*a
 return x/D,t*D*D,D
def sub(a,X,T):return a.subs({x:X,t:T},simultaneous=True)
def delta(a,f):
 H=f.subs(z,u)
 return S.diff(a,x)*S.diff(H,t)-S.diff(a,t)*S.diff(H,x)
fs=[z*z,z*z*(z+1),z*z*(z+1)**2,z**3-z*z,2*z**4+3*z*z+5]
for f in fs:
 a=S.diff(f,z).subs(z,u)
 check(S.diff(f,z).subs(z,0)==0,'source guard')
 zero(delta(u,f),'invariant derivation')
 zero(delta(x,f)-x*x*a,'x derivation')
 zero(delta(t,f)+2*x*t*a,'t derivation')
 zero(delta(delta(t,f),f)-2*x*x*t*a*a,'t second derivative')
 zero(delta(delta(delta(t,f),f),f),'t third derivative')
 p=x
 for n in range(1,5):
  p=delta(p,f)
  zero(p-S.factorial(n)*x**(n+1)*a**n,'x repeated derivative')
 X,T,D=flow(f,L)
 zero(X*X*T-u,'flow invariant')
 zero(S.diff(X,x)*S.diff(T,t)-S.diff(X,t)*S.diff(T,x)-1,'rational Jacobian')
 zero(D.subs(t,0)-1,'source denominator unit')
 zero(X.subs(t,0)-x,'fixed line x')
 zero(T.subs(t,0),'fixed line t')
 check(S.gcd(D,x)==1,'coprime numerator')
 zero(sub(x/(1+L*x*a),X,T)-x,'inverse x')
 zero(sub(t*(1+L*x*a)**2,X,T)-t,'inverse t')
 boundary={x:-1/s,t:-s*s}
 shift=s+L*S.diff(f,z).subs(z,-1)
 zero((X*T).subs(boundary,simultaneous=True)-shift,'boundary parameter')
 zero(X.subs(boundary,simultaneous=True)+1/shift,'boundary x')
 zero(T.subs(boundary,simultaneous=True)+shift*shift,'boundary t')
for f in fs:
 for g in fs:
  X,T,_=flow(f,L);A,B,_=flow(g,M);C,E,_=flow(L*f+M*g,1)
  zero(sub(A,X,T)-C,'composition x')
  zero(sub(B,X,T)-E,'composition t')
for g in [S.Integer(1),z+1,z*z+2*z+3]:
 f=z*z*(z+1)*g
 R=S.cancel(f.subs(z,u)/(t*t*q))
 zero(R-x**4*g.subs(z,u),'source carrier factor')
 check(all(2*j-i==-4 for (i,j),c in S.Poly(R,x,t).terms()),'outside Newton support')
 zero(S.diff(f,z).subs(z,-1)-g.subs(z,-1),'boundary velocity')
f=z*z*(z+1)**2
X,T,D=flow(f,L)
zero(D.subs({x:-1/s,t:-s*s},simultaneous=True)-1,'pointwise q denominator')
zero(X.subs({x:-1/s,t:-s*s},simultaneous=True)+1/s,'pointwise q x')
zero(T.subs({x:-1/s,t:-s*s},simultaneous=True)+s*s,'pointwise q t')
zero(flow(z*z-z*z,1)[0]-x,'identity cancellation x')
zero(flow(z*z-z*z,1)[1]-t,'identity cancellation t')
print(json.dumps({'status':'PASS','checks':checks,'polynomials':[str(f) for f in fs],
 'scope':'Exact finite algebraic controls; universal classification follows from the written proof.'},indent=2,sort_keys=True))
