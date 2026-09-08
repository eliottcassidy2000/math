#!/usr/bin/env python3
"""Exact interface checks for quartic pole closure with base-only Jacobian.

The unbounded pole contradictions are inherited proof dependencies, not
inferred from the finite Faber degrees checked here.
"""
from hashlib import sha256
import json
import sympy as S

x,z,w,r,b,t=S.symbols('x z w r b t')
p,q,c,dp,dq,dc=S.symbols('p q c dp dq dc')
V,beta,gamma,delta,epsilon,a,h=S.symbols('V beta gamma delta epsilon a h',nonzero=True)
gates=0

def check(value,label):
 global gates
 gates+=1
 if not value:raise RuntimeError(label)

def eq(left,right,label):check(S.cancel(left-right)==0,label)

def J(F,G,X=x,Z=z):return S.diff(F,X)*S.diff(G,Z)-S.diff(F,Z)*S.diff(G,X)

def root(P,Z,leading):
 P=S.Poly(P,Z)
 bb=P.coeff_monomial(Z**3);cc=P.coeff_monomial(Z**2)
 return leading*Z**2+bb*Z/(2*leading)+cc/(2*leading)-bb**2/(8*leading**3)

P=V*V*z**4+beta*z**3+gamma*z*z+delta*z+epsilon
H=root(P,z,V)
eq(H,V*z*z+beta*z/(2*V)+(4*gamma*V*V-beta*beta)/(8*V**3),'canonical rational root')
check(S.degree(S.cancel(P-H*H),z)<=1,'exact linear remainder')
Pchanged=S.expand(P.subs(z,a*w+h))
eq(root(Pchanged,w,V*a*a),H.subs(z,a*w+h),'all affine fibre root transport')

# Leading coefficient equation, without assuming a constant base Jacobian.
A=S.Function('A')(x);g=S.Function('g')(x)
for n in range(1,13):
 eq(J(A*z**4,g*z**n),
    (n*S.diff(A,x)*g-4*A*S.diff(g,x))*z**(n+3),'actual highest fibre coefficient')
 eq(S.diff(g**4/A**n,x),
    g**3/A**(n+1)*(4*A*S.diff(g,x)-n*S.diff(A,x)*g),'constant leading ratio')

# Exact Faber observables, with coefficients and their derivatives independent.
P0=z**4+p*z*z+q*z+c

def deriv(F):return S.diff(F,p)*dp+S.diff(F,q)*dq+S.diff(F,c)*dc

faber=[]
for n in range(1,13):
 alpha=S.Rational(n,4)
 ak=[S.Integer(1)]
 for k in range(n+3):
  value=0
  for j,v in [(2,p),(3,q),(4,c)]:
   if k-j+1>=0:value+=(alpha*j-(k-j+1))*v*ak[k-j+1]
  ak.append(S.expand(value/(k+1)))
 E=sum(ak[j]*z**(n-j) for j in range(n+1))
 Phi=4*ak[n+1];Psi=4*ak[n+2];Theta=4*ak[n+3]+p*ak[n+1]
 left=(dp*z*z+dq*z+dc)*S.diff(E,z)-S.diff(P0,z)*deriv(E)
 right=(z*z+p/4)*deriv(Phi)+z*deriv(Psi)+deriv(Theta)
 eq(left,right,'full Faber differential identity')
 eq(S.Poly(left,z).coeff_monomial(z*z),deriv(Phi),'first flux derivative is z squared coefficient')
 eq(S.Poly(left,z).coeff_monomial(z),deriv(Psi),'second flux derivative is z coefficient')
 faber.append((n,str(Phi),str(Psi)))

# Nonconstant base RHS, including genuine zeros where leading coefficients vanish.
for M,N,C in [(x,1,x**3),(x*x+1,x-2,(x-1)**4),(x**3,x*x,x*x+1)]:
 L=M*z+N
 for degree in [1,2]:
  if degree==1:
   FF=L**4+C;GG=L
  else:
   HH=L*L+C;FF=HH*HH+L;GG=-HH
  eq(J(FF,GG),M*S.diff(C,x),'base-only RHS with actual zero divisors')
  rr=S.cancel(root(FF,z,M*M))
  check(S.Poly(rr,x,z).total_degree()>=2,'root is an actual polynomial')
  check(S.degree(S.expand(FF-rr*rr),z)<=1,'positive canonical square prefix')
  for j in [1,2,3]:
   eq(J(FF,GG+FF**j),J(FF,GG),'raw mate degrees increase by removable global shears')

# Base-only hypothesis is necessary for this theorem.
Pbad=x*x*z**4+z*z
Hbad=root(Pbad,z,x)
eq(Hbad,x*z*z+1/(2*x),'finite-pole hostile canonical root')
eq(J(Pbad,x*z*z),-2*z**3,'fibre-dependent Jacobian hostile')
check(S.denom(S.cancel(Hbad)).subs(x,0)==0,'hostile pole is genuine')

# Exact DG chart and the parent-proposed global-first/root-nonglobal witness.
x2=1/r;t2=-r*r-r**4*b
F=t**4+2*x**4*t*t
H0=t*t+x**4
FB=S.expand(F.subs({x:x2,t:t2},simultaneous=True))
HB=S.expand(H0.subs({x:x2,t:t2},simultaneous=True))
eq(FB,r**8*(1+r*r*b)**4+2*(1+r*r*b)**2,'global quartic witness in full boundary chart')
eq(HB,r**4*(1+r*r*b)**2+r**(-4),'actual root pole order four')
eq(F-H0*H0,-x**8,'witness remainder degree zero and nonglobal')
eq(root(FB,b,r**8),HB,'boundary canonical root of witness')
eq(S.limit(r**4*HB,r,0),1,'nonremovable boundary principal coefficient')
eq(S.diff(F,x).subs(t,0),0,'witness critical line first derivative')
eq(S.diff(F,t).subs(t,0),0,'witness critical line second derivative')
coeff=S.Poly(FB,b);vv=r**8;be=coeff.coeff_monomial(b**3);ga=coeff.coeff_monomial(b*b)
eq(4*ga*vv*vv-be*be,8*r**20*(1+r**8),'exact boundary pole-congruence failure')
# Change-of-chart Jacobian follows from the literal Jacobian matrix.
eq(J(x2,t2,r,b),r*r,'actual Jacobian multiplier r squared')
for f0,g0 in [(x,t),(x*x+t,t*t+x),(F,H0)]:
 eq(J(f0.subs({x:x2,t:t2},simultaneous=True),
      g0.subs({x:x2,t:t2},simultaneous=True),r,b),
    r*r*J(f0,g0,x,t).subs({x:x2,t:t2},simultaneous=True),'chain-rule boundary bracket')
# Genuine global square-prefix controls; no constant-Jacobian mate is claimed.
for HH,LL in [(t*t+x*x*t,t),(t+x*x*t*t,x*t),(x*t*(1+x*x*t),1+x*x*t)]:
 ff=S.expand(HH*HH+LL)
 lead=S.Poly(HH,t).coeff_monomial(t*t)
 eq(root(ff,t,lead),HH,'global square-prefix canonical root')
 for item in [HH,LL,ff]:
  pulled=S.cancel(item.subs({x:x2,t:t2},simultaneous=True))
  check(S.denom(pulled)==1,'full boundary regularity of positive prefix')

blob=json.dumps(faber,separators=(',',':')).encode()
print('Quartic base-Jacobian and DG boundary interfaces PASS; unbounded proof audit required')
print('Faber universe: degrees 1..12 with independent coefficient derivatives')
print('Positive base-only Jacobian families: 6, with vanishing base RHS and unbounded-shear pattern')
print('Controls: finite-pole hostile, global quartic/root-pole hostile, exact chart and affine-root transport')
print('Always-active gates:',gates)
print('Semantic SHA256:',sha256(blob).hexdigest())
