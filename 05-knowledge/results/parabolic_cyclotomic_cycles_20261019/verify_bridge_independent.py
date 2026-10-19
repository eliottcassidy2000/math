#!/usr/bin/env python3
"""Independent generic dynatomic, resultant, and finite-set bridge checks."""
import sympy as S
import json
x,c,T,L,t,z=S.symbols('x c T L t z')
def require(ok,message):
 if not ok:raise ValueError(message)
def same(a,b):return S.cancel(a-b)==0
f=x*x+c
f2=S.expand(f.subs(x,f))
f3=S.expand(f.subs(x,f2))
phi3,rem=S.div(f3-x,f-x,x)
require(rem==0,'dynatomic division')
G=x**3-T*x*x-(T*T+2*T+3)*x+T**3+2*T*T+3*T+1
ct=-T*T-T-2
require(same(phi3.subs(c,ct),G*G.subs(T,-1-T)),'trace cubics factor dynatomic')
require(S.rem(G.subs(x,f.subs(c,ct)),G,x)==0,'map preserves cubic roots')
require(S.factor(S.discriminant(G,x))==(4*T*T+6*T+9)**2,'trace cubic discriminant')
mult=S.rem((8*x*f*f2).subs(c,ct),G,x)
require(same(mult,8*(ct*(T+1)+1)),'multiplier constant on trace cubic')
multpoly=S.resultant(T*T+T+c+2,L-8*(c*(T+1)+1),T)
expected=L*L-(8*c+16)*L+64*(c**3+2*c*c+c+1)
require(same(multpoly,expected),'multiplier resultant')
e=S.symbols('e',nonnegative=True)
for sign in (-1,1):
 value=1-4*e+sign*(14+8*e)*S.sqrt(e)
 require(S.simplify(expected.subs({c:-S.Rational(7,4)-e,L:value}))==0,'square root branches')
g=t**3+t*t-2*t-1
cyclo=sum(z**i for i in range(7))
require(S.factor(S.resultant(cyclo,t*z-z*z-1,z))==g*g,'cyclotomic resultant')
require(S.discriminant(g,t)==49 and g.subs(t,1)!=0 and g.subs(t,-1)!=0,'irreducible real cubic')
sigma=t*t-2
sigma2=S.rem(sigma.subs(t,sigma),g,t)
sigma3=S.rem(sigma.subs(t,sigma2),g,t)
require(sigma2==-t*t-t+1 and sigma3==t,'exact order-three Galois action')
require(S.rem(sigma-t,g,t)!=0,'nonidentity Galois action')
pointpoly=8*x**3+4*x*x-18*x-1
require(same(pointpoly,-8*g.subs(t,-x-S.Rational(1,2))),'affine root polynomial')
require(S.rem(((-t-S.Rational(1,2))**2-S.Rational(7,4))-(-sigma2-S.Rational(1,2)),g,t)==0,'orientation reversed')
require(same(phi3.subs(c,-S.Rational(7,4)),pointpoly**2/64),'full parabolic dynatomic factor')
require(S.gcd(pointpoly,x*x-S.Rational(7,4)-x)==1,'parabolic roots have exact period three')
require(S.rem((8*x*f*f2).subs(c,-S.Rational(7,4)),pointpoly,x)==1,'parabolic multiplier')
require(S.rem((8*x*f*f2).subs(c,-2),g.subs(t,x),x)==8,'cyclotomic multiplier')
# Recover the marked formulas from the parameter inverse t=x^2+x+c.
D=2*t*(t+1)
p0=(t**3+2*t*t+t+1)/D
C=S.cancel(t-p0*p0-p0)
p1=S.cancel(p0*p0+C);p2=S.cancel(p1*p1+C)
require(same(p2*p2+C,p0),'reconstructed marked cycle')
require(same(C+S.Rational(7,4),-g*g/(4*t*t*(t+1)**2)),'quadratic parameter unfolding')
lam=S.cancel(8*p0*p1*p2)
require(same(lam-1,-(t*t+t+1)**3*g/(t**3*(t+1)**3)),'multiplier unfolding')
require(S.gcd(g,S.diff(g,t))==1 and S.gcd(g,t*(t+1))==1,'simple roots away from poles')
require(S.gcd(g,t*t+t+1)==1,'nonzero unfolding coefficient')
control=[S.cancel(a.subs(t,S.Rational(5,4))) for a in (C,p0,p1,p2,lam)]
require(control[0]==-S.Rational(226801,129600) and control[-1]==S.Rational(5605019,5832000),'attracting control')
require(0<control[-1]<1,'strictly attracting')
# Exact neutral period-two and fixed-point boundary checks.
u=S.symbols('u');v=-1-u
require(S.factor(4*u*v-1)==-(2*u+1)**2,'period-two multiplier plus one degenerates')
require(S.discriminant(S.together(4*u*v+1),u)==32,'period-two multiplier minus one irrational')
require(S.Rational(1,2)-S.Rational(1,2)**2==S.Rational(1,4),'neutral plus fixed point')
require(-S.Rational(1,2)-S.Rational(1,2)**2==-S.Rational(3,4),'neutral minus fixed point')
print(json.dumps(dict(status='PASS',method='Independent generic dynatomic division, resultant identities, modular polynomial iteration, and reconstructed parameter inverse',source_cycle_multiplier=8,target_cycle_multiplier=1,affine_map_preserves_multiplier=False,control=[str(a) for a in control],scope='Symbolic claims verified independently; arbitrary-period neutral exclusion follows from the separately verified cited denominator theorem.'),indent=2))
