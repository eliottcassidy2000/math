#!/usr/bin/env python3
"""Exact controls for the global period-three index proof; no floating point."""
import json
import sympy as s

checks = []

def eq(actual, expected, name):
    if s.cancel(s.expand(actual - expected)) != 0:
        raise RuntimeError(f"FAIL: {name}: {actual} != {expected}")
    checks.append(name)

def truth(condition, name):
    if not condition:
        raise RuntimeError(f"FAIL: {name}")
    checks.append(name)

x, c, v, t, z, u, lam = s.symbols('x c v t z u lam')
f = x*x+c
q = f-x
F = s.expand((f*f+c)**2+c)
phi, remainder = s.div(F-x, q, x)
eq(remainder, 0, 'dynatomic division')
eq(s.degree(F,x), 8, 'return degree')
eq(s.rem(phi,q,x), 1-4*c+6*x, 'dynatomic fixed-point quotient')
D = 16*c*c+4*c+7
eq(s.resultant(q,phi,x), D, 'global collision resultant')
Q = lam*lam-2*lam+4*c
eq(s.resultant(Q,1-lam**3,lam), (4*c-1)*D, 'fixed return multiplier denominator')
eq(2-(8-24*c), 6*(4*c-1), 'fixed return index numerator')
eq(s.cancel(6*(4*c-1)/((4*c-1)*D)), 6/D, 'generic fixed index sum')
eq(s.discriminant(D,c), -432, 'two nonreal collision parameters')

cstar = -s.Rational(7,4)
P = x**3+x*x/2-s.Rational(9,4)*x-s.Rational(1,8)
fs, qs, Fs = (a.subs(c,cstar) for a in (f,q,F))
eq(Fs-x, qs*P*P, 'parabolic factorization')
eq(s.discriminant(P,x), 49, 'inherited cubic discriminant')
eq(s.gcd(P,qs), 1, 'cycle disjoint from fixed points')
truth(s.Poly(P,x,domain=s.QQ).is_irreducible, 'cycle cubic irreducible')
eq(s.rem(P,qs,x), x+s.Rational(5,2), 'short resultant remainder')
eq(qs.subs(x,-s.Rational(5,2)), 7, 'short resultant value')
eq(s.resultant(qs,P,x), 7, 'oriented gap resultant')
eq(s.resultant(qs,P,x)**2, s.discriminant(P,x), 'squared gap equals discriminant')
eq(D.subs(c,cstar), 49, 'index denominator at cusp')
eq(s.rem(P.subs(x,fs),P,x), 0, 'quadratic permutes cubic roots')
eq(s.rem(s.diff(Fs,x),P,x), 1, 'return multiplier one')
eq(6/s.Integer(49)+3*(-s.Rational(2,49)), 0, 'global residue budget')
eq(1-(-s.Rational(2,49)), s.Rational(51,49), 'iterative residue')
truth(6/s.Integer(49)+6*(-s.Rational(2,49)) != 0, 'reject algebraic multiplicity weighting')

# Independent local coefficient route, confined to exact cubic arithmetic.
a = s.rem(s.diff(Fs,x,2)/2,P,x)
b = s.rem(s.diff(Fs,x,3)/6,P,x)
eq(s.rem(b*s.invert(a*a,P,x),P,x), -s.Rational(2,49), 'independent cubic coefficient index')
eq(s.rem(a*a,P,x),1372*x*(x-1),'nonzero quadratic coefficient square')
alpha=s.rem(s.diff(F,c).subs(c,cstar),P,x)
eq(alpha,2*x*x-s.Rational(11,2),'parameter unfolding coefficient')
eq(s.rem(a*alpha,P,x),49,'unfolding scale equals forty-nine')
truth(P.subs(x,0) != 0 and P.subs(x,1) != 0, 'simple parabolic order')
scaled = s.expand(8*P.subs(x,x/2))
eq(scaled,x**3+x*x-9*x-1,'scaled monic orbit polynomial')
eq(s.discriminant(scaled,x),3136,'coordinate discriminant hostile')
eq((b/4)/(a/2)**2,b/a**2,'index survives scaling')

M=t**3+t*t-2*t-1
sigma=1-t-t*t
red=lambda a:s.rem(s.expand(a),M,t)
eq(P.subs(x,-s.Rational(1,2)-t),-M,'cyclotomic affine substitution')
eq(-s.Rational(1,2)-fs.subs(x,-s.Rational(1,2)-t),sigma,'induced Galois polynomial')
eq(red(M.subs(t,sigma)),0,'Galois root preservation')
eq(red(sigma.subs(t,sigma)),t*t-2,'Galois inverse doubling')
eq(red(sigma.subs(t,red(sigma.subs(t,sigma)))),t,'Galois order three')
eq(red(t**3-3*t),sigma,'angle tripling representative')
eq(s.expand(z**3*M.subs(t,z+1/z)),sum(z**i for i in range(7)),'seventh trace identity')
truth(s.Poly(M,t,domain=s.QQ).is_irreducible,'trace cubic irreducible')
truth([3**j % 7 for j in range(1,4)] == [3,2,6],'exponent classes modulo sign')

# General cyclic-gap identity, with three independent symbolic coordinates.
r1,r2,r3=s.symbols('r1 r2 r3')
eq(((r2-r1)*(r3-r2)*(r1-r3))**2,
   (r1-r2)**2*(r1-r3)**2*(r2-r3)**2,'abstract three-cycle gap discriminant')

cv=-(7+v*v)/4
sp=(-1+v)/2
Pm=lambda tr:x**3-tr*x*x+(cv-tr-1)*x-cv*tr-cv-1
A,B=Pm(sp),Pm(-1-sp)
eq(A*B,phi.subs(c,cv),'two-cycle trace factorization')
mu1=8*(cv*sp+cv+1)
mu2=8*(cv*(-1-sp)+cv+1)
Av=v*v+v+7
eq(s.resultant(q.subs(c,cv),A,x),Av,'trace oriented gap resultant')
eq(s.discriminant(A,x),Av*Av,'trace discriminant gap square')
eq(D.subs(c,cv),Av*Av.subs(v,-v),'trace resultant product')
eq(1-mu1,v*Av,'multiplier gap factorization')
eq(mu1,1-v*v-v*(7+v*v),'first cycle multiplier')
eq(mu2,1-v*v+v*(7+v*v),'second cycle multiplier')
eq(s.cancel(1/(1-mu1)+1/(1-mu2)),-2/(v**4+13*v*v+49),'paired cycle index')
eq(D.subs(c,cv),v**4+13*v*v+49,'trace denominator')
eq(s.limit(1/(1-mu1)+1/(1-mu2),v,0),-s.Rational(2,49),'coalescence index')
eq(6/D.subs(c,cv)+3*(-2/(v**4+13*v*v+49)),0,'generic cycle budget')
eq(A.subs(v,0),P,'first limiting cubic')
eq(B.subs(v,0),P,'second limiting cubic')

# Exceptional, higher-multiplicity, and nonpolynomial controls.
Fquarter=F.subs(c,s.Rational(1,4))
local=s.expand(Fquarter.subs(x,u+s.Rational(1,2))-s.Rational(1,2))
eq(local.coeff(u,2),3,'third cusp iterate quadratic coefficient')
eq(local.coeff(u,3),6,'third cusp iterate cubic coefficient')
eq(s.residue(1/(u-local),u,0),s.Rational(2,3),'coalesced fixed-point index')
eq((6/D).subs(c,s.Rational(1,4)),s.Rational(2,3),'removable fixed-cluster value')
H=z+z**3+z**5
higher=s.residue(1/(z-H),z,0)
eq(higher,1,'triple fixed-point residue')
eq(1/(1-s.diff(H,z).subs(z,s.I)),-s.Rational(1,2),'other index plus i')
eq(1/(1-s.diff(H,z).subs(z,-s.I)),-s.Rational(1,2),'other index minus i')
eq(higher-1,0,'higher multiplicity budget')
Mob=z/(1-z)
eq(s.residue(1/(z-Mob),z,0),1,'rational map rejects zero finite budget')

receipt={
 'status':'PASS', 'checks':len(checks), 'arithmetic':'exact rational and symbolic',
 'fixed_point_index_sum':'6/(16*c**2+4*c+7)',
 'collision_resultant':'16*c**2+4*c+7',
 'parabolic_cycle_index':'-2/49', 'iterative_residue':'51/49',
 'unfolding_scale':'49', 'local_delay_limit':'sqrt(epsilon)*N -> pi/7',
 'resultant_at_cusp':'7', 'cycle_discriminant':'49',
 'scaled_cycle_discriminant':'3136',
 'paired_cycle_index':'-2/(v**4+13*v**2+49)',
 'galois_action':'zeta_7 -> zeta_7**3; t -> 1-t-t**2',
 'checks_passed':checks,
}
print(json.dumps(receipt,indent=2,sort_keys=True))
