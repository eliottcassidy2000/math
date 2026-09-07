#!/usr/bin/env python3
"""Exact curve and embedded-resolution controls for the Nori application.

One declared curve, no coefficient census or numerical braid inference.
All gates remain active under optimized Python.
"""
import sympy as S

GATES = 0


def check(value, description):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(description)


def equal(left, right, description):
    check(S.cancel(left-right) == 0, description)


def jet(f, variable):
    num, den = S.fraction(S.cancel(f))
    if num == 0:
        raise ValueError('zero function has no leading jet')
    pn, pd = S.Poly(num, variable), S.Poly(den, variable)
    en = min(m[0] for m, c in pn.terms() if c)
    ed = min(m[0] for m, c in pd.terms() if c)
    cn = pn.nth(en)
    cd = pd.nth(ed)
    return en-ed, S.cancel(cn/cd)


t, s, p, q, u, v, A, w = S.symbols('t s p q u v A w')
U=t**4+t**2
V=t**6+t**5+t**2
N=S.cancel((U.subs(t,s)-U)/(s-t))
M=S.cancel((V.subs(t,s)-V)/(s-t))
T=S.diff(U,t).subs(t,s)*S.diff(V,t)-S.diff(U,t)*S.diff(V,t).subs(t,s)
equal(N, (s+t)*(s*s+t*t+1), 'first divided difference')
equal(M.subs(t,-s),s**4,'pair-sum-zero local contribution')
check(S.gcd(S.diff(U,t),S.diff(V,t)) == t,'unique critical parameter')
check(S.gcd(U,V) == t**2,'critical image has no other parameter')

H=-(p-1)*(p**4+2*p**3+4*p**2+8*p+1)
q0=(p*p+1)/2
pair=s*s-p*s+q0
equal(S.rem(4*M.subs(t,p-s)-H,pair,s),0,'pair-sum quotient')
equal(H.subs(p,0),1,'no pair-sum-zero off-diagonal points')
equal(S.discriminant(H,p),-11776000,'five simple pair sums')
equal(S.resultant(H,p*p+2,p),123,'all pairs have distinct parameters')
check(S.groebner([N,M,T],s,t,domain=S.QQ) ==
      S.groebner([s+t,t**4],s,t,domain=S.QQ),'transversality outside cusp')

J=(t*t-t+1)*(2*t**8+4*t**7+8*t**6+6*t**5+12*t**4+
             10*t**3+15*t*t+5*t+5)
equal(S.resultant(N,M,s),-t**4*J,'independent parameter projection')
equal(S.discriminant(J,t),-166571520000000000,'ten distinct off-diagonal parameters')
check(S.degree(J,t)==10 and J.subs(t,0)==5,'projection count and cusp separation')
equal(S.rem(U+1,t*t-t+1,t),0,'extra node first coordinate')
equal(S.rem(V-1,t*t-t+1,t),0,'extra node second coordinate')

# A triple common fibre forces this cubic to divide U-u.  All signs and
# coefficients are retained independently of the pair projection.
remainder=-t**3+(u+2)*t*t+u*t-u-v
equal(S.rem(V-v,U-u,t),remainder,'true common-fibre remainder')
factor=(t+A)*(t**3-A*t*t-u*t+u+v)
coefs=S.Poly(S.expand(factor-(U-u)),t).all_coeffs()
check(S.groebner(coefs+[A-u-2],A,u,v,domain=S.QQ) ==
      S.groebner([1],A,u,v,domain=S.QQ),'no cubic common divisor')
f=A*A+A-1
j=(A-2)*(A*A+1)
equal((A*A-A-1)*f-(A+2)*j,5,'short triple-fibre Bezout identity')

# Full implicit equation and all bad U-fibres, retained as a braid sidecar.
F=(u**6+2*u**5-10*u**4*v+10*u**4-2*u**3*v**2-22*u**3*v+
   11*u**3+10*u**2*v**2-28*u**2*v+5*u**2+6*u*v**3+
   13*u*v**2-10*u*v+v**4+4*v**3+5*v**2)
equal(S.resultant(U-u,V-v,t),F,'implicit resultant')
node_u=4*u**4+40*u**3+120*u*u+125*u+25
equal(S.discriminant(F,v),-16*u**5*(u+1)**2*(4*u+1)**2*node_u**2,
      'complete vertical discriminant')
equal(S.resultant(H,4*u+(p*p+1)**2,p),256*(u+1)*node_u,
      'node first coordinates')

# Projective degree and the unique point at infinity: homogeneous map
# [S:T] -> [S²T⁴+S⁴T²:T⁶+ST⁵+S⁴T²:S⁶] has no common zero.
X=(w*w+w**4)/(1+w+w**4)
Z=w**6/(1+w+w**4)
equal(X,U.subs(t,1/w)/V.subs(t,1/w),'infinity first chart')
equal(Z,1/V.subs(t,1/w),'infinity line equation')
check(jet(X,w)==(2,1) and jet(Z,w)==(6,1),'infinity leading orders')
check(jet(Z-X**3,w)==(7,2),'infinity (2,7) characteristic')

# Finite cusp: each row is the branch in the new blowup chart.  E is
# vertical in rows1,2; row3 lies at two axes; row4 meets only y=0.
x=U
y=S.expand(V-U+U**2)
equal(y,t**5+3*t**6+t**8,'valid triangular target change')
finite=[(x,y),(x,y/x),(x,y/x**2),(x**3/y,y/x**2),
        (x**5/y**2-1,y/x**2)]
expected_finite=[(2,5),(2,3),(2,1),(1,1),(1,1)]
for row,expected in zip(finite,expected_finite):
    check(tuple(jet(a,t)[0] for a in row)==expected,'finite chart orders')
check(jet(finite[3][0]/finite[3][1],t)==(0,1),'finite tangent slope')
check(jet(finite[4][1],t)==(1,1),'finite final exceptional transversality')
equal(finite[1][0]*finite[1][1],finite[0][1],'finite blowup1')
equal(finite[2][0]*finite[2][1],finite[1][1],'finite blowup2')
equal(finite[3][0]*finite[3][1],finite[2][0],'finite blowup3')
equal((finite[4][0]+1)*finite[4][1],finite[3][0],'finite blowup4')

# Infinity: L_infinity is y=0 initially and through row2; row3 makes
# its strict transform y=-1, so it misses the subsequent centres.
rho=Z/X**3-1
infinite=[(X,Z),(X,Z/X),(X,Z/X**2),(X,rho),
          (X/rho,rho),(X/rho**2-S.Rational(1,4),rho)]
expected_infinite=[(2,6),(2,4),(2,2),(2,1),(1,1),(1,1)]
for row,expected in zip(infinite,expected_infinite):
    check(tuple(jet(a,w)[0] for a in row)==expected,'infinity chart orders')
check(jet(Z/X**3,w)==(0,1),'infinity diagonal tangent')
check(jet(infinite[4][0]/infinite[4][1],w)==(0,S.Rational(1,4)),
      'infinity penultimate tangent slope')
check(jet(infinite[5][1],w)==(1,2),'infinity final transversality')
equal(infinite[1][0]*infinite[1][1],infinite[0][1],'infinity blowup1')
equal(infinite[2][0]*infinite[2][1],infinite[1][1],'infinity blowup2')
equal(infinite[3][0]*(infinite[3][1]+1),infinite[2][1],'infinity blowup3')
equal(infinite[4][0]*infinite[4][1],infinite[3][0],'infinity blowup4')
equal((infinite[5][0]+S.Rational(1,4))*infinite[5][1],infinite[4][0],
      'infinity blowup5')

finite_mult=tuple(min(row) for row in expected_finite[:-1])
infinite_mult=tuple(min(row) for row in expected_infinite[:-1])
check(finite_mult==(2,2,1,1),'finite multiplicity sequence')
check(infinite_mult==(2,2,2,1,1),'infinity multiplicity sequence')
loss=sum(a*a for a in finite_mult+infinite_mult)
check(loss==24 and 6**2-loss==12,'strict transform self intersection')
check(12>2*5,'strict Nori inequality')
check(not 12>2*6,'threshold equality is insufficient hostile')
check(2+3+5==(6-1)*(6-2)//2,'independent genus inventory')

print('FINITE-EXACT single target curve and complete embedded-resolution controls')
print('U=t^4+t^2; V=t^6+t^5+t^2')
print('finite singularities: one (2,5) cusp, five ordinary nodes; no triple fibre')
print('infinity: unique (2,7) cusp, L_infinity retained in the resolution')
print('H(p)='+str(H))
print('implicit F(u,v)='+str(F))
print('vertical discriminant='+str(S.factor(S.discriminant(F,v))))
print('finite chart (orders,leading coefficients)='+str(
    [[jet(a,t) for a in row] for row in finite]))
print('infinity chart (orders,leading coefficients)='+str(
    [[jet(a,w) for a in row] for row in infinite]))
print('centre multiplicities='+str((finite_mult,infinite_mult)))
print('D^2=36-10-14=12; 2r(D)=10; strict margin=2')
print('Nori theorem, complement identifications and monodromy consequence are analytic')
print('PASS always-active gates='+str(GATES))
