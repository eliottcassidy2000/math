#!/usr/bin/env python3
"""Independent polynomial-jet audit of the beta!=0 quartic/sextic argument.

No producer imports. This checks exact identities supporting the written
all-degree proof, not a finite search for Keller maps.
"""
import json
import sympy as s

checks = 0


def zero(expr, label):
    global checks
    checks += 1
    if s.cancel(s.expand(expr)) != 0:
        raise RuntimeError(label)


def require(value, label):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)


a,b,d,at,bt,dt,k,E,H,L,M,y,z=s.symbols('a b d at bt dt k E H L M y z')
D=lambda p: s.diff(p,a)*at+s.diff(p,b)*bt+s.diff(p,d)*dt
c=s.Rational(3,2)*a
B=(3*b+k)/2
e=s.Rational(3,2)*d+s.Rational(3,8)*a*a+E
h=s.Rational(3,4)*a*b+s.Rational(3,8)*k*a+H
f=s.Rational(3,4)*a*d+s.Rational(3,8)*b*b+s.Rational(3,8)*k*b-s.Rational(1,16)*a**3+E*a/2
P=y**4+a*y*y+b*y+d
Q=y**6+c*y**4+B*y**3+e*y*y+h*y+f
J=s.expand(D(P)*s.diff(Q,y)-s.diff(P,y)*D(Q))
I=3*b*d+3*k*d/2-3*a*a*b/4-3*k*a*a/16+2*E*b+H*a
K=3*d*d/2-3*a*a*d/4+2*E*d-3*a*b*b/4-3*k*a*b/8+H*b+3*a**4/32-E*a*a/2
zero(J-(D(I)*y*y+D(K)*y+h*dt-b*D(f)), 'whole Jacobian jet identity')
for j in range(3,8):
    zero(s.Poly(J,y).coeff_monomial(y**j),'top integrated row')

r=2*b+k
U=L+E*k
N=9*a*(k*k*a**3-16*b*r**3)+96*k*a*a*(H*a-U)+256*(H*a-U)**2+r*r*(384*H*b-256*E*E-384*M)
ds=a*a*(4*b+k)/(8*r)-2*E/3-2*H*a/(3*r)+2*U/(3*r)
zero((I-L).subs(d,ds),'solve first integral exactly')
zero(384*r*r*(K-M).subs(d,ds)-N,'independent denominator-cleared elimination')

# Constant-a chart, including a=0 and arbitrary constants E,H,L.
zero(ds-(a*a/4-2*E/3+(-k*a*a/8-2*H*a/3+2*U/3)/r),
     'constant plus constant-over-r form')
zero(s.diff(K,b)-(H-s.Rational(3,2)*a*b-s.Rational(3,8)*k*a),
     'constant-a,d branch derivative')
zero((h*dt-b*D(f)).subs({a:0,H:0,at:0,dt:0,d:-2*E/3})+
     s.Rational(3,8)*b*(2*b+k)*bt,'nonconstant-b scalar branch')

# The only candidate leading terms for positive degrees A,B are a^4 and ab^4.
np=s.Poly(s.expand(N),a,b)
require(np.coeff_monomial(a**4)==9*k*k,'a^4 coefficient')
require(np.coeff_monomial(a*b**4)==-1152,'ab^4 coefficient')
remaining=[monomial for monomial,coef in np.terms()
           if monomial not in ((4,0),(1,4))]
# Every remaining monomial is dominated componentwise by one of those two,
# with at least one strict coordinate: valid for all positive degrees A,B.
require(all((i<=4 and j==0 and i<4) or
            (i<=1 and j<=4 and (i,j)!=(1,4)) or
            (i<=3 and j==0) for i,j in remaining),
        'all-positive-degree dominance certificate')
correction=s.expand(N-9*a*(k*k*a**3-16*b*r**3))
require(max(4*i+3*j for (i,j),coef in s.Poly(correction,a,b).terms())==12,
        'correction weight at most 12m')

b0=k*k*z**3/16
a0=z*(2*b0+k)
d0=z*z*(2*b0+k)*(4*b0+k)/8
replacement={a:a0,b:b0,d:d0,E:0,H:0,L:0,M:0}
zero(I.subs(replacement),'comparison first integral')
zero(K.subs(replacement),'comparison second integral')
h0=h.subs(replacement)
f0=f.subs(replacement)
j0=s.factor(h0*s.diff(d0,z)-b0*s.diff(f0,z))
claimed=3*k**4*z*z*(k*z**3+8)*(10*k*k*z**6+61*k*z**3+64)/16384
zero(j0-claimed,'direct comparison determinant')
require(s.Poly(j0,z).degree()==11,'comparison derivative degree')
zero(s.Poly(j0,z).LC()-s.Rational(15,8192)*k**7,
     'comparison leading coefficient')

# Independent chain-rule controls include finite poles; no polynomial-chart
# assumption enters either identity.
x=s.symbols('x')
for zz in (x, x*x+1, (x*x+1)/x, (x**4+1)/(x+1)):
    aa,bb,dd=[v.subs({z:zz,k:2}) for v in (a0,b0,d0)]
    pp=P.subs({a:aa,b:bb,d:dd})
    qq=Q.subs({a:aa,b:bb,d:dd,k:2,E:0,H:0})
    actual=s.diff(pp,x)*s.diff(qq,y)-s.diff(pp,y)*s.diff(qq,x)
    expected=claimed.subs({z:zz,k:2})*s.diff(zz,x)
    zero(actual-expected,'full source determinant with rational chart')
    numer,denom=s.fraction(s.cancel(expected))
    m=s.degree(s.fraction(s.cancel(zz))[0],x)-s.degree(s.fraction(s.cancel(zz))[1],x)
    require(s.degree(numer,x)-s.degree(denom,x)==12*m-1,
            'rational degree-at-infinity control')

print(json.dumps({
    'status':'FINITE-EXACT supporting independent proof audit',
    'checks':checks,
    'sympy_version':s.__version__,
    'producer_imports':False,
    'universe':'formal polynomial jet ring; four exact rational comparison charts',
    'passed':True
},indent=2,sort_keys=True))
