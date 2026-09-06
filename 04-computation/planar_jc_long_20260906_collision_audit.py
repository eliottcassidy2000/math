#!/usr/bin/env python3
"""Independent sparse-polynomial/AlgebraicField audit; no producer import.

The collision is solved by full truncated-series fixed-Jacobian contraction,
not the producer's coefficient-at-a-time Fraction quotient implementation.
"""
import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.ring_series import rs_trunc, rs_series_inversion

CHECKS = 0
ORDER = 5


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


x, q, a, b, c, r = sp.symbols('x q a b c r')
P = x*x*(x*x-1)**2
Q = x**5+sp.Rational(9,2)*x**4-2*x**3-sp.Rational(27,4)*x*x+x-sp.Rational(3,4)+P*(a+b*x+c*x*x)
D = 1+x*x*q
B, C, E = (D-1)*(D+2)**2, x*D*(D+2), q*(D+3)
check(sp.expand(C*C*E-B*(B+4)) == 0, 'surface identity')
check(sp.expand(sp.diff(C,x)*sp.diff(E,q)-sp.diff(C,q)*sp.diff(E,x)-6*(B+2)) == 0,
      'compiler Jacobian identity')

# Derive the implicit Jacobian from the six actual collision equations.
zm, z0, zp, cs, es, h, s = sp.symbols('zm z0 zp cs es h s')
unknown = [zm, z0, zp, cs, es, h]
equations = []
for z in [zm, z0, zp]:
    value_q = Q.subs(x,z)+s+h*z
    equations.extend([C.subs({x:z,q:value_q}, simultaneous=True)-cs,
                      E.subs({x:z,q:value_q}, simultaneous=True)-es])
base = {zm:-1,z0:0,zp:1,cs:0,es:-3,h:0,s:0}
J = sp.Matrix(equations).jacobian(unknown).subs(base).applyfunc(sp.expand)
check(J.det() == -288, 'derived full implicit Jacobian determinant')
check(not any(v.free_symbols for v in J), 'Jacobian independent of all a,b,c')
JI = J.inv()
f4 = 72783360*r**4-77822208*r**3-28419741*r**2+7849770*r-1276420
g3 = 104569906176*r**3+87403912326*r*r-19653044202*r+4049817787
aa, bb, gg = sp.gcdex(f4,g3,r)
check(sp.expand(aa*f4+bb*g3) == 1 and gg == 1, 'exact all-root Bezout certificate')
check(sp.Poly(f4,r).is_irreducible, 'quartic coefficient domain is a field')
print('SYMBOLIC independently derived Jacobian=-288; compiler determinant; all-root Bezout PASS')


def audit_case(name, domain, abc, first, expected):
    R, t = ring('t', domain)
    one, zero = R.one, R.zero

    def trunc(v, cutoff=ORDER+1):
        return rs_trunc(v,t,cutoff)

    def mul(v,w):
        return trunc(v*w)

    def power(v,n):
        out = one
        for _ in range(n):
            out = mul(out,v)
        return out

    def const(value):
        return R.ground_new(domain.convert(value))

    coefficients = [domain.convert(sp.Poly(Q.subs({a:0,b:0,c:0}),x).nth(i))
                    for i in range(9)]
    for j, value in enumerate(abc):
        for degree, coefficient in [(2,1),(4,-2),(6,1)]:
            coefficients[j+degree] += domain.convert(coefficient)*value
    derivative_coefficients = [domain.convert(i)*coefficients[i]
                               for i in range(1,len(coefficients))]

    def evaluate(coefficients,z):
        out = zero
        for value in reversed(coefficients):
            out = mul(out,z)+R.ground_new(value)
        return out

    def evaluate_bivariate(expression,z,qz):
        out = zero
        for (ix,iq), coefficient in sp.Poly(expression,x,q).terms():
            out += mul(power(z,ix),power(qz,iq))*domain.convert(coefficient)
        return trunc(out)

    def compiler(z,chi):
        qz = evaluate(coefficients,z)+t+mul(chi,z)
        dz = one+mul(power(z,2),qz)
        values = [mul(dz-one,power(dz+2,2)),
                  mul(mul(z,dz),dz+2),mul(qz,dz+3)]
        qx = evaluate(derivative_coefficients,z)+chi
        tangent = [evaluate_bivariate(sp.diff(F,x),z,qz)+
                   mul(evaluate_bivariate(sp.diff(F,q),z,qz),qx)
                   for F in [C,E]]
        velocity_q = [evaluate_bivariate(sp.diff(F,q),z,qz) for F in [C,E]]
        return values,tangent,velocity_q

    def residual(state):
        values = [compiler(z,state[5])[0] for z in state[:3]]
        return [value[k]-state[k+2] for value in values for k in [1,2]]

    state = [const(value) for value in [-1,0,1,0,-3,0]]
    for stage in range(ORDER):
        old_residual = residual(state)
        correction = [sum((old_residual[j]*domain.convert(JI[i,j])
                           for j in range(6)),zero) for i in range(6)]
        state = [trunc(v-w) for v,w in zip(state,correction)]
        check(all(trunc(v,stage+2) == zero for v in residual(state)),
              name+' Hensel contraction gains one order')
    check(all(v == zero for v in residual(state)), name+' full retained collision')
    chi = state[5]
    check(min(k[0] for k in chi) == first, name+' first compensator order')
    check(chi[(first,)] == expected, name+' inherited leading coefficient')
    values = [compiler(z,chi) for z in state[:3]]
    check(values[0][0][0] == values[1][0][0] == values[2][0][0], name+' common B section')
    tangents = [v[1] for v in values]

    def det(v,w):
        return mul(v[0],w[1])-mul(v[1],w[0])

    raw = [det(tangents[1],tangents[2]),det(tangents[2],tangents[0]),
           det(tangents[0],tangents[1])]
    check([v.get((0,),domain.zero) for v in raw] ==
          [domain.convert(v) for v in [15,-54,39]], name+' raw orientation')
    inverse = rs_series_inversion(-raw[1],t,ORDER+1)
    check(mul(-raw[1],inverse) == one, name+' sparse ring inverse')
    lam = [mul(v,inverse) for v in raw]
    pi = sum(lam,zero)
    xi = sum((mul(l,z) for l,z in zip(lam,state[:3])),zero)
    check(xi[(0,)] == domain.convert(sp.Rational(4,9)), name+' first moment unit')
    check(trunc(pi+mul(chi.diff(t),xi),ORDER) == zero, name+' exact identity through known derivative order')
    check(min(k[0] for k in pi) == first-1, name+' actual period order')
    expected_period = expected*domain.convert(sp.Rational(-4*first,9))
    check(pi[(first-1,)] == expected_period, name+' actual leading period')
    for target_slot in range(3):
        period = zero
        for i,(z,v) in enumerate(zip(state[:3],values)):
            density = one+mul(chi.diff(t),z)
            normal_velocity = [mul(vq,density) for vq in v[2]]
            pullback = det(v[1],normal_velocity) if target_slot == 0 else v[1][target_slot-1]
            period += mul(lam[i],pullback)
        check(trunc(period,ORDER) == zero, name+' full target two-form slot '+str(target_slot))
    for z,v in zip(state[:3],values):
        density = one+mul(chi.diff(t),z)
        pulled = det(v[1],[mul(vq,density) for vq in v[2]])
        divisor = 6*(v[0][0]+2)
        check(trunc(pulled-mul(divisor,density),ORDER) == zero,
              name+' normalized regular form gives unit density')
    for clock in [t,t*t,t**3,t+t*t,2*t*t-t**3]:
        clock_order = min(k[0] for k in clock)
        clock_lead = clock[(clock_order,)]
        # Use the complete retained period, at sufficient precision for its
        # first live composite term. Ordinary polynomial composition is exact.
        composite = sum((const(v)*clock**k[0] for k,v in pi.items()),zero)
        live = (first-1)*clock_order
        check(min(k[0] for k in composite) == live, name+' full period clock order')
        check(composite[(live,)] == expected_period*clock_lead**(first-1),
              name+' full period clock leading coefficient')
    print('STRATUM',json.dumps({'name':name,'chi_order':first,'period_order':first-1,
                               'collision_through':ORDER,'period_through':ORDER-1},
                              separators=(',',':')))


QQ = sp.QQ
audit_case('generic',QQ,[QQ.zero]*3,2,QQ.convert(sp.Rational(-259,9)))
audit_case('plane',QQ,[QQ.convert(sp.Rational(-259,36)),QQ.zero,QQ.zero],3,
           QQ.convert(sp.Rational(16*5717,2187)))
p0 = sp.Rational(-5717,729)
audit_case('parabola',QQ,[QQ.convert(sp.Rational(-259,36)+p0),QQ.zero,QQ.convert(-p0)],4,
           QQ.convert(sp.Rational(-16*1276420,177147)))
K = QQ.algebraic_field((sp.Poly(f4,r),sp.Symbol('alpha')))
alpha = K.unit
pp = K.convert(sp.Rational(520,9))*alpha**2-K.convert(sp.Rational(1688,81))*alpha-K.convert(sp.Rational(5717,729))
gg = K.convert(104569906176)*alpha**3+K.convert(87403912326)*alpha**2-K.convert(19653044202)*alpha+K.convert(4049817787)
audit_case('quartic-field',K,[K.convert(sp.Rational(-259,36))+pp+4*alpha,-9*alpha,-pp],5,
           K.convert(sp.Rational(64,3**17))*gg)

# Formal inverse of the source coordinate escape, independently by Catalan
# reversion of X=x+h*x²/2. Its h-adic coefficients are polynomial in X.
X,H = sp.symbols('X H')
inverse = sum((-1)**n*sp.catalan(n)*H**n*X**(n+1)/2**n for n in range(7))
check(sp.series(inverse+H*inverse**2/2-X,H,0,7).removeO().expand() == 0,
      'source coordinate inverse through h6')
check(sp.diff(x+H*x*x/2,x) == 1+H*x, 'exact source coordinate Jacobian')
l1,l2,l3,z1,z2,z3 = sp.symbols('l1 l2 l3 z1 z2 z3')
check(sp.expand(sum(l*(1+H*z) for l,z in zip([l1,l2,l3],[z1,z2,z3]))-
                (l1+l2+l3)-H*(l1*z1+l2*z2+l3*z3)) == 0,
      'transported relation sum')
print('SOURCE_ESCAPE exact Jacobian and relation covariance; Catalan inverse through h6 PASS')
print('PASS checks='+str(CHECKS))
print('SOURCE_SHA256 '+hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
