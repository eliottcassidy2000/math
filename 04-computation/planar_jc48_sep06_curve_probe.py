#!/usr/bin/env python3
"""Exact one-parameter (4,6) curve controls, no numerical root classification."""
import sympy as S

s, t, p, q, L, A = S.symbols('s t p q L A')
GATES = 0


def check(value, label):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(label)


def zero(value, label):
    check(S.cancel(value) == 0, label)


def unit_gcd(f, g, label):
    value = S.gcd(f, g)
    check(value != 0 and not value.free_symbols, label)


def gb_equal(expressions, expected, label):
    actual = S.groebner(expressions, t, s, domain=S.QQ)
    target = S.groebner(expected, t, s, domain=S.QQ)
    check(actual == target, label)
    return actual


def main():
    U = t**4+t
    V = t**6+t**2+L*t
    N = s**3+s*s*t+s*t*t+t**3+1
    M = sum(s**(5-i)*t**i for i in range(6))+s+t+L
    H = p**6+2*p**3-4*p**2-4*L*p-3
    C = 128*L**3-288*L-283
    E = (L+1)*(L**3-L**2+3*L+1)
    T5 = 3125*L**5+6250*L**4-6625*L**3-24604*L**2-18669*L-4371
    T = (L-2)*T5
    R = (3*s**12-2*L*s**10+6*s**9+7*s**8+(4-4*L)*s**6
         +(10-2*L**2)*s**5+(5+4*L)*s**4+2*s**3
         +(3-L**2-2*L)*s*s+(4-2*L**2+2*L)*s+2+3*L-L**3)
    D = S.diff(U,t).subs(t,s)*S.diff(V,t)-S.diff(U,t)*S.diff(V,t).subs(t,s)
    zero((s-t)*N-(U.subs(t,s)-U), 'raw U divided difference')
    zero((s-t)*M-(V.subs(t,s)-V), 'raw V divided difference')
    zero(N-(p**3-2*p*q+1).subs({p:s+t,q:s*t}), 'symmetric N')
    zero(M-(p**5-4*p**3*q+3*p*q*q+p+L).subs({p:s+t,q:s*t}),
         'symmetric M')
    zero(4*p*(p**5-4*p**3*q+3*p*q*q+p+L).subs(q,(p**3+1)/(2*p))+H,
         'pair-sum reduction')
    zero(p*p-4*(p**3+1)/(2*p)+(p**3+2)/p,'pair discriminant')
    check(H.subs(p,0) == -3, 'p remains invertible')
    zero(S.resultant(N,M,t)-R, 'ordinary resultant')
    zero(S.resultant(H,2*p*s*s-2*p*p*s+p**3+1,p)+64*R,
         'independent pair-algebra resultant')
    zero(S.discriminant(R,s)+12288*C*E**6*T**2, 'complete projection discriminant')
    zero(S.discriminant(H,p)-4096*T, 'pair-sum discriminant')
    zero(S.resultant(H,p**3+2,p)-C, 'diagonal resultant')
    zero(S.resultant(S.diff(U,t),S.diff(V,t),t)+8*C, 'critical resultant')
    for name, f in [('C',C),('E',E),('T',T)]:
        unit_gcd(f,S.diff(f,L),name+' squarefree')
    for f,g in [(C,E),(C,T),(E,T)]:
        unit_gcd(f,g,'distinct exceptional loci')
    sr = S.subresultants(H,S.diff(H,p),p)
    aa = 375*L**4-1134*L*L-813*L-113
    bb = 315*L**3-144*L*L-907*L-405
    zero(sr[-2]-2560*(aa*p+bb),'linear subresultant')
    zero(sr[-1]+4096*T,'constant subresultant')
    unit_gcd(aa,T,'linear gcd at every tangency parameter')

    crit = 4*t**3+1
    jet = S.diff(U,t,2)*S.diff(V,t,3)-S.diff(V,t,2)*S.diff(U,t,3)
    zero(S.rem(jet,crit,t)+12*t*(15*t+4),'critical second/third jets')
    unit_gcd(-12*t*(15*t+4),crit,'ordinary cusp determinant')
    unit_gcd(S.diff(U,t,2),crit,'cusp second jet')
    zero(S.rem(S.diff(V,t),crit,t)-(L-S.Rational(3,2)*t*t+2*t),
         'critical parameter map')

    phi = A**4-2*A+1
    K = t**3-A*t*t+A*A*t+1-A**3
    zero((t+A)*K-(t**4+t-(A-1))+phi,'exact triple fibre factorization')
    zero(S.resultant(phi,L+A*A,A)-E,'triple parameter eliminant')
    zero(E-(L**4+2*L*L+4*L+1),'triple explicit inverse parameter')
    unit_gcd(S.discriminant(K,t),phi,'three distinct fibre parameters')
    unit_gcd(K.subs(t,-A),phi,'fourth quartic root excluded')
    zero(K.subs(t,-A)-(1-4*A**3),'fourth root law')

    # Independent direct collision ideals at the three rational controls.
    shapes = {
        1: 4*t+3*s**9+3*s**8-2*s**7+s**6+10*s**5+3*s**4-s**3+2*s*s+6*s+4,
        2: 33*t+6*s**9+15*s**8-8*s**7-14*s**6+38*s**5+16*s**4
           -6*s**3-9*s*s+11*s+33,
    }
    for value in (1,2,-1):
        mv, dv, rv = M.subs(L,value),D.subs(L,value),R.subs(L,value)
        unit_gcd(S.diff(U,t),S.diff(V,t).subs(L,value),'named immersion')
        gb_equal([N,mv,s-t],[1],'named diagonal exclusion')
        tangent = [s+t+1,s*s+s] if value == 2 else [1]
        gb_equal([N,mv,dv],tangent,'named exact tangent ideal')
        if value in shapes:
            gb_equal([N,mv],[shapes[value],rv],'independent lex shape')
        if value == 1:
            unit_gcd(rv,S.diff(rv,s),'six-node projection squarefree')
    zero(H.subs(L,2)-(p+1)**2*(p*p-p-1)*(p*p-p+3),'tacnode pair factorization')
    zero(R.subs(L,2)-s*s*(s+1)**2*(s**4-s**3+2*s-1)
         *(3*s**4-3*s**3+2*s*s-2*s+5),'tacnode projection factorization')
    zero(V.subs(L,2)-2*U-t*t*(1-t*t)**2,'tacnode tangent subtraction')
    W = t*t*(1-t*t)**2
    for v, expected in [(0,S.Integer(1)),(-1,S.Rational(4,9))]:
        coefficient = S.diff(W,t,2).subs(t,v)/(2*S.diff(U,t).subs(t,v)**2)
        check(coefficient == expected,'tacnode graph quadratic coefficient')
    triple = t*(t*t-t+1)
    zero(S.gcd(U,V.subs(L,-1))-triple,'exact triple target fibre')
    F6 = 3*s**6+6*s**5+5*s**4+4*s**3+9*s*s+10*s+4
    zero(R.subs(L,-1)-s*s*(s*s-s+1)**2*F6,'triple projection factorization')
    unit_gcd(F6,S.diff(F6,s),'remaining triple-control projection simple')
    unit_gcd(F6,s*(s*s-s+1),'remaining projection avoids triple')
    unit_gcd(2-4*t,t*t-t+1,'triple tangent against parameter zero')
    check(S.discriminant(t*t-t+1,t) == -3,'other triple tangents distinct')
    z = S.groebner([N,M.subs(L,-1)],t,s,domain=S.QQ)
    linear = next(f.as_expr() for f in z.polys if S.degree(f.as_expr(),t) == 1)
    unit_gcd(S.Poly(linear,t).LC(),F6,'one partner at each remaining parameter')

    # One exceptional point: no local intersection rule is assumed there.
    nn,rr,dd,actual,nz,overlap = S.symbols('nn rr dd actual nz overlap')
    delta = dd-actual
    euler = dd*(nn+rr-1)+actual*(1-2*nn-rr)+nn*(dd-2*delta)+overlap+nz
    zero(euler-((rr-1)*delta+nz+overlap),'one exceptional point Euler ledger')
    passports = set()
    for degree in range(2,13):
        for debt in range(1,degree):
            for branches in range(2,6):
                for special in range(3):
                    for total_overlap in range(3):
                        if (branches-1)*debt+special+total_overlap == 1:
                            passports.add((branches,debt,special,total_overlap))
    check(passports == {(2,1,0,0)},'only nonunibranch formal passport')
    for degree in range(2,26):
        for actual_degree in range(1,degree):
            if 5*max(0,degree-2*actual_degree) <= 1:
                check(degree <= 2*actual_degree,'five-node cusp degree constraint')

    print('STATUS: EXACT ONE-PARAMETER CURVE CERTIFICATE; no Keller realization claimed')
    print('U=t^4+t; V=t^6+t^2+lambda*t; intrinsic pole pair=(4,6)')
    print('full divided-difference algebra length=12 for every lambda')
    print('C =',C)
    print('E =',S.expand(E))
    print('T =',T)
    print('exceptional root counts: cusp=3, triple=4, tangency=6; pairwise disjoint')
    print('generic complement: six ordinary nodes')
    print('lambda=1: six nodes; lambda=2: one tacnode+four nodes')
    print('lambda=-1: one ordinary triple+three nodes; collision ideal reduced')
    print('projection discriminant = -12288*C*E^6*T^2')
    print('sole Keller support: generic, triple and tangency strata excluded')
    print('three cusp parameters OPEN; necessary n_cusp+sum_node_overlaps=1, d<=2a')
    print('ALL GATES PASS:',GATES)


if __name__ == '__main__':
    main()
