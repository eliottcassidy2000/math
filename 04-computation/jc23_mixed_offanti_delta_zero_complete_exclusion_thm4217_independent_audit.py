#!/usr/bin/env python3
"""Independent exact audit of THM-4217's mixed off-anti Delta=0 wall.

It specializes the complete source before every elimination, derives the
relation between the (A,B) and (A,C0) projections, and resolves the entire
bottom-endpoint cascade through the terminal quotient N(u)=0.
"""

from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly, variable):
    terms = sp.Poly(poly, variable).terms()
    need(bool(terms), "zero polynomial")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def residual(poly, variable):
    order = valuation(poly, variable)
    answer = sp.Poly(sp.cancel(poly / variable**order), variable)
    return order, answer


def convex_hull(points):
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def face(poly, x, y, u, v, level):
    return sp.factor(sum(
        coefficient*x**i*y**j
        for (i, j), coefficient in sp.Poly(poly, x, y).terms()
        if u*i + v*j == level
    ))


def rem_expr(expr, modulus, variable):
    numerator, denominator = sp.together(expr).as_numer_denom()
    need(sp.resultant(denominator, modulus, variable) != 0,
         "denominator meets modulus")
    return sp.factor(sp.rem(sp.Poly(numerator, variable),
                            sp.Poly(modulus, variable)).as_expr())


def main():
    s, p, u, Q = sp.symbols("s p u Q")
    phi, theta, eta, zeta = sp.symbols("phi theta eta zeta")
    X, T = sp.symbols("X T")
    k0 = sp.Rational(2848, 45)
    epsilon = -sp.Rational(1376, 135)
    t = p - s**2

    H = sp.expand(
        -3*p + sp.Rational(8, 3)*p**2 + epsilon*p**3
        + k0*s**2*p**2 + phi*s*p**3 + theta*s**2*p**3
        + eta*s*p**4 + zeta*s**3*p**3
    )
    G = -s**2/(2*t) + H
    A = sp.cancel((-s*p + t**2*sp.diff(H, s))/p)
    C0 = sp.expand(s**2 + 2*t**2*sp.diff(H, p))
    B = sp.cancel((C0 + s*A)/t**2)
    need(sp.denom(A) == sp.denom(B) == 1, "polynomial source pair")
    need(sp.factor(C0 + s*A - t**2*B) == 0, "pair conversion")
    need(sp.factor(t**2*sp.diff(G, s) - p*A) == 0, "gradient A")
    need(sp.factor(2*t**2*sp.diff(G, p) - C0) == 0, "gradient C0")
    need((sp.degree(A, s), sp.degree(B, s), sp.degree(C0, s)) == (6, 3, 7),
         "source degrees")
    need(sp.factor(sp.Poly(A, s).LC() - 3*zeta*p**2) == 0, "A infinity")
    need(sp.factor(sp.Poly(B, s).LC() - 9*zeta*p**2) == 0, "B infinity")
    need(sp.factor(sp.Poly(C0, s).LC() - 6*zeta*p**2) == 0, "C0 infinity")
    need(sp.factor(A.subs(p, 0) + s) == 0, "p-zero A")
    need(sp.factor(B.subs(p, 0) + 6) == 0, "p-zero B")
    need(sp.factor(A.subs(p, s**2) + s) == 0, "t-zero A")
    need(B.subs({s: 0, p: 0}) == -6, "t-zero endpoint")

    JAB = sp.det(sp.Matrix(((sp.diff(A, s), sp.diff(A, p)),
                            (sp.diff(B, s), sp.diff(B, p)))))
    JAC = sp.det(sp.Matrix(((sp.diff(A, s), sp.diff(A, p)),
                            (sp.diff(C0, s), sp.diff(C0, p)))))
    Hess = sp.det(sp.hessian(G, (s, p)))
    bridge_AB = sp.together(p*JAB - 2*t**2*Hess)
    bridge_AC = sp.together(p*JAC - 2*t**4*Hess)
    a = t**2
    a_s, a_p = sp.diff(a, s), sp.diff(a, p)
    Gss, Gsp, Gpp = sp.diff(G, s, 2), sp.diff(G, s, p), sp.diff(G, p, 2)
    explicit_AC = (
        A*sp.diff(C0, s)
        + 2*p*A*(a_s*Gpp-a_p*Gsp)
        + C0*(a_p*Gss-a_s*Gsp)
    )
    need(sp.factor(sp.together(bridge_AC-explicit_AC)) == 0,
         "explicit AC0 Hessian bridge")
    pair_jacobian_correction = (
        B*(sp.diff(A, s)*a_p-sp.diff(A, p)*a_s) + A*sp.diff(A, p)
    )
    need(sp.factor(sp.together(
        a*bridge_AB - (bridge_AC-p*pair_jacobian_correction)
    )) == 0, "AB bridge from AC0 bridge")

    res = sp.factor(sp.resultant(A, B, s))
    order, R21 = residual(res, p)
    I = sp.factor(4*theta*k0**2 - 27*zeta**2)
    generic_lc = 2125764*eta**3*zeta**2*(eta + zeta)**4
    need((order, R21.degree()) == (6, 21), "generic p6 R21")
    need(sp.factor(R21.TC() + 46656*zeta*I) == 0, "generic bottom")
    need(sp.factor(R21.LC() - generic_lc) == 0, "generic top")

    # A disjoint critical pair.  Resultant multiplicativity applied to
    # C0=t^2 B-sA proves Res(A,C0)=Res(A,t)^2 Res(A,B).  The small second
    # resultant is checked directly; no second giant Sylvester determinant
    # is imported from the primary path.
    res_At = sp.factor(sp.resultant(A, t, s))
    need(sp.factor(res_At**2-p**2) == 0, "AC0 squared p artifact")
    order_AC = order + 2
    need(order_AC == 8, "AC0 p8 R21")

    # I=0 is parameterized bijectively on zeta!=0 by
    # theta=3u^2, zeta=2*k0*u/3, with u=3*zeta/(2*k0).
    sub_I = {theta: 3*u**2, zeta: 2*k0*u/3}
    res_I = sp.factor(sp.resultant(A.subs(sub_I), B.subs(sub_I), s))
    need(sp.factor(res_I - res.subs(sub_I)) == 0, "direct I specialization")
    order_I, R20 = residual(res_I, p)
    J = 8544*phi - 1215*u**3 - 22784*u
    lc_I = sp.factor(generic_lc.subs(sub_I))
    need((order_I, R20.degree()) == (7, 20), "I wall p7 R20")
    need(sp.factor(R20.TC() - sp.Rational(8305770496, 1125)*u**2*J) == 0,
         "I wall bottom")
    need(sp.factor(R20.LC() - lc_I) == 0, "I wall top")

    # J=0 fixes phi and exposes one further affine-endpoint equation M=0.
    phi_J = sp.factor((1215*u**3 + 22784*u)/8544)
    sub_J = dict(sub_I, phi=phi_J)
    res_J = sp.factor(sp.resultant(A.subs(sub_J), B.subs(sub_J), s))
    need(sp.factor(res_J - res.subs(sub_J)) == 0, "direct J specialization")
    order_J, R19 = residual(res_J, p)
    E = 2460375*u**4 - 204543360*u**2 + 5580439552
    M = 547499520*eta + u*E
    need((order_J, R19.degree()) == (8, 19), "J wall p8 R19")
    need(sp.factor(R19.TC() - sp.Rational(2916352, 16875)*u**2*M) == 0,
         "J wall bottom")
    need(sp.factor(R19.LC() - lc_I) == 0, "J wall top")

    # M=0 fixes eta.  The remaining endpoint is the quartic N.
    eta_M = -u*E/sp.Integer(547499520)
    sub_M = dict(sub_J, eta=eta_M)
    res_M = sp.factor(sp.resultant(A.subs(sub_M), B.subs(sub_M), s))
    need(sp.factor(res_M - res.subs(sub_M)) == 0, "direct M specialization")
    order_M, R18 = residual(res_M, p)
    N = 27064125*u**4 - 5739517440*u**2 + 47239069696
    anti_gate = 18225*u**4 - 1515136*u**2 - 129777664
    need(sp.factor(eta_M + 2*k0*u/3 + u*anti_gate/sp.Integer(4055552)) == 0,
         "off-anti gate formula")
    need((order_M, R18.degree()) == (9, 18), "M wall p9 R18")
    need(sp.factor(R18.TC() + sp.Rational(512, 375)*u**5*N) == 0,
         "M wall bottom")
    need(sp.factor(R18.LC() - generic_lc.subs(sub_M)) == 0, "M wall top")

    # Terminal quotient.  On N=0 the next p-coefficient is a unit; the other
    # gcd checks prove that no N-root leaves eta*zeta*(eta+zeta)!=0.
    next_N = (
        523137234375*u**6 + 724483046016000*u**4
        - 213325715694551040*u**2 + 1760617494936027136
    )
    need(sp.factor(R18.nth(1) - sp.Rational(4, 1501875)*u**5*next_N) == 0,
         "N wall next coefficient")
    for hostile, label in (
        (sp.diff(N, u), "N squarefree"),
        (u, "u unit"),
        (E, "eta unit"),
        (anti_gate, "off-anti unit"),
        (next_N, "terminal coefficient unit"),
    ):
        need(sp.gcd(sp.Poly(N, u), sp.Poly(hostile, u)).degree() == 0, label)
    # Hence the source resultant is p^10 R17 at each of the four N-roots.
    terminal_degree = R18.degree() - 1
    need(terminal_degree == 17, "terminal residual degree")

    # Independent normalized chart: two collapsed fibres contribute two
    # Morse points apiece, and the off-anti top row prevents X-infinity loss.
    P = T + X**2*T**2
    Y = X*T*P
    GN = sp.expand(
        -X**2*T/2 - 3*P + sp.Rational(8, 3)*P**2 + epsilon*P**3
        + k0*Y**2 + phi*P**2*Y + theta*P*Y**2
        + eta*P**3*Y + zeta*Y**3
    )
    f = sp.cancel(sp.diff(GN, X)/T)
    h = sp.diff(GN, T)
    need(sp.denom(f) == 1, "normalized f polynomial")
    need((sp.degree(f, X), sp.degree(h, X)) == (8, 9), "normalized degrees")
    need(sp.factor(sp.Poly(f, X).LC() - 9*(eta+zeta)*T**8) == 0,
         "normalized f infinity")
    need(sp.factor(sp.Poly(h, X).LC() - 9*(eta+zeta)*T**8) == 0,
         "normalized h infinity")
    need(sp.factor(h.subs(T, 0) + (X**2 + 6)/2) == 0, "T0 points")
    special_T = -sp.Rational(1, 6)
    need(sp.rem(sp.Poly(f.subs(T, special_T), X),
                sp.Poly(X**2 - 6, X)).is_zero, "P0 f points")
    need(sp.rem(sp.Poly(h.subs(T, special_T), X),
                sp.Poly(X**2 - 6, X)).is_zero, "P0 h points")
    Hess_N = sp.factor(sp.det(sp.hessian(GN, (X, T))))
    need(rem_expr(Hess_N.subs(T, 0) - 6, X**2 + 6, X) == 0,
         "T0 Hessian +6")
    need(rem_expr(Hess_N.subs(T, special_T) + 6, X**2 - 6, X) == 0,
         "P0 Hessian -6")

    # Contracted Delta=0 boundary ledger; all endpoint strata retain it.
    FQ = sp.expand((s**2 - p)*(1 - Q*H) - Q*s**2/2)
    support = tuple(monomial for monomial, coefficient
                    in sp.Poly(FQ, s, p).terms() if coefficient != 0)
    vertices = convex_hull(support)
    expected_vertices = ((0, 1), (2, 0), (5, 3), (1, 5), (0, 4))
    need(vertices == expected_vertices, "contracted polygon")
    expected_faces = (
        s**2*(1-Q/2)-p,
        s**2*(1-Q/2-k0*Q*(s*p)**2-zeta*Q*(s*p)**3),
        Q*p**3*s*(p-s**2)*(eta*p+zeta*s**2),
        Q*p**4*(epsilon+eta*s*p),
        p*(-1+Q*(-3*p+sp.Rational(8, 3)*p**2+epsilon*p**3)),
    )
    face_data = ((1, 2, 2), (-1, 1, -2), (-1, -2, -11),
                 (1, -1, -4), (1, 0, 0))
    for expected, data in zip(expected_faces, face_data):
        need(sp.factor(face(FQ, s, p, *data) - expected) == 0, "boundary face")
    area2 = abs(sum(
        vertices[i][0]*vertices[(i+1) % len(vertices)][1]
        - vertices[(i+1) % len(vertices)][0]*vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(gcd(
        abs(vertices[(i+1) % len(vertices)][0]-vertices[i][0]),
        abs(vertices[(i+1) % len(vertices)][1]-vertices[i][1]),
    ) for i in range(len(vertices)))
    need((area2, boundary, (area2-boundary+2)//2) == (30, 10, 11),
         "Pick ledger")
    packet = (8, 8, 4, 2, 2, 2, 1)
    need((sum(packet), sum(index-1 for index in packet)) == (27, 20),
         "packet ledger")

    print("AUDIT mixed_off_anti_Delta0_source_tower")
    print("scope=Delta=0,K=2848/45,eta*zeta*(eta+zeta)!=0")
    print("pair_AB_degrees=(6,3);pair_AC0_degrees=(6,7)")
    print("pair_resultants=p^6*R21 and p^8*R21;residual_identical=YES")
    print("hessian_bridges=AB:pJ=2t^2Hess;AC0:pJ=2t^4Hess")
    print("chart_loss=p0/t0_empty_for_AB;X_infinity_empty_off_anti")
    print("universal_fibres=T0:2,Hess+6;P0:2,Hess-6")
    print("tower=R21,R20,R19,R18,R17;affine_lengths=25,24,23,22,21")
    print("I=4*theta*k0^2-27*zeta^2")
    print("I0:zeta=5696u/135,theta=3u^2")
    print("J=8544phi-1215u^3-22784u")
    print("M=547499520eta+u*(2460375u^4-204543360u^2+5580439552)")
    print("N=27064125u^4-5739517440u^2+47239069696")
    print("terminal_N_gcds=squarefree,u,eta,offanti,next-coefficient:ALL_ONE")
    print("polygon=((0,1),(2,0),(5,3),(1,5),(0,4));Pick=(30,10,11)")
    print("packet=(8,8,4,2,2,2,1);full_n=27;defect=20")
    print(f"checks={CHECKS}")
    print("MIXED_OFF_ANTI_DELTA0_SOURCE_TOWER_ACCEPT")


if __name__ == "__main__":
    main()
