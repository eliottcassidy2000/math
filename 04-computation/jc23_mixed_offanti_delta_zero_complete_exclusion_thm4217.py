#!/usr/bin/env python3
"""Primary exact certificate for THM-4217.

The certificate specializes every source wall before elimination, proves the
complete D/J/S/T0 critical-length cascade, restores the two universal pairs,
and reconstructs every primitive face of the contracted Newton polygon.
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
    need(bool(terms), "zero polynomial has no valuation")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


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


def polygon_ledger(points):
    vertices = convex_hull(points)
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]))
        for i in range(len(vertices))
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return vertices, area2, boundary, (area2 - boundary + 2) // 2


def face(poly, s, p, u, v, level):
    return sp.factor(sum(
        coefficient * s**i * p**j
        for (i, j), coefficient in sp.Poly(poly, s, p).terms()
        if u*i + v*j == level
    ))


def main():
    s, p, Phi, Theta, eta, zeta, Q, u = sp.symbols(
        "s p Phi Theta eta zeta Q u"
    )
    X, T = sp.symbols("X T")
    k0 = sp.Rational(2848, 45)
    epsilon = -sp.Rational(1376, 135)
    t = p - s**2
    H = (
        -3*p + sp.Rational(8, 3)*p**2 + epsilon*p**3
        + k0*s**2*p**2 + Phi*s*p**3 + Theta*s**2*p**3
        + eta*s*p**4 + zeta*s**3*p**3
    )
    G = -s**2/(2*t) + H
    A = sp.cancel((-s*p + t**2*sp.diff(H, s))/p)
    C0 = sp.expand(s**2 + 2*t**2*sp.diff(H, p))
    B = sp.cancel((C0 + s*A)/t**2)

    # Polynomial source chart, loss ledger, and Morse bridge.
    need(sp.denom(A) == sp.denom(B) == 1, "source pair is polynomial")
    need((sp.degree(A, s), sp.degree(B, s)) == (6, 3),
         "source s-degrees")
    need(sp.factor(sp.Poly(A, s).LC() - 3*zeta*p**2) == 0,
         "A leading row")
    need(sp.factor(sp.Poly(B, s).LC() - 9*zeta*p**2) == 0,
         "B leading row")
    need(sp.factor(t**2*sp.diff(G, s) - p*A) == 0,
         "first gradient identity")
    need(sp.factor(2*t**2*sp.diff(G, p) - (t**2*B-s*A)) == 0,
         "second gradient identity")
    jacobian = sp.det(sp.Matrix(((sp.diff(A, s), sp.diff(A, p)),
                                 (sp.diff(B, s), sp.diff(B, p)))))
    hessian = sp.det(sp.hessian(G, (s, p)))
    bridge = sp.together(p*jacobian - 2*t**2*hessian).as_numer_denom()[0]
    need(sp.factor(sp.reduced(bridge, [A, B], s, p)[1]) == 0,
         "Hessian bridge modulo source ideal")
    need(sp.factor(A.subs(p, 0) + s) == 0, "p-zero A row")
    need(sp.factor(B.subs(p, 0) + 6) == 0, "p-zero B row")
    need(sp.factor(A.subs(p, s**2) + s) == 0, "t-zero A row")
    need(B.subs({p: s**2, s: 0}) == -6, "t-zero endpoint")

    # General mixed row.
    resultant = sp.factor(sp.resultant(A, B, s))
    generic_v = valuation(resultant, p)
    generic = sp.Poly(sp.cancel(resultant/p**generic_v), p)
    D = 4*k0**2*Theta - 27*zeta**2
    need((generic_v, generic.degree()) == (6, 21), "generic p^6 R21")
    need(sp.factor(generic.TC() + 46656*zeta*D) == 0,
         "generic bottom endpoint")
    need(sp.factor(generic.LC()
                   - 2**2*3**12*eta**3*zeta**2*(eta+zeta)**4) == 0,
         "generic top endpoint")

    # D=0 is covered bijectively by zeta=2*k0*u/3, Theta=3*u^2.
    d_sub = {zeta: 2*k0*u/3, Theta: 3*u**2}
    need(sp.factor(D.subs(d_sub)) == 0, "D-wall chart lands")
    d_direct = sp.factor(sp.resultant(A.subs(d_sub), B.subs(d_sub), s))
    need(sp.factor(d_direct-resultant.subs(d_sub)) == 0,
         "D-wall specialize-before-eliminate")
    d_v = valuation(d_direct, p)
    d_residual = sp.Poly(sp.cancel(d_direct/p**d_v), p)
    J = 8544*Phi - 1215*u**3 - 22784*u
    need((d_v, d_residual.degree()) == (7, 20), "D-wall p^7 R20")
    need(sp.factor(d_residual.TC()
                   - sp.Rational(8305770496, 1125)*u**2*J) == 0,
         "D-wall bottom endpoint J")
    need(sp.factor(d_residual.LC()
                   - sp.Rational(129777664, 11390625)
                   * eta**3*u**2*(135*eta+5696*u)**4) == 0,
         "D-wall top endpoint")

    # J=0 fixes Phi uniquely.
    phi_j = u*(1215*u**2+22784)/8544
    j_direct = sp.factor(sp.resultant(
        A.subs(d_sub).subs(Phi, phi_j),
        B.subs(d_sub).subs(Phi, phi_j), s
    ))
    need(sp.factor(j_direct-d_direct.subs(Phi, phi_j)) == 0,
         "J-wall specialize-before-eliminate")
    j_v = valuation(j_direct, p)
    j_residual = sp.Poly(sp.cancel(j_direct/p**j_v), p)
    E = 2460375*u**4 - 204543360*u**2 + 5580439552
    S = 547499520*eta + u*E
    need((j_v, j_residual.degree()) == (8, 19), "J-wall p^8 R19")
    need(sp.factor(j_residual.TC()
                   - sp.Rational(2916352, 16875)*u**2*S) == 0,
         "J-wall bottom endpoint S")
    need(sp.factor(j_residual.LC()-d_residual.LC()) == 0,
         "J-wall top endpoint")

    # S=0 fixes eta uniquely.  The last exceptional polynomial is T0.
    eta_s = -u*E/sp.Integer(547499520)
    s_direct = sp.factor(sp.resultant(
        A.subs(d_sub).subs({Phi: phi_j, eta: eta_s}),
        B.subs(d_sub).subs({Phi: phi_j, eta: eta_s}), s
    ))
    need(sp.factor(s_direct-j_direct.subs(eta, eta_s)) == 0,
         "S-wall specialize-before-eliminate")
    s_v = valuation(s_direct, p)
    s_residual = sp.Poly(sp.cancel(s_direct/p**s_v), p)
    T0 = 27064125*u**4 - 5739517440*u**2 + 47239069696
    A0 = 18225*u**4 - 1515136*u**2 - 129777664
    need((s_v, s_residual.degree()) == (9, 18), "S-wall p^9 R18")
    need(sp.factor(s_residual.TC()
                   + sp.Rational(512, 375)*u**5*T0) == 0,
         "S-wall bottom endpoint T0")
    need(sp.factor(s_residual.LC()
                   + u**9*A0**4*E**3
                   / sp.Integer(11731773052904797284477266998693409587200000))
         == 0, "S-wall top endpoint")

    # T0=0 is an all-roots finite algebra, not a sample algebraic root.
    # Coprimality makes the next p-row and every scope unit invertible.
    next_row = sp.Poly(s_residual.nth(1), u)
    terminal_modulus = sp.Poly(T0, u)
    need(sp.gcd(terminal_modulus, sp.diff(terminal_modulus)).degree() == 0,
         "T0 squarefree")
    for name, scope_unit in (
        ("u", u), ("eta", E), ("eta+zeta", A0),
        ("next p-row", next_row.as_expr()),
    ):
        need(sp.gcd(terminal_modulus, sp.Poly(scope_unit, u)).degree() == 0,
             "T0 coprime to " + name)
    need(s_residual.degree()-1 == 17, "terminal residual degree R17")

    # The two excluded coordinate fibres are four universal Morse points.
    P = T + X**2*T**2
    Y = X*T*P
    GN = (
        -X**2*T/2 - 3*P + sp.Rational(8, 3)*P**2 + epsilon*P**3
        + k0*Y**2 + Phi*P**2*Y + Theta*P*Y**2
        + eta*P**3*Y + zeta*Y**3
    )
    f = sp.cancel(sp.diff(GN, X)/T)
    h = sp.diff(GN, T)
    need(sp.denom(f) == 1, "normalized first derivative")
    need(sp.factor(f.subs(T, 0) + X) == 0, "T=0 f row")
    need(sp.factor(h.subs(T, 0) + (X**2+6)/2) == 0, "T=0 h row")
    special_T = -sp.Rational(1, 6)
    for equation in (f, h):
        need(sp.rem(sp.Poly(equation.subs(T, special_T), X),
                    sp.Poly(X**2-6, X)).is_zero,
             "T=-1/6 universal pair")
    hess_n = sp.factor(sp.det(sp.hessian(GN, (X, T))))
    need(sp.rem(sp.Poly(hess_n.subs(T, 0)-6, X),
                sp.Poly(X**2+6, X)).is_zero, "T=0 Hessian +6")
    need(sp.rem(sp.Poly(hess_n.subs(T, special_T)+6, X),
                sp.Poly(X**2-6, X)).is_zero, "T=-1/6 Hessian -6")
    lengths = tuple(degree+4 for degree in (21, 20, 19, 18, 17))
    need(lengths == (25, 24, 23, 22, 21), "critical-length cascade")

    # Contracted Newton polygon and every primitive nonvertical face.
    FQ = sp.expand((s**2-p)*(1-Q*H)-Q*s**2/2)
    ledger = polygon_ledger(sp.Poly(FQ, s, p).monoms())
    vertices = ((0, 1), (2, 0), (5, 3), (1, 5), (0, 4))
    need(ledger == (vertices, 30, 10, 11), "Newton/Pick ledger")
    expected_faces = (
        ((1, 2, 2), s**2*(1-Q/2)-p),
        ((-1, 1, -2),
         s**2*(1-Q/2-k0*Q*(s*p)**2-zeta*Q*(s*p)**3)),
        ((-1, -2, -11), Q*s*p**3*(p-s**2)*(eta*p+zeta*s**2)),
        ((1, -1, -4), Q*p**4*(epsilon+eta*s*p)),
    )
    for (normal_u, normal_v, level), expected in expected_faces:
        need(sp.factor(face(FQ, s, p, normal_u, normal_v, level)-expected)
             == 0, "primitive face")

    packet = []
    for start, end in zip(vertices, vertices[1:]+vertices[:1]):
        dx, dy = end[0]-start[0], end[1]-start[1]
        length = gcd(abs(dx), abs(dy))
        normal_u, normal_v = -dy//length, dx//length
        level = normal_u*start[0] + normal_v*start[1]
        index = normal_u+normal_v-level
        if not (start[0] == end[0] == 0):
            packet.extend([index]*length)
    packet = tuple(sorted(packet, reverse=True))
    need(packet == (8, 8, 4, 2, 2, 2, 1), "complete candidate packet")
    need(sum(index-1 for index in packet) == 20 == 2*11-2,
         "packet defect equals 2g-2")

    # Split top roots, replacement diagonal, and cubic carrier are live
    # precisely under eta*zeta*(eta+zeta)!=0.
    a, W, q = sp.symbols("a W q")
    top = a*(1-a)**3*(eta*(1-a)+zeta)
    need(sp.factor(sp.resultant(a, eta*(1-a)+zeta, a)-(eta+zeta)) == 0,
         "two top factors collide exactly on eta+zeta=0")
    need(sp.diff(epsilon+eta*W, W) == eta, "diagonal root is simple")
    carrier = zeta*W**3+k0*W**2-(q-sp.Rational(1, 2))
    carrier_disc = sp.factor(sp.discriminant(carrier, W))
    need(sp.factor(carrier_disc - (q-sp.Rational(1, 2))
                   * (4*k0**3-27*zeta**2*(q-sp.Rational(1, 2)))) == 0,
         "cubic carrier separability")
    need(top != 0, "top face polynomial")

    # Full and finite response inequalities.
    full_n, full_defect = 27, 20
    full_caps = tuple(2*(full_n-L) for L in lengths)
    need(full_caps == (4, 6, 8, 10, 12), "full caps")
    need(all(cap < full_defect for cap in full_caps), "full contradictions")
    finite_n, carriers, finite_origin_index = 21, 3, 17
    need(2*finite_n-lengths[0] < finite_n-carriers,
         "L25 finite support-union contradiction")
    finite_caps = tuple(2*(finite_n+carriers-L)+carriers
                        for L in lengths[1:])
    need(finite_caps == (3, 5, 7, 9), "finite carrier-orbit caps")
    need(all(cap < finite_origin_index for cap in finite_caps),
         "finite contradictions")

    print("AUDIT complete_mixed_offanti_Delta_zero_exact_M9")
    print("scope=Delta=0,K=2848/45,Phi/Theta arbitrary,eta*zeta*(eta+zeta)!=0")
    print("source_resultants=p^6R21,p^7R20,p^8R19,p^9R18,p^10R17")
    print("endpoint_cascade=D,J,S,T0;lengths=25,24,23,22,21")
    print("terminal_T0=squarefree_and_coprime_to_scope_units_and_next_row")
    print("universal_coordinate_points=2+2;Hessian=+6,-6")
    print("polygon=((0,1),(2,0),(5,3),(1,5),(0,4));Pick=(30,10,11)")
    print("packet=(8,8,4,2,2,2,1);defect=20;carrier=prime_cubic")
    print("full_caps=4,6,8,10,12<20;finite_caps=3,5,7,9<17")
    print(f"checks={CHECKS}")
    print("COMPLETE_MIXED_OFFANTI_DELTA_ZERO_ACCEPT")


if __name__ == "__main__":
    main()
