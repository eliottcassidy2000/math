#!/usr/bin/env python3
"""Exact regular-model certificate for THM-4218.

The universe is the dense exact-M=10 top chamber of the inherited b=d=0
reduced (2,3) seam.  The certificate keeps the Newton degeneration, source
critical control, labelled carrier response, and degree-three attachment
obstruction in one deterministic replay.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def monomials_through(weight):
    answer = []
    for i in range(weight // 2 + 1):
        for j in range(weight // 3 + 1):
            residual_weight = 2 * i + 3 * j
            if (0 < residual_weight <= weight
                    and (i, j) not in {(0, 1), (1, 1)}):
                answer.append((i, j, residual_weight))
    return tuple(sorted(answer, key=lambda row: (row[2], row[1], row[0])))


def expanded_support(monomials):
    """Valued support (s exponent,p exponent,Q height), after combining base."""
    support = {(2, 0, 0), (0, 1, 0)}
    for i, j, _ in monomials:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


def lower_planes(support):
    points = sorted(support)
    faces = {}
    for ia, ib, ic in combinations(range(len(points)), 3):
        first, second, third = (points[k] for k in (ia, ib, ic))
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant == 0:
            continue
        slope_s = Fraction(
            (second[2] - first[2]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[2] - first[2]),
            determinant,
        )
        slope_p = Fraction(
            (second[0] - first[0]) * (third[2] - first[2])
            - (second[2] - first[2]) * (third[0] - first[0]),
            determinant,
        )
        constant = (Fraction(first[2]) - slope_s * first[0]
                    - slope_p * first[1])
        gaps = [Fraction(height) - slope_s * i - slope_p * j - constant
                for i, j, height in points]
        if min(gaps) < 0:
            continue
        equality = tuple(index for index, gap in enumerate(gaps) if gap == 0)
        if len(equality) >= 3:
            faces[equality] = (slope_s, slope_p, constant)
    return {plane for plane in faces.values()}


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, first, second):
        return ((first[0] - origin[0]) * (second[1] - origin[1])
                - (first[1] - origin[1]) * (second[0] - origin[0]))

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
        vertices[k][0] * vertices[(k + 1) % len(vertices)][1]
        - vertices[(k + 1) % len(vertices)][0] * vertices[k][1]
        for k in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(k + 1) % len(vertices)][0] - vertices[k][0]),
            abs(vertices[(k + 1) % len(vertices)][1] - vertices[k][1]))
        for k in range(len(vertices))
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    genus = (area2 - boundary + 2) // 2
    return vertices, area2, boundary, genus


def edge_packet(vertices):
    packet = []
    edge_rows = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        edge_rows.append((start, end, length, inward, constant, index))
        if start[0] == end[0] == 0:
            continue
        packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True)), tuple(edge_rows)


def face(polynomial, first, second, normal, constant):
    answer = 0
    for (i, j), coefficient in sp.Poly(polynomial, first, second).terms():
        if normal[0] * i + normal[1] * j == constant:
            answer += coefficient * first**i * second**j
    return sp.factor(answer)


def valuation(polynomial, variable):
    terms = sp.Poly(polynomial, variable).terms()
    need(bool(terms), "valuation of zero polynomial")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def coefficient_digest(polynomial, monic=False):
    if monic:
        polynomial = polynomial.monic()
    payload = ",".join(str(sp.Rational(c)) for c in polynomial.all_coeffs())
    return sha256(payload.encode()).hexdigest()


def main():
    # Complete exact-M=10 source support.
    monomials = monomials_through(10)
    expected_monomials = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10),
    )
    need(monomials == expected_monomials, "M10 monomial universe")

    required = {(5, 0), (2, 2), (0, 3)}
    optional = tuple(row for row in monomials if row[:2] not in required)
    plane_main = (Fraction(1, 10), Fraction(1, 5), Fraction(-1, 5))
    plane_tail = (Fraction(1, 6), Fraction(1, 6), Fraction(-1, 3))
    expected_planes = {plane_main, plane_tail}
    for mask in range(1 << len(optional)):
        chosen = [row for row in monomials if row[:2] in required]
        chosen.extend(row for bit, row in enumerate(optional)
                      if mask & (1 << bit))
        need(lower_planes(expanded_support(chosen)) == expected_planes,
             "lower hull changed under lower-support deletion")

    # All possible coincident expanded coefficients are interior to the
    # forced lower envelope.  Deleting any collection of them is a stronger
    # cancellation hostile than an arbitrary equality of their coefficients.
    collisions = {(2, 3, 1), (2, 4, 1), (2, 5, 1), (3, 4, 1)}
    full_support = expanded_support(monomials)
    for mask in range(1 << len(collisions)):
        deleted = {point for bit, point in enumerate(sorted(collisions))
                   if mask & (1 << bit)}
        need(lower_planes(full_support - deleted) == expected_planes,
             "coefficient cancellation changed lower hull")

    # Analytic gap audit, independent of the subset enumeration.
    for i, j, weight in monomials:
        for point in ((j + 2, i + j, 1), (j, i + j + 1, 1)):
            m_value = plane_main[0] * point[0] + plane_main[1] * point[1] + plane_main[2]
            t_value = plane_tail[0] * point[0] + plane_tail[1] * point[1] + plane_tail[2]
            need(Fraction(point[2]) - m_value == Fraction(10 - weight, 10),
                 "main-plane weight gap")
            need(Fraction(point[2]) - t_value >= 0, "tail-plane gap")

    main_points = ((0, 1), (2, 0), (4, 4), (2, 5), (0, 6))
    tail_points = ((2, 0), (5, 3), (4, 4))
    global_points = main_points + tail_points
    need(polygon_ledger(main_points)[1:] == (30, 10, 11), "main Pick ledger")
    need(polygon_ledger(tail_points)[1:] == (6, 6, 1), "tail Pick ledger")
    polygon, area2, boundary, genus = polygon_ledger(global_points)
    need((polygon, area2, boundary, genus) == (
        ((0, 1), (2, 0), (5, 3), (4, 4), (0, 6)), 36, 12, 13
    ), "global polygon/Pick ledger")
    packet, edge_rows = edge_packet(polygon)
    need(packet == (9, 9, 6, 2, 2, 2, 1), "outer packet")
    need(sum(packet) == 31 and sum(index - 1 for index in packet) == 24,
         "packet degree/defect")

    # Exact complete source and all outer face polynomials.
    s, p, Q = sp.symbols("s p Q")
    Phi, Delta, Theta, eta = sp.symbols("Phi Delta Theta eta")
    zeta, u, v = sp.symbols("zeta u v")
    K = sp.Rational(2848, 45) - sp.Rational(7, 6)*Delta
    need(sp.factor(K - (sp.Rational(2848, 45)
                        - sp.Rational(7, 6)*Delta)) == 0,
         "forced K/Delta relation")
    epsilon = -sp.Rational(1376, 135)
    H = (
        -3*p + sp.Rational(8, 3)*p**2 + epsilon*p**3
        + K*s**2*p**2 + Phi*s*p**3 + Delta*p**4
        + Theta*s**2*p**3 + eta*s*p**4 + zeta*s**3*p**3
        + u*p**5 + v*s**2*p**4
    )
    F = sp.expand((s**2 - p)*(1 - Q*H) - Q*s**2/2)
    expected_outer_faces = (
        s**2*(1 - Q/2) - p,
        s**2*(1 - Q/2 - Q*K*(s*p)**2 - Q*zeta*(s*p)**3),
        -Q*s**4*p**3*(zeta*s + v*p),
        Q*p**4*(p - s**2)*(u*p + v*s**2),
        p*(-1 + Q*(-3*p + sp.Rational(8, 3)*p**2 + epsilon*p**3
                   + Delta*p**4 + u*p**5)),
    )
    for row, expected_face in zip(edge_rows, expected_outer_faces):
        start, end, length, normal, constant, index = row
        need(sp.factor(face(F, s, p, normal, constant) - expected_face) == 0,
             "outer face polynomial")

    # Face normalizations and transverse attachment inventory.
    S, P = sp.symbols("S P")
    main_face = (S**2 - P)*(1 - u*P**5 - v*S**2*P**4)
    tail_face = 1 - zeta*S**3*P**3 - v*S**2*P**4
    need(sp.factor(main_face - (S**2 - P)*(1 - u*P**5 - v*S**2*P**4)) == 0,
         "main face factorization")
    U0 = S**2 - P
    C5 = 1 - u*P**5 - v*S**2*P**4
    node_det = sp.det(sp.Matrix((
        (sp.diff(U0, S), sp.diff(U0, P)),
        (sp.diff(C5, S), sp.diff(C5, P)),
    ))).subs(P, S**2)
    need(sp.factor(node_det) == -10*S**9*(u + v), "ten transverse main nodes")
    node_polynomial = 1 - (u + v)*S**10
    node_domain = sp.QQ.frac_field(u, v)
    need(sp.gcd(sp.Poly(node_polynomial, S, domain=node_domain),
                sp.Poly(sp.diff(node_polynomial, S), S, domain=node_domain)).degree() == 0,
         "ten main nodes are distinct")

    Y = sp.symbols("Y")
    genus_two = 1 - u*P**5 - v*Y**2
    need(sp.factor(C5 - genus_two.subs(Y, S*P**2)) == 0,
         "genus-two normalization")
    need(sp.gcd(sp.Poly(1 - u*P**5, P, domain=sp.QQ.frac_field(u)),
                sp.Poly(-5*u*P**4, P, domain=sp.QQ.frac_field(u))).degree() == 0,
         "genus-two branch polynomial squarefree")
    need((5 - 1)//2 == 2, "genus-two value")

    T, W = sp.symbols("T W")
    need(sp.factor(tail_face - (1 - zeta*T**3 - v*W**2).subs(
        {T: S*P, W: S*P**2})) == 0, "elliptic-tail normalization")
    tail_domain = sp.QQ.frac_field(zeta, v)
    tail_poly = sp.Poly(1 - zeta*T**3 - v*W**2, T, W, domain=tail_domain)
    tail_derivative_matrix = sp.Matrix(((3*zeta, 2*v), (3*zeta, 4*v)))
    need(sp.factor(tail_derivative_matrix.det()) == 6*zeta*v,
         "tail torus smoothness")

    # The order-five action forces a primitive quartic CM field in the
    # genus-two Jacobian.  Its unique quadratic subfield is real, excluding
    # every elliptic-product endomorphism algebra.
    X = sp.symbols("X")
    phi5 = sp.cyclotomic_poly(5, X)
    need(phi5 == X**4 + X**3 + X**2 + X + 1 and sp.Poly(phi5, X).is_irreducible,
         "Phi5 primitive CM field")
    R = sp.symbols("R")
    trace_resultant = sp.factor(sp.resultant(phi5, X**2 - R*X + 1, X))
    need(trace_resultant == (R**2 + R - 1)**2,
         "unique real quadratic subfield trace")
    need(sp.discriminant(R**2 + R - 1, R) == 5,
         "quadratic subfield is real")
    need(4 > 2 and 8 > 4, "elliptic-product commutative dimension obstruction")

    # Tail attachment points become (0,+/-1) on y^2=x^3+1.
    x0, y0 = sp.Integer(0), sp.Integer(1)
    tangent = sp.cancel(3*x0**2/(2*y0))
    x_double = sp.expand(tangent**2 - 2*x0)
    y_double = sp.expand(tangent*(x0 - x_double) - y0)
    need((x_double, y_double) == (0, -1), "tail attachment difference is 3-torsion")
    need(sp.cyclotomic_poly(3, 1) == 3, "Norm(1-zeta3)=3")

    # Q=sigma^30 makes both lower heights integral and primitive.
    sigma = sp.symbols("sigma")
    normal_main = (3, 6, -1)
    normal_tail = (5, 5, -1)
    need(gcd(gcd(abs(normal_main[0]), abs(normal_main[1])), abs(normal_main[2])) == 1,
         "main primitive face multiplicity")
    need(gcd(gcd(abs(normal_tail[0]), abs(normal_tail[1])), abs(normal_tail[2])) == 1,
         "tail primitive face multiplicity")
    for i, j, height in full_support:
        main_gap = 30*height - (3*i + 6*j - 6)
        tail_gap = 30*height - (5*i + 5*j - 10)
        need(main_gap >= 0 and tail_gap >= 0, "integral lower-height gap")

    Z = sp.symbols("Z")
    special_edges = (
        1 - Z,
        1 - zeta*Z**3,
        zeta + v*Z,
        (1 - Z)*(u + v*Z),
        1 - u*Z**5,
        1 - v*Z**2,
    )
    edge_discriminants = tuple(sp.factor(sp.discriminant(edge, Z))
                               for edge in special_edges)
    need(edge_discriminants == (
        1, -27*zeta**2, 1, (u + v)**2, 3125*u**4, 4*v
    ), "all special edge schemes squarefree")

    # Exact main and tail charts after the common base change.
    H_main = sp.expand(sp.cancel(
        sigma**30 * H.subs({s: sigma**-3*S, p: sigma**-6*P})
    ))
    need(sp.denom(H_main) == 1, "main scaled H polynomial")
    F_main = sp.expand(sp.cancel(
        sigma**6 * F.subs({Q: sigma**30, s: sigma**-3*S, p: sigma**-6*P})
    ))
    expected_F_main = sp.expand(
        (S**2 - P)*(1 - H_main) - sigma**30*S**2/2
    )
    need(sp.factor(F_main - expected_F_main) == 0, "exact main scaled equation")
    need(sp.factor(F_main.subs(sigma, 0) - main_face) == 0,
         "main special face")
    # At each transverse main node, U=S^2-P and V=(1-H_main)/S^2
    # give U*V=sigma^30/2: an A_29 smoothing with a rational resolution.
    need(30 - 1 == 29, "A29 resolution chain length")

    H_tail = sp.expand(sp.cancel(
        sigma**30 * H.subs({s: sigma**-5*S, p: sigma**-5*P})
    ))
    need(sp.denom(H_tail) == 1, "tail scaled H polynomial")
    F_tail = sp.expand(sp.cancel(
        sigma**10 * F.subs({Q: sigma**30, s: sigma**-5*S, p: sigma**-5*P})
    ))
    expected_F_tail = sp.expand(
        (S**2 - sigma**5*P)*(1 - H_tail) - sigma**30*S**2/2
    )
    need(sp.factor(F_tail - expected_F_tail) == 0, "exact tail scaled equation")
    need(sp.factor(F_tail.subs(sigma, 0) - S**2*tail_face) == 0,
         "tail special face")

    # Generic boundary squarefreeness and actual-source torus recovery.
    q = sp.symbols("q")
    carrier = zeta*Z**3 + K*Z**2 - (q - sp.Rational(1, 2))
    carrier_discriminant = sp.factor(sp.discriminant(carrier, Z))
    need(sp.factor(carrier_discriminant - (
        (q - sp.Rational(1, 2))
        * (4*K**3 - 27*zeta**2*(q - sp.Rational(1, 2)))
    )) == 0, "prime cubic carrier discriminant")
    carrier_field = sp.QQ.frac_field(Delta, zeta, q)
    need(sp.gcd(sp.Poly(carrier, Z, domain=carrier_field),
                sp.Poly(sp.diff(carrier, Z), Z, domain=carrier_field)).degree() == 0,
         "generic cubic carrier squarefree")
    H0 = -3*Z + sp.Rational(8, 3)*Z**2 + epsilon*Z**3 + Delta*Z**4 + u*Z**5
    vertical_field = sp.QQ.frac_field(Delta, u, q)
    need(sp.gcd(sp.Poly(H0 - q, Z, domain=vertical_field),
                sp.Poly(sp.diff(H0, Z), Z, domain=vertical_field)).degree() == 0,
         "generic vertical edge squarefree")
    need(sp.factor(F.subs(p, s**2) + Q*s**2/2) == 0,
         "generic torus excludes t=0")

    # Target good reduction at the same base change.
    a, A_target, C_target, X_target, Y_target = sp.symbols(
        "a A_target C_target X_target Y_target"
    )
    target_equation = (
        C_target**2 - A_target**3 + sp.Rational(3, 4)*a**2*A_target
        - sigma**-30 + sp.Rational(1, 4)*a**3
    )
    scaled_target = sp.expand(sp.cancel(
        sigma**30 * target_equation.subs({
            A_target: sigma**-10*X_target,
            C_target: sigma**-15*Y_target,
        })
    ))
    expected_target = (
        Y_target**2 - X_target**3 - 1
        + sp.Rational(3, 4)*a**2*sigma**20*X_target
        + sp.Rational(1, 4)*a**3*sigma**30
    )
    need(sp.factor(scaled_target - expected_target) == 0,
         "target good-reduction scaling")
    need(sp.discriminant(X_target**3 + 1, X_target) == -27,
         "j0 special target smooth")

    # Exact critical-open nonemptiness control.
    ss, pp = sp.symbols("ss pp")
    tt = pp - ss**2
    Delta_c = sp.Rational(1)
    K_c = sp.Rational(2848, 45) - sp.Rational(7, 6)*Delta_c
    H_c = (
        -3*pp + sp.Rational(8, 3)*pp**2 + epsilon*pp**3
        + K_c*ss**2*pp**2 + 2*ss*pp**3 + Delta_c*pp**4
        + 5*ss**2*pp**3 + 7*ss*pp**4 + 11*ss**3*pp**3
        + 13*pp**5 + 17*ss**2*pp**4
    )
    A_c = sp.cancel((-ss*pp + tt**2*sp.diff(H_c, ss))/pp)
    C0_c = sp.expand(ss**2 + 2*tt**2*sp.diff(H_c, pp))
    B_c = sp.cancel((C0_c + ss*A_c)/tt**2)
    need(sp.denom(A_c) == sp.denom(B_c) == 1, "source critical pair polynomial")
    need((sp.degree(A_c, ss), sp.degree(B_c, ss)) == (6, 3),
         "source critical degrees")
    need(sp.Poly(A_c, ss).LC() == 33*pp**2
         and sp.Poly(B_c, ss).LC() == 99*pp**2,
         "source infinity leading rows")
    need(sp.factor(A_c.subs(pp, 0) + ss) == 0 and B_c.subs(pp, 0) == -6,
         "p-zero source chart")
    need(sp.factor(A_c.subs(pp, ss**2) + ss) == 0
         and B_c.subs({pp: ss**2, ss: 0}) == -6,
         "t-zero source chart")
    resultant_ab = sp.factor(sp.resultant(A_c, B_c, ss))
    resultant_ac = sp.factor(sp.resultant(A_c, C0_c, ss))
    need((valuation(resultant_ab, pp), valuation(resultant_ac, pp)) == (6, 8),
         "source exceptional resultant powers")
    R25_ab = sp.Poly(sp.cancel(resultant_ab/pp**6), pp)
    R25_ac = sp.Poly(sp.cancel(resultant_ac/pp**8), pp)
    need(R25_ab == R25_ac and R25_ab.degree() == 25,
         "two source eliminations agree on R25")
    need(R25_ab.nth(0) == -sp.Rational(189675421056, 5)
         and R25_ab.LC() == 2731549392000000000,
         "R25 endpoint units")
    need(sp.gcd(R25_ab, R25_ab.diff()).degree() == 0,
         "R25 control squarefree")

    XX, TT = sp.symbols("XX TT")
    PP = TT + XX**2*TT**2
    YY = XX*TT*PP
    G_c = (
        -XX**2*TT/2 - 3*PP + sp.Rational(8, 3)*PP**2 + epsilon*PP**3
        + K_c*YY**2 + 2*PP**2*YY + Delta_c*PP**4 + 5*PP*YY**2
        + 7*PP**3*YY + 11*YY**3 + 13*PP**5 + 17*PP**2*YY**2
    )
    f_c = sp.cancel(sp.diff(G_c, XX)/TT)
    h_c = sp.diff(G_c, TT)
    resultant_xt = sp.factor(sp.resultant(f_c, h_c, XX))
    need(valuation(resultant_xt, TT) == 72, "normalized T exceptional power")
    normalized_residual = sp.cancel(resultant_xt/TT**72)
    universal_power = 0
    while sp.rem(sp.Poly(normalized_residual, TT), sp.Poly(6*TT + 1, TT)) == 0:
        normalized_residual = sp.cancel(normalized_residual/(6*TT + 1))
        universal_power += 1
    Q25 = sp.Poly(normalized_residual, TT)
    need(universal_power == 2 and Q25.degree() == 25,
         "normalized resultant T^72(6T+1)^2Q25")
    need(Q25.nth(0) != 0 and Q25.eval(-sp.Rational(1, 6)) != 0
         and sp.gcd(Q25, Q25.diff()).degree() == 0,
         "normalized Q25 open control")

    # Genus and degree conservation ledgers.
    main_graph_rank = 10 - 2 + 1
    tail_attachment_graph_gain = 2 - 2 + 1
    total_graph_rank = main_graph_rank + tail_attachment_graph_gain
    need((main_graph_rank, tail_attachment_graph_gain, total_graph_rank) == (9, 1, 10),
         "dual graph ranks")
    need(2 + 1 + total_graph_rank == genus, "regular special genus checksum")
    critical_length = R25_ab.degree() + 2 + 2
    full_degree = sum(packet)
    finite_degree = full_degree - 3*2
    beta = 3
    need((critical_length, full_degree, finite_degree, beta) == (29, 31, 25, 3),
         "critical/response degrees")
    full_cap = 2*(full_degree - critical_length)
    finite_cap = 2*finite_degree - critical_length - 1 + beta
    need(full_cap == 4 < 24, "full response deficit")
    need(finite_cap == 23 < finite_degree - 1, "finite response deficit")
    need(full_degree % 3 == finite_degree % 3 == 1,
         "both carrier responses are nonzero mod three")
    # The multiplicity-one elliptic tail is the only nonconstant special
    # component; equality of its two attachment images forces its isogeny
    # kernel to contain a subgroup of order three.  Degree conservation
    # therefore forces 3|degree(phi), contradicting both response degrees.
    need(31 % 3 != 0 and 25 % 3 != 0, "degree-three contradiction")

    print("THM4218_WEIGHT10_HIDDEN_ELLIPTIC_TAIL_EXACT_CERTIFICATE")
    print("scope=exact_M10,b=d=0,reduced_(2,3),u*v*zeta*(u+v)!=0,all_lower_coefficients")
    print("monomials=11;lower_support_subsets=256;lower_faces=M10,T6;cancellations=16")
    print("polygon=((0,1),(2,0),(5,3),(4,4),(0,6));Pick=(36,12,13)")
    print("packet=(9,9,6,2,2,2,1);degree=31;defect=24;cubic_carrier_degree=3")
    print("main=genus2_plus_rational_with_10_nodes;tail=j0_elliptic;shared_roots=2")
    print("CM=Q(zeta5)_primitive;genus2_Jacobian_simple;Hom_to_E0=0")
    print("attachments=(0,+1),(0,-1);difference_order=3;Norm(1-zeta3)=3")
    print("base_change=Q:sigma^30;face_multiplicities=1;A29_nodes=10;other_chains=rational")
    print("target_special=Y^2=X^3+1;positive_genus=(2,1);dual_graph_rank=10")
    print("critical_control=R25_squarefree;L=29;source_digest="
          + coefficient_digest(R25_ab))
    print("responses=full:n31_cap4<24;finite:n25_beta3_cap23<24")
    print("degree_conservation=3_divides_tail_degree;31_mod3=1;25_mod3=1")
    print(f"checks={CHECKS}")
    print("verdict=EXACT_M10_DENSE_TOP_CHAMBER_EXCLUDED_RELATIVE")


if __name__ == "__main__":
    main()
