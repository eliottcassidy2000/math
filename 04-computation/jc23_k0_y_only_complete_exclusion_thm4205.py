#!/usr/bin/env python3
"""Exact source/owner certificate for THM-4205's Y-only K=0 closure.

This certificate specializes the complete exact-weight-nine source before
elimination.  Its decisive object is the two-generator source critical ideal,
not a squarefree coordinate projection.
"""

from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(poly, variable):
    terms = sp.Poly(poly, variable).terms()
    need(bool(terms), "zero polynomial")
    return min(monomial[0] for monomial, coefficient in terms if coefficient != 0)


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, left, right):
        return (
            (left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0])
        )

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


def zero_mod_x2(expression, variable, square):
    return sp.rem(sp.Poly(expression, variable), sp.Poly(variable**2 - square, variable)).is_zero


def main():
    s, p, Phi, Theta, zeta, Q = sp.symbols("s p Phi Theta zeta Q")
    t = p - s**2
    Delta = sp.Rational(5696, 105)
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    need(K == 0, "K=0/Delta relation")

    H = (
        -3 * p
        + sp.Rational(8, 3) * p**2
        - sp.Rational(1376, 135) * p**3
        + Phi * s * p**3
        + Delta * p**4
        + Theta * s**2 * p**3
        + zeta * s**3 * p**3
    )
    G = -s**2 / (2 * t) + H
    A = sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p)
    C0 = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    B = sp.cancel((C0 + s * A) / t**2)
    need(sp.denom(A) == sp.denom(B) == 1, "polynomial source pair")
    need((sp.degree(A, s), sp.degree(B, s)) == (6, 3), "source degrees")
    need(sp.factor(sp.Poly(A, s).LC() - 3 * zeta * p**2) == 0, "A infinity row")
    need(sp.factor(sp.Poly(B, s).LC() - 9 * zeta * p**2) == 0, "B infinity row")
    need(sp.factor(t**2 * sp.diff(G, s) - p * A) == 0, "first gradient identity")
    need(
        sp.factor(2 * t**2 * sp.diff(G, p) - (t**2 * B - s * A)) == 0,
        "second gradient identity",
    )
    jacobian = sp.det(
        sp.Matrix(((sp.diff(A, s), sp.diff(A, p)), (sp.diff(B, s), sp.diff(B, p))))
    )
    hessian = sp.det(sp.hessian(G, (s, p)))
    bridge = sp.together(p * jacobian - 2 * t**2 * hessian).as_numer_denom()[0]
    need(sp.factor(sp.reduced(bridge, [A, B], s, p)[1]) == 0, "Hessian bridge")
    need(sp.factor(A.subs(p, 0) + s) == 0, "p=0 A row")
    need(sp.factor(B.subs(p, 0) + 6) == 0, "p=0 B row")
    need(sp.factor(A.subs(p, s**2) + s) == 0, "t=0 A row")
    need(B.subs({s: 0, p: 0}) == -6, "t=0 endpoint")

    resultant = sp.factor(sp.resultant(A, B, s))
    p_valuation = valuation(resultant, p)
    residual = sp.Poly(sp.cancel(resultant / p**p_valuation), p)
    need((p_valuation, residual.degree()) == (6, 20), "p^6 R20")
    need(sp.factor(residual.nth(0) - 1259712 * zeta**3) == 0, "R20 constant")
    need(
        sp.factor(residual.LC() - sp.Rational(40870620168192, 1225) * zeta**7)
        == 0,
        "R20 leading coefficient",
    )

    # Restore the two coordinate fibres lost by the rational source chart.
    X, T = sp.symbols("X T")
    P = T + X**2 * T**2
    Y = X * T * P
    GN = (
        -X**2 * T / 2
        - 3 * P
        + sp.Rational(8, 3) * P**2
        - sp.Rational(1376, 135) * P**3
        + Phi * P**2 * Y
        + Delta * P**4
        + Theta * P * Y**2
        + zeta * Y**3
    )
    f = sp.cancel(sp.diff(GN, X) / T)
    h = sp.diff(GN, T)
    need(sp.denom(f) == 1, "normalized derivative")
    need((sp.degree(f, X), sp.degree(h, X)) == (8, 9), "normalized degrees")
    need(sp.factor(sp.Poly(f, X).LC() - 9 * zeta * T**8) == 0, "normalized f infinity")
    need(sp.factor(sp.Poly(h, X).LC() - 9 * zeta * T**8) == 0, "normalized h infinity")
    hess_xt = sp.det(sp.hessian(GN, (X, T)))
    for t_value, x_square, g_value, hess_value in (
        (0, -6, 0, 6),
        (-sp.Rational(1, 6), 6, sp.Rational(1, 2), -6),
    ):
        substitutions = {T: t_value}
        need(zero_mod_x2(sp.expand(sp.diff(GN, X).subs(substitutions)), X, x_square), "universal G_X")
        need(zero_mod_x2(sp.expand(sp.diff(GN, T).subs(substitutions)), X, x_square), "universal G_T")
        need(zero_mod_x2(sp.expand((GN - g_value).subs(substitutions)), X, x_square), "universal value")
        need(zero_mod_x2(sp.expand((hess_xt - hess_value).subs(substitutions)), X, x_square), "universal Hessian")

    # Literal Y-only Newton polygon and every primitive face.
    fibre = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
    support = {
        monomial: sp.factor(coefficient)
        for monomial, coefficient in sp.Poly(fibre, s, p).terms()
        if coefficient != 0
    }
    vertices = convex_hull(support)
    expected_vertices = ((0, 1), (2, 0), (5, 3), (3, 4), (0, 5))
    need(vertices == expected_vertices, "Y-only Newton polygon")
    expected_faces = (
        s**2 * (1 - Q / 2) - p,
        s**2 * (1 - Q / 2 - Q * zeta * s**3 * p**3),
        Q * p**3 * s**3 * zeta * (p - s**2),
        Q * p**4 * (Delta * p + zeta * s**3),
    )
    packet = []
    boundary = 0
    for edge_index, (start, stop) in enumerate(zip(vertices, vertices[1:] + vertices[:1])):
        dx, dy = stop[0] - start[0], stop[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        u, v = -dy // length, dx // length
        constant = u * start[0] + v * start[1]
        index = u + v - constant
        boundary += length
        if not (start[0] == stop[0] == 0):
            packet.extend([index] * length)
        if edge_index < 4:
            face = sum(
                coefficient * s**i * p**j
                for (i, j), coefficient in support.items()
                if u * i + v * j == constant
            )
            need(sp.factor(face - expected_faces[edge_index]) == 0, "primitive face")
    area2 = abs(
        sum(
            vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
            - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
            for index in range(len(vertices))
        )
    )
    genus = (area2 - boundary + 2) // 2
    packet = tuple(sorted(packet, reverse=True))
    need((area2, boundary, genus) == (30, 10, 11), "Pick ledger")
    need(packet == (11, 8, 2, 2, 2, 1), "Y-only packet")

    # The prime cubic is one owner.  Its three geometric conjugates cannot be
    # assigned independently; THM-4147's transport keeps beta=3 as sidecar.
    W, q = sp.symbols("W q")
    carrier = sp.Poly(zeta * W**3 - (q - sp.Rational(1, 2)), W)
    need(carrier.degree() == 3, "cubic carrier degree")
    need(sp.resultant(carrier.as_expr(), sp.diff(carrier.as_expr(), W), W) != 0, "generic carrier separability")
    need(sp.Poly(q - sp.Rational(1, 2), q).degree() == 1, "carrier valuation-one prime")

    critical_length = residual.degree() + 4
    defect = sum(index - 1 for index in packet)
    full_degree = sum(packet)
    finite_degree, beta = full_degree - 6, 3
    need((critical_length, defect, full_degree, finite_degree, beta) == (24, 20, 26, 20, 3), "response ledger")
    need(2 * finite_degree - critical_length - 1 + beta == 18 < finite_degree - 1, "finite response gate")
    need(beta < finite_degree - 1, "identity-handle gate")
    need(2 * (full_degree - critical_length) == 4 < defect, "full response gate")

    # Exact M=9 means (eta,zeta)!=(0,0).  These four rows are disjoint and
    # exhaustive and are owned respectively by THM-4192, this audit,
    # THM-4209, and THM-4176.
    seen_rows = set()
    for eta_value in range(-2, 3):
        for zeta_value in range(-2, 3):
            if (eta_value, zeta_value) == (0, 0):
                continue
            matches = {
                "P": zeta_value == 0 and eta_value != 0,
                "Y": eta_value == 0 and zeta_value != 0,
                "G": eta_value * zeta_value * (eta_value + zeta_value) != 0,
                "A": eta_value * zeta_value != 0 and eta_value + zeta_value == 0,
            }
            need(sum(matches.values()) == 1, "coefficient partition")
            seen_rows.add(next(row for row, present in matches.items() if present))
    need(seen_rows == {"P", "Y", "G", "A"}, "all partition rows witnessed")

    print("JC23_K0_Y_ONLY_COMPLETE_SOURCE")
    print("scope=K=0,Delta=5696/105,eta=0,zeta!=0,Phi/Theta arbitrary")
    print("source_pair=(A,B);degrees=(6,3);infinity=(3*zeta*p^2,9*zeta*p^2)")
    print("hessian_bridge=p*detD(A,B)=2*t^2*detHess mod(A,B)")
    print("resultant=p^6*R20;R0=1259712*zeta^3")
    print("RLC=(40870620168192/1225)*zeta^7")
    print("universal=T0:2,T-1/6:2;critical_length=24")
    print("polygon=((0,1),(2,0),(5,3),(3,4),(0,5));Pick=(30,10,11)")
    print("packet=(11,8,2,2,2,1);full_n=26;defect=20")
    print("carrier=zeta*W^3=q-1/2;prime_owner_degree=3;finite=(n,beta)=(20,3)")
    print("responses=finite:18<19,identity:3<19,full:4<20")
    print("partition=P:4192,Y:4205,G:4209,A:4176")
    print(f"checks={CHECKS}")
    print("K0_Y_ONLY_COMPLETE_ACCEPT")


if __name__ == "__main__":
    main()
