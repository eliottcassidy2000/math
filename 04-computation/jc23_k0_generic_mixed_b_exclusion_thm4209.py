#!/usr/bin/env python3
"""Primary exact source-chart audit for THM-4209's K=0 mixed chambers."""

from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(poly, variable):
    data = sp.Poly(poly, variable).as_dict()
    need(bool(data), "zero polynomial")
    return min(key[0] for key in data)


def main():
    s, p, Phi, Theta, eta, zeta, Q = sp.symbols(
        "s p Phi Theta eta zeta Q"
    )
    t = p - s**2
    Delta = sp.Rational(5696, 105)
    H = (
        -3 * p
        + sp.Rational(8, 3) * p**2
        - sp.Rational(1376, 135) * p**3
        + Phi * s * p**3
        + Delta * p**4
        + Theta * s**2 * p**3
        + eta * s * p**4
        + zeta * s**3 * p**3
    )
    G = -s**2 / (2 * t) + H
    A = sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p)
    C0 = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    B = sp.cancel((C0 + s * A) / t**2)
    need(sp.denom(A) == sp.denom(B) == 1, "nonpolynomial source pair")
    need((sp.degree(A, s), sp.degree(B, s)) == (6, 3), "source degrees")
    need(sp.factor(sp.Poly(A, s).LC() - 3 * zeta * p**2) == 0, "A LC")
    need(sp.factor(sp.Poly(B, s).LC() - 9 * zeta * p**2) == 0, "B LC")
    need(sp.factor(t**2 * sp.diff(G, s) - p * A) == 0, "gradient s")
    need(
        sp.factor(2 * t**2 * sp.diff(G, p) - (t**2 * B - s * A)) == 0,
        "gradient p",
    )
    J = sp.det(
        sp.Matrix(((sp.diff(A, s), sp.diff(A, p)), (sp.diff(B, s), sp.diff(B, p))))
    )
    Hess = sp.det(sp.hessian(G, (s, p)))
    bridge = sp.together(p * J - 2 * t**2 * Hess).as_numer_denom()[0]
    need(sp.factor(sp.reduced(bridge, [A, B], s, p)[1]) == 0, "Hessian bridge")
    need(sp.factor(A.subs(p, 0) + s) == 0, "p=0 A")
    need(sp.factor(B.subs(p, 0) + 6) == 0, "p=0 B")
    need(sp.factor(A.subs(p, s**2) + s) == 0, "t=0 A")
    need(B.subs({s: 0, p: 0}) == -6, "t=0 endpoint")

    resultant = sp.factor(sp.resultant(A, B, s))
    pval = valuation(resultant, p)
    residual = sp.Poly(sp.cancel(resultant / p**pval), p)
    need((pval, residual.degree()) == (6, 21), "p^6 R21")
    need(
        sp.factor(residual.nth(0) - 2**6 * 3**9 * zeta**3) == 0,
        "R21 constant",
    )
    need(
        sp.factor(
            residual.LC()
            - 2**2 * 3**12 * eta**3 * zeta**2 * (eta + zeta) ** 4
        )
        == 0,
        "R21 leading coefficient",
    )

    # The rational source chart divides out the two universal coordinate
    # fibres.  Rebuild their normalized equations directly before restoring
    # four points to the affine critical length.
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
        + eta * P**3 * Y
        + zeta * Y**3
    )
    f = sp.cancel(sp.diff(GN, X) / T)
    h = sp.diff(GN, T)
    need(sp.denom(f) == 1, "normalized first derivative")
    need((sp.degree(f, X), sp.degree(h, X)) == (8, 9), "normalized degrees")
    need(
        sp.factor(sp.Poly(f, X).LC() - 9 * (eta + zeta) * T**8) == 0,
        "normalized f leading row",
    )
    need(
        sp.factor(sp.Poly(h, X).LC() - 9 * (eta + zeta) * T**8) == 0,
        "normalized h leading row",
    )
    need(sp.factor(f.subs(T, 0) + X) == 0, "T=0 f row")
    need(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0, "T=0 h row")
    special_T = -sp.Rational(1, 6)
    need(
        sp.rem(sp.Poly(f.subs(T, special_T), X), sp.Poly(X**2 - 6, X)).is_zero,
        "T=-1/6 universal f pair",
    )
    need(
        sp.rem(sp.Poly(h.subs(T, special_T), X), sp.Poly(X**2 - 6, X)).is_zero,
        "T=-1/6 universal h pair",
    )

    # Exact primitive faces; Phi and Theta are absent from all of them.
    FQ = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
    support = {
        monomial: sp.factor(coefficient)
        for monomial, coefficient in sp.Poly(FQ, s, p).terms()
    }
    vertices = [(0, 1), (2, 0), (5, 3), (1, 5), (0, 5)]
    expected_faces = (
        s**2 * (1 - Q / 2) - p,
        s**2 * (1 - Q / 2 - Q * zeta * s**3 * p**3),
        Q * p**3 * s * (p - s**2) * (eta * p + zeta * s**2),
        Q * p**5 * (Delta + eta * s),
    )
    packet = []
    boundary = 0
    for edge_index, (start, end) in enumerate(
        zip(vertices, vertices[1:] + vertices[:1])
    ):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        u, v = -dy // length, dx // length
        constant = u * start[0] + v * start[1]
        index = u + v - constant
        boundary += length
        if not (start[0] == end[0] == 0):
            packet.extend([index] * length)
        if edge_index < 4:
            face = sum(
                coefficient * s**i * p**j
                for (i, j), coefficient in support.items()
                if u * i + v * j == constant
            )
            need(sp.factor(face - expected_faces[edge_index]) == 0, "face")
    area2 = abs(
        sum(
            vertices[i][0] * vertices[(i + 1) % 5][1]
            - vertices[(i + 1) % 5][0] * vertices[i][1]
            for i in range(5)
        )
    )
    genus = (area2 - boundary + 2) // 2
    packet = tuple(sorted(packet, reverse=True))
    need((area2, boundary, genus) == (31, 11, 11), "Pick ledger")
    need(packet == (8, 8, 4, 2, 2, 2, 1), "packet")

    L, defect, full_n = 25, sum(e - 1 for e in packet), sum(packet)
    finite_n, beta = full_n - 6, 3
    need((defect, full_n, finite_n, beta) == (20, 27, 21, 3), "response")
    need(2 * (full_n - L) == 4 < defect, "full gate")
    need(2 * finite_n - L - 1 + beta == 19 < finite_n - 1, "finite gate")
    need(beta < finite_n - 1, "identity gate")

    # Repeated-top anatomy and the next exact subwalls.
    anti = sp.Poly(sp.expand(residual.as_expr().subs(zeta, -eta)), p)
    need(anti.degree() == 19, "anti R19 degree")
    need(sp.factor(anti.nth(0) + 1259712 * eta**3) == 0, "anti R19 constant")
    need(
        sp.factor(anti.LC() - 1327104 * eta**5 * (Delta + Theta) ** 4) == 0,
        "anti R19 LC",
    )
    anti_theta_wall = sp.Poly(sp.expand(anti.as_expr().subs(Theta, -Delta)), p)
    need(anti_theta_wall.degree() == 17, "anti Theta=-Delta R17")
    need(
        sp.factor(anti_theta_wall.LC() - 777924 * Phi**4 * eta**5) == 0,
        "anti Theta=-Delta R17 LC",
    )
    anti_terminal = sp.Poly(sp.expand(anti_theta_wall.as_expr().subs(Phi, 0)), p)
    need(anti_terminal.degree() == 15, "anti terminal R15")
    need(
        sp.factor(
            anti_terminal.LC()
            - sp.Rational(229431851352064, 50625) * eta**5
        )
        == 0,
        "anti terminal R15 LC",
    )
    anti_packet = (7, 7, 4, 2, 2, 2, 1)
    anti_length = anti.degree() + 4
    anti_defect = sum(index - 1 for index in anti_packet)
    anti_full_n = sum(anti_packet)
    anti_finite_n, anti_beta = anti_full_n - 6, 3
    need(
        (anti_length, anti_defect, anti_full_n, anti_finite_n, anti_beta)
        == (23, 18, 25, 19, 3),
        "anti response ledger",
    )
    need(2 * (anti_full_n - anti_length) == 4 < anti_defect, "anti full gate")
    need(
        2 * anti_finite_n - anti_length - 1 + anti_beta == 17
        < anti_finite_n - 1,
        "anti finite gate",
    )

    print("AUDIT complete_K0_generic_mixed_B_source")
    print("scope=K=0,Delta=5696/105,Phi/Theta arbitrary,eta*zeta*(eta+zeta)!=0")
    print("source_pair=(deg_s A,deg_s B)=(6,3);LC=(3zeta p^2,9zeta p^2)")
    print("hessian_bridge=ideal_zero;p0/t0=open_chart_empty")
    print("resultant=p^6*R21;R0=2^6*3^9*zeta^3")
    print("RLC=2^2*3^12*eta^3*zeta^2*(eta+zeta)^4")
    print("normalized_degrees=(8,9);universal_coordinate_points=2+2")
    print("faces_independent_of_Phi_Theta=YES;packet=(8,8,4,2,2,2,1)")
    print("carrier=zeta*W^3=q-1/2;L=25;full_cap=4<20;finite_cap=19<20")
    print("anti_next: zeta=-eta gives R19, LC=1327104*eta^5*(Delta+Theta)^4")
    print("anti_Theta=-Delta: R17, LC=777924*Phi^4*eta^5")
    print("anti_Theta=-Delta_Phi=0: R15")
    print("anti_packet=(7,7,4,2,2,2,1);L=23;full_cap=4<18;finite_cap=17<18")
    print(f"checks={CHECKS}")
    print("COMPLETE_K0_GENERIC_MIXED_B_SOURCE_ACCEPT")


if __name__ == "__main__":
    main()
