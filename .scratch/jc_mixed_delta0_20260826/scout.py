#!/usr/bin/env python3
"""Exact scout for the mixed off-antidiagonal Delta=0 exact-M=9 wall.

This is session scratch, not a canonical certificate.  It specializes the
source before elimination, checks chart loss and Hessian transport, and
reconstructs the contracted Newton polygon and all primitive faces.
"""

from math import gcd

import sympy as sp


def need(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly, variable):
    terms = sp.Poly(poly, variable).terms()
    need(terms, "zero polynomial")
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
    genus = (area2 - boundary + 2) // 2
    return vertices, area2, boundary, genus


def primitive_faces(poly, s, p):
    vertices = convex_hull(sp.Poly(poly, s, p).monoms())
    data = []
    packet = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        u, v = -dy // length, dx // length
        level = u * start[0] + v * start[1]
        index = u + v - level
        face = sum(
            coefficient * s**i * p**j
            for (i, j), coefficient in sp.Poly(poly, s, p).terms()
            if u * i + v * j == level
        )
        data.append((start, end, (u, v), length, level, index,
                     sp.factor(face)))
        if not (start[0] == end[0] == 0):
            packet.extend([index] * length)
    return vertices, tuple(data), tuple(sorted(packet, reverse=True))


def main():
    s, p, Phi, Theta, eta, zeta, Q = sp.symbols(
        "s p Phi Theta eta zeta Q"
    )
    epsilon = -sp.Rational(1376, 135)
    k0 = sp.Rational(2848, 45)
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
    need(sp.denom(A) == sp.denom(B) == 1, "nonpolynomial source pair")
    need(sp.factor(t**2*sp.diff(G, s) - p*A) == 0,
         "first gradient identity")
    need(sp.factor(2*t**2*sp.diff(G, p) - (t**2*B-s*A)) == 0,
         "second gradient identity")
    J = sp.det(sp.Matrix(((sp.diff(A, s), sp.diff(A, p)),
                          (sp.diff(B, s), sp.diff(B, p)))))
    Hess = sp.det(sp.hessian(G, (s, p)))
    bridge = sp.together(p*J - 2*t**2*Hess).as_numer_denom()[0]
    need(sp.factor(sp.reduced(bridge, [A, B], s, p)[1]) == 0,
         "Hessian bridge")
    need(sp.factor(A.subs(p, 0) + s) == 0, "p=0 A")
    need(sp.factor(B.subs(p, 0) + 6) == 0, "p=0 B")
    need(sp.factor(A.subs(p, s**2) + s) == 0, "t=0 A")

    resultant = sp.factor(sp.resultant(A, B, s))
    pval = valuation(resultant, p)
    residual = sp.Poly(sp.cancel(resultant/p**pval), p)
    print("source_s_degrees", sp.degree(A, s), sp.degree(B, s))
    print("source_lc", sp.factor(sp.Poly(A, s).LC()),
          sp.factor(sp.Poly(B, s).LC()))
    print("resultant", "pval", pval, "degree", residual.degree())
    print("resultant_TC", sp.factor(residual.TC()))
    print("resultant_LC", sp.factor(residual.LC()))

    # Inner endpoint wall.  The parametrization is bijective for zeta!=0:
    # zeta=(2*k0/3)u and Theta=3u^2.  First scout successive p-adic
    # strict transforms from the already exact general resultant; the
    # canonical certificate must additionally recompute each wall directly.
    u = sp.symbols("u")
    inner_sub = {zeta: sp.Rational(5696, 135)*u, Theta: 3*u**2}
    inner_resultant = sp.factor(resultant.subs(inner_sub))
    inner_pval = valuation(inner_resultant, p)
    inner = sp.Poly(sp.cancel(inner_resultant/p**inner_pval), p)
    print("inner_resultant", "pval", inner_pval, "degree", inner.degree())
    print("inner_TC", sp.factor(inner.TC()))
    print("inner_LC", sp.factor(inner.LC()))

    J_inner = 8544*Phi - 22784*u - 1215*u**3
    phi_J = u*(1215*u**2+22784)/8544
    j_resultant = sp.factor(inner_resultant.subs(Phi, phi_J))
    j_pval = valuation(j_resultant, p)
    j_residual = sp.Poly(sp.cancel(j_resultant/p**j_pval), p)
    print("J_resultant", "pval", j_pval, "degree", j_residual.degree())
    print("J_TC", sp.factor(j_residual.TC()))
    print("J_LC", sp.factor(j_residual.LC()))
    S_inner = 2460375*u**4 - 204543360*u**2 + 5580439552
    M_inner = 547499520*eta + u*S_inner
    print("inner_J", J_inner)
    print("inner_S", S_inner)
    print("inner_M", M_inner)
    eta_M = -u*S_inner/sp.Integer(547499520)
    m_resultant = sp.factor(j_resultant.subs(eta, eta_M))
    m_pval = valuation(m_resultant, p)
    m_residual = sp.Poly(sp.cancel(m_resultant/p**m_pval), p)
    print("M_resultant", "pval", m_pval, "degree", m_residual.degree())
    print("M_TC", sp.factor(m_residual.TC()))
    print("M_LC", sp.factor(m_residual.LC()))
    print("M_eta", sp.factor(eta_M))
    print("M_eta_plus_zeta", sp.factor(
        eta_M + sp.Rational(5696, 135)*u))
    N_inner = 27064125*u**4 - 5739517440*u**2 + 47239069696
    A_offanti = 18225*u**4 - 1515136*u**2 - 129777664
    next_coefficient = sp.factor(m_residual.nth(1))
    next_mod_N = sp.factor(sp.rem(sp.Poly(next_coefficient, u),
                                  sp.Poly(N_inner, u)).as_expr())
    lc_mod_N = sp.factor(sp.rem(sp.Poly(m_residual.LC(), u),
                                sp.Poly(N_inner, u)).as_expr())
    print("inner_N", N_inner)
    print("N_next_mod", next_mod_N)
    print("N_LC_mod", lc_mod_N)
    for name, hostile in (("u", u), ("eta", S_inner),
                          ("eta_plus_zeta", A_offanti),
                          ("next", next_mod_N), ("LC", lc_mod_N)):
        print("N_gcd_" + name,
              sp.gcd(sp.Poly(N_inner, u), sp.Poly(hostile, u)).degree())

    FQ = sp.expand((s**2-p)*(1-Q*H)-Q*s**2/2)
    vertices, faces, packet = primitive_faces(FQ, s, p)
    print("polygon", polygon_ledger(sp.Poly(FQ, s, p).monoms()))
    for item in faces:
        print("face", item)
    print("packet", packet, "defect", sum(index-1 for index in packet),
          "degree", sum(packet))

    # Hostile coefficient controls must not change the lower hull or faces
    # claimed independent of Phi and Theta.
    for control in ({Phi: 0, Theta: 0, eta: 2, zeta: 3},
                    {Phi: 17, Theta: -19, eta: 2, zeta: -5}):
        specialized = sp.Poly(FQ.subs(control), s, p)
        print("control", control, polygon_ledger(specialized.monoms()))


if __name__ == "__main__":
    main()
