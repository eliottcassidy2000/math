#!/usr/bin/env python3
"""Independent exact audit of the repeated-top-wall critical schemes.

This file deliberately does not import the primary implementation.
Finite critical schemes are counted by Groebner quotient algebras and certified
reduced/separated by a primitive linear element, not by either producer
coordinate resultant.  The named row-D coordinate resultants are then rebuilt
directly only to audit their symbolic projection portfolio.
"""

from hashlib import sha256
from pathlib import Path

import sympy as sp


X, T, U = sp.symbols("X T U")
P = T + X**2 * T**2
Y = X * T * P


def source(delta, theta, phi, eta):
    kappa = sp.Rational(2848, 45) - sp.Rational(7, 6) * delta
    return sp.expand(
        -X**2*T/2 - 3*P + sp.Rational(8, 3)*P**2
        - sp.Rational(1376, 135)*P**3 + kappa*Y**2
        + phi*P**2*Y + delta*P**4 + theta*P*Y**2
        + eta*P**3*Y - eta*Y**3
    )


CONTROLS = (
    ("A", sp.Rational(1), sp.Rational(19, 11), sp.Rational(11, 7), sp.Rational(23, 13), 23),
    ("B", sp.Rational(1), -sp.Rational(1), sp.Rational(11, 7), sp.Rational(23, 13), 21),
    ("C", sp.Rational(1), -sp.Rational(1), sp.Rational(0), sp.Rational(23, 13), 19),
    ("D", sp.Rational(2048, 45), -sp.Rational(2048, 45), sp.Rational(0), sp.Rational(23, 13), 16),
)


def standard_monomials(gb):
    lms = [tuple(poly.LM(order=gb.order).exponents) for poly in gb.polys]
    pure_x = min(a for a, b in lms if b == 0)
    pure_t = min(b for a, b in lms if a == 0)
    result = []
    for i in range(pure_x):
        for j in range(pure_t):
            if not any(i >= a and j >= b for a, b in lms):
                result.append((i, j))
    return tuple(result), tuple(lms)


def quotient_critical_audit(label, delta, theta, phi, eta, expected):
    g = source(delta, theta, phi, eta)
    gx, gt = sp.diff(g, X), sp.diff(g, T)
    print(f"START_QUOTIENT {label}", flush=True)
    gb = sp.groebner((gx, gt), X, T, order="grevlex", domain=sp.QQ,
                     method="f5b")
    basis, lms = standard_monomials(gb)
    assert len(basis) == expected, (label, len(basis), expected, lms)
    index = {exponent: i for i, exponent in enumerate(basis)}
    # Try small integral slopes until a primitive/separating element is found.
    chosen = None
    for slope in (1, 2, 3, 5, 7, 11):
        matrix = sp.zeros(expected)
        for col, (i, j) in enumerate(basis):
            _, remainder = gb.reduce(sp.expand(X**i * T**j * (X + slope*T)))
            rpoly = sp.Poly(remainder, X, T, domain=sp.QQ)
            for exponent, coefficient in rpoly.terms():
                matrix[index[exponent], col] = coefficient
        cp = sp.Poly(matrix.charpoly(U).as_expr(), U, domain=sp.QQ)
        if sp.gcd(cp, cp.diff()).degree() == 0:
            chosen = slope, cp
            break
    assert chosen is not None, label
    slope, cp = chosen
    # Universal pairs: T=0, X^2=-6 and T=-1/6, X^2=6.
    universal_zero = sp.Poly(U**2 + 6, U, domain=sp.QQ)
    universal_half = sp.Poly((U + sp.Rational(slope, 6))**2 - 6, U, domain=sp.QQ)
    universal = universal_zero * universal_half
    residual, rem = sp.div(cp, universal)
    assert rem.is_zero
    assert residual.degree() == expected - 4
    assert sp.gcd(residual, universal).degree() == 0
    assert sp.gcd(residual, residual.diff()).degree() == 0
    payload = ";".join(str(q) for q in cp.monic().all_coeffs())
    digest = sha256(payload.encode()).hexdigest()
    return {
        "label": label,
        "length": len(basis),
        "residual": residual.degree(),
        "slope": slope,
        "primitive_sha256": digest,
        "leading_monomials": lms,
    }


def monic_digest(poly):
    p = sp.Poly(poly).monic()
    payload = f"degree={p.degree()};" + ",".join(
        f"{sp.Rational(c).p}/{sp.Rational(c).q}" for c in p.all_coeffs()
    )
    return sha256(payload.encode()).hexdigest()


def symbolic_monic_digest(poly, variable):
    p = sp.Poly(poly, variable).monic()
    payload = f"degree={p.degree()};" + ",".join(str(c) for c in p.all_coeffs())
    return sha256(payload.encode()).hexdigest()


def valuation(expr, variable):
    return min(key[0] for key in sp.Poly(expr, variable).as_dict())


def deep_portfolio():
    eta = sp.symbols("eta")
    delta = sp.Rational(2048, 45)
    theta = -delta
    kappa = sp.Rational(1376, 135)
    g = source(delta, theta, sp.Integer(0), eta)
    f = sp.cancel(sp.diff(g, X) / T)
    h = sp.diff(g, T)
    raw_q = sp.resultant(f, h, X)
    assert valuation(raw_q, T) == 35
    q_expr = sp.cancel(raw_q / (T**35 * (6*T + 1)**2))
    q12 = sp.Poly(q_expr, T)
    assert q12.degree() == 12
    J = 22143375*eta**2 + 15510536192
    assert sp.factor(q12.nth(0) + sp.Rational(3**11, 2**5)*eta**5) == 0
    assert sp.factor(q12.LC() - eta**3*J**2/sp.Integer(3690562500)) == 0

    s, p = sp.symbols("s p")
    t = p - s**2
    H = sp.expand(
        -3*p + sp.Rational(8, 3)*p**2 - sp.Rational(1376, 135)*p**3
        + kappa*s**2*p**2 + delta*p**4 + theta*s**2*p**3
        + eta*s*p**3*(p-s**2)
    )
    # Clear the rational chart derivatives independently from the direct G(s,p).
    Gsp = -s**2/(2*t) + H
    A = sp.cancel(t**2 * sp.diff(Gsp, s) / p)
    C = sp.cancel(2*t**2 * sp.diff(Gsp, p))
    assert sp.denom(A) == 1 and sp.denom(C) == 1
    raw_r = sp.resultant(A, C, s)
    assert valuation(raw_r, p) == 8
    r12 = sp.Poly(sp.cancel(raw_r / p**8), p)
    assert r12.degree() == 12
    assert sp.factor(r12.nth(0) + sp.Rational(64, 1125)*eta*J) == 0
    assert sp.factor(r12.LC() + 531441*eta**7) == 0

    q_disc = sp.factor_list(sp.discriminant(q12.as_expr(), T), eta)[1]
    r_disc = sp.factor_list(sp.discriminant(r12.as_expr(), p), eta)[1]
    q_by_degree = {sp.degree(f, eta): (m, sp.Poly(f, eta)) for f, m in q_disc}
    r_by_degree = {sp.degree(f, eta): (m, sp.Poly(f, eta)) for f, m in r_disc}
    assert sorted((d, m) for d, (m, _) in q_by_degree.items()) == [(1, 32), (30, 1), (36, 2)]
    assert sorted((d, m) for d, (m, _) in r_by_degree.items()) == [(1, 40), (24, 2), (30, 1)]
    h30q = q_by_degree[30][1].monic()
    h30r = r_by_degree[30][1].monic()
    assert h30q == h30r
    h36 = q_by_degree[36][1]
    h24 = r_by_degree[24][1]
    q_univ_num = sp.together(q12.eval(-sp.Rational(1, 6))).as_numer_denom()[0]
    gcd_24_36 = sp.gcd(h24, h36).degree()
    gcd_24_univ = sp.gcd(h24, sp.Poly(q_univ_num, eta)).degree()
    gcd_j_30 = sp.gcd(sp.Poly(J, eta), h30q).degree()
    assert (gcd_24_36, gcd_24_univ, gcd_j_30) == (0, 0, 0)

    # Audit the formerly excluded endpoint wall J(eta)=0 inside the exact
    # quadratic field.  Reduce coefficients first, then factor over either
    # conjugate embedding of Q[eta]/(J).
    jpoly = sp.Poly(J, eta, domain=sp.QQ)
    qj_expr = sp.Integer(0)
    for i in range(q12.degree() + 1):
        coefficient = sp.rem(sp.Poly(q12.nth(i), eta, domain=sp.QQ), jpoly).as_expr()
        qj_expr += coefficient*T**i
    qj = sp.Poly(qj_expr, T, domain=sp.QQ.frac_field(eta))
    assert qj.degree() == 11
    eta0 = sp.sqrt(-sp.Rational(15510536192, 22143375))
    qj_ext = sp.Poly(sp.expand(qj.as_expr().subs(eta, eta0)), T,
                     extension=eta0)
    assert sp.gcd(qj_ext, qj_ext.diff()).degree() == 0
    qj_factors = sp.factor_list(qj_ext.as_expr(), T, extension=eta0)[1]
    factor_degrees = [
        (sp.Poly(factor, T, extension=eta0).degree(), multiplicity)
        for factor, multiplicity in qj_factors
    ]
    assert factor_degrees == [(11, 1)]
    qj_at_universal = sp.factor(qj.eval(-sp.Rational(1, 6)))
    assert qj.nth(0) != 0 and qj_at_universal != 0 and qj.LC() != 0

    pf = sp.Poly(f, X)
    ph = sp.Poly(h, X)
    assert pf.degree() == 6 and sp.factor(pf.LC() - 7*T**7*eta) == 0
    assert ph.degree() == 7 and sp.factor(ph.LC() - 8*T**7*eta) == 0
    hessian = sp.det(sp.hessian(g, (X, T)))
    for t_value, x_square, g_value, determinant in (
        (-sp.Rational(1, 6), sp.Rational(6), sp.Rational(1, 2), -sp.Rational(6)),
        (sp.Rational(0), -sp.Rational(6), sp.Rational(0), sp.Rational(6)),
    ):
        modulus = sp.Poly(X**2-x_square, X)
        for expression in (sp.diff(g, X), sp.diff(g, T), g-g_value,
                           hessian-determinant):
            remainder = sp.rem(sp.Poly(sp.expand(expression.subs(T, t_value)), X),
                               modulus)
            assert remainder.is_zero

    jwall = {
        "degree": qj.degree(),
        "constant": str(sp.factor(qj.nth(0))),
        "leading": str(sp.factor(qj.LC())),
        "at_minus_one_sixth": str(qj_at_universal),
        "gcd_degree": sp.gcd(qj_ext, qj_ext.diff()).degree(),
        "factor_degrees": factor_degrees,
        "QJ_sha256": symbolic_monic_digest(qj.as_expr(), T),
        "f_X": (pf.degree(), str(sp.factor(pf.LC()))),
        "h_X": (ph.degree(), str(sp.factor(ph.LC()))),
        "critical_length": 11+2+2,
        "finite": {
            "n": 12, "beta": 3, "rational_origin_defect": 9,
            "both_nonidentity_capacity": 2*12-(11+2+2)-2+3,
            "transitive_threshold": 11,
            "identity_handle_carrier_cap": 3,
        },
        "full": {
            "n": 18, "overlap": 18-(11+2+2),
            "commutator_cap": 2*(18-(11+2+2)),
            "origin_defect": 12,
        },
    }

    # Formal Boolean audit of the claimed switch: under eta,J,H30 nonzero,
    # R works whenever H24!=0; otherwise the two gcds force H36 and Q(-1/6)
    # nonzero, so Q works.
    return {
        "J": str(J),
        "Q_constant": str(sp.factor(q12.nth(0))),
        "Q_leading": str(sp.factor(q12.LC())),
        "R_constant": str(sp.factor(r12.nth(0))),
        "R_leading": str(sp.factor(r12.LC())),
        "Q_factor_degrees": sorted((d, m) for d, (m, _) in q_by_degree.items()),
        "R_factor_degrees": sorted((d, m) for d, (m, _) in r_by_degree.items()),
        "H30_sha256": monic_digest(h30q),
        "H36_sha256": monic_digest(h36),
        "H24_sha256": monic_digest(h24),
        "gcd_degrees": (gcd_24_36, gcd_24_univ, gcd_j_30),
        "Q_universal_degree": sp.degree(q_univ_num, eta),
        "J_wall": jwall,
    }


def main():
    rows = []
    for control in reversed(CONTROLS):
        row = quotient_critical_audit(*control)
        rows.append(row)
        print("QUOTIENT", row, flush=True)
    print("START_DEEP", flush=True)
    portfolio = deep_portfolio()
    print("DEEP", portfolio)


if __name__ == "__main__":
    main()
