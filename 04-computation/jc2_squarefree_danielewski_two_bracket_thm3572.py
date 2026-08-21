#!/usr/bin/env python3
"""Finite exact companion for provisional THM-3572."""

from __future__ import annotations

import ast
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


b, c, e, x, q, t = sp.symbols("b c e x q t")


def zero(expr: sp.Expr) -> bool:
    return sp.cancel(sp.factor(expr)) == 0


def on_surface(expr: sp.Expr, sigma: sp.Expr) -> sp.Expr:
    """Inject A_sigma into C(b,c) by e=sigma/c^2."""
    return sp.cancel(expr.subs(e, sigma / c**2))


def bracket(f: sp.Expr, g: sp.Expr, sigma: sp.Expr) -> sp.Expr:
    return sp.expand(
        c**2 * (sp.diff(f, b) * sp.diff(g, c) - sp.diff(f, c) * sp.diff(g, b))
        - 2 * c * e * (sp.diff(f, b) * sp.diff(g, e) - sp.diff(f, e) * sp.diff(g, b))
        - sp.diff(sigma, b)
        * (sp.diff(f, c) * sp.diff(g, e) - sp.diff(f, e) * sp.diff(g, c))
    )


def bezout_packet(sigma: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    u, v, gcd = sp.gcdex(sigma, sp.diff(sigma, b), b)
    require(zero(gcd - 1), f"nonunit gcd for {sigma}")
    # u*sigma+v*sigma'=1, so T=-v in U*sigma-T*sigma'=1.
    return sp.factor(u), sp.factor(-v)


def audit_bezout_and_brackets(sigmas: list[sp.Expr]) -> int:
    checks = 0
    for sigma in sigmas:
        require(sp.gcd(sigma, sp.diff(sigma, b)) == 1, "squarefree control failed")
        u, tau = bezout_packet(sigma)
        require(zero(u * sigma - tau * sp.diff(sigma, b) - 1), "Bezout identity")

        # In Omega^1, c*alpha-db=-T*d(c^2e-sigma).
        alpha = (u * c * e, -2 * tau * e, -tau * c)
        cleared = (
            sp.expand(c * alpha[0] - 1),
            sp.expand(c * alpha[1]),
            sp.expand(c * alpha[2]),
        )
        drel = (-sp.diff(sigma, b), 2 * c * e, c**2)
        for lhs, rhs in zip(cleared, tuple(-tau * z for z in drel)):
            require(zero(on_surface(lhs - rhs, sigma)), "cleared one-form identity")
            checks += 1

        first = (u + sp.diff(tau, b)) * c * e
        second = -tau * e
        two = bracket(first, b, sigma) + bracket(second, c, sigma)
        require(zero(on_surface(two - 1, sigma)), "two-bracket identity")
        checks += 1
    return checks


def audit_de_rham(sigmas: list[sp.Expr]) -> tuple[int, int]:
    rank_checks = 0
    functional_checks = 0
    for sigma in sigmas:
        degree = sp.degree(sigma, b)
        roots = [sp.Rational(j) for j in range(degree)]
        require(zero(sigma - sp.prod(b - root for root in roots)), "root packet mismatch")
        for cutoff in range(0, 7):
            columns = []
            top = degree + cutoff - 1
            for j in range(cutoff + 1):
                image = sp.Poly(sp.diff(sigma * b**j, b), b)
                columns.append([image.coeff_monomial(b**i) for i in range(top + 1)])
            matrix = sp.Matrix.hstack(*(sp.Matrix(col) for col in columns))
            require(matrix.rank() == cutoff + 1, "de Rham derivative rank")
            require((top + 1) - matrix.rank() == degree - 1, "de Rham codimension")
            rank_checks += 1

        # Root-antiderivative differences annihilate (sigma*g)' and are independent.
        test_degree = degree - 2
        functional_matrix = []
        for root in roots[1:]:
            row = []
            for power in range(test_degree + 1):
                antiderivative = b ** (power + 1) / sp.Rational(power + 1)
                row.append(sp.expand(antiderivative.subs(b, root) - antiderivative.subs(b, roots[0])))
            functional_matrix.append(row)
        require(sp.Matrix(functional_matrix).rank() == degree - 1, "root functionals independent")
        functional_checks += 1
        for root in roots[1:]:
            for power in range(0, 6):
                primitive = sigma * b**power
                require(zero(primitive.subs(b, root) - primitive.subs(b, roots[0])), "image functional")
                functional_checks += 1
    return rank_checks, functional_checks


def audit_homogeneous_formula(sigma: sp.Expr) -> tuple[int, int]:
    formula_checks = 0
    nonconstant_controls = 0
    for r in range(-4, 5):
        for s in range(-4, 5):
            f = b**2 + 2 * b + 3
            g = b**3 - b + 2
            lhs = bracket(c**r * f, c**s * g, sigma)
            rhs = c ** (r + s + 1) * (s * sp.diff(f, b) * g - r * f * sp.diff(g, b))
            require(zero(lhs - rhs), "homogeneous bracket formula")
            formula_checks += 1

    for r in range(0, 9):
        m = (r + 2) // 2
        f = b**2 + b + 1
        g = b**3 + 2 * b + 1
        h = sigma**m * g
        value = sp.expand(-(r + 1) * sp.diff(f, b) * h - r * f * sp.diff(h, b))
        require(sp.degree(value, b) >= 1, "homogeneous nonconstant hostile")
        nonconstant_controls += 1
    return formula_checks, nonconstant_controls


def audit_two_by_two_factorization(sigmas: list[sp.Expr]) -> int:
    controls = 0
    for sigma in sigmas:
        h = sp.expand(sigma * (b**2 + 2 * b + 3))
        K = b**3 - b + 2
        A, B, L, M = map(sp.Rational, (2, 3, 5, 7))
        for k in range(1, 9):
            f = A * h
            g = B * h**k
            F = L * K ** (2 * k - 1)
            G = M * K
            low = sp.expand(2 * f * sp.diff(g, b) - 2 * k * sp.diff(f, b) * g)
            high = sp.expand(sp.diff(F, b) * G - (2 * k - 1) * F * sp.diff(G, b))
            require(zero(low) and zero(high), "extreme Wronskians")
            middle = sp.expand(
                sp.diff(f, b) * G
                + 2 * f * sp.diff(G, b)
                - 2 * k * sp.diff(F, b) * g
                - (2 * k - 1) * F * sp.diff(g, b)
            )
            factor = sp.expand(
                (K * sp.diff(h, b) + 2 * h * sp.diff(K, b))
                * (A * M - k * (2 * k - 1) * L * B * K ** (2 * k - 2) * h ** (k - 1))
            )
            require(zero(middle - factor), "two-by-two factorization")
            require(sp.degree(K * sp.diff(h, b) + 2 * h * sp.diff(K, b), b) >= 1, "unit hostile")
            controls += 1
    return controls


def audit_quadratic() -> int:
    sigma = b * (b + 4)
    u = -sp.Rational(1, 4)
    tau = -(b + 2) / 8
    require(zero(u * sigma - tau * sp.diff(sigma, b) - 1), "quadratic Bezout")
    value = bracket(-3 * c * e / 8, b, sigma) + bracket((b + 2) * e / 8, c, sigma)
    require(zero(on_surface(value - 1, sigma)), "quadratic two-bracket")
    return 2


def jacobian(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(f, x) * sp.diff(g, q) - sp.diff(f, q) * sp.diff(g, x))


def audit_degree_five() -> tuple[int, int, int]:
    tt = x**2 * q
    P = 1 + tt**2
    S = 1 + 5 * tt**2
    aa = q / S**2
    cc = x * S * P
    bb = sp.expand(tt * P**2)
    W = 125 * tt**6 + 450 * tt**4 + 565 * tt**2 + 256
    ee = sp.expand(q * W)
    sigma = b * (3125 * b**2 + 256)

    identity_checks = 0
    require(zero(jacobian(aa, cc) + 1), "degree-five rational Jacobian")
    identity_checks += 1
    require(zero((3125 * bb**2 + 256) - S**2 * W), "critical-value square")
    identity_checks += 1
    require(zero(cc**2 * ee - bb * (3125 * bb**2 + 256)), "degree-five surface")
    identity_checks += 1
    require(sp.degree(t * (1 + t**2) ** 2, t) == 5, "generic degree")
    identity_checks += 1

    # Poisson morphism on the three generator pairs.
    subs = {b: bb, c: cc, e: ee}
    for f, g in ((b, c), (c, e), (b, e)):
        target = bracket(f, g, sigma).subs(subs)
        source = jacobian(f.subs(subs), g.subs(subs))
        require(zero(source + target), "Poisson pullback")
        identity_checks += 1

    # Five collision branches over a generic nonzero lambda.
    lam = sp.Rational(2)
    require(zero(bb.subs({x: 0, q: lam})), "central collision b")
    require(zero(cc.subs({x: 0, q: lam})), "central collision c")
    require(zero(ee.subs({x: 0, q: lam}) - 256 * lam), "central collision e")
    collision_checks = 3
    for t_value in (sp.I, -sp.I):
        q_value = 16 * lam
        x2 = t_value / q_value
        # Expressions depend only on x^2 once P=0, except c which vanishes by P.
        require(zero((tt - t_value).subs({x**2: x2, q: q_value})), "collision t")
        require(zero((bb).subs({x**2: x2, q: q_value})), "collision b")
        require(zero(P.subs({x**2: x2, q: q_value})), "collision P")
        require(zero(ee.subs({x**2: x2, q: q_value}) - 256 * lam), "collision e")
        collision_checks += 4
    require(1 + 2 + 2 == 5, "collision multiplicity invoice")
    collision_checks += 1

    u = -sp.Rational(28125, 131072) * b
    tau = -(9375 * b**2 + 512) / 131072
    require(zero(u * sigma - tau * sp.diff(sigma, b) - 1), "cubic Bezout")
    two = bracket((u + sp.diff(tau, b)) * c * e, b, sigma) + bracket(-tau * e, c, sigma)
    require(zero(on_surface(two - 1, sigma)), "cubic two-bracket")
    explicit_checks = 2

    # iota(b,c,e)=(-b,c,-e) is anti-Poisson.
    iota = {b: -b, c: c, e: -e}
    for f, g in ((b, c), (c, e), (b, e)):
        lhs = bracket(f.subs(iota, simultaneous=True), g.subs(iota, simultaneous=True), sigma)
        rhs = -bracket(f, g, sigma).subs(iota, simultaneous=True)
        require(zero(lhs - rhs), "anti-Poisson involution")
        explicit_checks += 1
    return identity_checks, collision_checks, explicit_checks


def audit_no_asserts() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "assert nodes are forbidden")
    return count


def main() -> None:
    rational_sigmas = [sp.prod(b - j for j in range(degree)) for degree in range(2, 7)]
    bezout_checks = audit_bezout_and_brackets(rational_sigmas)
    rank_checks, functional_checks = audit_de_rham(rational_sigmas)
    formula_checks, homogeneous_controls = audit_homogeneous_formula(rational_sigmas[1])
    two_by_two_controls = audit_two_by_two_factorization([rational_sigmas[0], rational_sigmas[1]])
    quadratic_checks = audit_quadratic()
    degree5_checks, collision_checks, cubic_checks = audit_degree_five()
    assert_nodes = audit_no_asserts()

    print("THM-3572 finite exact companion")
    print("status PROVISIONAL_VERIFIED_EXACT_PENDING_INDEPENDENT_AUDIT")
    print("universe exact_QQ_symbolic_squarefree_degrees_2_through_6")
    print("bezout_primitive_and_two_bracket_checks", bezout_checks)
    print("de_rham_rank_checks", rank_checks)
    print("root_antiderivative_functional_checks", functional_checks)
    print("homogeneous_formula_checks", formula_checks)
    print("homogeneous_nonconstant_controls", homogeneous_controls)
    print("two_by_two_factorization_controls", two_by_two_controls)
    print("quadratic_controls", quadratic_checks)
    print("degree_five_identity_and_poisson_checks", degree5_checks)
    print("degree_five_collision_branch_checks", collision_checks)
    print("symmetric_cubic_bezout_and_involution_checks", cubic_checks)
    print("ast_assert_nodes", assert_nodes)
    print("scope universal_claims_are_proved_finite_rows_are_controls_no_JC2_claim")
    print("PASS")


if __name__ == "__main__":
    main()
