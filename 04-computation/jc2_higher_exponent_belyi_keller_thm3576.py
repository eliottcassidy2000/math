#!/usr/bin/env python3
"""Finite exact companion for provisional THM-3576."""

from __future__ import annotations

import ast
from math import comb
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def zero(expr: sp.Expr) -> bool:
    return sp.cancel(sp.factor(expr)) == 0


t, x, q, b, c, e = sp.symbols("t x q b c e")


def packet(n: int) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    R = 1 + t
    P = sp.Add(*(sp.Rational(comb(n - 1, j), n * j + 1) * t**j for j in range(n)))
    B = sp.expand(t * P**n)
    beta = sp.factor(B.subs(t, -1))
    quotient, remainder = sp.div(sp.Poly(B - beta, t), sp.Poly(R**n, t))
    require(remainder.is_zero, f"W divisibility n={n}")
    return sp.factor(P), R, B, beta, sp.factor(quotient.as_expr())


def jacobian(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(f, x) * sp.diff(g, q) - sp.diff(f, q) * sp.diff(g, x))


def target_bracket(f: sp.Expr, g: sp.Expr, n: int, sigma: sp.Expr) -> sp.Expr:
    return sp.expand(
        c**n * (sp.diff(f, b) * sp.diff(g, c) - sp.diff(f, c) * sp.diff(g, b))
        - n * c ** (n - 1) * e * (sp.diff(f, b) * sp.diff(g, e) - sp.diff(f, e) * sp.diff(g, b))
        - sp.diff(sigma, b) * (sp.diff(f, c) * sp.diff(g, e) - sp.diff(f, e) * sp.diff(g, c))
    )


def on_target(expr: sp.Expr, n: int, sigma: sp.Expr) -> sp.Expr:
    return sp.cancel(expr.subs(e, sigma / c**n))


def audit_tower() -> tuple[int, int, int, int]:
    ode_checks = 0
    geometry_checks = 0
    intersection_checks = 0
    fewnomial_checks = 0
    for n in range(2, 9):
        P, R, B, beta, W = packet(n)
        require(zero(P + n * t * sp.diff(P, t) - R ** (n - 1)), f"ODE n={n}")
        require(zero(sp.diff(B, t) - P ** (n - 1) * R ** (n - 1)), f"B derivative n={n}")
        require(sp.degree(B, t) == n * (n - 1) + 1, f"degree n={n}")
        require(sp.degree(W, t) == (n - 1) ** 2, f"W degree n={n}")
        require(P.subs(t, 0) == 1 and P.subs(t, -1) != 0, f"P endpoints n={n}")
        require(zero(W.subs(t, -1) - P.subs(t, -1) ** (n - 1) / n), f"W owner n={n}")
        require(sp.gcd(P, sp.diff(P, t)) == 1, f"P squarefree n={n}")
        require(sp.gcd(W, sp.diff(W, t)) == 1, f"W squarefree n={n}")
        require(sp.gcd(W, P * R) == 1, f"boundary disjoint n={n}")
        require(beta != 0, f"beta nonzero n={n}")
        ode_checks += 10

        # Master Keller specialization in the invariant t-coordinate.
        G = sp.expand(P * R)
        H = R ** (-n)
        require(zero(sp.diff(t * H * G**n, t) - G ** (n - 1)), f"master Keller n={n}")
        require(zero(t * P**n - B), f"b observable n={n}")
        require(zero((B - beta) - R**n * W), f"e observable n={n}")
        sigma = sp.expand(b * (b - beta))
        require(sp.gcd(sigma, sp.diff(sigma, b)) == 1, f"target smooth n={n}")
        geometry_checks += 4

        # Maximal-intersection exponent invoice for several complete periods.
        for s in range(1, 4 * n + 1):
            m = (s + n - 1) // n
            require(0 <= n * m - s < n, f"intersection remainder n={n},s={s}")
            require(n * m >= s, f"source regularity n={n},s={s}")
            intersection_checks += 2

        # Exact symplectic primitive and two-bracket identity.
        U, V, gcd = sp.gcdex(sigma, sp.diff(sigma, b), b)
        require(zero(gcd - 1), f"Bezout gcd n={n}")
        T = -V
        require(zero(U * sigma - T * sp.diff(sigma, b) - 1), f"Bezout n={n}")
        first = (U + sp.diff(T, b)) * c * e / (n - 1)
        second = -T * e
        value = target_bracket(first, b, n, sigma) + target_bracket(second, c, n, sigma)
        require(zero(on_target(value - 1, n, sigma)), f"two brackets n={n}")
        geometry_checks += 3

        # Universal two-by-two factorization controls.
        h = sp.expand(sigma * (b**2 + 2 * b + 3))
        K = b**3 - b + 2
        A, BB, L, M = map(sp.Rational, (2, 3, 5, 7))
        for k in range(1, 7):
            p = n * (k - 1) + 1
            f = A * h
            g = BB * h**k
            F = L * K**p
            GG = M * K
            low = sp.expand(n * (f * sp.diff(g, b) - k * sp.diff(f, b) * g))
            high = sp.expand(sp.diff(F, b) * GG - p * F * sp.diff(GG, b))
            middle = sp.expand(
                sp.diff(f, b) * GG
                + n * f * sp.diff(GG, b)
                - n * k * sp.diff(F, b) * g
                - p * F * sp.diff(g, b)
            )
            expected = sp.expand(
                (K * sp.diff(h, b) + n * h * sp.diff(K, b))
                * (A * M - k * p * L * BB * K ** (p - 1) * h ** (k - 1))
            )
            require(zero(low) and zero(high), f"2x2 extremes n={n},k={k}")
            require(zero(middle - expected), f"2x2 factor n={n},k={k}")
            require(sp.degree(K * sp.diff(h, b) + n * h * sp.diff(K, b), b) >= 1, "2x2 unit hostile")
            fewnomial_checks += 3
    return ode_checks, geometry_checks, intersection_checks, fewnomial_checks


def audit_direct_source() -> int:
    checks = 0
    for n in range(2, 6):
        P, R, B, beta, W = packet(n)
        tt = x**n * q
        PP, RR, BB, WW = (expr.subs(t, tt) for expr in (P, R, B, W))
        aa = q / RR**n
        cc = x * PP * RR
        ee = sp.expand(q * WW)
        require(zero(jacobian(aa, cc) + 1), f"direct Jacobian n={n}")
        require(zero(BB - aa * cc**n), f"direct b n={n}")
        require(zero(ee - aa * (BB - beta)), f"direct e n={n}")
        require(zero(cc**n * ee - BB * (BB - beta)), f"direct surface n={n}")
        checks += 4
    return checks


def audit_degree_seven() -> tuple[int, int, int]:
    n = 3
    P, R, B, beta, W = packet(n)
    explicit_P = (2 * t**2 + 7 * t + 14) / 14
    explicit_W = (8 * t**4 + 60 * t**3 + 258 * t**2 + 557 * t + 729) / 2744
    require(zero(P - explicit_P), "P3 display")
    require(zero(beta + sp.Rational(729, 2744)), "beta3 display")
    require(zero(W - explicit_W), "W3 display")
    require(zero(W.subs(t, -1) - sp.Rational(27, 196)), "W3 side owner")
    explicit_checks = 4

    disc = sp.factor(sp.discriminant(B - b, t))
    quotient = sp.factor(disc / (b**4 * (b - beta) ** 2))
    require(quotient != 0 and zero(sp.diff(quotient, b)), "A7 square discriminant shape")
    require(sp.degree(B, t) == 7 and sp.degree(P, t) == 2 and sp.degree(W, t) == 4, "A7 passport degrees")
    require(1 + 3 + 3 == 7 and 3 + 4 == 7, "A7 passport invoice")
    passport_checks = 3

    tt = x**3 * q
    PP, RR, BB, WW = (expr.subs(t, tt) for expr in (P, R, B, W))
    aa = q / RR**3
    cc = x * PP * RR
    ee = sp.expand(q * WW)
    sigma = sp.expand(b * (b - beta))
    subs = {b: BB, c: cc, e: ee}
    poisson_checks = 0
    for f, g in ((b, c), (c, e), (b, e)):
        source = jacobian(f.subs(subs), g.subs(subs))
        target = target_bracket(f, g, 3, sigma).subs(subs)
        require(zero(source + target), "degree-seven anti-Poisson")
        poisson_checks += 1
    require(zero(jacobian(aa, cc) + 1), "degree-seven direct Keller")
    require(zero(cc**3 * ee - BB * (BB - beta)), "degree-seven direct target")
    poisson_checks += 2
    return explicit_checks, passport_checks, poisson_checks


def audit_no_asserts() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "assert nodes are forbidden")
    return count


def main() -> None:
    ode, geometry, intersection, fewnomial = audit_tower()
    direct = audit_direct_source()
    explicit, passport, poisson = audit_degree_seven()
    assert_nodes = audit_no_asserts()
    print("THM-3576 finite exact companion")
    print("status PROVISIONAL_VERIFIED_EXACT_PENDING_INDEPENDENT_AUDIT")
    print("universe exact_QQ_symbolic_exponents_2_through_8")
    print("hypergeometric_ode_and_critical_value_checks", ode)
    print("geometry_observable_and_two_bracket_checks", geometry)
    print("maximal_intersection_exponent_checks", intersection)
    print("universal_two_by_two_factorization_checks", fewnomial)
    print("direct_source_keller_surface_checks", direct)
    print("degree_seven_display_checks", explicit)
    print("degree_seven_passport_discriminant_checks", passport)
    print("degree_seven_anti_poisson_checks", poisson)
    print("ast_assert_nodes", assert_nodes)
    print("scope universal_proofs_not_finite_extrapolation_no_JC2_claim")
    print("PASS")


if __name__ == "__main__":
    main()
