#!/usr/bin/env python3
"""Exact referee for THM-2063 (one-fiber-linear planar Keller pairs)."""

from __future__ import annotations

import random

import sympy as sp


x, y, u, v, t = sp.symbols("x y u v t")


def require(condition: bool) -> None:
    if not condition:
        raise AssertionError("exact referee check failed")


def jacobian(p: sp.Expr, q: sp.Expr, xx: sp.Symbol, yy: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.diff(p, xx) * sp.diff(q, yy) - sp.diff(p, yy) * sp.diff(q, xx))


def random_poly(rng: random.Random, variable: sp.Symbol, degree: int) -> sp.Expr:
    coeffs = [rng.randint(-4, 4) for _ in range(degree + 1)]
    if degree and coeffs[-1] == 0:
        coeffs[-1] = rng.choice((-4, -3, -2, -1, 1, 2, 3, 4))
    return sp.expand(sum(c * variable**j for j, c in enumerate(coeffs)))


def check_leading_coefficient() -> int:
    checks = 0
    a0, a1, a2, a3, b0, b1, b2 = sp.symbols("a0:4 b0:3")
    A = a0 + a1 * x + a2 * x**2 + a3 * x**3
    B = b0 + b1 * x + b2 * x**2
    for n in range(1, 7):
        lower = sum((j + 2 * x + x**2) * y**j for j in range(n))
        qn = (n + 1) - 3 * x + 2 * x**2 + x**3
        Q = sp.expand(lower + qn * y**n)
        coefficient = sp.Poly(jacobian(A * y + B, Q, x, y), y).coeff_monomial(y**n)
        expected = sp.expand(n * sp.diff(A, x) * qn - A * sp.diff(qn, x))
        require(sp.expand(coefficient - expected) == 0)
        checks += 1
    return checks


def check_nonzero_a_forms(rng: random.Random, trials: int, inverse_trials: int) -> tuple[int, int]:
    for _ in range(trials):
        a = sp.Integer(rng.choice([i for i in range(-5, 6) if i]))
        kappa = sp.Integer(rng.choice([i for i in range(-7, 8) if i]))
        b = sp.Integer(rng.randint(-5, 5))
        B = random_poly(rng, x, rng.randint(0, 5))
        Ht = random_poly(rng, t, rng.randint(0, 5))
        P = sp.expand(a * y + B)
        H_of_P = sp.expand(Ht.subs(t, P))
        Q = sp.expand(H_of_P - kappa * x / a + b)
        require(jacobian(P, Q, x, y) == kappa)

    # Full two-sided composition is substantially larger than the determinant
    # identity. Keep a separate exact low-degree sample instead of turning the
    # referee into an expansion benchmark.
    for _ in range(inverse_trials):
        a = sp.Integer(rng.choice([i for i in range(-5, 6) if i]))
        kappa = sp.Integer(rng.choice([i for i in range(-7, 8) if i]))
        b = sp.Integer(rng.randint(-5, 5))
        B = random_poly(rng, x, rng.randint(0, 2))
        Ht = random_poly(rng, t, rng.randint(0, 2))
        P = sp.expand(a * y + B)
        Q = sp.expand(Ht.subs(t, P) - kappa * x / a + b)
        H_of_u = sp.expand(Ht.subs(t, u))
        X = sp.expand(a * (H_of_u + b - v) / kappa)
        Y = sp.expand((u - B.subs(x, X)) / a)
        require(sp.expand(P.subs({x: X, y: Y}, simultaneous=True) - u) == 0)
        require(sp.expand(Q.subs({x: X, y: Y}, simultaneous=True) - v) == 0)
        require(sp.expand(X.subs({u: P, v: Q}, simultaneous=True) - x) == 0)
        require(sp.expand(Y.subs({u: P, v: Q}, simultaneous=True) - y) == 0)
    return trials, inverse_trials


def check_zero_a_forms(rng: random.Random, trials: int) -> int:
    for _ in range(trials):
        alpha = sp.Integer(rng.choice([i for i in range(-5, 6) if i]))
        kappa = sp.Integer(rng.choice([i for i in range(-7, 8) if i]))
        beta = sp.Integer(rng.randint(-5, 5))
        Ht = random_poly(rng, t, rng.randint(0, 6))
        P = alpha * x + beta
        Q = sp.expand(kappa * y / alpha + Ht.subs(t, P))
        require(jacobian(P, Q, x, y) == kappa)
        X = (u - beta) / alpha
        Y = sp.expand(alpha * (v - Ht.subs(t, u)) / kappa)
        require(sp.expand(P.subs({x: X, y: Y}, simultaneous=True) - u) == 0)
        require(sp.expand(Q.subs({x: X, y: Y}, simultaneous=True) - v) == 0)
    return trials


def check_hostile_top_layer() -> None:
    A = x**2 + 1
    n = 5
    qn = 3 * A**n
    require(sp.expand(n * sp.diff(A, x) * qn - A * sp.diff(qn, x)) == 0)
    # Passing the top layer is not enough: the final equation -A*R'=kappa
    # cannot hold with nonzero constant kappa when A is nonconstant.
    r0, r1, r2, kappa = sp.symbols("r0 r1 r2 kappa")
    remainder = sp.Poly(sp.expand(-A * sp.diff(r0 + r1 * x + r2 * x**2, x) - kappa), x)
    equations = remainder.all_coeffs()
    solution = sp.solve(equations + [sp.Symbol("w") * kappa - 1], [r0, r1, r2, kappa, sp.Symbol("w")], dict=True)
    require(solution == [])


def main() -> None:
    rng = random.Random(2063)
    print("JC2 ONE-FIBER-LINEAR AUDIT -- exact coefficient descent")
    print(f"leading-coefficient identities: {check_leading_coefficient()} PASS")
    determinants, inverses = check_nonzero_a_forms(rng, 24, 8)
    print(f"nonzero-A normal-form determinants: {determinants} PASS")
    print(f"nonzero-A two-sided inverses: {inverses} PASS")
    print(f"zero-A normal forms and inverses: {check_zero_a_forms(rng, 16)} PASS")
    check_hostile_top_layer()
    print("hostile nonconstant-A top-layer trap: PASS")
    print("PASS")


if __name__ == "__main__":
    main()
