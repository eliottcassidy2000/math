"""Exact controls for the deformed Hamiltonian torsion ladder.

This is diagnostic support for
the matching canonical theorem and reflections.  The proof there is
uniform; these checks exercise exact polynomial identities and hostile
nonlinear jets over QQ.
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path
import sys

import sympy as sp


x, z, T = sp.symbols("x z T")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_polynomial(expr: sp.Expr) -> bool:
    numerator, denominator = sp.fraction(sp.cancel(expr))
    if denominator != 1:
        return False
    try:
        sp.Poly(numerator, x, z)
    except sp.PolynomialError:
        return False
    return True


def truncated_inverse(poly: sp.Expr, order: int) -> sp.Expr:
    """Return the representative of poly^(-1) modulo x^order."""

    return sp.expand(sp.series(1 / poly, x, 0, order).removeO())


def truncated_bridge(phi: sp.Expr, order: int) -> sp.Expr:
    """Find F(T), deg F < order, from F(phi)=(order)*u/phi'."""

    u = sp.cancel(phi / x)
    coeffs = sp.symbols(f"c0:{order}")
    candidate = sum(coeffs[j] * T**j for j in range(order))
    residue = sp.series(
        candidate.subs(T, phi) - order * u / sp.diff(phi, x),
        x,
        0,
        order,
    ).removeO().expand()
    equations = [sp.Eq(residue.coeff(x, j), 0) for j in range(order)]
    solution = sp.solve(equations, coeffs, dict=True)
    if len(solution) != 1:
        raise RuntimeError(f"bridge interpolation was not unique: {solution}")
    return sp.expand(candidate.subs(solution[0]))


def exact_case(r: int, lam: sp.Rational, phi: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    if r < 2 or lam == 0 or phi.subs(x, 0) != 0 or sp.diff(phi, x).subs(x, 0) == 0:
        raise ValueError("the audited regime is r>=2, lambda!=0, phi(0)=0, phi'(0)!=0")

    P = sp.expand(phi + lam * x**r * z)
    Px = sp.diff(P, x)
    Pz = sp.diff(P, z)
    A = truncated_inverse(Px, r)
    B = sp.cancel((1 - A * Px) / Pz)
    m = sp.expand(sp.diff(A, x) + sp.diff(B, z))
    h = -A / Pz
    Q0 = x ** (1 - r) / (lam * (r - 1))
    U = sp.cancel(P / x)
    F = truncated_bridge(phi, r - 1)
    g = sp.cancel(P * h + F.subs(T, P) * Q0)

    D = lambda q: sp.expand(Px * sp.diff(q, z) - Pz * sp.diff(q, x))

    checks = {
        "bezout": sp.cancel(A * Px + B * Pz - 1),
        "mu_primitive": sp.cancel(D(h) - m),
        "unit_primitive": sp.cancel(D(Q0) - 1),
        "theta_kill": sp.cancel(P ** (r - 1) * Q0),
        "mu_kill": sp.cancel(P**r * h),
        "g_polynomial_denominator": sp.denom(g),
        "bridge": sp.cancel(D(g) - P * m - F.subs(T, P)),
    }
    require(checks["bezout"] == 0, f"Bezout identity failed at r={r}")
    require(checks["mu_primitive"] == 0, f"mu primitive failed at r={r}")
    require(checks["unit_primitive"] == 0, f"unit primitive failed at r={r}")
    require(is_polynomial(checks["theta_kill"]), f"theta kill failed at r={r}")
    require(is_polynomial(checks["mu_kill"]), f"mu kill failed at r={r}")
    require(
        not is_polynomial(sp.cancel(P ** (r - 2) * Q0)),
        f"theta lower power lost its pole at r={r}",
    )
    require(
        not is_polynomial(sp.cancel(P ** (r - 1) * h)),
        f"mu lower power lost its pole at r={r}",
    )
    require(checks["g_polynomial_denominator"] == 1, f"bridge pole at r={r}")
    require(checks["bridge"] == 0, f"bridge identity failed at r={r}")
    require(F.subs(T, 0) == r - 1, f"bridge constant failed at r={r}")

    return sp.factor(F), sp.factor(m)


def main() -> None:
    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    tree = ast.parse(source.decode("utf-8"))
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assertion_nodes == 0, "source contains assertion nodes")
    require(float_literals == 0, "source contains floating literals")

    print("Linear-in-z jet torsion and rational-primitive trichotomy audit")
    print("scope=GENERIC-SYMBOLIC_plus_FINITE-EXACT_HOSTILES;not_JC2_closure")
    print(f"python={sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
    print(f"sympy={sp.__version__}")
    print(f"assertion_nodes={assertion_nodes};float_literals={float_literals}")

    a = sp.symbols("a")
    symbolic = [
        (2, sp.Rational(1), x + a * x**2),
        (3, sp.Rational(1), x + a * x**2),
        (4, sp.Rational(1), x + a * x**2),
        (5, sp.Rational(1), x + a * x**3),
    ]
    for r, lam, phi in symbolic:
        F, m = exact_case(r, lam, phi)
        print(f"symbolic r={r} phi={phi}: F={F}; mu_rep={m}")

    hostile = [
        (2, sp.Rational(-3), 2 * x - 3 * x**2 + 5 * x**4),
        (3, sp.Rational(2), -x + 4 * x**2 - 2 * x**3),
        (4, sp.Rational(-2), 3 * x + x**2 - 5 * x**3 + 2 * x**6),
        (5, sp.Rational(3), 2 * x - x**2 + 3 * x**3 - 4 * x**4),
        (6, sp.Rational(-1), -2 * x + 5 * x**2 + 7 * x**5 + x**8),
        (7, sp.Rational(4), x - 2 * x**3 + 3 * x**4 - x**6),
    ]
    for r, lam, phi in hostile:
        F, _ = exact_case(r, lam, phi)
        print(f"hostile r={r} lambda={lam} phi={phi}: F={F}")

    for r, s, lam, aa in [
        (2, 2, sp.Rational(1), sp.Rational(1)),
        (2, 3, sp.Rational(-2), sp.Rational(3)),
        (4, 2, sp.Rational(5), sp.Rational(-1)),
    ]:
        g = sp.expand(lam * x**r * (1 + aa * x) ** s)
        P = x + g * z
        A = 1 - sp.diff(g, x) * z
        B = sp.cancel(sp.diff(g, x) ** 2 * z**2 / g)
        require(is_polynomial(B), f"multi-root B is not polynomial: {(r, s, lam, aa)}")
        require(
            sp.expand(A * sp.diff(P, x) + B * sp.diff(P, z) - 1) == 0,
            f"multi-root Bezout identity failed: {(r, s, lam, aa)}",
        )
        m = sp.factor(sp.diff(A, x) + sp.diff(B, z))
        formula = sp.factor(z * (2 * sp.diff(g, x) ** 2 / g - sp.diff(g, x, 2)))
        require(sp.simplify(m - formula) == 0, "multi-root divergence formula failed")
        residue_zero = sp.simplify(sp.residue(1 / g, x, 0))
        residue_second = sp.simplify(sp.residue(1 / g, x, -1 / aa))
        require(residue_zero != 0 and residue_second != 0, "residue hostile vanished")
        print(
            f"multi-root r={r} s={s} lambda={lam} a={aa}: mu_rep={m};"
            f"residues=({residue_zero},{residue_second})"
        )

    print(f"source_lf_sha256={hashlib.sha256(source).hexdigest()}")
    print("PASS: 4 symbolic families, 6 exact hostile jets, and 3 multi-root rows")


if __name__ == "__main__":
    main()
