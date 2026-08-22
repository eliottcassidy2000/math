#!/usr/bin/env python3
"""Exact audit for the canonical divergence class of P=x+lambda*x^r*z.

The uniform proof is written in the matching reflection.  This companion
checks its symbolic identities, sign conventions, exact torsion witnesses,
Kummer equations, and a bounded hostile linear-algebra bank.  It performs no
floating-point computation and uses explicit runtime gates rather than Python
assertions, so ``python`` and ``python -O`` execute the same checks.
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path
import sys

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def exact_zero(expr: sp.Expr) -> bool:
    return sp.simplify(sp.cancel(sp.together(expr))) == 0


x, z = sp.symbols("x z")


def hamiltonian(P: sp.Expr, f: sp.Expr) -> sp.Expr:
    """D_P(f)=P_x f_z-P_z f_x, the repository's Jacobian convention."""

    return sp.expand(sp.diff(P, x) * sp.diff(f, z) - sp.diff(P, z) * sp.diff(f, x))


def x_boundary_valuation(expr: sp.Expr) -> int:
    """Order along x=0 in QQ(x,z), after exact cancellation."""

    numerator, denominator = sp.fraction(sp.cancel(expr))
    numerator_terms = sp.Poly(numerator, x, z, domain=sp.QQ).terms()
    denominator_terms = sp.Poly(denominator, x, z, domain=sp.QQ).terms()
    require(bool(numerator_terms), "zero numerator has no boundary valuation")
    require(bool(denominator_terms), "zero denominator")
    numerator_order = min(monomial[0] for monomial, _ in numerator_terms)
    denominator_order = min(monomial[0] for monomial, _ in denominator_terms)
    return numerator_order - denominator_order


def is_polynomial(expr: sp.Expr) -> bool:
    numerator, denominator = sp.fraction(sp.cancel(expr))
    if denominator != 1:
        return False
    try:
        sp.Poly(numerator, x, z, domain=sp.QQ)
    except sp.PolynomialError:
        return False
    return True


def bounded_rank_audit(r: int) -> tuple[tuple[int, int], int, int, int]:
    """Hostile finite bank: divergence is outside a rectangular D_P image.

    This is deliberately only a finite control.  The all-degree exclusion is
    the localization-and-pole proof in the reflection.
    """

    P = x + x**r * z
    x_cap = 2 * r + 2
    z_cap = 4
    domain = [x**i * z**j for j in range(z_cap + 1) for i in range(x_cap + 1)]
    images = [hamiltonian(P, monomial) for monomial in domain]
    negative_target = r * (r + 1) * x ** (r - 2) * z
    witness = x**2 * z**2 + 3 * x * z + x
    positive_target = hamiltonian(P, witness)

    support: set[tuple[int, int]] = set()
    for polynomial in images + [negative_target, positive_target]:
        support.update(sp.Poly(polynomial, x, z, domain=sp.QQ).monoms())
    ordered_support = sorted(support)

    def vector(polynomial: sp.Expr) -> sp.Matrix:
        terms = dict(sp.Poly(polynomial, x, z, domain=sp.QQ).terms())
        return sp.Matrix([terms.get(monomial, 0) for monomial in ordered_support])

    matrix = sp.Matrix.hstack(*(vector(image) for image in images))
    negative_vector = vector(negative_target)
    positive_vector = vector(positive_target)
    rank = matrix.rank()
    negative_rank = matrix.row_join(negative_vector).rank()
    positive_rank = matrix.row_join(positive_vector).rank()
    require(negative_rank == rank + 1, f"bounded hostile unexpectedly solved for r={r}")
    require(positive_rank == rank, f"positive image control failed for r={r}")
    return matrix.shape, rank, negative_rank, positive_rank


def generic_symbolic_audit() -> None:
    r = sp.symbols("r", integer=True, positive=True)
    lam = sp.symbols("lambda", nonzero=True)
    P = x + lam * x**r * z
    Px = sp.diff(P, x)
    Pz = sp.diff(P, z)
    A = 1 - r * lam * x ** (r - 1) * z
    B = r**2 * lam * x ** (r - 2) * z**2
    mu_representative = lam * r * (r + 1) * x ** (r - 2) * z
    correction = r * z / x - x ** (-r) / lam
    local_mate = x ** (1 - r) / (lam * (r - 1))
    bridge_witness = (r - 1) * z + r * lam * x ** (r - 1) * z**2

    require(exact_zero(A * Px + B * Pz - 1), "generic Bezout identity")
    require(
        exact_zero(sp.diff(A, x) + sp.diff(B, z) - mu_representative),
        "generic divergence formula",
    )
    require(exact_zero(hamiltonian(P, P)), "D_P(P)=0")
    require(
        exact_zero(hamiltonian(P, correction) - mu_representative),
        "localized divergence correction",
    )
    require(exact_zero(A + correction * Pz), "adjusted first row is Q_z")
    require(
        exact_zero(B - correction * Px + sp.diff(local_mate, x)),
        "adjusted second row is -Q_x",
    )
    require(exact_zero(hamiltonian(P, local_mate) - 1), "local mate has response +1")
    opposite = sp.diff(local_mate, x) * Pz - sp.diff(local_mate, z) * Px
    require(exact_zero(opposite + 1), "opposite bracket convention has response -1")
    require(
        exact_zero(P * correction + (r - 1) * local_mate - bridge_witness),
        "polynomial bridge primitive",
    )
    require(
        exact_zero(
            hamiltonian(P, bridge_witness) - (P * mu_representative + (r - 1))
        ),
        "P*mu=-(r-1)*theta bridge",
    )
    require(
        exact_zero(
            hamiltonian(P, P**r * correction) - P**r * mu_representative
        ),
        "P^r kills the divergence class",
    )
    require(
        exact_zero(
            hamiltonian(P, P ** (r - 1) * local_mate) - P ** (r - 1)
        ),
        "P^(r-1) kills the unit-response class",
    )


def finite_row(r: int, lam_value: int) -> tuple[str, ...]:
    lam = sp.Integer(lam_value)
    P = x + lam * x**r * z
    Px = sp.diff(P, x)
    Pz = sp.diff(P, z)
    A = 1 - r * lam * x ** (r - 1) * z
    B = r**2 * lam * x ** (r - 2) * z**2
    mu_representative = lam * r * (r + 1) * x ** (r - 2) * z
    correction = r * z / x - x ** (-r) / lam
    local_mate = x ** (1 - r) / (lam * (r - 1))
    bridge_witness = (r - 1) * z + r * lam * x ** (r - 1) * z**2

    checks = (
        exact_zero(A * Px + B * Pz - 1),
        exact_zero(sp.diff(A, x) + sp.diff(B, z) - mu_representative),
        exact_zero(hamiltonian(P, correction) - mu_representative),
        exact_zero(A + correction * Pz),
        exact_zero(B - correction * Px + sp.diff(local_mate, x)),
        exact_zero(hamiltonian(P, local_mate) - 1),
        exact_zero(P * correction + (r - 1) * local_mate - bridge_witness),
        exact_zero(
            hamiltonian(P, bridge_witness) - (P * mu_representative + (r - 1))
        ),
    )
    require(all(checks), f"finite identity row failed for r={r}, lambda={lam}")

    theta_lower = P ** (r - 2) * local_mate
    theta_kill = P ** (r - 1) * local_mate
    mu_lower = P ** (r - 1) * correction
    mu_kill = P**r * correction
    theta_valuations = (
        x_boundary_valuation(theta_lower),
        x_boundary_valuation(theta_kill),
    )
    mu_valuations = (
        x_boundary_valuation(mu_lower),
        x_boundary_valuation(mu_kill),
    )
    require(theta_valuations == (-1, 0), f"theta pole ladder failed for r={r}")
    require(mu_valuations == (-1, 0), f"mu pole ladder failed for r={r}")
    require(not is_polynomial(theta_lower), f"theta lower power became polynomial for r={r}")
    require(is_polynomial(theta_kill), f"theta killing power not polynomial for r={r}")
    require(not is_polynomial(mu_lower), f"mu lower power became polynomial for r={r}")
    require(is_polynomial(mu_kill), f"mu killing power not polynomial for r={r}")

    w, q = sp.symbols("w q")
    deck = w ** (r - 1) - lam * (r - 1) * q
    factorization = sp.factor_list(
        sp.Poly(deck, w, domain=sp.QQ.frac_field(q))
    )[1]
    factor_degrees = tuple(factor.degree() for factor, _ in factorization)
    require(factor_degrees == (r - 1,), f"Kummer polynomial factored for r={r}")
    deck_on_chart = deck.subs({w: 1 / x, q: local_mate})
    require(exact_zero(deck_on_chart), f"Kummer chart equation failed for r={r}")

    return (
        f"r={r}",
        f"lambda={lam_value}",
        "identities=PASS",
        f"theta_boundary_orders={theta_valuations}",
        f"mu_boundary_orders={mu_valuations}",
        f"deck_degree={r - 1}",
        f"deck_factor_degrees={factor_degrees}",
    )


def main() -> None:
    source_path = Path(__file__)
    source = source_path.read_bytes().replace(b"\r\n", b"\n")
    tree = ast.parse(source.decode("utf-8"))
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assertion_nodes == 0, "source contains assertion nodes")
    require(float_literals == 0, "source contains floating literals")

    generic_symbolic_audit()
    rows: list[str] = []
    for r in range(2, 10):
        for lam in (1, -1, 2, -3):
            rows.append(";".join(finite_row(r, lam)))

    rank_rows: list[str] = []
    for r in range(2, 9):
        shape, rank, negative_rank, positive_rank = bounded_rank_audit(r)
        rank_rows.append(
            f"r={r};matrix={shape};rank={rank};"
            f"rank_with_mu={negative_rank};rank_with_image_control={positive_rank}"
        )

    semantic_material = "\n".join(rows + rank_rows).encode("utf-8")
    semantic_sha256 = hashlib.sha256(semantic_material).hexdigest()
    script_sha256 = hashlib.sha256(source).hexdigest()

    print("jacobian_divergence_family_agent_exact")
    print("status=GENERIC-SYMBOLIC + FINITE-EXACT HOSTILE BANK")
    print(f"python={sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
    print(f"sympy={sp.__version__}")
    print(f"assertion_nodes={assertion_nodes}")
    print(f"float_literals={float_literals}")
    print("convention=D_P(f)=P_x*f_z-P_z*f_x=Jac(P,f)")
    print("generic_bezout=A*P_x+B*P_z=1")
    print("generic_mu=lambda*r*(r+1)*x^(r-2)*z")
    print("generic_local_correction=r*z/x-lambda^(-1)*x^(-r)")
    print("generic_local_mate=x^(1-r)/(lambda*(r-1))")
    print("annihilator_theta=(P^(r-1));annihilator_mu=(P^r)")
    print("torsion_bridge=P*mu=-(r-1)*theta")
    print("generic_symbolic_checks=PASS")
    print("finite_identity_bank_begin")
    for row in rows:
        print(row)
    print("finite_identity_bank_end")
    print("bounded_rank_bank_begin")
    for row in rank_rows:
        print(row)
    print("bounded_rank_bank_end")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={script_sha256}")
    print("PASS")


if __name__ == "__main__":
    main()
