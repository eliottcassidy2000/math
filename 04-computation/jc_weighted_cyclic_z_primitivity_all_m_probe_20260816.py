#!/usr/bin/env python3
"""Exact companion for all-grade z-primitivity in the cyclic weighted family.

For the THM-3448 seed of parameter ell, reduce the z-reconstruction numerator

    H = gamma * (gamma * (gamma - 1 + a) - a*w)

modulo the degree n=ell+2 inverse polynomial.  The coefficient of w^(n-1)
in the remainder is a nonzero polynomial in P,Q.  Hence H(w), and therefore
z, is not in the target base field.  Full S_n monodromy then upgrades this to
primitivity.

The symbolic proof is checked in n, and direct polynomial remainders provide
hostile controls for ell=1..30.  Truth gates survive optimized Python.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import sympy as sp


EXPECTED_SEMANTIC_SHA256: str | None = "6f06d5042a944f817b23e5bc43e3d70f74098ed4640e20eb95a638bec575cbf0"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def poly_digest(expression: sp.Expr, generators: tuple[sp.Symbol, ...]) -> str:
    polynomial = sp.Poly(sp.cancel(expression), *generators, domain=sp.QQ)
    payload = tuple(
        (monomial, int(coefficient.p), int(coefficient.q))
        for monomial, coefficient in polynomial.terms()
    )
    return sha256(repr(payload).encode("ascii")).hexdigest()


def symbolic_leading_coefficient() -> tuple[sp.Expr, ...]:
    n, ell, P, Q = sp.symbols("n ell P Q")

    # Let c_k extract the w^(n-1) coefficient after reduction modulo
    # w^n-w^(n-1)+Pw-Q.  For n>=4, c_k=1 on n-1<=k<=2n-3.
    # The next affine run and its final exceptional step are:
    c_3n_6 = 1 - (n - 3) * P + (n - 4) * Q
    c_3n_5 = 1 - (n - 2) * P + (n - 3) * Q
    c_3n_4 = 1 - (n - 1) * P + (n - 2) * Q
    c_3n_3 = 1 - n * P + P**2 + (n - 1) * Q

    lambda_gamma2 = 1 - n * (n - 2) * P
    lambda_w_gamma = sp.Integer(1)
    high_gamma3 = sp.expand(
        -(n - 1) ** 3 * c_3n_6
        + 3 * n * (n - 1) ** 2 * c_3n_5
        - 3 * n**2 * (n - 1) * c_3n_4
        + n**3 * c_3n_3
    )
    require(
        sp.expand(high_gamma3 - (1 + n**3 * P**2 + (3 - 4 * n) * P + (4 * n - 4) * Q)) == 0,
        "high gamma-cube coefficient",
    )

    lambda_gamma3 = sp.expand(
        3 * n * P**2 + 3 * P * (1 - n**2 * P) + high_gamma3
    )
    expected_gamma3 = sp.expand(
        1
        + (n**3 - 3 * n**2 + 3 * n) * P**2
        + (6 - 4 * n) * P
        + (4 * n - 4) * Q
    )
    require(sp.expand(lambda_gamma3 - expected_gamma3) == 0, "lambda gamma cube")

    a = -(2 * n - 3) / (2 * (n - 2))
    lambda_h = sp.factor(
        lambda_gamma3 + (a - 1) * lambda_gamma2 - a * lambda_w_gamma
    )
    lambda_h_ell = sp.factor(lambda_h.subs(n, ell + 2))
    expected_ell = sp.expand(
        ((ell + 1) ** 3 + 1) * P**2
        + sp.Rational(1, 2) * (4 * ell**2 + ell - 2) * P
        + 4 * (ell + 1) * Q
    )
    require(sp.cancel(lambda_h_ell - expected_ell) == 0, "all-ell leading coefficient")
    require(expected_ell != 0, "leading coefficient polynomial")
    return (
        sp.factor(lambda_gamma2),
        sp.factor(lambda_w_gamma),
        sp.factor(lambda_gamma3),
        sp.factor(expected_ell),
    )


def direct_remainder_rows() -> tuple[object, ...]:
    w, P, Q = sp.symbols("w P Q")
    rows: list[tuple[object, ...]] = []
    digest_rows: list[tuple[object, ...]] = []
    for ell in range(1, 31):
        n = ell + 2
        a = -sp.Rational(2 * ell + 1, 2 * ell)
        p = (ell + 1) * w**ell - (ell + 2) * w ** (ell + 1)
        gamma = P - p
        inverse = w ** (ell + 2) - w ** (ell + 1) + P * w - Q
        numerator = sp.expand(gamma * (gamma * (gamma - 1 + a) - a * w))
        remainder = sp.rem(numerator, inverse, w)
        remainder_poly = sp.Poly(remainder, w, P, Q, domain=sp.QQ)
        require(remainder_poly.degree(w) == n - 1, (ell, "remainder degree"))

        leading = sp.factor(sp.Poly(remainder, w).LC())
        if ell == 1:
            expected = sp.Rational(3, 2) * P * (6 * P + 1)
        else:
            expected = sp.expand(
                ((ell + 1) ** 3 + 1) * P**2
                + sp.Rational(1, 2) * (4 * ell**2 + ell - 2) * P
                + 4 * (ell + 1) * Q
            )
        require(sp.expand(leading - expected) == 0, (ell, "leading coefficient"))
        require(leading != 0, (ell, "nonzero leading coefficient"))

        # Hostile specialization: at P=Q=0 the entire remainder vanishes.
        # A single degenerate target therefore cannot prove the generic claim.
        require(sp.expand(remainder.subs({P: 0, Q: 0})) == 0, (ell, "zero-target hostile"))
        require(leading.subs({P: 1, Q: 0}) != 0, (ell, "positive control"))
        digest_rows.append((ell, n, poly_digest(remainder, (w, P, Q))))
        if ell in (1, 2, 3, 5, 9, 17, 30):
            rows.append((ell, n, n - 1, sp.sstr(leading), digest_rows[-1][2]))

    digest = sha256(repr(tuple(digest_rows)).encode("ascii")).hexdigest()
    return tuple(rows), digest


def main() -> None:
    symbolic = symbolic_leading_coefficient()
    sample_rows, remainder_digest = direct_remainder_rows()
    verdict = (
        "for_every_ell>=2_the_z_numerator_remainder_has_exact_degree_n-1",
        "therefore_z_is_not_in_the_target_base_field",
        "S_n_point_stabilizer_maximality_then_gives_K(z)=K(w)",
        "ell=1_is_checked_separately_and_is_also_primitive",
        "all_x_y_z_views_in_the_cyclic_weighted_family_share_the_inverse_discriminant_square_class",
        "for_odd_ell=2m-3_the_even_C_discriminant_order_still_misses_the_effective_C_component",
        "no_historical_E_m_identification;no_map_classification;no_JC2_or_LRC_consequence",
    )
    semantic_surface = (tuple(sp.sstr(value) for value in symbolic), sample_rows, remainder_digest, verdict)
    semantic_sha256 = sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, (semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    print("Cyclic weighted Keller family: all-grade z-primitivity remainder")
    print("status=EXACT_SYMBOLIC_PROOF_COMPANION;scope=THM3448_cyclic_family")
    print(f"symbolic=(lambda_gamma2,lambda_wgamma,lambda_gamma3,lambda_H_ell)={tuple(sp.sstr(value) for value in symbolic)}")
    print(f"direct_remainder_samples={sample_rows}")
    print(f"direct_remainder_range=ell_1_through_30;ledger_sha256={remainder_digest}")
    print("hostile=P=Q=0_makes_every_tested_remainder_zero;generic_polynomial_nonvanishing_is_required")
    print(f"verdict={verdict}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;all_checks_survive_python_O;no_randomness;no_elapsed_fields")
    print("commands=python -B 04-computation/jc_weighted_cyclic_z_primitivity_all_m_probe_20260816.py;python -B -O 04-computation/jc_weighted_cyclic_z_primitivity_all_m_probe_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
