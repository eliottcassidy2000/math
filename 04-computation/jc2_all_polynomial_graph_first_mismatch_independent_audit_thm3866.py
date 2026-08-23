#!/usr/bin/env python3
"""Independent exact hostile checks for provisional THM-3866.

The canonical checker is intentionally neither imported nor executed here.
This script rebuilds the formal quotient by monic coefficient division,
checks completion-to-polynomial descent on nonconstant marked graphs, attacks
all degree boundaries, and freezes a repeated-total-branch companion.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")


A, C, Z = sp.symbols("A C Z")
ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "01-canon" / "theorems" / (
    "THM-3866-all-polynomial-graph-branches-force-projective-companion.md"
)
CANON_PROOF_BODY_SHA256 = "a60cbf21b92bf844f9cbdfe695f7d28157f3cb38d9e4bac834be5224516ae3ed"
MAX_N = 4
EXPECTED_SEMANTIC_SHA256 = "44c2dcd32765728b70180debeb973e63c6745145782ca955006460b1e815438b"
GATES = 0


def require(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def equal(left: sp.Expr, right: sp.Expr, label: object) -> None:
    require(sp.expand(left - right) == 0, label)


def nonzero(value: sp.Expr, label: object) -> None:
    require(sp.expand(value) != 0, label)


def coeff_a(value: sp.Expr, exponent: int) -> sp.Expr:
    return sp.expand(value).coeff(A, exponent)


def degree_c(value: sp.Expr) -> int:
    value = sp.expand(value)
    nonzero(value, "degree of zero")
    return int(sp.Poly(value, C).degree())


def leading_c(value: sp.Expr) -> sp.Expr:
    return sp.expand(sp.Poly(sp.expand(value), C).LC())


def beta(index: int) -> sp.Rational:
    return sp.factor(
        sp.binomial(sp.Rational(3, 2), index + 2)
        * sp.Rational(2, 3) ** (index + 2)
    )


def delta(profile: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * profile**2
        + 8 * A * C**3
        - 54 * A * C * profile
        + 9 * C**2
        - 54 * profile
    )


def homogenize_c(value: sp.Expr, degree: int) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(value), C)
    return sp.expand(
        sum(
            coefficient * C**exponent * Z ** (degree - exponent)
            for (exponent,), coefficient in polynomial.terms()
        )
    )


def formal_quotient(profile_s: sp.Expr, through: int) -> tuple[
    sp.Expr, sp.Expr, sp.Expr, tuple[sp.Expr, ...]
]:
    """Construct Q_* coefficients without using the canonical companion."""
    profile_s = sp.expand(profile_s)
    graph = sp.expand(6 * profile_s * (1 + A * profile_s))
    base = sp.expand(2 * profile_s**2 * (3 + 4 * A * profile_s))
    factor = sp.expand(C - graph)
    marked = sp.expand(1 + 2 * A * profile_s)

    # Direct marked-root controls.
    equal(1 + sp.Rational(2, 3) * A * graph, marked**2,
          ("marked square", profile_s))
    equal(1 + A * graph + A**2 * base, marked**3,
          ("marked cube", profile_s))
    equal(delta(base).subs(C, graph), 0, ("graph divides branch", profile_s))

    graph_zero = graph.subs(A, 0)
    quotients: list[sp.Expr] = []
    for index in range(through + 1):
        numerator = beta(index) * C ** (index + 2) - coeff_a(base, index)
        numerator += sum(
            coeff_a(graph, shift) * quotients[index - shift]
            for shift in range(1, index + 1)
        )
        quotient, remainder = sp.div(
            sp.Poly(sp.expand(numerator), C),
            sp.Poly(C - graph_zero, C),
        )
        equal(remainder.as_expr(), 0,
              ("monic formal remainder", profile_s, index))
        item = sp.expand(quotient.as_expr())
        require(degree_c(item) == index + 1,
                ("q degree", profile_s, index))
        equal(leading_c(item), beta(index),
              ("q leading coefficient", profile_s, index))
        quotients.append(item)

    # Coefficientwise multiplication and graph evaluation are separate paths.
    quotient_truncation = sp.expand(
        sum(A**index * quotients[index] for index in range(through + 1))
    )
    binomial_truncation = sp.expand(
        sum(beta(index) * A**index * C ** (index + 2)
            for index in range(through + 1))
    )
    product_error = sp.expand(
        factor * quotient_truncation - (binomial_truncation - base)
    )
    graph_error = sp.expand(binomial_truncation.subs(C, graph) - base)
    for exponent in range(through + 1):
        equal(coeff_a(product_error, exponent), 0,
              ("formal product coefficient", profile_s, exponent))
        equal(coeff_a(graph_error, exponent), 0,
              ("formal graph coefficient", profile_s, exponent))

    return factor, graph, base, tuple(quotients)


def completion_descent(
    profile_name: str,
    factor: sp.Expr,
    base: sp.Expr,
    quotients: tuple[sp.Expr, ...],
) -> tuple[tuple[object, ...], ...]:
    """Check F-divisibility first, then coefficientwise A^N descent."""
    rows = []
    for index in range(MAX_N + 1):
        truncation = sp.expand(
            sum(A**j * quotients[j] for j in range(index))
        )
        b_n = sp.expand(base + factor * truncation)
        quotient, remainder = sp.div(sp.Poly(delta(b_n), C), sp.Poly(factor, C))
        equal(remainder.as_expr(), 0,
              ("polynomial F divisibility", profile_name, index))
        ell = sp.expand(quotient.as_expr())
        for exponent in range(index):
            equal(coeff_a(ell, exponent), 0,
                  ("completion coefficient descends", profile_name, index, exponent))
        rational_r = sp.cancel(ell / A**index)
        equal(sp.denom(rational_r), 1,
              ("A-adic quotient polynomial", profile_name, index))
        r_n = sp.expand(rational_r)
        equal(r_n.subs(A, 0), 54 * quotients[index],
              ("exact first A coefficient", profile_name, index))
        if index == 0:
            require(degree_c(r_n) == 2, ("N0 generic degree", profile_name))
            require(degree_c(r_n.subs(A, 0)) == 1,
                    ("N0 special degree", profile_name))
            equal(leading_c(r_n), 8 * A, ("N0 leading", profile_name))
        else:
            require(degree_c(r_n) == 2 * index + 1,
                    ("RN generic degree", profile_name, index))
            require(degree_c(r_n.subs(A, 0)) == index + 1,
                    ("RN special degree", profile_name, index))
            equal(
                leading_c(r_n),
                -27 * beta(index - 1) ** 2 * A**index,
                ("RN leading", profile_name, index),
            )
        rows.append(
            (
                profile_name,
                index,
                degree_c(r_n),
                degree_c(r_n.subs(A, 0)),
                str(sp.factor(leading_c(r_n))),
            )
        )
    return tuple(rows)


def build_base(
    factor: sp.Expr, base: sp.Expr, quotients: tuple[sp.Expr, ...], index: int
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    truncation = sp.expand(sum(A**j * quotients[j] for j in range(index)))
    b_n = sp.expand(base + factor * truncation)
    r_n = sp.cancel(delta(b_n) / (A**index * factor))
    equal(sp.denom(r_n), 1, ("base residual denominator", index))
    return b_n, sp.expand(1 + A * C + A**2 * b_n), sp.expand(r_n)


def mismatch_audit(
    factor: sp.Expr, base: sp.Expr, quotients: tuple[sp.Expr, ...]
) -> tuple[tuple[object, ...], ...]:
    cases: tuple[tuple[str, int, sp.Expr], ...] = (
        ("N0_T0", 0, sp.Integer(0)),
        ("N0_d0", 0, 2 + A),
        ("N0_d1", 0, (1 + A) * C + 3),
        ("N1_lower", 1, 3 + A),
        ("N1_resonant_fibre_cancel", 1, -beta(0) * C + 2),
        ("N1_higher", 1, (1 + A) * C**3 + C),
        ("N2_T0", 2, sp.Integer(0)),
        ("N2_lower", 2, (A + 2) * C + 1),
        ("N2_resonant", 2, -beta(1) * C**2 + C + 1),
        ("N2_higher", 2, (A + 1) * C**3 + 2),
        ("N3_special_top_cancel", 3, quotients[3] + 1 + A * C**5),
        ("N4_resonant_fibre_cancel", 4, -beta(3) * C**4 + C + 2),
    )
    rows = []
    for label, index, transverse in cases:
        b_n, u_n, r_n = build_base(factor, base, quotients, index)
        t_zero = sp.expand(transverse.subs(A, 0))
        nonzero(quotients[index] - t_zero, ("first mismatch", label))
        residual = sp.expand(
            r_n
            - 54 * u_n * transverse
            - 27 * A ** (index + 2) * factor * transverse**2
        )
        b = sp.expand(b_n + A**index * factor * transverse)
        equal(delta(b), A**index * factor * residual,
              ("nonlinear response factorization", label))
        special = sp.expand(residual.subs(A, 0))
        equal(special, 54 * (quotients[index] - t_zero),
              ("special first mismatch", label))
        nonzero(special, ("special residual nonzero", label))
        generic_degree = degree_c(residual)
        special_degree = degree_c(special)
        require(generic_degree > special_degree,
                ("strict degree drop", label, generic_degree, special_degree))

        if transverse == 0:
            expected_degree = 2 if index == 0 else 2 * index + 1
            expected_leading = 8 * A if index == 0 else (
                -27 * beta(index - 1) ** 2 * A**index
            )
            degree_label = "T=0"
        else:
            d = degree_c(transverse)
            tau = leading_c(transverse)
            if index == 0 and d == 0:
                expected_degree = 2
                expected_leading = 8 * A
                degree_label = "N=0,d=0"
            elif index == 0 or d > index:
                expected_degree = 2 * d + 1
                expected_leading = -27 * A ** (index + 2) * tau**2
                degree_label = "d>N"
            elif d < index:
                expected_degree = 2 * index + 1
                expected_leading = -27 * beta(index - 1) ** 2 * A**index
                degree_label = "d<N"
            else:
                expected_degree = 2 * index + 1
                expected_leading = (
                    -27 * A**index * (beta(index - 1) + A * tau) ** 2
                )
                degree_label = "d=N"
            if degree_label == "d=N" and sp.expand(tau + beta(index - 1)) == 0:
                equal(expected_leading.subs(A, 1), 0,
                      ("resonant finite-fibre cancellation", label))
                nonzero(expected_leading, ("resonant generic noncancellation", label))
        require(generic_degree == expected_degree,
                ("generic degree law", label, generic_degree, expected_degree))
        equal(leading_c(residual), expected_leading,
              ("leading coefficient law", label))
        rows.append(
            (
                label,
                index,
                None if transverse == 0 else degree_c(transverse),
                degree_label,
                generic_degree,
                special_degree,
                str(sp.factor(leading_c(residual))),
            )
        )
    return tuple(rows)


def repeated_graph_control() -> tuple[object, ...]:
    """A total F^2 branch whose selected infinity component is not F."""
    profile_s = sp.Integer(1)
    graph = 6 * (1 + A)
    factor = sp.expand(C - graph)
    base = 2 * (3 + 4 * A)
    b = sp.expand(base + 2 * factor)
    companion = sp.expand(12 * A**2 - 8 * A * C + 12 * A - 9)
    branch = delta(b)
    equal(branch, -factor**2 * companion, "repeated total graph factorization")
    residual = sp.expand(-factor * companion)
    equal(branch, factor * residual, "residual includes graph")

    factor_h = homogenize_c(factor, 1)
    companion_h = homogenize_c(companion, 1)
    residual_h = homogenize_c(residual, 2)
    equal(residual_h, -factor_h * companion_h,
          "individual-degree homogenization product")
    require(factor_h.subs({A: 0, C: 1, Z: 0}) == 1,
            "marked graph misses finite-base infinity")
    require(companion_h.subs({A: 0, C: 1, Z: 0}) == 0,
            "distinct companion carries finite-base infinity")
    nonzero(companion.subs(A, 0), "companion nonvertical")
    require(sp.factor(companion) == companion, "primitive linear companion")
    require(sp.expand(factor - companion) != 0, "companion distinct from graph")

    # This relation makes A a unit on the companion and exhibits G_m.
    equal(companion, A * (12 * A + 12 - 8 * C) - 9,
          "Laurent inverse relation")
    parameter_c = (12 * A**2 + 12 * A - 9) / (8 * A)
    equal(sp.together(companion.subs(C, parameter_c)), 0,
          "Laurent parametrization")
    return (
        "s=1,Q=2",
        "Delta=-F^2*G",
        str(factor),
        str(companion),
        "A^{-1}=(12A+12-8C)/9",
        "normalization=G_m=P1-{0,infinity}",
    )


def main() -> None:
    theorem_body = CANON.read_text(encoding="utf-8").split("---", 2)[2]
    proof_body = theorem_body.split("## 7. Exact replay", 1)[0]
    require(hashlib.sha256(proof_body.encode()).hexdigest() == CANON_PROOF_BODY_SHA256,
            "canonical theorem proof-body hash")

    # Algebraic UFD factor identity used in the classification, checked before
    # any graph examples.  The proof of uniqueness is recorded in REPORT.md.
    z = sp.symbols("z")
    equal(z**3 - 1 - sp.Rational(3, 2) * (z**2 - 1),
          sp.Rational(1, 2) * (z - 1) ** 2 * (2 * z + 1),
          "UFD marked-root factor identity")

    profiles = (
        ("zero", sp.Integer(0)),
        ("constant_one", sp.Integer(1)),
        ("nonconstant", 1 + 2 * A - 3 * A**2),
    )
    built = {}
    quotient_rows = []
    completion_rows = []
    for name, profile_s in profiles:
        factor, graph, base, quotients = formal_quotient(profile_s, MAX_N)
        built[name] = (factor, graph, base, quotients)
        quotient_rows.append(
            (
                name,
                tuple(degree_c(item) for item in quotients),
                tuple(str(beta(index)) for index in range(MAX_N + 1)),
            )
        )
        completion_rows.extend(
            completion_descent(name, factor, base, quotients)
        )

    # Attack degree cases on a genuinely nonconstant graph.
    factor, _, base, quotients = built["nonconstant"]
    mismatch_rows = mismatch_audit(factor, base, quotients)
    repeated_control = repeated_graph_control()

    semantic = (
        CANON_PROOF_BODY_SHA256,
        tuple(quotient_rows),
        tuple(completion_rows),
        mismatch_rows,
        repeated_control,
        GATES,
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "frozen independent semantic packet")

    print("THM3866_INDEPENDENT_HOSTILE_AUDIT_20260823")
    print("verdict=PASS;scope=all_polynomial_graphs;no_JC_consequence")
    print(f"canonical_theorem_proof_body_sha256={CANON_PROOF_BODY_SHA256}")
    print("formal_quotient_rows=" + repr(tuple(quotient_rows)).replace(" ", ""))
    print("completion_rows=" + repr(tuple(completion_rows)).replace(" ", ""))
    print("mismatch_rows=" + repr(mismatch_rows).replace(" ", ""))
    print("repeated_graph_control=" + repr(repeated_control).replace(" ", ""))
    print(f"semantic_sha256={semantic_sha256}")
    print(f"gates={GATES}")
    print(
        "interpretation=unique_formal_comparator_and_first_mismatch_survive;"
        "static_degree_drop_selects_nonvertical_reduced_companion;"
        "repeated_total_graph_does_not_absorb_infinity_component;"
        "two_missing_places_have_A_images_0_and_infinity"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
