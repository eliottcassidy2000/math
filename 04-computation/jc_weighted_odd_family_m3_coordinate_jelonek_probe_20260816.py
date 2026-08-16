#!/usr/bin/env python3
"""Exact first beyond-cubic test for the explicit odd weighted Keller family.

THM-1605 preserves only a historical, mangled description of an ``E_m``
family.  The literal family currently available in canon is the cyclic
subfamily of THM-3448.  Reindexing its parameter by ``ell = 2*m - 3`` gives
generic degree ``2*m - 1`` and recovers the quadratic seed at ``m = 2``.

This companion performs the first lawful beyond-cubic calculation, ``m = 3``.
It reconstructs the determinant-one degree-five map, its inverse quintic, the
three actual source-coordinate resultants, their discriminants, and the two
Jelonek components.  In particular it checks exactly that

    Disc(E_x) / Disc(T), Disc(E_y) / Disc(T), Disc(E_z) / Disc(T)

are nonzero squares in Q(P,Q,C), while the pullback of ``Disc(T)`` is
``C^4*L_5``.  Thus all three primitive coordinate views have square class
``[L_5]`` on the target chart, but that sign class misses the genuine
``C = 0`` Jelonek component, whose generic inertia is a 3-cycle.

All truth gates use ``require`` and survive optimized Python.  There is no
randomness and no elapsed-time field in the transcript.
"""

from __future__ import annotations

from hashlib import sha256
from math import isqrt
from pathlib import Path

import sympy as sp


EXPECTED_SEMANTIC_SHA256: str | None = "ed6845e743f8554327653521f243817264b08d1ca864c8513c0b2af7ce17ac81"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def polynomial_digest(expression: sp.Expr, generators: tuple[sp.Symbol, ...]) -> str:
    polynomial = sp.Poly(sp.cancel(expression), *generators, domain=sp.QQ)
    payload = tuple(
        (monomial, int(coefficient.p), int(coefficient.q))
        for monomial, coefficient in polynomial.terms()
    )
    return sha256(repr(payload).encode("ascii")).hexdigest()


def polynomial_stats(expression: sp.Expr, generators: tuple[sp.Symbol, ...]) -> tuple[object, ...]:
    polynomial = sp.Poly(sp.cancel(expression), *generators, domain=sp.QQ)
    return (
        len(polynomial.terms()),
        polynomial.total_degree(),
        tuple(polynomial.degree(generator) for generator in generators),
        polynomial_digest(expression, generators),
    )


def rational_square_root(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    require(value >= 0, ("negative rational square", value))
    numerator_root = isqrt(int(value.p))
    denominator_root = isqrt(int(value.q))
    require(
        numerator_root * numerator_root == value.p
        and denominator_root * denominator_root == value.q,
        ("nonsquare rational", value),
    )
    return sp.Rational(numerator_root, denominator_root)


def perfect_square_root(
    expression: sp.Expr, generators: tuple[sp.Symbol, ...]
) -> sp.Expr:
    expression = sp.cancel(expression)
    numerator, denominator = expression.as_numer_denom()
    numerator_coefficient, numerator_factors = sp.factor_list(numerator, *generators)
    denominator_coefficient, denominator_factors = sp.factor_list(denominator, *generators)
    root = rational_square_root(
        sp.Rational(numerator_coefficient) / sp.Rational(denominator_coefficient)
    )
    for factor, exponent in numerator_factors:
        require(exponent % 2 == 0, ("odd numerator exponent", exponent))
        root *= factor ** (exponent // 2)
    for factor, exponent in denominator_factors:
        require(exponent % 2 == 0, ("odd denominator exponent", exponent))
        root /= factor ** (exponent // 2)
    root = sp.factor(root)
    require(sp.cancel(root * root - expression) == 0, "square-root reconstruction")
    return root


x, y, z = sp.symbols("x y z")
w, P, Q, C = sp.symbols("w P Q C")
X, Y, Z = sp.symbols("X Y Z")
A_t, B_t, C_t = sp.symbols("A B C_t")


def cyclic_seed(m: int) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Rational, int]:
    require(m >= 2, m)
    ell = 2 * m - 3
    p = sp.expand((ell + 1) * w**ell - (ell + 2) * w ** (ell + 1))
    R = sp.expand(w ** (ell + 1) * (1 - w))
    q = sp.expand(ell * w ** (ell + 1) - (ell + 1) * w ** (ell + 2))
    a = -sp.Rational(2 * ell + 1, 2 * ell)
    require(sp.diff(R, w) == p, (m, "R prime"))
    require(sp.diff(q, w) == sp.expand(w * sp.diff(p, w)), (m, "q prime"))
    require(p.subs(w, 0) == 0 and p.subs(w, 1) == -1, (m, "endpoints"))
    require(R.subs(w, 1) == 0, (m, "integral endpoint"))
    kappa = sp.diff(p, w).subs(w, 1)
    require(a == sp.factor(-(1 + kappa) / (2 + kappa)), (m, "a"))
    return p, R, q, a, ell


def source_map_m3() -> tuple[sp.Expr, ...]:
    p, _, q, a, ell = cyclic_seed(3)
    require((p, q, a, ell) == (4 * w**3 - 5 * w**4, 3 * w**4 - 4 * w**5, -sp.Rational(7, 6), 3), "m3 seed")
    u = 1 + x * y
    gamma = 1 + a * x * y + x**2 * z
    alpha_numerator = sp.expand(u + 3 * u**4 * gamma**2 - 4 * u**5 * gamma**3)
    beta_numerator = sp.expand(1 + 4 * u**3 * gamma**2 - 5 * u**4 * gamma**3)
    alpha = sp.cancel(alpha_numerator / x**2)
    beta = sp.cancel(beta_numerator / x)
    require(sp.denom(alpha) == 1 and sp.denom(beta) == 1, "apparent quotients cancel")
    F = (sp.expand(alpha), sp.expand(beta), sp.expand(x * gamma))
    degrees = tuple(sp.Poly(component, x, y, z).total_degree() for component in F)
    z_degrees = tuple(sp.Poly(component, x, y, z).degree(z) for component in F)
    require(degrees == (17, 16, 4) and z_degrees == (3, 3, 1), (degrees, z_degrees))
    jacobian = sp.factor(sp.Matrix(F).jacobian((x, y, z)).det())
    require(jacobian == 1, ("Jacobian", jacobian))

    source_w = sp.expand(u * gamma)
    source_P = sp.expand(F[1] * F[2])
    source_Q = sp.expand(F[0] * F[2] ** 2)
    inverse_identity = sp.factor(source_w**5 - source_w**4 + source_P * source_w - source_Q)
    require(inverse_identity == 0, ("inverse identity", inverse_identity))
    return F + (gamma, source_w, source_P, source_Q)


def coordinate_row(
    label: str,
    coordinate: sp.Symbol,
    relation: sp.Expr,
    inverse_equation: sp.Expr,
    branch_discriminant: sp.Expr,
    expected_c_root_exponent: int,
) -> tuple[object, ...]:
    eliminant = sp.cancel(sp.resultant(inverse_equation, relation, w))
    eliminant_polynomial = sp.Poly(eliminant, coordinate, P, Q, C, domain=sp.QQ)
    require(eliminant_polynomial.degree(coordinate) == 5, (label, "coordinate degree"))
    factor_coefficient, factors = sp.factor_list(eliminant, coordinate, P, Q, C)
    require(
        factor_coefficient != 0 and len(factors) == 1 and factors[0][1] == 1,
        (label, "irreducibility", factor_coefficient, len(factors)),
    )

    coordinate_discriminant = sp.factor(sp.discriminant(eliminant, coordinate))
    require(coordinate_discriminant != 0, (label, "zero discriminant"))
    ratio = sp.cancel(coordinate_discriminant / branch_discriminant)
    index_root = perfect_square_root(ratio, (P, Q, C))
    index_core = sp.cancel(index_root / C**expected_c_root_exponent)
    require(not index_core.has(C), (label, "unexpected C in index core"))
    require(sp.cancel(index_root**2 - ratio) == 0, (label, "index square"))
    require(
        sp.gcd(
            sp.Poly(branch_discriminant, P, Q, domain=sp.QQ),
            sp.Poly(index_core, P, Q, domain=sp.QQ),
        ).total_degree()
        == 0,
        (label, "branch/index common factor"),
    )
    return (
        label,
        sp.Rational(factor_coefficient),
        polynomial_stats(eliminant, (coordinate, P, Q, C)),
        polynomial_stats(index_core, (P, Q)),
        expected_c_root_exponent,
        polynomial_digest(coordinate_discriminant, (P, Q, C)),
        sp.sstr(eliminant) if label == "x" else "hash-pinned",
    )


def m3_coordinate_atlas() -> tuple[object, ...]:
    p, R, _, a, ell = cyclic_seed(3)
    require(ell == 3, ell)
    inverse_equation = sp.expand(-R + P * w - Q)
    require(inverse_equation == w**5 - w**4 + P * w - Q, inverse_equation)
    gamma = sp.expand(P - p)
    require(sp.diff(inverse_equation, w) == gamma, "T prime equals gamma")

    branch_discriminant = sp.factor(sp.discriminant(inverse_equation, w))
    resultant_discriminant = sp.factor(
        sp.resultant(inverse_equation, sp.diff(inverse_equation, w), w)
    )
    require(branch_discriminant == resultant_discriminant, "odd-degree resultant sign")
    branch_factors = sp.factor_list(branch_discriminant, P, Q)[1]
    require(len(branch_factors) == 1 and branch_factors[0][1] == 1, "branch irreducible and reduced")

    rows = (
        coordinate_row("x", X, X * gamma - C, inverse_equation, branch_discriminant, 10),
        coordinate_row("y", Y, C * Y - w + gamma, inverse_equation, branch_discriminant, 10),
        coordinate_row(
            "z",
            Z,
            C**2 * Z - gamma * (gamma * (gamma - 1 + a) - a * w),
            inverse_equation,
            branch_discriminant,
            20,
        ),
    )

    flat_eliminant = sp.factor(sp.resultant(inverse_equation, X - P, w))
    require(sp.expand(flat_eliminant - (X - P) ** 5) == 0, flat_eliminant)
    require(sp.discriminant(flat_eliminant, X) == 0, "flat hostile must have zero discriminant")
    return inverse_equation, branch_discriminant, rows, "flat_view=(X-P)^5_has_zero_discriminant"


def target_boundary(
    inverse_equation: sp.Expr, branch_discriminant: sp.Expr
) -> tuple[object, ...]:
    pullback = sp.factor(branch_discriminant.subs({P: B_t * C_t, Q: A_t * C_t**2}))
    L5 = sp.expand(
        3125 * A_t**4 * C_t**4
        - 2500 * A_t**3 * B_t * C_t**3
        + 256 * A_t**3 * C_t**2
        - 50 * A_t**2 * B_t**2 * C_t**2
        - 36 * A_t * B_t**3 * C_t
        + 256 * B_t**5 * C_t
        - 27 * B_t**4
    )
    require(pullback == C_t**4 * L5, ("pullback", pullback))
    require(L5.subs(C_t, 0) == -27 * B_t**4, "component separation at C=0")
    L5_factors = sp.factor_list(L5, A_t, B_t, C_t)[1]
    require(len(L5_factors) == 1 and L5_factors[0][1] == 1, "L5 irreducible and reduced")

    s, W = sp.symbols("s W")
    scaled = sp.expand(inverse_equation.subs({w: s * W, P: B_t * s**3, Q: A_t * s**6}) / s**4)
    require(scaled == s * W**5 - W**4 + B_t * W - A_t * s**2, scaled)
    require(
        sp.expand(scaled.subs(s, 0) - W * (B_t - W**3)) == 0,
        "C3 Newton residual",
    )

    # The m=2 control has the same even chart-index factor, but C=0 is not a
    # separate Jelonek component; THM-3448 proves the effectivity statement.
    _, R2, _, _, ell2 = cyclic_seed(2)
    require(ell2 == 1, ell2)
    T3 = sp.expand(-R2 + P * w - Q)
    D3 = sp.factor(sp.discriminant(T3, w))
    pullback3 = sp.factor(D3.subs({P: B_t * C_t, Q: A_t * C_t**2}))
    L3 = sp.factor(pullback3 / C_t**2)
    require(L3.subs(C_t, 0) == B_t**2 - 4 * A_t, L3.subs(C_t, 0))

    return (
        polynomial_stats(L5, (A_t, B_t, C_t)),
        sp.sstr(L5),
        "m3_Jelonek=V(C_t)_union_V(L5)",
        "C_t_inertia=C3_even;raw_discriminant_order=4;order_index=1;three_sheets_escape",
        "L5_inertia=transposition_odd;raw_discriminant_order=1;two_sheets_escape",
        "m3_square_class_after_pullback=[L5]_and_misses_C_t",
        polynomial_stats(L3, (A_t, B_t, C_t)),
        "m2_Jelonek=V(L3)_only;C_t^2_is_index_without_separate_effective_component",
    )


def family_rows() -> tuple[tuple[int, ...], ...]:
    rows = []
    for m in range(2, 11):
        _, _, _, a, ell = cyclic_seed(m)
        generic_degree = ell + 2
        first_coordinate_degree = 5 * ell + 2
        second_coordinate_degree = 5 * ell + 1
        c_discriminant_order = ell + 1
        c_inertia_sign = (-1) ** (ell - 1)
        require(generic_degree == 2 * m - 1, (m, generic_degree))
        require(a == -sp.Rational(4 * m - 5, 4 * m - 6), (m, a))
        require(c_discriminant_order == 2 * m - 2 and c_inertia_sign == 1, (m, c_discriminant_order))
        rows.append(
            (
                m,
                ell,
                generic_degree,
                first_coordinate_degree,
                second_coordinate_degree,
                4,
                c_discriminant_order,
                c_inertia_sign,
            )
        )
    return tuple(rows)


def main() -> None:
    F = source_map_m3()
    inverse_equation, branch_discriminant, coordinate_rows, flat_hostile = m3_coordinate_atlas()
    boundary = target_boundary(inverse_equation, branch_discriminant)
    rows = family_rows()

    map_stats = tuple(polynomial_stats(component, (x, y, z)) for component in F[:3])
    branch_stats = polynomial_stats(branch_discriminant, (P, Q))
    verdict = (
        "THM1605_historical_E_m_formula_absent;explicit_object_here_is_THM3448_cyclic_weighted_subfamily_ell=2m-3",
        "m3_map_det=1;ordinary_degrees=(17,16,4);generic_degree=5;global_monodromy=S5_by_THM3438",
        "x_y_z_resultants_are_irreducible_quintics_and_each_discriminant_is_D5_times_a_square",
        "common_generic_square_class=[D5];after_P=BC,Q=AC2_it_is_[L5]",
        "literal_fixed_map_-4_does_not_persist_as_a_normalization_constant",
        "primitive_coordinate_common_class_persists_by_trace_index_squares",
        "full_Jelonek_detection_fails:V(C)_is_real_nonproperness_but_C3_has_even_sign_and_C4_is_square",
        "m2_is_exceptional:C_factor_is_not_a_separate_Jelonek_component",
        "all_explicit_odd_subfamily_C_cycles_have_odd_length_and_even_sign;the_sign_resolvent_misses_C_for_every_m>=3",
        "no_identification_with_an_unstored_outside_formula;no_arbitrary_map_classification;no_JC2_or_LRC_consequence",
    )
    semantic_surface = (
        map_stats,
        sp.sstr(inverse_equation),
        branch_stats,
        sp.sstr(branch_discriminant),
        coordinate_rows,
        flat_hostile,
        boundary,
        rows,
        verdict,
    )
    semantic_sha256 = sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, (semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    print("Explicit odd weighted Keller family: first m=3 coordinate/Jelonek test")
    print("status=EXACT_DETERMINISTIC_COMPANION;family_scope=THM3448_cyclic_weighted_subfamily_ell=2m-3")
    print("historical_boundary=THM1605_contains_no_literal_E_m_definition;no_identity_with_an_unstored_family_is_claimed")
    print(f"finite_family_rows=(m,ell,generic_degree,degA,degB,degC,C_disc_order,C_inertia_sign)={rows}")
    print(f"m3_map_stats={map_stats};det=1;inverse_T={inverse_equation}")
    print(f"m3_branch_stats={branch_stats};D5={branch_discriminant}")
    for row in coordinate_rows:
        print(f"m3_coordinate_row={row}")
    print(f"hostile={flat_hostile}")
    print(f"target_boundary={boundary}")
    print(f"verdict={verdict}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;all_checks_survive_python_O;no_randomness;no_elapsed_fields")
    print("commands=python -B 04-computation/jc_weighted_odd_family_m3_coordinate_jelonek_probe_20260816.py;python -B -O 04-computation/jc_weighted_odd_family_m3_coordinate_jelonek_probe_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
