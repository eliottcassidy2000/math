#!/usr/bin/env python3
"""Exact companion for THM-2943.

For every normalized four-slot support

    (0, a, b, M),  1 <= a < b < M,  M in {6, 7, 8},

this script reconstructs the denominator-cleared factorial moment forms
Q,C,F of THM-2925.  It interpolates the original and stable-mutated
36-by-36 Macaulay minors, verifies the universal THM-2942 resultant
factorization, and proves that their exact polynomial gcd has
nonnegative coefficients and nonzero constant term.  The two residual
flag cofactors are coprime.  Hence the two charts cannot vanish
simultaneously at any integer depth n >= 0.

Width six is retained as an inherited control.  The new theorem
consequence is exactly widths seven and eight.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from flint import fmpq_poly, fmpz_mat


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


SOURCE = Path(__file__).with_name(
    "gmc_macaulay_extraneous_flag_pluecker_thm2942.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
SOURCE_SHA256 = (
    "85ce40de2aa777b8af091dfec934de68beeb42676bea3b48dbf815990a51e0e9"
)
require(
    sha256(SOURCE_BYTES).hexdigest() == SOURCE_SHA256,
    "THM-2942 factorization dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2942_exact", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2942")
thm2942 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(thm2942)

thm2925 = thm2942.thm2925
t = thm2942.t
X = fmpq_poly([0, 1])

BASE_Q = tuple(range(20))
BASE_C = tuple(range(21, 30)) + (35,)
F_GLOBAL = tuple(range(36, 46))
ORIGINAL_F = (0, 1, 2, 3, 4, 5)
MUTATED_F = (0, 3, 4, 5, 6, 7)
PLUECKER_CHARTS = (
    ORIGINAL_F,
    MUTATED_F,
    (0, 1, 3, 4, 5, 6),
    (0, 2, 3, 4, 5, 7),
    (0, 1, 3, 4, 5, 7),
    (0, 2, 3, 4, 5, 6),
)

require(
    tuple(t.SELECTED_ROWS)
    == BASE_Q + BASE_C + tuple(F_GLOBAL[index] for index in ORIGINAL_F),
    "original selected-row chart changed",
)


def lf_digest(records: list[str]) -> str:
    return sha256(("\n".join(records) + "\n").encode()).hexdigest()


def normalized_positive_leading(polynomial):
    require(polynomial != 0, "zero polynomial cannot be normalized")
    if polynomial.leading_coefficient() < 0:
        return -polynomial
    return polynomial


def primitive_positive_leading(polynomial):
    polynomial = normalized_positive_leading(polynomial)
    content = polynomial.content()
    require(content > 0, "polynomial content is not positive")
    return polynomial // content


def polynomial_digest(polynomial) -> str:
    return sha256(
        ",".join(map(str, polynomial.coeffs())).encode()
    ).hexdigest()


def polynomial_forms(
    width: int,
    interior: tuple[int, int, int],
):
    return tuple(
        {
            monomial: thm2925.scaled_coefficient(
                width,
                monomial,
                interior,
            )
            for monomial in t.MONOMIALS[order]
        }
        for order in t.ORDERS
    )


def evaluate_rows(rows, depth: int) -> list[list[int]]:
    return [
        [
            int(coefficient(depth))
            if callable(coefficient)
            else int(coefficient)
            for coefficient in row
        ]
        for row in rows
    ]


def selection(chart: tuple[int, ...]) -> tuple[int, ...]:
    return BASE_Q + BASE_C + tuple(F_GLOBAL[index] for index in chart)


def determinant(
    numeric_rows: list[list[int]],
    chart: tuple[int, ...],
) -> int:
    selected = selection(chart)
    return int(fmpz_mat([numeric_rows[index] for index in selected]).det())


def interpolate(values: list[int]):
    differences: list[int] = []
    row = values
    while row:
        differences.append(row[0])
        row = [
            row[index + 1] - row[index]
            for index in range(len(row) - 1)
        ]
    polynomial = fmpq_poly(0)
    basis = fmpq_poly(1)
    for index, coefficient in enumerate(differences):
        if coefficient:
            polynomial += coefficient * basis
        basis *= (X - index) / (index + 1)
    require(
        polynomial.denom() == 1,
        "interpolated Macaulay minor is not integral",
    )
    return polynomial.numer()


def interpolate_pair(rows, degree_bound: int):
    original_values = []
    mutated_values = []
    for depth in range(degree_bound + 1):
        numeric_rows = evaluate_rows(rows, depth)
        original_values.append(determinant(numeric_rows, ORIGINAL_F))
        mutated_values.append(determinant(numeric_rows, MUTATED_F))
    return interpolate(original_values), interpolate(mutated_values)


def flag_polynomials(forms):
    quadratic, cubic, _quartic = forms
    q200 = quadratic[(2, 0, 0)]
    q110 = quadratic[(1, 1, 0)]
    q020 = quadratic[(0, 2, 0)]
    q011 = quadratic[(0, 1, 1)]
    q002 = quadratic[(0, 0, 2)]
    c300 = cubic[(3, 0, 0)]
    c210 = cubic[(2, 1, 0)]
    c120 = cubic[(1, 2, 0)]
    c021 = cubic[(0, 2, 1)]
    c012 = cubic[(0, 1, 2)]
    curvature = (
        c120 * q200**2
        - c210 * q110 * q200
        - c300 * q020 * q200
        + c300 * q110**2
    )
    alternate = (
        c012 * q020 * q200**2
        - c021 * q011 * q200**2
        - c210 * q002 * q020 * q200
        + c210 * q011**2 * q200
        + c300 * q002 * q020 * q110
        - c300 * q011**2 * q110
    )
    return q200, c300, curvature, alternate


def audit_family(
    width: int,
    first: int,
    second: int,
) -> tuple[str, bool, int, int]:
    interior = (0, first, second)
    forms = polynomial_forms(width, interior)
    rows, _metadata = thm2942.macaulay_rows(forms)
    degree_bound = 58 * width - 36
    original, mutated = interpolate_pair(rows, degree_bound)
    require(
        original.degree() <= degree_bound
        and mutated.degree() <= degree_bound,
        "THM-2925 determinant degree invoice failed",
    )

    q200, c300, curvature, alternate = flag_polynomials(forms)
    original_flag = q200**6 * c300 * curvature
    mutated_flag = q200**5 * c300 * alternate
    require(
        original % original_flag == 0
        and mutated % mutated_flag == 0,
        "THM-2942 flag factor failed exact division",
    )
    resultant_original = original // original_flag
    resultant_mutated = mutated // mutated_flag
    require(
        resultant_original == resultant_mutated,
        "the two charts produced different resultant quotients",
    )
    resultant = resultant_original

    flag_common = normalized_positive_leading(curvature.gcd(alternate))
    flag_coefficients = tuple(flag_common.coeffs())
    require(
        flag_common(0) > 0
        and all(coefficient >= 0 for coefficient in flag_coefficients),
        "common flag content lost one-sign positivity",
    )
    primitive_original = curvature // flag_common
    primitive_mutated = alternate // flag_common
    require(
        primitive_original.gcd(primitive_mutated).degree() == 0,
        "primitive flag cofactors retained a common root",
    )
    require(
        (q200 * primitive_original).gcd(primitive_mutated).degree() == 0,
        "mutated primitive cofactor acquired a q200 common root",
    )

    raw_common = normalized_positive_leading(original.gcd(mutated))
    expected_common = normalized_positive_leading(
        q200**5 * c300 * flag_common * resultant
    )
    require(
        primitive_positive_leading(raw_common)
        == primitive_positive_leading(expected_common),
        "raw determinant gcd is not associated to resultant times content",
    )
    association_numerator = raw_common.content()
    association_denominator = expected_common.content()
    raw_coefficients = tuple(raw_common.coeffs())
    require(
        raw_common(0) > 0
        and all(coefficient >= 0 for coefficient in raw_coefficients),
        "raw determinant gcd lost one-sign positivity",
    )
    reduced_original = original // raw_common
    reduced_mutated = mutated // raw_common
    require(
        reduced_original.gcd(reduced_mutated).degree() == 0,
        "raw reduced chart cofactors retained a common root",
    )

    direct_controls = 0
    for depth in (
        degree_bound + 1,
        degree_bound + 2,
        degree_bound + 3,
    ):
        numeric_rows = evaluate_rows(rows, depth)
        require(
            determinant(numeric_rows, ORIGINAL_F) == original(depth)
            and determinant(numeric_rows, MUTATED_F) == mutated(depth),
            "outside-grid direct determinant control failed",
        )
        direct_controls += 1

    pluecker_depth = degree_bound + 4
    numeric_rows = evaluate_rows(rows, pluecker_depth)
    p12, p67, p16, p27, p17, p26 = (
        determinant(numeric_rows, chart) for chart in PLUECKER_CHARTS
    )
    require(
        p12 * p67 - p16 * p27 + p17 * p26 == 0,
        "three-term Pluecker exchange failed",
    )

    original_coefficients = tuple(reduced_original.coeffs())
    original_mixed = (
        any(coefficient > 0 for coefficient in original_coefficients)
        and any(coefficient < 0 for coefficient in original_coefficients)
    )
    record = ":".join(
        (
            str(width),
            str(first),
            str(second),
            str(original.degree()),
            str(mutated.degree()),
            str(resultant.degree()),
            str(flag_common.degree()),
            str(raw_common.degree()),
            str(reduced_original.degree()),
            str(reduced_mutated.degree()),
            str(int(original_mixed)),
            str(association_numerator),
            str(association_denominator),
            polynomial_digest(flag_common),
            polynomial_digest(raw_common),
            polynomial_digest(reduced_original),
            polynomial_digest(reduced_mutated),
        )
    )
    return record, original_mixed, direct_controls, 1


def main() -> None:
    records = []
    family_counts = {}
    mixed_counts = {}
    direct_controls = 0
    pluecker_checks = 0
    for width in (6, 7, 8):
        family_counts[width] = 0
        mixed_counts[width] = 0
        for first, second in combinations(range(1, width), 2):
            record, mixed, controls, pluecker = audit_family(
                width,
                first,
                second,
            )
            records.append(record)
            family_counts[width] += 1
            mixed_counts[width] += int(mixed)
            direct_controls += controls
            pluecker_checks += pluecker

    require(
        family_counts == {6: 10, 7: 15, 8: 21},
        "normalized support census changed",
    )
    require(
        mixed_counts == {6: 7, 7: 12, 8: 18},
        "original mixed-pivot census changed",
    )
    require(
        direct_controls == 138 and pluecker_checks == 46,
        "independent control census changed",
    )

    print("THM-2943 WIDTH-SEVEN/EIGHT TWO-CHART RESULTANT CLOSURE")
    print(f"factorization_dependency_sha256={SOURCE_SHA256}")
    print(
        "normalized_supports=(0,a,b,M);"
        "widths=6:10,7:15,8:21;total=46"
    )
    print(
        "original_rows=Q:0..19,C:21..29+35,F:36..41;"
        "mutated_local_F=0,3,4,5,6,7;"
        "mutated_global_F=36,39,40,41,42,43"
    )
    print(
        "factorization=Delta0=q200^6*c300*K*Res;"
        "Delta1=q200^5*c300*P_alt*Res"
    )
    print(
        "raw_gcd~q200^5*c300*Res*gcd(K,P_alt)_over_Q[n];"
        "raw_gcd_coefficients=NONNEGATIVE;"
        "raw_gcd_constant=POSITIVE"
    )
    print(
        "primitive_flag_cofactors=COPRIME;"
        "raw_reduced_chart_cofactors=COPRIME"
    )
    print(
        "original_mixed_pivots=width6:7,width7:12,width8:18,total:37;"
        "two_chart_exchange_removes_all_common_walls"
    )
    print(
        f"outside_interpolation_direct_controls={direct_controls};"
        f"three_term_Pluecker_checks={pluecker_checks}"
    )
    print(f"family_census_digest_sha256={lf_digest(records)}")
    print(
        "consequence=first_window_SFC4_exact_widths_7_and_8;"
        "width6=independent_inherited_control"
    )
    print(
        "scope=integer_depth_n>=0;four_slot_first_window_only;"
        "no width>=9 or arbitrary_width claim"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
