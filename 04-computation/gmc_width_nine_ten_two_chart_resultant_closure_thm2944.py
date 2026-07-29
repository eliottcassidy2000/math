#!/usr/bin/env python3
"""Exact companion for THM-2944.

This is the width-nine/ten continuation of THM-2943.  It interpolates
only the original raw Macaulay determinant; the audited universal
factorization then reconstructs the common resultant and the mutated
determinant.  Direct outside-grid determinants audit both charts.

For every support (0,a,b,M), M in {9,10}, the script proves:

* gcd(q200*K/g, P_alt/g)=1 for g=gcd(K,P_alt);
* the full common factor q200^5*c300*g*Res has nonnegative
  coefficients and positive constant term;
* the original/mutated determinants therefore have no common
  nonnegative integral-depth zero; and
* the radical repeated divisor of g*Res is the explicit negative-root
  product depending only on (M,a), not on b.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from itertools import combinations
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


SOURCE = Path(__file__).with_name(
    "gmc_width_seven_eight_two_chart_resultant_closure_thm2943.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
SOURCE_SHA256 = (
    "d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba"
)
require(
    sha256(SOURCE_BYTES).hexdigest() == SOURCE_SHA256,
    "THM-2943 two-chart dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2943_exact", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2943")
thm2943 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(thm2943)

thm2942 = thm2943.thm2942
thm2925 = thm2943.thm2925
t = thm2943.t
X_Z = t.POLY_X


def interpolate_original(rows, degree_bound: int):
    values = []
    for depth in range(degree_bound + 1):
        numeric_rows = thm2943.evaluate_rows(rows, depth)
        values.append(
            thm2943.determinant(numeric_rows, thm2943.ORIGINAL_F)
        )
    return thm2943.interpolate(values)


def squarefree_part(polynomial):
    polynomial = thm2943.primitive_positive_leading(polynomial)
    derivative = polynomial.derivative()
    require(derivative != 0, "constant repeated-divisor polynomial")
    return thm2943.primitive_positive_leading(
        polynomial // polynomial.gcd(derivative)
    )


def expected_radical(width: int, first: int):
    answer = type(X_Z)(1)
    for index in range(1, first + 1):
        answer *= 2 * X_Z + 2 * index - 1
    for shift in range(first, width + 1):
        answer *= X_Z + shift
    return thm2943.primitive_positive_leading(answer)


def audit_family(
    width: int,
    first: int,
    second: int,
) -> tuple[str, bool, int, int]:
    interior = (0, first, second)
    forms = thm2943.polynomial_forms(width, interior)
    rows, _metadata = thm2942.macaulay_rows(forms)
    degree_bound = 58 * width - 36
    original = interpolate_original(rows, degree_bound)
    require(
        original.degree() <= degree_bound,
        "THM-2925 determinant degree invoice failed",
    )

    q200, c300, curvature, alternate = thm2943.flag_polynomials(forms)
    original_flag = q200**6 * c300 * curvature
    require(
        original % original_flag == 0,
        "THM-2942 original flag factor failed exact division",
    )
    resultant = original // original_flag
    mutated = q200**5 * c300 * alternate * resultant

    flag_common = thm2943.normalized_positive_leading(
        curvature.gcd(alternate)
    )
    require(
        flag_common(0) > 0
        and all(coefficient >= 0 for coefficient in flag_common.coeffs()),
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

    common_factor = thm2943.normalized_positive_leading(
        q200**5 * c300 * flag_common * resultant
    )
    require(
        common_factor(0) > 0
        and all(
            coefficient >= 0 for coefficient in common_factor.coeffs()
        ),
        "full common resultant content lost one-sign positivity",
    )
    require(
        thm2943.primitive_positive_leading(original.gcd(mutated))
        == thm2943.primitive_positive_leading(common_factor),
        "reconstructed common factor differs from the raw chart gcd",
    )
    reduced_original = original // common_factor
    reduced_mutated = mutated // common_factor
    require(
        reduced_original.gcd(reduced_mutated).degree() == 0,
        "raw reduced chart cofactors retained a common root",
    )

    discriminant_carrier = thm2943.primitive_positive_leading(
        flag_common * resultant
    )
    repeated_divisor = thm2943.primitive_positive_leading(
        discriminant_carrier.gcd(discriminant_carrier.derivative())
    )
    radical = squarefree_part(repeated_divisor)
    expected = expected_radical(width, first)
    require(
        radical == expected,
        "g*Res repeated-divisor radical changed",
    )
    require(
        all(root_shift > 0 for root_shift in range(first, width + 1)),
        "negative-root invoice changed",
    )

    direct_controls = 0
    for depth in (
        degree_bound + 1,
        degree_bound + 2,
        degree_bound + 3,
    ):
        numeric_rows = thm2943.evaluate_rows(rows, depth)
        require(
            thm2943.determinant(
                numeric_rows,
                thm2943.ORIGINAL_F,
            )
            == original(depth)
            and thm2943.determinant(
                numeric_rows,
                thm2943.MUTATED_F,
            )
            == mutated(depth),
            "outside-grid direct determinant control failed",
        )
        direct_controls += 1

    pluecker_depth = degree_bound + 4
    numeric_rows = thm2943.evaluate_rows(rows, pluecker_depth)
    p12, p67, p16, p27, p17, p26 = (
        thm2943.determinant(numeric_rows, chart)
        for chart in thm2943.PLUECKER_CHARTS
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
            str(common_factor.degree()),
            str(reduced_original.degree()),
            str(reduced_mutated.degree()),
            str(int(original_mixed)),
            thm2943.polynomial_digest(flag_common),
            thm2943.polynomial_digest(common_factor),
            thm2943.polynomial_digest(reduced_original),
            thm2943.polynomial_digest(reduced_mutated),
            thm2943.polynomial_digest(repeated_divisor),
            thm2943.polynomial_digest(radical),
        )
    )
    return record, original_mixed, direct_controls, 1


def main() -> None:
    records = []
    width_records = {9: [], 10: []}
    family_counts = {}
    mixed_counts = {}
    direct_controls = 0
    pluecker_checks = 0
    for width in (9, 10):
        family_counts[width] = 0
        mixed_counts[width] = 0
        for first, second in combinations(range(1, width), 2):
            record, mixed, controls, pluecker = audit_family(
                width,
                first,
                second,
            )
            records.append(record)
            width_records[width].append(record)
            family_counts[width] += 1
            mixed_counts[width] += int(mixed)
            direct_controls += controls
            pluecker_checks += pluecker

    require(
        family_counts == {9: 28, 10: 36},
        "normalized support census changed",
    )
    require(
        direct_controls == 192 and pluecker_checks == 64,
        "independent control census changed",
    )

    print("THM-2944 WIDTH-NINE/TEN TWO-CHART RESULTANT CLOSURE")
    print(f"two_chart_dependency_sha256={SOURCE_SHA256}")
    print(
        "normalized_supports=(0,a,b,M);"
        "widths=9:28,10:36;total=64"
    )
    print(
        "original_local_F=0,1,2,3,4,5;"
        "mutated_local_F=0,3,4,5,6,7"
    )
    print(
        "common_factor~q200^5*c300*Res*gcd(K,P_alt);"
        "coefficients=NONNEGATIVE;constant=POSITIVE"
    )
    print(
        "gcd(q200*K0,P0)=1;"
        "reconstructed_raw_reduced_cofactors=COPRIME"
    )
    print(
        "rad_gcd(g*Res,(g*Res)')="
        "prod_j=1..a(2n+2j-1)*prod_r=a..M(n+r);"
        "independent_of_b;roots=NEGATIVE"
    )
    print(
        f"original_mixed_pivots=width9:{mixed_counts[9]},"
        f"width10:{mixed_counts[10]},"
        f"total:{mixed_counts[9]+mixed_counts[10]}"
    )
    print(
        f"outside_interpolation_direct_controls={direct_controls};"
        f"three_term_Pluecker_checks={pluecker_checks}"
    )
    print(
        f"width9_digest_sha256={thm2943.lf_digest(width_records[9])}"
    )
    print(
        f"width10_digest_sha256={thm2943.lf_digest(width_records[10])}"
    )
    print(f"combined_digest_sha256={thm2943.lf_digest(records)}")
    print("consequence=first_window_SFC4_exact_widths_9_and_10")
    print(
        "scope=integer_depth_n>=0;four_slot_first_window_only;"
        "no width>=11 or arbitrary_width claim"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
