#!/usr/bin/env python3
"""Exact companion for THM-2948.

For every normalized four-slot support

    (0,a,b,M),  1 <= a < b < M,  6 <= M <= 10,

this script reconstructs the denominator-cleared first-window ternary
forms Q,C,F and one fixed degree-seven Macaulay determinant.  The
PROVED THM-2942 factorization

    Delta_0 = q200^6*c300*K*Res(Q,C,F)

is then used as an exact calibration: all selected-chart and common
Pluecker flag factors are divided out before any derivative-gcd test.

For the resulting bare resultant R the script proves, on all 110
families, that

    rad gcd(R,R')
      = prod_(j=1)^a (2*n+2*j-1) prod_(r=a)^M (n+r).

After every full seam multiplicity is removed, the residual core is
squarefree and coprime to the seam.  Thus no nonlinear repeated factor
is hidden by the two-chart common flag content.
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
    "THM-2943 calibrated constructor dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2943_exact", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2943")
thm2943 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(thm2943)

thm2942 = thm2943.thm2942
X = thm2943.t.POLY_X


def interpolate_original(rows, degree_bound: int):
    values = []
    for depth in range(degree_bound + 1):
        numeric_rows = thm2943.evaluate_rows(rows, depth)
        values.append(
            thm2943.determinant(
                numeric_rows,
                thm2943.ORIGINAL_F,
            )
        )
    return thm2943.interpolate(values)


def normalized(polynomial):
    return thm2943.primitive_positive_leading(polynomial)


def squarefree_part(polynomial):
    polynomial = normalized(polynomial)
    if polynomial.degree() <= 0:
        return polynomial
    return normalized(
        polynomial // polynomial.gcd(polynomial.derivative())
    )


def expected_seam(width: int, first: int):
    answer = type(X)(1)
    for index in range(1, first + 1):
        answer *= 2 * X + 2 * index - 1
    for shift in range(first, width + 1):
        answer *= X + shift
    return normalized(answer)


def remove_full_factor(polynomial, factor):
    exponent = 0
    while polynomial % factor == 0:
        polynomial //= factor
        exponent += 1
    return polynomial, exponent


def coefficient_variations(polynomial) -> int:
    signs = [
        1 if coefficient > 0 else -1
        for coefficient in polynomial.coeffs()
        if coefficient
    ]
    return sum(
        signs[index] != signs[index - 1]
        for index in range(1, len(signs))
    )


def family_record(width: int, first: int, second: int):
    forms = thm2943.polynomial_forms(
        width,
        (0, first, second),
    )
    rows, _metadata = thm2942.macaulay_rows(forms)
    degree_bound = 58 * width - 36
    determinant = interpolate_original(rows, degree_bound)
    require(
        determinant.degree() <= degree_bound,
        "THM-2925 determinant degree invoice failed",
    )

    q200, c300, curvature, _alternate = thm2943.flag_polynomials(forms)
    selected_flag = q200**6 * c300 * curvature
    require(
        selected_flag != 0 and determinant % selected_flag == 0,
        "THM-2942 selected-chart calibration failed exact division",
    )
    resultant = normalized(determinant // selected_flag)
    require(
        resultant.degree() == 46 * width - 26,
        "bare resultant degree changed",
    )
    require(
        resultant(0) > 0,
        "bare resultant lost its positive endpoint orientation",
    )

    repeated = normalized(resultant.gcd(resultant.derivative()))
    repeated_radical = squarefree_part(repeated)
    seam = expected_seam(width, first)
    require(
        repeated_radical == seam,
        "bare-resultant repeated-root radical left the Gamma seam",
    )

    core = resultant
    seam_exponents = []
    for index in range(1, first + 1):
        factor = 2 * X + 2 * index - 1
        core, exponent = remove_full_factor(core, factor)
        require(exponent >= 2, "half-integral seam stopped being repeated")
        seam_exponents.append(exponent)
    for shift in range(first, width + 1):
        factor = X + shift
        core, exponent = remove_full_factor(core, factor)
        require(exponent >= 2, "integral seam stopped being repeated")
        seam_exponents.append(exponent)
    core = normalized(core)
    require(
        core.gcd(core.derivative()).degree() == 0,
        "wall-removed bare-resultant core is not squarefree",
    )
    require(
        core.gcd(seam).degree() == 0,
        "wall-removed core retained a Gamma seam factor",
    )

    # The finite atlas has a second, independent positivity feature.  It
    # is recorded here but is not extrapolated beyond width ten.
    variations = coefficient_variations(core)
    require(
        variations == 0 and core(0) > 0,
        "wall-removed core lost coefficientwise positivity",
    )

    control_depth = degree_bound + 1
    numeric_rows = thm2943.evaluate_rows(rows, control_depth)
    require(
        thm2943.determinant(
            numeric_rows,
            thm2943.ORIGINAL_F,
        )
        == determinant(control_depth),
        "outside-grid determinant control failed",
    )

    record = ":".join(
        (
            str(width),
            str(first),
            str(second),
            str(resultant.degree()),
            str(repeated.degree()),
            str(core.degree()),
            ",".join(map(str, seam_exponents)),
            str(variations),
            thm2943.polynomial_digest(resultant),
            thm2943.polynomial_digest(repeated),
            thm2943.polynomial_digest(repeated_radical),
            thm2943.polynomial_digest(core),
        )
    )
    return record, repeated.degree(), core.degree()


def digest(records: list[str]) -> str:
    return sha256(("\n".join(records) + "\n").encode()).hexdigest()


def main() -> None:
    all_records = []
    width_records = {}
    family_counts = {}
    repeated_ranges = {}
    core_ranges = {}

    for width in range(6, 11):
        records = []
        repeated_degrees = []
        core_degrees = []
        for first, second in combinations(range(1, width), 2):
            record, repeated_degree, core_degree = family_record(
                width,
                first,
                second,
            )
            records.append(record)
            repeated_degrees.append(repeated_degree)
            core_degrees.append(core_degree)
        width_records[width] = records
        all_records.extend(records)
        family_counts[width] = len(records)
        repeated_ranges[width] = (
            min(repeated_degrees),
            max(repeated_degrees),
        )
        core_ranges[width] = (min(core_degrees), max(core_degrees))

    require(
        family_counts == {6: 10, 7: 15, 8: 21, 9: 28, 10: 36},
        "normalized support census changed",
    )
    require(len(all_records) == 110, "full finite atlas count changed")

    print("THM-2948 PURE RESULTANT NEGATIVE SEAM ATLAS")
    print(f"calibrated_constructor_sha256={SOURCE_SHA256}")
    print(
        "universe=(0,a,b,M),1<=a<b<M,6<=M<=10;"
        "counts=10,15,21,28,36;total=110"
    )
    print(
        "load_bearing_calibration="
        "R=Delta0/(q200^6*c300*K);"
        "no_gcd(K,P_alt)_factor_retained"
    )
    print("bare_resultant_degree=46*M-26")
    print(
        "rad_gcd(R,Rprime)="
        "prod_j=1..a(2n+2j-1)*prod_r=a..M(n+r);"
        "independent_of_b"
    )
    print(
        "all_seam_multiplicities>=2;"
        "wall_removed_core=SQUAREFREE;"
        "core_coprime_to_seam=YES"
    )
    print(
        "finite_atlas_core_coefficients=ONE_SIGN_POSITIVE;"
        "R(0)=POSITIVE"
    )
    for width in range(6, 11):
        print(
            f"width{width}_families={family_counts[width]};"
            f"repeated_degree_range={repeated_ranges[width][0]}"
            f"..{repeated_ranges[width][1]};"
            f"core_degree_range={core_ranges[width][0]}"
            f"..{core_ranges[width][1]};"
            f"digest_sha256={digest(width_records[width])}"
        )
    print(f"combined_digest_sha256={digest(all_records)}")
    print("outside_interpolation_direct_controls=110;PASS")
    print(
        "consequence=with_THM2945_no_nonnegative_depth_resultant_zero_"
        "for_widths6_to10"
    )
    print(
        "scope=finite_exact_first_window_atlas;"
        "all_width_seam_formula_and_transversality_OPEN"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
