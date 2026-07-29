#!/usr/bin/env python3
"""Exact fixed-cofactor atlas for THM-2949.

For every four-slot support (0,a,b,M) in the finite width range below,
this companion interpolates the same 35-by-35 degree-seven Macaulay
minor.  It checks its exact degree, three outside-grid values, ordinary
coefficient sign structure, and the first shifted Gregory--Newton base.

The width-three sidecar reconstructs all 216 distinguished cofactors by
an exact adjugate calculation and computes their primitive common gcd.
Only integer/rational FLINT arithmetic is truth-bearing.
"""

from __future__ import annotations

import importlib.util
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path

from flint import fmpz_mat


WIDTH_MIN = 3
WIDTH_MAX = 12
MAX_NEWTON_BASE = 300
OUTSIDE_OFFSETS = (1, 2)

HERE = Path(__file__).resolve().parent
CONSTRUCTOR = HERE / "gmc_width_seven_eight_two_chart_resultant_closure_thm2943.py"
PARITY_GATE = HERE / "gmc_conjugate_pair_corank_parity_thm2947.py"
EXPECTED_CONSTRUCTOR_SHA256 = (
    "d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba"
)
EXPECTED_PARITY_GATE_SHA256 = (
    "d1bd09ff20925183f5488fcd8850469867f1dfad2bdb808504fc896708605744"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_digest(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


require(
    lf_digest(CONSTRUCTOR) == EXPECTED_CONSTRUCTOR_SHA256,
    "THM-2943 constructor dependency changed",
)
require(
    lf_digest(PARITY_GATE) == EXPECTED_PARITY_GATE_SHA256,
    "THM-2947 parity-gate dependency changed",
)

SPEC = importlib.util.spec_from_file_location("thm2943_for_2949", CONSTRUCTOR)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2943")
thm2943 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(thm2943)
thm2942 = thm2943.thm2942
t = thm2943.t

SELECTED = tuple(t.SELECTED_ROWS)
DELETE_ROW_POSITION = 30
DELETE_GLOBAL_ROW = SELECTED[DELETE_ROW_POSITION]
DELETE_COLUMN = 0
DELETE_TARGET = t.TARGET_MONOMIALS[DELETE_COLUMN]
COFACTOR_ROWS = tuple(
    row
    for position, row in enumerate(SELECTED)
    if position != DELETE_ROW_POSITION
)
COFACTOR_COLUMNS = tuple(column for column in range(36) if column != DELETE_COLUMN)

require(len(SELECTED) == 36, "selected chart size changed")
require(SELECTED[30:] == tuple(range(36, 42)), "selected F-row bank changed")
require(DELETE_GLOBAL_ROW == 36, "fixed deleted row changed")
require(DELETE_TARGET == (0, 0, 7), "fixed deleted target changed")


def normalized_positive_leading(polynomial):
    require(polynomial != 0, "zero polynomial cannot be normalized")
    return -polynomial if polynomial.leading_coefficient() < 0 else polynomial


def polynomial_digest(polynomial) -> str:
    return sha256(",".join(map(str, polynomial.coeffs())).encode()).hexdigest()


def sequence_digest(records: list[str]) -> str:
    return sha256(("\n".join(records) + "\n").encode()).hexdigest()


def newton_coefficients(values: list[int]) -> list[int]:
    answer: list[int] = []
    row = values
    while row:
        answer.append(row[0])
        row = [
            row[index + 1] - row[index]
            for index in range(len(row) - 1)
        ]
    return answer


def shift_newton(coefficients: list[int]) -> list[int]:
    return [
        coefficients[index] + coefficients[index + 1]
        for index in range(len(coefficients) - 1)
    ] + [coefficients[-1]]


def first_positive_newton_base(values: list[int]):
    coefficients = newton_coefficients(values)
    for base in range(MAX_NEWTON_BASE + 1):
        if all(value > 0 for value in coefficients):
            return base, coefficients
        coefficients = shift_newton(coefficients)
    raise RuntimeError("positive Newton base exceeds the declared search")


def ordinary_sign_profile(polynomial):
    signs = [
        1 if coefficient > 0 else -1
        for coefficient in polynomial.coeffs()
        if coefficient
    ]
    variations = sum(
        signs[index] != signs[index - 1]
        for index in range(1, len(signs))
    )
    negative_to_positive_one_cut = (
        variations == 1
        and all(
            signs[index] <= signs[index + 1]
            for index in range(len(signs) - 1)
        )
    )
    return variations, negative_to_positive_one_cut


def forms_and_rows(width: int, first: int, second: int):
    forms = thm2943.polynomial_forms(width, (0, first, second))
    rows, metadata = thm2942.macaulay_rows(forms)
    require(
        metadata[DELETE_GLOBAL_ROW] == (4, (0, 0, 3)),
        "deleted F-row metadata changed",
    )
    return rows


def fixed_cofactor_value(rows, depth: int) -> int:
    numeric = thm2943.evaluate_rows(rows, depth)
    return int(
        fmpz_mat(
            [
                [numeric[row][column] for column in COFACTOR_COLUMNS]
                for row in COFACTOR_ROWS
            ]
        ).det()
    )


def audit_width(width: int):
    degree_bound = 55 * width - 35
    records = []
    coefficient_payload = []
    newton_payload = []
    continuous_hostile = None
    for first, second in combinations(range(1, width), 2):
        rows = forms_and_rows(width, first, second)
        values = [
            fixed_cofactor_value(rows, depth)
            for depth in range(degree_bound + 1)
        ]
        polynomial = thm2943.interpolate(values)
        polynomial = normalized_positive_leading(polynomial)
        require(
            polynomial.degree() == degree_bound,
            f"cofactor degree dropped at {(width, first, second)}",
        )
        for offset in OUTSIDE_OFFSETS:
            depth = degree_bound + offset
            require(
                polynomial(depth) == fixed_cofactor_value(rows, depth),
                f"outside-grid mismatch at {(width, first, second, depth)}",
            )
        far_depth = 2 * degree_bound + 3
        require(
            polynomial(far_depth) == fixed_cofactor_value(rows, far_depth),
            f"far outside-grid mismatch at {(width, first, second)}",
        )

        oriented_values = [
            int(polynomial(depth))
            for depth in range(degree_bound + 1)
        ]
        base, newton = first_positive_newton_base(oriented_values)
        require(polynomial(base) > 0, "positive Newton base has bad value")
        require(
            all(polynomial(depth) != 0 for depth in range(base)),
            f"integer prefix zero at {(width, first, second)}",
        )
        variations, one_cut = ordinary_sign_profile(polynomial)
        if variations == 1:
            require(one_cut, "one variation has the wrong orientation")
            require(
                polynomial(0) < 0 < polynomial(base),
                "one-cut crossing is not confined below its Newton base",
            )
        if variations == 0:
            require(
                all(coefficient > 0 for coefficient in polynomial.coeffs()),
                "zero-variation polynomial is not coefficientwise positive",
            )
        require(
            (polynomial(0) > 0) == (variations % 2 == 0),
            "origin sign does not match the ordinary variation parity",
        )

        coefficient_payload.append(
            f"{first},{second}:"
            + ",".join(map(str, polynomial.coeffs()))
        )
        newton_payload.append(
            f"{first},{second}@{base}:" + ",".join(map(str, newton))
        )
        records.append(
            {
                "support": (first, second),
                "base": base,
                "variations": variations,
                "origin_sign": 1 if polynomial(0) > 0 else -1,
            }
        )
        if (width, first, second) == (11, 1, 2):
            continuous_hostile = tuple(
                1 if polynomial(depth) > 0 else -1
                for depth in (11, 12, 19, 20)
            )

    shifted = tuple(
        (record["support"], record["base"])
        for record in records
        if record["base"]
    )
    negative_origin = tuple(
        record["support"]
        for record in records
        if record["origin_sign"] < 0
    )
    variation_counts = tuple(
        (
            variation,
            sum(record["variations"] == variation for record in records),
        )
        for variation in sorted(
            {record["variations"] for record in records}
        )
    )
    require(
        len(records) == comb(width - 1, 2),
        "support count changed",
    )
    if width == 11:
        require(
            continuous_hostile == (1, -1, -1, 1),
            "width-eleven continuous-sign hostile changed",
        )
    return {
        "width": width,
        "families": len(records),
        "degree": degree_bound,
        "base0": sum(record["base"] == 0 for record in records),
        "max_base": max(record["base"] for record in records),
        "shifted": shifted,
        "negative_origin": negative_origin,
        "variation_counts": variation_counts,
        "coefficient_digest": sequence_digest(coefficient_payload),
        "newton_digest": sequence_digest(newton_payload),
        "continuous_hostile": continuous_hostile,
    }


def width_thirteen_sign_boundary():
    """Exact first width-thirteen control in the higher-cut regime."""
    width, first, second = 13, 1, 4
    degree_bound = 55 * width - 35
    rows = forms_and_rows(width, first, second)
    values = [
        fixed_cofactor_value(rows, depth)
        for depth in range(degree_bound + 1)
    ]
    polynomial = normalized_positive_leading(thm2943.interpolate(values))
    require(polynomial.degree() == degree_bound, "boundary degree changed")
    for depth in (
        degree_bound + 1,
        degree_bound + 2,
        2 * degree_bound + 3,
    ):
        require(
            polynomial(depth) == fixed_cofactor_value(rows, depth),
            "width-thirteen boundary outside-grid mismatch",
        )
    base, newton = first_positive_newton_base(
        [int(polynomial(depth)) for depth in range(degree_bound + 1)]
    )
    variations, one_cut = ordinary_sign_profile(polynomial)
    require((base, variations, one_cut) == (21, 3, False), "boundary changed")
    require(
        all(polynomial(depth) != 0 for depth in range(base)),
        "boundary integer prefix acquired a zero",
    )
    return {
        "support": (0, 1, 4, 13),
        "degree": degree_bound,
        "base": base,
        "variations": variations,
        "coefficient_digest": polynomial_digest(polynomial),
        "newton_digest": sequence_digest([",".join(map(str, newton))]),
    }


def all_distinguished_cofactor_values(rows, depth: int):
    numeric = thm2943.evaluate_rows(rows, depth)
    matrix = fmpz_mat(
        [
            [numeric[row][column] for column in range(36)]
            for row in SELECTED
        ]
    )
    determinant = matrix.det()
    require(determinant != 0, f"width-three chart singular at depth {depth}")
    inverse = matrix.inv()
    answer = {}
    for row_position in range(30, 36):
        for column in range(36):
            signed_cofactor = determinant * inverse[column, row_position]
            value = int(signed_cofactor)
            require(signed_cofactor == value, "nonintegral adjugate entry")
            if (row_position + column) % 2:
                value = -value
            answer[(row_position, column)] = value
    return answer


def primitive_positive(polynomial):
    polynomial = normalized_positive_leading(polynomial)
    return polynomial // polynomial.content()


def width_three_full_atlas():
    width = 3
    degree_bound = 55 * width - 35
    rows = forms_and_rows(width, 1, 2)
    addresses = tuple(
        (row_position, column)
        for row_position in range(30, 36)
        for column in range(36)
    )
    values = {address: [] for address in addresses}
    for depth in range(degree_bound + 1):
        level = all_distinguished_cofactor_values(rows, depth)
        for address in addresses:
            values[address].append(level[address])
    polynomials = {
        address: thm2943.interpolate(value_list)
        for address, value_list in values.items()
    }
    for polynomial in polynomials.values():
        require(polynomial != 0, "distinguished cofactor vanished identically")
        require(
            polynomial.degree() == degree_bound,
            "distinguished cofactor degree changed",
        )
    for depth in (
        degree_bound + 1,
        degree_bound + 2,
        2 * degree_bound + 3,
    ):
        level = all_distinguished_cofactor_values(rows, depth)
        for address, polynomial in polynomials.items():
            require(
                polynomial(depth) == level[address],
                f"216-bank outside-grid mismatch at {(address, depth)}",
            )

    common = None
    for polynomial in polynomials.values():
        primitive = primitive_positive(polynomial)
        common = primitive if common is None else common.gcd(primitive)
    require(common is not None, "common cofactor gcd missing")
    common = primitive_positive(common)
    require(common.degree() == 60, "width-three common seam degree changed")

    expected_quintic = (
        8788 * thm2943.X**5
        + 54873 * thm2943.X**4
        + 126718 * thm2943.X**3
        + 132729 * thm2943.X**2
        + 61288 * thm2943.X
        + 9732
    )
    expected_common = (
        (thm2943.X + 1) ** 2
        * (2 * thm2943.X + 1) ** 5
        * (49 * thm2943.X**2 + 99 * thm2943.X + 38) ** 5
        * (thm2943.X + 2) ** 19
        * (thm2943.X + 3) ** 19
        * expected_quintic
    )
    require(common == expected_common, "width-three common seam changed")

    fixed = primitive_positive(
        polynomials[(DELETE_ROW_POSITION, DELETE_COLUMN)]
    )
    quotient, remainder = divmod(fixed, common)
    require(remainder == 0, "common seam does not divide the fixed cofactor")
    quotient_unit, quotient_factors = quotient.factor()
    require(quotient_unit == 1, "fixed quotient unit changed")
    factor_profile = tuple(
        (factor.degree(), exponent)
        for factor, exponent in quotient_factors
    )
    require(
        factor_profile == ((66, 1), (1, 4))
        or factor_profile == ((1, 4), (66, 1)),
        "fixed quotient factor profile changed",
    )
    residual66 = next(
        factor
        for factor, exponent in quotient_factors
        if factor.degree() == 66 and exponent == 1
    )
    require(
        all(coefficient > 0 for coefficient in residual66.coeffs()),
        "degree-66 residual lost coefficient positivity",
    )

    bank_payload = [
        f"{row},{column}:"
        + ",".join(map(str, primitive_positive(polynomial).coeffs()))
        for (row, column), polynomial in sorted(polynomials.items())
    ]
    return {
        "addresses": len(addresses),
        "nonzero": sum(polynomial != 0 for polynomial in polynomials.values()),
        "degree": degree_bound,
        "common_degree": common.degree(),
        "common_digest": polynomial_digest(common),
        "bank_digest": sequence_digest(bank_payload),
        "fixed_quotient_profile": tuple(sorted(factor_profile)),
        "residual66_digest": polynomial_digest(residual66),
    }


def format_support(support) -> str:
    return f"{support[0]}:{support[1]}"


def main() -> None:
    widths = tuple(range(WIDTH_MIN, WIDTH_MAX + 1))
    work_widths = tuple(reversed(widths))
    with ProcessPoolExecutor(max_workers=min(4, len(widths))) as executor:
        summaries = sorted(
            executor.map(audit_width, work_widths),
            key=lambda summary: summary["width"],
        )
    boundary = width_thirteen_sign_boundary()
    atlas = width_three_full_atlas()

    require(
        sum(summary["families"] for summary in summaries)
        == sum(comb(width - 1, 2) for width in widths),
        "global support count changed",
    )

    print("THM-2949 FIXED RANK-THIRTY-FIVE COFACTOR NEWTON ATLAS")
    print(
        f"widths={WIDTH_MIN}..{WIDTH_MAX};"
        f"families={sum(summary['families'] for summary in summaries)};"
        f"fixed_delete=rowpos{DELETE_ROW_POSITION}/global{DELETE_GLOBAL_ROW}"
        f"/column{DELETE_COLUMN}/target{DELETE_TARGET}"
    )
    print(
        "outside_grid=degree+1,degree+2,2*degree+3;"
        f"max_Newton_search={MAX_NEWTON_BASE}"
    )
    for summary in summaries:
        shifted = ",".join(
            f"{format_support(support)}@{base}"
            for support, base in summary["shifted"]
        ) or "none"
        negative = ",".join(
            format_support(support)
            for support in summary["negative_origin"]
        ) or "none"
        variations = ",".join(
            f"{variation}:{count}"
            for variation, count in summary["variation_counts"]
        )
        print(
            f"M={summary['width']};families={summary['families']};"
            f"degree={summary['degree']};base0={summary['base0']};"
            f"max_base={summary['max_base']};"
            f"ordinary_variations={variations};negative_origin={negative};"
            f"shifted={shifted};"
            f"coeff_digest={summary['coefficient_digest']};"
            f"Newton_digest={summary['newton_digest']}"
        )
    print(
        "width3_all216="
        f"addresses:{atlas['addresses']};nonzero:{atlas['nonzero']};"
        f"degree:{atlas['degree']};common_degree:{atlas['common_degree']};"
        f"common_digest:{atlas['common_digest']};"
        f"bank_digest:{atlas['bank_digest']}"
    )
    print(
        "width3_fixed_quotient="
        f"profile:{atlas['fixed_quotient_profile']};"
        f"residual66_positive:true;"
        f"residual66_digest:{atlas['residual66_digest']}"
    )
    width_eleven = next(
        summary for summary in summaries if summary["width"] == 11
    )
    print(
        "width11_continuous_sign_hostile="
        "support:0,1,2,11;"
        "signs_at_11,12,19,20:"
        f"{width_eleven['continuous_hostile']};"
        "two_positive_real_crossings_but_no_integer_zero"
    )
    print(
        "width13_sign_boundary="
        f"support:{boundary['support']};degree:{boundary['degree']};"
        f"first_positive_Newton_base:{boundary['base']};"
        f"ordinary_variations:{boundary['variations']};"
        f"coeff_digest:{boundary['coefficient_digest']};"
        f"Newton_digest:{boundary['newton_digest']};"
        "scope:single_control_not_width13_closure"
    )
    print(
        "consequence=the fixed 35-minor is nonzero at every integer depth "
        f"for all supports of widths {WIDTH_MIN} through {WIDTH_MAX}"
    )
    print(
        "scope=finite exact atlas; shared degree-60 seam is chart-dependent; "
        "no arbitrary-width or real-ray resultant-positivity claim"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
