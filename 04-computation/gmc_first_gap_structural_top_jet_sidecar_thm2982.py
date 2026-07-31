#!/usr/bin/env python3
"""Durable structural top-jet sidecar for THM-2982.

This companion packages three exact checks without importing any discovery
probe: the closed factorial top-three formula, the edge-local pure-resultant
Hessian, and the finite curvature-sign census.  In scratch it imports the
candidate THM-2982 atlas companion; after promotion it imports the identical
canonical byte surface.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

from flint import fmpq_mat


ROOT = Path(__file__).resolve().parents[1]
CANONICAL_ATLAS = (
    ROOT
    / "04-computation"
    / "gmc_first_gap_wall_stripped_norm_core_strict_ulc_through_thirty_four_thm2982.py"
)
SCRATCH_ATLAS = (
    ROOT
    / ".scratch"
    / "gmc_first_gap_wall_stripped_norm_core_strict_ulc_through_thirty_four_thm2982.py"
)
ATLAS_PATH = CANONICAL_ATLAS if CANONICAL_ATLAS.exists() else SCRATCH_ATLAS
EXPECTED_ATLAS_SHA256 = (
    "645353f04d2143f91b7b213e7223761c65c72a81c0ebb522184643fdd97ca24b"
)
EXPECTED_RECORD_DIGEST = (
    "c3b116e42666d32eb6712091413c550db6d0ed920cb9d513b439ff1a0cbc2b5e"
)
DEGREES = (2, 3, 4)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


require(lf_sha256(ATLAS_PATH) == EXPECTED_ATLAS_SHA256, "THM-2982 atlas hash changed")
SPEC = importlib.util.spec_from_file_location("thm2982_structural_atlas", ATLAS_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load THM-2982 atlas companion")
ATLAS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ATLAS
SPEC.loader.exec_module(ATLAS)
BASE = ATLAS.load_base("structural_sidecar")


def form_top_slices(width: int):
    forms = BASE.thm2943.polynomial_forms(width, (0, 1, 2))
    slices = [[], [], []]
    for form in forms:
        degree = max(poly.degree() for poly in form.values())
        for drop in range(3):
            slices[drop].append(
                {
                    monomial: int(poly[degree - drop])
                    for monomial, poly in form.items()
                }
            )
    return forms, slices


def chart_log_jets(slices):
    m0, m1, m2 = [
        fmpq_mat(ATLAS.selected_rows(BASE, forms_at_drop))
        for forms_at_drop in slices
    ]
    inverse = m0.inv()
    first = inverse * m1
    return (
        ATLAS.as_fraction(ATLAS.matrix_trace(first)),
        2 * ATLAS.as_fraction(ATLAS.matrix_trace(inverse * m2))
        - ATLAS.as_fraction(ATLAS.matrix_trace(first * first)),
    )


def power_sum_1(value: int) -> Fraction:
    return Fraction(value * (value + 1), 2)


def power_sum_2(value: int) -> Fraction:
    return Fraction(value * (value + 1) * (2 * value + 1), 6)


def tensor_top_three(order: int, directions: tuple[int, ...]):
    total = sum(directions)
    log1 = power_sum_1(total) / order - sum(
        (power_sum_1(value) for value in directions), start=Fraction(0)
    )
    log2 = -power_sum_2(total) / (2 * order * order) + sum(
        (power_sum_2(value) / 2 for value in directions), start=Fraction(0)
    )
    leading = Fraction(order**total)
    return leading, leading * log1, leading * (log1 * log1 / 2 + log2)


def denominator_top_three(width: int, order: int):
    first = (
        (order - 1) * power_sum_1(width - 1) + (order - 2) * width
    )
    log2 = -(
        (order - 1) * power_sum_2(width - 1) + (order - 2) * width * width
    )
    return Fraction(1), first, (first * first + log2) / 2


def coefficient_top_three(width: int, order: int, monomial):
    directions = tuple(
        offset
        for offset, count in enumerate(monomial)
        for _ in range(count)
    )
    signed = [Fraction(0), Fraction(0), Fraction(0)]
    for mask in range(1 << order):
        selected = tuple(
            width if mask & (1 << position) else direction
            for position, direction in enumerate(directions)
        )
        sign = -1 if mask.bit_count() % 2 else 1
        row = tensor_top_three(order, selected)
        for index in range(3):
            signed[index] += sign * row[index]
    multinomial = math.factorial(order)
    for count in monomial:
        multinomial //= math.factorial(count)
    signed = [multinomial * value for value in signed]
    denominator = denominator_top_three(width, order)
    return (
        signed[0],
        denominator[1] * signed[0] + signed[1],
        denominator[2] * signed[0]
        + denominator[1] * signed[1]
        + signed[2],
    )


def add_exp(left, right):
    return tuple(left[index] + right[index] for index in range(3))


def multiply(left, right):
    answer = {}
    for a, ca in left.items():
        for b, cb in right.items():
            exponent = add_exp(a, b)
            answer[exponent] = answer.get(exponent, Fraction(0)) + ca * cb
    return {exponent: value for exponent, value in answer.items() if value}


def linear_power(coefficients, exponent: int):
    answer = {(0, 0, 0): Fraction(1)}
    linear = {
        tuple(int(index == axis) for index in range(3)): coefficients[axis]
        for axis in range(3)
        if coefficients[axis]
    }
    for _ in range(exponent):
        answer = multiply(answer, linear)
    return answer


def transform_form(form, inverse):
    powers = {
        (row, exponent): linear_power(
            tuple(ATLAS.as_fraction(inverse[row, col]) for col in range(3)), exponent
        )
        for row in range(3)
        for exponent in range(5)
    }
    answer = {}
    for exponent, coefficient in form.items():
        term = {(0, 0, 0): Fraction(int(coefficient))}
        for row, power in enumerate(exponent):
            term = multiply(term, powers[(row, power)])
        for monomial, value in term.items():
            answer[monomial] = answer.get(monomial, Fraction(0)) + value
    return {monomial: value for monomial, value in answer.items() if value}


def transformed_slices(width: int, slices):
    response = fmpq_mat(
        [
            [order**offset - order**width for offset in (0, 1, 2)]
            for order in DEGREES
        ]
    )
    inverse = response.inv()
    transformed = [
        [transform_form(slices[drop][index], inverse) for drop in range(3)]
        for index in range(3)
    ]
    for index, degree in enumerate(DEGREES):
        base = tuple(degree if axis == index else 0 for axis in range(3))
        require(transformed[index][0] == {base: Fraction(1)}, (width, index))
    return transformed


def local_resultant_data(transformed):
    first = Fraction(0)
    second = Fraction(0)
    self_curvatures = []
    for index, degree in enumerate(DEGREES):
        base = tuple(degree if axis == index else 0 for axis in range(3))
        alpha = transformed[index][1].get(base, Fraction(0))
        beta = transformed[index][2].get(base, Fraction(0))
        curvature = 2 * beta - alpha * alpha
        self_curvatures.append(curvature)
        multiplicity = math.prod(DEGREES) // degree
        first += multiplicity * alpha
        second += multiplicity * curvature

    products = []
    zero_labels = []
    for left in range(3):
        for right in range(left + 1, 3):
            remaining = 3 - left - right
            for transfer in range(1, min(DEGREES[left], DEGREES[right]) + 1):
                left_exp = [0, 0, 0]
                left_exp[left] = DEGREES[left] - transfer
                left_exp[right] = transfer
                right_exp = [0, 0, 0]
                right_exp[right] = DEGREES[right] - transfer
                right_exp[left] = transfer
                c_left = transformed[left][1].get(tuple(left_exp), Fraction(0))
                c_right = transformed[right][1].get(tuple(right_exp), Fraction(0))
                product = c_left * c_right
                products.append(product)
                second -= 2 * DEGREES[remaining] * transfer * product
                if not product:
                    zero_labels.append(f"{DEGREES[left]}-{DEGREES[right]}:a{transfer}")
    return first, second, self_curvatures, products, zero_labels


def trace_resultant_log_jets(forms, slices):
    chart_u, chart_ell2 = chart_log_jets(slices)
    q200, c300, curvature, _alternate = BASE.thm2943.flag_polynomials(forms)
    q_u, q_ell2 = ATLAS.polynomial_log_jets(q200)
    c_u, c_ell2 = ATLAS.polynomial_log_jets(c300)
    k_u, k_ell2 = ATLAS.polynomial_log_jets(curvature)
    return (
        chart_u - 6 * q_u - c_u - k_u,
        chart_ell2 - 6 * q_ell2 - c_ell2 - k_ell2,
    )


def audit_width(width: int):
    forms, slices = form_top_slices(width)
    top_match = width <= 34
    if top_match:
        for form_index, order in enumerate(DEGREES):
            for monomial in BASE.thm2943.t.MONOMIALS[order]:
                closed = coefficient_top_three(width, order, monomial)
                actual = tuple(
                    Fraction(slices[drop][form_index][monomial]) for drop in range(3)
                )
                require(closed == actual, (width, order, monomial))
                require(all(value.denominator == 1 for value in closed), (width, monomial))

    transformed = transformed_slices(width, slices)
    local_u, local_ell2, self_curvatures, products, zero_labels = local_resultant_data(
        transformed
    )
    if top_match:
        trace_u, trace_ell2 = trace_resultant_log_jets(forms, slices)
        require(local_u == trace_u, f"first jet mismatch at M={width}")
        require(local_ell2 == trace_ell2, f"second jet mismatch at M={width}")

    require(all(value < 0 for value in self_curvatures), f"self sign at M={width}")
    require(sum(value > 0 for value in products) == 6, f"edge sign at M={width}")
    require(sum(value == 0 for value in products) == 1, f"edge zero at M={width}")
    require(zero_labels == ["3-4:a3"], f"edge zero label at M={width}")
    require(local_ell2 < 0, f"resultant curvature at M={width}")
    return {
        "M": width,
        "top_slice_closed_match": True if top_match else None,
        "hessian_trace_match": True if top_match else None,
        "negative_self_curvatures": 3,
        "positive_edge_products": 6,
        "zero_edge_label": "3-4:a3",
        "pure_resultant_ell2_negative": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = [audit_width(width) for width in range(6, 51)]
    lines = [json.dumps(row, sort_keys=True, separators=(",", ":")) for row in records]
    digest = sha256("\n".join(lines).encode()).hexdigest()
    require(digest == EXPECTED_RECORD_DIGEST, "structural record digest changed")
    transcript = "\n".join(
        [
            "FIRST-GAP STRUCTURAL TOP-JET SIDECAR",
            f"atlas_dependency_sha256={lf_sha256(ATLAS_PATH)}",
            "closed_top_slices=M6..34;coefficient_records=899;scalar_equalities=2697",
            "resultant_Hessian_trace_matches=M6..34;jet_equalities=58",
            "curvature_sign_scope=M6..50;negative_self=135;positive_edges=270;structural_zeros=45",
            "structural_zero_label=3-4:a3;wall_square_sum=only_remaining_k1_sign_competitor",
            *lines,
            f"record_digest={digest}",
            "all_width_sign=NOT_ASSERTED",
            "all_exact_checks=PASS",
        ]
    ) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
