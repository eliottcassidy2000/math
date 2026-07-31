#!/usr/bin/env python3
"""Exact all-width self-curvature certificate for THM-2988."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

import sympy as sp


M, U, V = sp.symbols("M U V")
DEGREES = (2, 3, 4)
POWERS = {2: U, 3: V, 4: U**2}
ROOT = Path(__file__).resolve().parents[1]
STRUCTURAL_PATH = (
    ROOT / "04-computation/gmc_first_gap_structural_top_jet_sidecar_thm2982.py"
)
EXPECTED_STRUCTURAL_SHA256 = (
    "ca8d03d4ff3b647797623346eee285374b1595942ae66295a19b050ac8bf4547"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


def p1(value):
    return sp.Rational(1, 2) * value * (value + 1)


def p2(value):
    return sp.Rational(1, 6) * value * (value + 1) * (2 * value + 1)


def tensor_top_three(order: int, directions):
    count_m = sum(value == M for value in directions)
    fixed_sum = sum((value for value in directions if value != M), start=0)
    total = count_m * M + fixed_sum
    leading = POWERS[order] ** count_m * order**fixed_sum
    log1 = p1(total) / order - sum((p1(value) for value in directions), start=0)
    log2 = -p2(total) / (2 * order**2) + sum(
        (p2(value) / 2 for value in directions), start=0
    )
    return leading, leading * log1, leading * (log1**2 / 2 + log2)


def denominator_top_three(order: int):
    first = (order - 1) * p1(M - 1) + (order - 2) * M
    log2 = -((order - 1) * p2(M - 1) + (order - 2) * M**2)
    return sp.Integer(1), first, (first**2 + log2) / 2


def monomials(degree: int):
    return tuple(
        (a, b, degree - a - b)
        for a in range(degree + 1)
        for b in range(degree - a + 1)
    )


def coefficient_top_three(order: int, monomial):
    directions = tuple(
        offset for offset, count in enumerate(monomial) for _ in range(count)
    )
    signed = [sp.Integer(0), sp.Integer(0), sp.Integer(0)]
    for mask in range(1 << order):
        selected = tuple(
            M if mask & (1 << position) else direction
            for position, direction in enumerate(directions)
        )
        sign = -1 if mask.bit_count() % 2 else 1
        row = tensor_top_three(order, selected)
        for index in range(3):
            signed[index] += sign * row[index]
    multinomial = math.factorial(order)
    for count in monomial:
        multinomial //= math.factorial(count)
    signed = [sp.expand(multinomial * value) for value in signed]
    denominator = denominator_top_three(order)
    return (
        signed[0],
        sp.expand(denominator[1] * signed[0] + signed[1]),
        sp.expand(
            denominator[2] * signed[0]
            + denominator[1] * signed[1]
            + signed[2]
        ),
    )


def evaluate_form(form, point):
    answer = sp.Integer(0)
    for exponent, coefficient in form.items():
        term = coefficient
        for axis, power in enumerate(exponent):
            term *= point[axis] ** power
        answer += term
    return sp.expand(answer)


def symbolic_self_curvatures():
    rows = []
    raw_forms = []
    for order in DEGREES:
        rows.append([order**offset - POWERS[order] for offset in (0, 1, 2)])
        triples = {
            monomial: coefficient_top_three(order, monomial)
            for monomial in monomials(order)
        }
        raw_forms.append(
            [
                {monomial: triples[monomial][drop] for monomial in triples}
                for drop in range(3)
            ]
        )
    response = sp.Matrix(rows)
    determinant = sp.factor(response.det())
    adjugate = response.adjugate()
    curvatures = {}
    for index, degree in enumerate(DEGREES):
        column = tuple(adjugate[row, index] for row in range(3))
        alpha = evaluate_form(raw_forms[index][1], column)
        beta = evaluate_form(raw_forms[index][2], column)
        # alpha/det^degree and beta/det^degree are the two base coefficients.
        # Thus curvature=(2 beta det^degree-alpha^2)/det^(2 degree).
        curvatures[degree] = sp.expand(
            2 * beta * determinant**degree - alpha**2
        )
    return determinant, curvatures


def grouped_representation(expression):
    polynomial = sp.Poly(sp.expand(expression), U, V, M)
    grouped = {}
    for (u_power, v_power, m_power), coefficient in polynomial.terms():
        grouped[u_power, v_power] = (
            grouped.get((u_power, v_power), sp.Integer(0))
            + coefficient * M**m_power
        )
    return polynomial, grouped


def representation_digest(grouped) -> str:
    rows = []
    for pair in sorted(grouped, reverse=True):
        polynomial = sp.Poly(grouped[pair], M)
        rows.append(
            json.dumps(
                {
                    "u": pair[0],
                    "v": pair[1],
                    "degree": polynomial.degree(),
                    "coefficients": [str(value) for value in polynomial.all_coeffs()],
                },
                sort_keys=True,
                separators=(",", ":"),
            )
        )
    return sha256("\n".join(rows).encode()).hexdigest()


def self_curvature_record(degree: int, curvature_numerator):
    # Prove -curvature_numerator positive by one dominant exponential.
    expression, grouped = grouped_representation(-curvature_numerator)
    ordered = sorted(
        grouped,
        key=lambda pair: 2 ** pair[0] * 3 ** pair[1],
        reverse=True,
    )
    dominant = ordered[0]
    runner_up = ordered[1]
    dominant_poly = sp.Poly(grouped[dominant], M)
    require(dominant_poly.LC() > 0, (degree, "dominant coefficient"))
    negative_tariff = sum(
        (
            -coefficient
            for coefficient in dominant_poly.all_coeffs()[1:]
            if coefficient < 0
        ),
        sp.Integer(0),
    )
    lower_polys = [sp.Poly(grouped[pair], M) for pair in ordered[1:]]
    lower_degree = max(poly.degree() for poly in lower_polys)
    delta = max(0, lower_degree - dominant_poly.degree())
    absolute_tariff = sum(
        (
            sum(abs(coefficient) for coefficient in poly.all_coeffs())
            for poly in lower_polys
        ),
        sp.Integer(0),
    )
    dominant_base = 2 ** dominant[0] * 3 ** dominant[1]
    runner_up_base = 2 ** runner_up[0] * 3 ** runner_up[1]
    base_ratio = sp.Rational(runner_up_base, dominant_base)
    require(base_ratio == sp.Rational(3, 4), (degree, base_ratio))

    tail_start = 6
    while dominant_poly.LC() - negative_tariff / tail_start <= 0:
        tail_start += 1
    dominant_floor = dominant_poly.LC() - negative_tariff / tail_start
    tariff = sp.cancel(absolute_tariff / dominant_floor)
    while (
        tariff * tail_start**delta * base_ratio**tail_start >= 1
        or sp.Rational(tail_start + 1, tail_start) ** delta * base_ratio >= 1
    ):
        tail_start += 1
        dominant_floor = dominant_poly.LC() - negative_tariff / tail_start
        tariff = sp.cancel(absolute_tariff / dominant_floor)
    tail_margin = sp.cancel(
        1 - tariff * tail_start**delta * base_ratio**tail_start
    )
    require(tail_margin > 0, (degree, "tail margin"))

    prefix_rows = []
    for width in range(6, tail_start):
        value = sp.Integer(
            expression.as_expr().subs({M: width, U: 2**width, V: 3**width})
        )
        require(value > 0, (degree, width, "finite prefix"))
        prefix_rows.append(f"{width}:{value}")
    return {
        "degree": degree,
        "denominator_power": 2 * degree,
        "expression_terms": len(expression.terms()),
        "group_count": len(grouped),
        "grouped_digest": representation_digest(grouped),
        "dominant_exponents": dominant,
        "dominant_base": dominant_base,
        "dominant_coefficient": str(sp.factor(grouped[dominant])),
        "dominant_degree": dominant_poly.degree(),
        "dominant_negative_tariff": str(negative_tariff),
        "runner_up_exponents": runner_up,
        "runner_up_base": runner_up_base,
        "base_ratio": str(base_ratio),
        "lower_degree": lower_degree,
        "coarse_degree_delta": delta,
        "absolute_lower_tariff": str(absolute_tariff),
        "tail_start": tail_start,
        "dominant_floor": str(dominant_floor),
        "tail_margin": str(tail_margin),
        "finite_prefix_count": len(prefix_rows),
        "finite_prefix_digest": sha256("\n".join(prefix_rows).encode()).hexdigest(),
    }


def as_fraction(value) -> Fraction:
    return Fraction(int(value.p), int(value.q))


def load_structural():
    require(
        lf_sha256(STRUCTURAL_PATH) == EXPECTED_STRUCTURAL_SHA256,
        "THM-2982 structural sidecar hash changed",
    )
    spec = importlib.util.spec_from_file_location(
        "thm2988_direct_structural_check", STRUCTURAL_PATH
    )
    require(spec is not None and spec.loader is not None, "cannot load sidecar")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def direct_equality_digest(determinant, curvatures):
    structural = load_structural()
    rows = []
    for width in range(6, 51):
        _forms, slices = structural.form_top_slices(width)
        transformed = structural.transformed_slices(width, slices)
        substitution = {M: width, U: 2**width, V: 3**width}
        determinant_value = int(determinant.subs(substitution))
        for index, degree in enumerate(DEGREES):
            base = tuple(degree if axis == index else 0 for axis in range(3))
            alpha = Fraction(transformed[index][1].get(base, 0))
            beta = Fraction(transformed[index][2].get(base, 0))
            direct = 2 * beta - alpha * alpha
            symbolic = Fraction(
                int(curvatures[degree].subs(substitution)),
                determinant_value ** (2 * degree),
            )
            require(direct == symbolic, ("direct equality", width, degree))
            rows.append(f"{width}:{degree}:{direct.numerator}/{direct.denominator}")
    return len(rows), sha256("\n".join(rows).encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    determinant, curvatures = symbolic_self_curvatures()
    expected = -2 * (U**2 + 3 * U - 3 * V - 1)
    require(sp.expand(determinant - expected) == 0, "response determinant")
    records = [
        self_curvature_record(degree, curvatures[degree]) for degree in DEGREES
    ]
    direct_count, direct_digest = direct_equality_digest(determinant, curvatures)
    record_lines = [
        json.dumps(record, sort_keys=True, separators=(",", ":"))
        for record in records
    ]
    record_digest = sha256("\n".join(record_lines).encode()).hexdigest()
    transcript = "\n".join(
        [
            "FIRST-GAP SELF-CURVATURE ALL-WIDTH CERTIFICATE",
            f"response_determinant={determinant}",
            *record_lines,
            f"direct_transformed_slice_equalities={direct_count}",
            f"direct_equality_digest={direct_digest}",
            f"record_digest={record_digest}",
            "self_curvatures=STRICTLY_NEGATIVE_FOR_ALL_M_GE_6",
            "all_exact_checks=PASS",
        ]
    ) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
