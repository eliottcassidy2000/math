#!/usr/bin/env python3
"""Exact conditional all-width leading-edge certificate for THM-2989."""

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
THM2982 = ROOT / "04-computation/gmc_first_gap_wall_stripped_norm_core_strict_ulc_through_thirty_four_thm2982.py"
EXPECTED_THM2982_SHA256 = (
    "645353f04d2143f91b7b213e7223761c65c72a81c0ebb522184643fdd97ca24b"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


def polynomial_group_digest(grouped, variable) -> str:
    rows = []
    for pair in sorted(grouped, reverse=True):
        polynomial = sp.Poly(grouped[pair], variable)
        rows.append(
            json.dumps(
                {
                    "first_exponent": pair[0],
                    "second_exponent": pair[1],
                    "degree": polynomial.degree(),
                    "coefficients": [str(value) for value in polynomial.all_coeffs()],
                },
                sort_keys=True,
                separators=(",", ":"),
            )
        )
    return sha256("\n".join(rows).encode()).hexdigest()


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
    directions = tuple(offset for offset, count in enumerate(monomial) for _ in range(count))
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
        sp.expand(denominator[2] * signed[0] + denominator[1] * signed[1] + signed[2]),
    )


def line_coefficient(form, base, direction, transfer: int):
    """Coefficient of y^(degree-transfer) z^transfer in F(base*y+direction*z)."""
    answer = sp.Integer(0)
    for exponent, coefficient in form.items():
        choices = [(0,)]
        # Tiny degree (<=4): enumerate how many z factors each variable supplies.
        choices = [()]
        for power in exponent:
            choices = [old + (take,) for old in choices for take in range(power + 1)]
        for taken in choices:
            if sum(taken) != transfer:
                continue
            term = coefficient
            for axis, (power, take) in enumerate(zip(exponent, taken)):
                term *= (
                    math.comb(power, take)
                    * base[axis] ** (power - take)
                    * direction[axis] ** take
                )
            answer += term
    return sp.expand(answer)


def resultant_log_jets_sparse():
    rows = []
    raw_forms = []
    for order in DEGREES:
        rows.append([order**offset - POWERS[order] for offset in (0, 1, 2)])
        triples = {monomial: coefficient_top_three(order, monomial) for monomial in monomials(order)}
        raw_forms.append(
            [{monomial: triples[monomial][drop] for monomial in triples} for drop in range(3)]
        )
    response = sp.Matrix(rows)
    determinant = sp.factor(response.det())
    adjugate = response.adjugate()
    columns = [tuple(adjugate[row, column] for row in range(3)) for column in range(3)]

    # A degree-j form evaluated on an adjugate column must be divided by
    # determinant^j.  Use common denominators D^4 and D^8 and never ask
    # SymPy to cancel the large intermediate rational functions.
    first_numerators = {}
    second_numerators = {}
    self_curvature_numerators = {}
    transfers = {}
    for index, degree in enumerate(DEGREES):
        alpha = line_coefficient(raw_forms[index][1], columns[index], columns[index], 0)
        beta = line_coefficient(raw_forms[index][2], columns[index], columns[index], 0)
        multiplicity = math.prod(DEGREES) // degree
        first_numerators[index] = (degree, multiplicity, alpha)
        second_numerators[index] = (degree, multiplicity, alpha, beta)
        self_curvature_numerators[degree] = sp.expand(2 * beta * determinant**degree - alpha**2)
        for other in range(3):
            if other == index:
                continue
            for transfer in range(1, min(degree, DEGREES[other]) + 1):
                transfers[index, other, transfer] = line_coefficient(
                    raw_forms[index][1], columns[index], columns[other], transfer
                )
    first_common = sp.Integer(0)
    second_common = sp.Integer(0)
    for degree, multiplicity, alpha in first_numerators.values():
        first_common += multiplicity * alpha * determinant ** (4 - degree)
    for degree, multiplicity, alpha, beta in second_numerators.values():
        second_common += multiplicity * (
            2 * beta * determinant ** (8 - degree)
            - alpha**2 * determinant ** (8 - 2 * degree)
        )
    for left in range(3):
        for right in range(left + 1, 3):
            remaining = 3 - left - right
            for transfer in range(1, min(DEGREES[left], DEGREES[right]) + 1):
                second_common -= (
                    2
                    * DEGREES[remaining]
                    * transfer
                    * transfers[left, right, transfer]
                    * transfers[right, left, transfer]
                    * determinant ** (8 - DEGREES[left] - DEGREES[right])
                )
    return (
        determinant,
        sp.expand(first_common),
        sp.expand(second_common),
        self_curvature_numerators,
    )


def sum1(upper):
    return sp.Rational(1, 2) * upper * (upper + 1)


def sum2(upper):
    return sp.Rational(1, 6) * upper * (upper + 1) * (2 * upper + 1)


def interval_sum(power: int, low, high):
    if power == 0:
        return high - low + 1
    function = sum1 if power == 1 else sum2
    return function(high) - function(low - 1)


def wall_stat_residue(residue: int, power: int):
    T = sp.symbols("T", integer=True, nonnegative=True)
    m = 30 * T + residue
    a = 10 * T + residue // 3
    b = 15 * T + residue // 2
    c = 20 * T + (2 * residue) // 3
    value = sp.Integer(0)
    value += 6 * (1**power)
    value += 21 * (2**power)
    value += 26 * interval_sum(power, 3, a)
    value += 24 * interval_sum(power, a + 1, b)
    value += 20 * interval_sum(power, b + 1, c)
    value += 19 * interval_sum(power, c + 1, m)
    value += 3 * interval_sum(power, 2, b)
    value += 4 * interval_sum(power, b + 1, m - 1)
    value += 2 * m**power
    value -= interval_sum(power, 3, b)
    value -= m**power
    if residue % 10 == 1:
        quartic = (4 * m + 1) / 5
        value += quartic**power
    if power == 0:
        value += 6
    elif power == 1:
        value += 3
    else:
        value += sp.Rational(3, 2)
    return T, sp.factor(value)


def residue_record(determinant, resultant_u_num, resultant_ell2_num, residue: int, dump: Path | None):
    T, wall_degree = wall_stat_residue(residue, 0)
    _T, wall_u = wall_stat_residue(residue, 1)
    _T, wall_square = wall_stat_residue(residue, 2)
    x, y = sp.symbols("X Y")
    m = 30 * T + residue
    substitution = {M: m, U: 2**residue * x, V: 3**residue * y}
    degree = sp.factor(46 * m - 26 - wall_degree)
    determinant_sub = determinant.subs(substitution)
    u_num = resultant_u_num.subs(substitution) - wall_u * determinant_sub**4
    ell2_num = resultant_ell2_num.subs(substitution) + wall_square * determinant_sub**8
    numerator = sp.Poly(sp.expand(-u_num**2 - degree * ell2_num), x, y, T)
    denominator = determinant_sub**8
    grouped = {}
    for (x_power, y_power, t_power), coefficient in numerator.terms():
        grouped[x_power, y_power] = grouped.get((x_power, y_power), 0) + coefficient * T**t_power
    ordered_bases = sorted(grouped, key=lambda pair: 2 ** pair[0] * 3 ** pair[1], reverse=True)
    dominant = ordered_bases[0]
    runner_up = ordered_bases[1]
    dominant_poly = sp.Poly(grouped[dominant], T)
    dominant_degree = dominant_poly.degree()
    dominant_leading = dominant_poly.LC()
    require(
        dominant_leading > 0,
        (residue, "dominant leading coefficient is not positive"),
    )
    dominant_coefficients = dominant_poly.all_coeffs()
    dominant_negative_tariff = sum(
        (-coefficient for coefficient in dominant_coefficients[1:] if coefficient < 0),
        sp.Integer(0),
    )
    lower_polys = [sp.Poly(grouped[pair], T) for pair in ordered_bases[1:]]
    lower_degree = max(poly.degree() for poly in lower_polys)
    degree_delta = max(0, lower_degree - dominant_degree)
    absolute_tariff = sum(
        (sum(abs(coefficient) for coefficient in poly.all_coeffs()) for poly in lower_polys),
        sp.Integer(0),
    )
    base_ratio = sp.Rational(2 ** runner_up[0] * 3 ** runner_up[1], 2 ** dominant[0] * 3 ** dominant[1])
    minimum_t = 2 if residue < 4 else 1
    dominance_start = minimum_t
    while dominant_leading - dominant_negative_tariff / dominance_start <= 0:
        dominance_start += 1
    tail_start = dominance_start
    dominant_floor = dominant_leading - dominant_negative_tariff / tail_start
    tariff = sp.cancel(absolute_tariff / dominant_floor)
    while tariff * tail_start**degree_delta * base_ratio ** (30 * tail_start) >= 1:
        tail_start += 1
        dominant_floor = dominant_leading - dominant_negative_tariff / tail_start
        tariff = sp.cancel(absolute_tariff / dominant_floor)
    # The tail majorant decreases from tail_start onward.  The factor
    # ((T+1)/T)^delta is at most 2^delta and q^30 is already tiny.
    step_ratio_bound = sp.cancel(2**degree_delta * base_ratio**30)
    require(
        step_ratio_bound < 1,
        (residue, "coarse tail is not monotone"),
    )
    tail_ratio = sp.cancel(
        tariff * tail_start**degree_delta * base_ratio ** (30 * tail_start)
    )
    tail_margin = sp.cancel(1 - tail_ratio)
    require(tail_margin > 0, (residue, "tail margin"))
    finite_values = []
    finite_rows = []
    for value in range(minimum_t, tail_start):
        evaluated = numerator.as_expr().subs({T: value, x: 2 ** (30 * value), y: 3 ** (30 * value)})
        require(
            evaluated > 0,
            (residue, value, "finite edge numerator is not positive"),
        )
        finite_values.append(value)
        finite_rows.append(f"{value}:{evaluated}")
    if dump is not None:
        dump.write_text(
            "numerator=" + str(numerator.as_expr()) + "\n" + "denominator=" + str(sp.factor(denominator)) + "\n",
            encoding="utf-8",
            newline="\n",
        )
    return {
        "residue": residue,
        "degree": str(degree),
        "numerator_terms": len(numerator.terms()),
        "numerator_group_count": len(grouped),
        "grouped_digest": polynomial_group_digest(grouped, T),
        "numerator_total_degree": numerator.total_degree(),
        "numerator_ops": sp.count_ops(numerator.as_expr()),
        "denominator_ops": sp.count_ops(denominator),
        "denominator": str(sp.factor(denominator)),
        "dominant_exponents": dominant,
        "dominant_base": 2 ** dominant[0] * 3 ** dominant[1],
        "dominant_coefficient": str(sp.factor(grouped[dominant])),
        "dominant_degree": dominant_degree,
        "dominant_leading_coefficient": str(dominant_leading),
        "dominant_negative_tariff": str(dominant_negative_tariff),
        "dominant_floor": str(dominant_floor),
        "dominance_start": dominance_start,
        "runner_up_exponents": runner_up,
        "runner_up_base": 2 ** runner_up[0] * 3 ** runner_up[1],
        "base_ratio": str(base_ratio),
        "absolute_lower_tariff": str(absolute_tariff),
        "lower_degree": lower_degree,
        "coarse_tail_degree_delta": degree_delta,
        "coarse_tail_start": tail_start,
        "tail_ratio": str(tail_ratio),
        "tail_margin": str(tail_margin),
        "tail_step_ratio_bound": str(step_ratio_bound),
        "minimum_T": minimum_t,
        "minimum_width": 30 * minimum_t + residue,
        "finite_prefix_count": len(finite_values),
        "finite_prefix_digest": sha256("\n".join(finite_rows).encode()).hexdigest(),
        "coarse_tail_ratio_num_bits": int(sp.numer(tariff * tail_start**degree_delta * base_ratio ** (30 * tail_start))).bit_length(),
        "coarse_tail_ratio_den_bits": int(sp.denom(tariff * tail_start**degree_delta * base_ratio ** (30 * tail_start))).bit_length(),
        "wall_degree": str(wall_degree),
        "wall_u": str(wall_u),
        "wall_square": str(wall_square),
    }


def as_fraction(value) -> Fraction:
    value = sp.Rational(value)
    return Fraction(int(value.p), int(value.q))


def validate_m34(determinant, resultant_u_num, resultant_ell2_num) -> Fraction:
    require(
        lf_sha256(THM2982) == EXPECTED_THM2982_SHA256,
        "THM-2982 companion hash changed",
    )
    spec = importlib.util.spec_from_file_location("thm2982_edge_expected", THM2982)
    require(spec is not None and spec.loader is not None, "cannot load THM-2982 companion")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    substitution = {M: 34, U: 2**34, V: 3**34}
    determinant_value = int(determinant.subs(substitution))
    resultant_u = Fraction(int(resultant_u_num.subs(substitution)), determinant_value**4)
    resultant_ell2 = Fraction(int(resultant_ell2_num.subs(substitution)), determinant_value**8)
    T, wall_degree = wall_stat_residue(4, 0)
    _T, wall_u = wall_stat_residue(4, 1)
    _T, wall_square = wall_stat_residue(4, 2)
    wall_degree = int(wall_degree.subs(T, 1))
    wall_u = as_fraction(wall_u.subs(T, 1))
    wall_square = as_fraction(wall_square.subs(T, 1))
    degree = 46 * 34 - 26 - wall_degree
    core_u = resultant_u - wall_u
    core_ell2 = resultant_ell2 + wall_square
    ratio = Fraction(degree - 1, degree) * core_u * core_u / (core_ell2 + core_u * core_u)
    require(ratio == module.EXPECTED_M34_RATIO, "sparse M34 edge ratio mismatch")
    return ratio


def self_curvature_record(determinant, degree: int, curvature_numerator):
    """Prove -curvature_numerator positive by one dominant exponential."""
    expression = sp.Poly(sp.expand(-curvature_numerator), U, V, M)
    grouped = {}
    for (u_power, v_power, m_power), coefficient in expression.terms():
        grouped[u_power, v_power] = grouped.get((u_power, v_power), 0) + coefficient * M**m_power
    ordered = sorted(grouped, key=lambda pair: 2 ** pair[0] * 3 ** pair[1], reverse=True)
    dominant = ordered[0]
    runner_up = ordered[1]
    dominant_poly = sp.Poly(grouped[dominant], M)
    if dominant_poly.LC() <= 0:
        raise RuntimeError((degree, "negative dominant self-curvature coefficient"))
    negative_tariff = sum(
        (-coefficient for coefficient in dominant_poly.all_coeffs()[1:] if coefficient < 0),
        sp.Integer(0),
    )
    lower_polys = [sp.Poly(grouped[pair], M) for pair in ordered[1:]]
    lower_degree = max(poly.degree() for poly in lower_polys)
    delta = max(0, lower_degree - dominant_poly.degree())
    absolute_tariff = sum(
        (sum(abs(coefficient) for coefficient in poly.all_coeffs()) for poly in lower_polys),
        sp.Integer(0),
    )
    base_ratio = sp.Rational(2 ** runner_up[0] * 3 ** runner_up[1], 2 ** dominant[0] * 3 ** dominant[1])
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
    for width in range(6, tail_start):
        value = expression.as_expr().subs({M: width, U: 2**width, V: 3**width})
        if value <= 0:
            raise RuntimeError((degree, width, "self-curvature finite prefix failed"))
    return {
        "degree": degree,
        "denominator_power": 2 * degree,
        "expression_terms": len(expression.terms()),
        "dominant_exponents": dominant,
        "dominant_base": 2 ** dominant[0] * 3 ** dominant[1],
        "dominant_coefficient": str(sp.factor(grouped[dominant])),
        "runner_up_exponents": runner_up,
        "runner_up_base": 2 ** runner_up[0] * 3 ** runner_up[1],
        "tail_start": tail_start,
        "finite_prefix_count": tail_start - 6,
        "coarse_degree_delta": delta,
        "dominant_floor": str(dominant_floor),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    determinant, resultant_u_num, resultant_ell2_num, _self_numerators = (
        resultant_log_jets_sparse()
    )
    expected = -2 * (3 * U - 3 * V + U**2 - 1)
    require(sp.expand(determinant - expected) == 0, "response determinant changed")
    m34_ratio = validate_m34(determinant, resultant_u_num, resultant_ell2_num)
    records = [
        residue_record(
            determinant,
            resultant_u_num,
            resultant_ell2_num,
            residue,
            None,
        )
        for residue in range(30)
    ]
    record_lines = [
        json.dumps(record, sort_keys=True, separators=(",", ":"))
        for record in records
    ]
    record_digest = sha256("\n".join(record_lines).encode()).hexdigest()
    wall_lines = [
        json.dumps(
            {
                "residue": record["residue"],
                "degree": record["wall_degree"],
                "u": record["wall_u"],
                "square": record["wall_square"],
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        for record in records
    ]
    wall_digest = sha256("\n".join(wall_lines).encode()).hexdigest()
    transcript = "\n".join(
        [
            "FIRST-GAP WALL-STRIPPED ALL-WIDTH LEADING-EDGE CERTIFICATE",
            f"response_determinant={determinant}",
            "determinant_nonvanishing=PROVED_FOR_ALL_M_GE_6",
            (
                "M34_sparse_edge_ratio="
                f"{m34_ratio.numerator}/{m34_ratio.denominator};EXACT_MATCH"
            ),
            "substitution=X=2^(30T),Y=3^(30T);residues=0..29",
            *record_lines,
            f"wall_invoice_digest={wall_digest}",
            f"record_digest={record_digest}",
            "encoded_continuation_scope=M_GE_34",
            "leading_edge_numerator=-u_core^2-degree*ell2_core",
            "encoded_continuation_sign=STRICTLY_POSITIVE",
            "actual_core_scope=M6..32_PROVED;M33_PLUS_REQUIRES_STATED_CONTINUATION",
            "full_ULC_no_return_arbitrary_radial=NOT_ASSERTED",
            "all_exact_checks=PASS",
        ]
    ) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
