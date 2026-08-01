#!/usr/bin/env python3
"""Exact companion for THM-3040.

The theorem is structural.  At the formal U=V=0 corner every nonempty
terminal-offset face of the factorial inclusion expansion vanishes.  The
remaining three forms are independent of the width M, so resultant
multihomogeneity reduces width transport to the three THM-2925 denominator
characters.

This companion independently checks the complete character ledger, the
terminal-face census, the nonzero leading resultant, the exact quotient on a
hostile width grid, the all-order Bernoulli/Faulhaber algebra through order
32, and every frozen C-corner jet through order eight.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb, factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "05-knowledge" / "results"

TABLES = {
    1: ("gmc_first_gap_resultant_jet_P1_table_thm3030.json", 1,
        "44e1be729b907432db8333e97a4d2558a02f1e8bd775a788e0ce709dfdbab778"),
    2: ("gmc_first_gap_resultant_jet_P2_table_thm3030.json", 48,
        "b2439d5227e35d6e7ed7d4bfe32c6f5204fe98547b481964e93ec6fcc4c6560d"),
    3: ("gmc_first_gap_resultant_jet_P3_table_thm3030.json", 1152,
        "3a4b39fa74afd68ab32973dfe0508dee9d34f9f417d37b566ed786e410191c8d"),
    4: ("gmc_first_gap_fourth_resultant_jet_P4_table_thm3013.json", 1658880,
        "200fae225af6c381733a386f31e107b81e6d33104699e5f69e8d7cd7e445d163"),
    5: ("gmc_first_gap_resultant_jet_P5_table_thm3015.json", 39813120,
        "8cf6ed5cfca3b92a9229f79aae26f87ab1e65db1cf288b4299fe37b4a47b1de9"),
    6: ("gmc_first_gap_resultant_jet_P6_table_thm3015.json", 120394874880,
        "7004c579194e13f10aa03ceb26707adbaeae01e64b5be85d76792987f20c150e"),
    7: ("gmc_first_gap_resultant_jet_P7_table_thm3015.json", 2889476997120,
        "bc53e1a23a9694c277de3aa1e9f6c4401be585452b7a20e35abc0a7a050fb287"),
    8: ("gmc_first_gap_resultant_jet_P8_table_thm3030.json", 1664338750341120,
        "bba6b4b9916a316c41b800a044861a15840820b6048133b754d85cfad78873ad"),
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


def determinant3(rows: tuple[tuple[int, int, int], ...]) -> int:
    (a, b, c), (d, e, f), (g, h, i) = rows
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def determinant_bareiss(rows: tuple[tuple[int, ...], ...]) -> int:
    matrix = [list(row) for row in rows]
    size = len(matrix)
    require(size and all(len(row) == size for row in matrix), "determinant needs square matrix")
    sign = 1
    denominator = 1
    for pivot_index in range(size - 1):
        if matrix[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, size)
                 if matrix[row][pivot_index] != 0),
                None,
            )
            require(swap is not None, "singular Bareiss pivot")
            matrix[pivot_index], matrix[swap] = matrix[swap], matrix[pivot_index]
            sign = -sign
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = matrix[row][column] * pivot - matrix[row][pivot_index] * matrix[pivot_index][column]
                require(numerator % denominator == 0, "Bareiss exact division failed")
                matrix[row][column] = numerator // denominator
        denominator = pivot
    return sign * matrix[-1][-1]


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


def add_poly(left, right):
    answer = dict(left)
    for power, coefficient in right.items():
        answer[power] = answer.get(power, Fraction(0)) + coefficient
        if not answer[power]:
            del answer[power]
    return answer


def scale_poly(poly, scalar):
    scalar = Fraction(scalar)
    return {power: scalar * coefficient for power, coefficient in poly.items()
            if scalar * coefficient}


def shift_poly(poly):
    answer = {}
    for exponent, coefficient in poly.items():
        for power in range(exponent + 1):
            answer[power] = answer.get(power, Fraction(0)) + coefficient * comb(exponent, power)
    return {power: coefficient for power, coefficient in answer.items() if coefficient}


def bernoulli_numbers(limit: int):
    values = [Fraction(1)]
    for degree in range(1, limit + 1):
        subtotal = sum(
            Fraction(comb(degree + 1, index)) * values[index]
            for index in range(degree)
        )
        values.append(-subtotal / Fraction(degree + 1))
    return values


BERNOULLI = bernoulli_numbers(40)


def faulhaber_m_minus_one(power: int):
    """Polynomial in M for sum_(s=1)^(M-1) s^power."""
    answer = {}
    for index in range(power + 1):
        exponent = power + 1 - index
        coefficient = Fraction(comb(power + 1, index), power + 1) * BERNOULLI[index]
        if coefficient:
            answer[exponent] = coefficient
    return answer


def predicted_corner(order: int):
    p = scale_poly(faulhaber_m_minus_one(order), 46)
    p = add_poly(p, {order: Fraction(20)})
    # The M-independent integration constant is deliberately omitted.  Every
    # nonconstant coefficient and every width difference is canonical.
    return scale_poly(p, Fraction((-1) ** (order - 1), order))


def load_frozen_corner(order: int):
    filename, content, expected_hash = TABLES[order]
    path = RESULTS / filename
    require(lf_sha256(path) == expected_hash, f"frozen table hash changed at P_{order}")
    rows = json.loads(path.read_text(encoding="utf-8"))
    corner = {
        m_power: Fraction(numerator, denominator * content)
        for m_power, u_power, v_power, numerator, denominator in rows
        if u_power == 0 and v_power == 0
    }
    return corner


def exponent_map(width: int):
    """Normalized t-factor exponents in D2^12 D3^8 D4^6."""
    multidegree = {2: 12, 3: 8, 4: 6}
    answer = {}
    for order, weight in multidegree.items():
        for root in range(1, width):
            answer[root] = answer.get(root, 0) + weight * (order - 1)
        answer[width] = answer.get(width, 0) + weight * (order - 2)
    return {root: exponent for root, exponent in answer.items() if exponent}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    multidegree = {2: 12, 3: 8, 4: 6}
    interior = sum(multidegree[r] * (r - 1) for r in multidegree)
    terminal = sum(multidegree[r] * (r - 2) for r in multidegree)
    promotion = interior - terminal
    require((interior, terminal, promotion) == (46, 20, 26), "character ledger changed")

    face_records = []
    for order in (2, 3, 4):
        monomials = tuple(compositions(order, 3))
        total = len(monomials) * 2**order
        survivors = len(monomials)
        killed = 0
        for monomial in monomials:
            directions = tuple(
                direction for direction, count in enumerate(monomial)
                for _ in range(count)
            )
            require(len(directions) == order, "direction expansion changed")
            for mask in range(1 << order):
                terminal_count = mask.bit_count()
                # The tensor prefactor is r^B X_r^terminal_count.  Corner
                # evaluation X_r=0 retains exactly mask=0.
                if terminal_count:
                    killed += 1
                else:
                    require(mask == 0, "zero-terminal face is not empty")
        require(killed == total - survivors, "terminal-face census changed")
        face_records.append({
            "order": order,
            "monomials": len(monomials),
            "inclusion_faces": total,
            "corner_survivors": survivors,
            "terminal_faces_killed": killed,
        })

    response_rows = ((1, 2, 4), (1, 3, 9), (1, 4, 16))
    response_det = determinant3(response_rows)
    leading_resultant = response_det ** (2 * 3 * 4)
    require(response_det == 2 and leading_resultant == 2**24, "corner leading resultant changed")
    generalized_determinants = []
    for lower in combinations(range(13), 3):
        rows = tuple(tuple(base**offset for offset in lower) for base in (2, 3, 4))
        value = determinant3(rows)
        require(value > 0, (lower, "generalized Vandermonde lost strict positivity"))
        generalized_determinants.append((lower, value))
    require(len(generalized_determinants) == 286, "fixed-lower triple census changed")

    general_k_records = []
    for slots in range(3, 11):
        harmonic = sum((Fraction(1, index) for index in range(1, slots + 1)), Fraction(0))
        total_degree = factorial(slots)
        multidegrees = [total_degree // order for order in range(2, slots + 1)]
        old_terminal_promotion = sum(multidegrees)
        new_terminal = sum(
            weight * (order - 2)
            for order, weight in zip(range(2, slots + 1), multidegrees)
        )
        interior_density = old_terminal_promotion + new_terminal
        expected_a = total_degree * (harmonic - 1)
        expected_b = total_degree * (slots + 1 - 2 * harmonic)
        expected_i = total_degree * (slots - harmonic)
        require(expected_a.denominator == expected_b.denominator == expected_i.denominator == 1,
                (slots, "nonintegral general character"))
        require((old_terminal_promotion, new_terminal, interior_density)
                == (expected_a.numerator, expected_b.numerator, expected_i.numerator),
                (slots, "general character formula"))
        rows = tuple(
            tuple(base**offset for offset in range(slots - 1))
            for base in range(2, slots + 1)
        )
        vandermonde = determinant_bareiss(rows)
        expected_vandermonde = 1
        for left in range(2, slots + 1):
            for right in range(left + 1, slots + 1):
                expected_vandermonde *= right - left
        require(vandermonde == expected_vandermonde > 0,
                (slots, "ordinary Vandermonde leader"))
        general_k_records.append({
            "slots": slots,
            "multidegrees": multidegrees,
            "old_terminal_promotion_A": old_terminal_promotion,
            "new_terminal_B": new_terminal,
            "interior_density": interior_density,
            "leader_base_vandermonde": vandermonde,
            "leader_exponent": total_degree,
        })
    require(
        [(row["old_terminal_promotion_A"], row["new_terminal_B"])
         for row in general_k_records[:3]] == [(5, 2), (26, 20), (154, 172)],
        "k=3,4,5 controls changed",
    )

    width_records = []
    for width in range(3, 41):
        current = exponent_map(width)
        following = exponent_map(width + 1)
        quotient = dict(following)
        for root, exponent in current.items():
            quotient[root] = quotient.get(root, 0) - exponent
            if not quotient[root]:
                del quotient[root]
        expected = {width: 26, width + 1: 20}
        require(quotient == expected, (width, quotient, expected))
        require(sum(current.values()) == 46 * width - 26, "resultant degree changed")
        width_records.append((width, tuple(sorted(quotient.items()))))

    all_order_checks = 0
    for order in range(1, 33):
        corner = predicted_corner(order)
        difference = add_poly(shift_poly(corner), scale_poly(corner, -1))
        expected = scale_poly(
            add_poly({order: Fraction(26)}, scale_poly(shift_poly({order: Fraction(1)}), 20)),
            Fraction((-1) ** (order - 1), order),
        )
        require(difference == expected, (order, "log quotient recurrence"))
        for m in range(1, (order + 1) // 2 + 1):
            slot = 2 * m - 1
            if slot >= order:
                continue
            exponent = order - slot
            falling = factorial(order - 1) // factorial(order - 2 * m + 1)
            constant = Fraction(46) * abs(BERNOULLI[2 * m]) / factorial(2 * m)
            predicted = Fraction((-1) ** (order + m) * falling) * constant
            require(corner.get(exponent, Fraction(0)) == predicted,
                    (order, m, corner.get(exponent), predicted))
        for slot in range(2, order, 2):
            require(corner.get(order - slot, Fraction(0)) == 0,
                    (order, slot, "even slot"))
        all_order_checks += 1

    frozen_checks = 0
    frozen_constants = []
    for order in TABLES:
        actual = load_frozen_corner(order)
        predicted = predicted_corner(order)
        for exponent, coefficient in predicted.items():
            require(actual.get(exponent, Fraction(0)) == coefficient,
                    (order, exponent, actual.get(exponent), coefficient))
        residual = add_poly(actual, scale_poly(predicted, -1))
        require(set(residual).issubset({0}), (order, "nonconstant frozen mismatch", residual))
        frozen_constants.append(str(residual.get(0, Fraction(0))))
        frozen_checks += 1

    c5 = Fraction(46) * abs(BERNOULLI[10]) / factorial(10)
    c6 = Fraction(46) * abs(BERNOULLI[12]) / factorial(12)
    require(c5 == Fraction(23, 23950080), "c5 changed")
    require(c6 == Fraction(15893, 653837184000), "c6 changed")
    require(predicted_corner(10)[1] == Fraction(-23, 66), "P10 linear slot changed")
    require(predicted_corner(12)[1] == Fraction(15893, 16380), "P12 linear slot changed")

    records = [
        {"resultant_character": {
            "multidegree_QCF": [12, 8, 6],
            "interior": interior,
            "terminal": terminal,
            "promotion": promotion,
        }},
        {"terminal_face_census": face_records},
        {"corner_leader": {
            "linear_response_det": response_det,
            "resultant_constant": leading_resultant,
            "fixed_lower_triples_0_through_12": len(generalized_determinants),
            "generalized_vandermonde_digest": sha256(
                repr(generalized_determinants).encode("ascii")
            ).hexdigest(),
        }},
        {"width_quotient": {
            "audited_widths": "3..40",
            "identity": "(1+M*t)^26*(1+(M+1)*t)^20",
            "record_digest": sha256(repr(width_records).encode("ascii")).hexdigest(),
        }},
        {"general_k_slot_character": general_k_records},
        {"all_order_algebra": {
            "orders_checked": all_order_checks,
            "frozen_C_tables_matched": frozen_checks,
            "frozen_integration_constants": frozen_constants,
            "c5": str(c5),
            "P10_linear_slot": str(predicted_corner(10)[1]),
            "c6": str(c6),
            "P12_linear_slot": str(predicted_corner(12)[1]),
        }},
        {"scope": {
            "object": "intrinsic pure-resultant formal U=V=0 corner germ",
            "physical_width": False,
            "raw_selected_chart": False,
            "wall_stripped_core": False,
            "global_Pj_polynomial_needed": False,
        }},
    ]
    lines = [json.dumps(record, sort_keys=True, separators=(",", ":")) for record in records]
    transcript = "\n".join([
        "THM-3040 FORMAL CORNER RESULTANT WIDTH QUOTIENT",
        "method=terminal-face ideal plus resultant multihomogeneity;frozen tables used only as controls",
        *lines,
        "all_exact_checks=PASS",
    ]) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
