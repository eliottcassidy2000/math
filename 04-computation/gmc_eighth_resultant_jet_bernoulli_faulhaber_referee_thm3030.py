#!/usr/bin/env python3
"""Independent table-level referee for THM-3030.

This companion does not import the primary verifier and does not claim to
replay the missing interpolation-grid build.  It reads the eight frozen jet
tables directly, checks their immutable hashes and support, and audits every
corner slot over its full visible range.

Its main new observation is exact on j <= 8.  If

    C_j(M) = Q_j(M,0,0),  p_j(M) = (-1)^(j-1) j C_j(M),

then

    p_j(M) - 46 sum_(s=1)^(M-1) s^j - 20 M^j

is constant in M for every frozen j.  Faulhaber's formula therefore explains
all visible nonterminal C-slice slots and identifies

    c_m^C = 46 |B_(2m)|/(2m)!,  m=1,2,3,4.

The all-order continuation is a prediction, not a theorem.  It first gives a
new nonterminal datum at P_10 (m=5); P_9 is only the exceptional terminal slot.
At P_12 it predicts the first Bernoulli-numerator break of the observed
constant-23 numerator, through B_12=-691/2730.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from hashlib import sha256
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

EXPECTED_RECORD_DIGEST = "39db53bc1efd476701186ea5f31493a56f783367293a579a55abb04d26ab7ac2"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


def bernoulli_numbers(limit: int):
    """Return B_0,...,B_limit with the convention B_1=-1/2."""
    values = [Fraction(1)]
    for degree in range(1, limit + 1):
        subtotal = sum(
            Fraction(comb(degree + 1, index)) * values[index]
            for index in range(degree)
        )
        values.append(-subtotal / Fraction(degree + 1))
    return values


BERNOULLI = bernoulli_numbers(16)


def add_polynomials(left, right):
    answer = dict(left)
    for exponent, coefficient in right.items():
        answer[exponent] = answer.get(exponent, Fraction(0)) + coefficient
        if answer[exponent] == 0:
            del answer[exponent]
    return answer


def scale_polynomial(polynomial, scalar):
    scalar = Fraction(scalar)
    return {
        exponent: scalar * coefficient
        for exponent, coefficient in polynomial.items()
        if scalar * coefficient
    }


def shift_polynomial(polynomial):
    """Return p(M+1) from p(M), in the monomial basis."""
    answer = {}
    for exponent, coefficient in polynomial.items():
        for power in range(exponent + 1):
            answer[power] = answer.get(power, Fraction(0)) + coefficient * comb(exponent, power)
    return {power: value for power, value in answer.items() if value}


def faulhaber_to_m_minus_one(power: int):
    """Polynomial for sum_(s=1)^(M-1) s^power."""
    answer = {}
    for index in range(power + 1):
        exponent = power + 1 - index
        coefficient = Fraction(comb(power + 1, index), power + 1) * BERNOULLI[index]
        if coefficient:
            answer[exponent] = coefficient
    return answer


def falling(start: int, length: int) -> int:
    answer = 1
    for offset in range(length):
        answer *= start - offset
    return answer


def load_table(order: int):
    filename, content, expected_hash = TABLES[order]
    path = RESULTS / filename
    require(lf_sha256(path) == expected_hash, f"table hash changed at P_{order}")
    rows = json.loads(path.read_text(encoding="utf-8"))
    require(all(len(row) == 5 and row[4] > 0 for row in rows), f"bad row at P_{order}")
    keys = [(row[0], row[1], row[2]) for row in rows]
    require(len(keys) == len(set(keys)), f"duplicate exponent at P_{order}")
    polynomial = {
        (a, b, c): Fraction(numerator, denominator * content)
        for a, b, c, numerator, denominator in rows
    }
    return rows, polynomial


def slices(polynomial, order: int):
    return {
        "A": {a: value for (a, b, c), value in polynomial.items() if b == 4 * order and c == 0},
        "E": {a: value / Fraction(9) ** order for (a, b, c), value in polynomial.items()
              if b == 0 and c == 2 * order},
        "C": {a: value for (a, b, c), value in polynomial.items() if b == 0 and c == 0},
    }


def shape_record(tables):
    term_counts = []
    for order, (rows, polynomial) in tables.items():
        expected = (2 * order + 1) ** 3 - 3 * (2 * order - 2 - order // 2) - int(order == 3)
        require(len(rows) == expected, f"term count at P_{order}")
        require(all(
            0 <= a <= 2 * order and 0 <= b and 0 <= c and b + 2 * c <= 4 * order
            for a, b, c in polynomial
        ), f"support at P_{order}")
        term_counts.append(len(rows))
    return {"orders": 8, "term_counts": term_counts, "all_table_hashes": "PASS"}


def faulhaber_record(tables):
    constants = []
    recurrence_checks = 0
    for order, (_rows, polynomial) in tables.items():
        corner = slices(polynomial, order)["C"]
        power_sum = scale_polynomial(corner, (-1) ** (order - 1) * order)
        bulk = scale_polynomial(faulhaber_to_m_minus_one(order), 46)
        bulk = add_polynomials(bulk, {order: Fraction(20)})
        residual = add_polynomials(power_sum, scale_polynomial(bulk, -1))
        require(set(residual).issubset({0}), (order, "nonconstant Faulhaber residual", residual))
        constants.append(str(residual.get(0, Fraction(0))))

        difference = add_polynomials(shift_polynomial(power_sum), scale_polynomial(power_sum, -1))
        expected = add_polynomials(
            {order: Fraction(26)},
            scale_polynomial(shift_polynomial({order: Fraction(1)}), 20),
        )
        require(difference == expected, (order, "width recurrence"))
        recurrence_checks += 1
    return {
        "orders": 8,
        "constant_residuals": constants,
        "width_recurrence_checks": recurrence_checks,
        "identity": "p_j=46*sum_1^(M-1)s^j+20*M^j+K_j",
    }


def slot_record(tables):
    visible_counts = {name: 0 for name in "AEC"}
    constants = {name: {} for name in "AEC"}
    even_zero_checks = 0
    for order, (_rows, polynomial) in tables.items():
        current = slices(polynomial, order)
        for slot in range(2, order, 2):
            exponent = order - slot
            for name in "AEC":
                require(current[name].get(exponent, Fraction(0)) == 0,
                        (order, name, slot, "even slot"))
                even_zero_checks += 1
        for m in range(1, order // 2 + 1):
            slot = 2 * m - 1
            exponent = order - slot
            fall = falling(order - 1, 2 * m - 2)
            for name in "AEC":
                coefficient = current[name].get(exponent, Fraction(0))
                value = coefficient / Fraction((-1) ** (order + m) * fall)
                prior = constants[name].setdefault(m, value)
                require(value == prior, (order, name, m, value, prior))
                visible_counts[name] += 1

    bernoulli_constants = {}
    denominators = []
    for m in range(1, 5):
        predicted = Fraction(46) * abs(BERNOULLI[2 * m]) / factorial(2 * m)
        require(constants["C"][m] == predicted, (m, constants["C"][m], predicted))
        bernoulli_constants[str(m)] = str(predicted)
        denominators.append(predicted.denominator)
        d_c = predicted.denominator
        require(constants["A"][m] == Fraction(3 + 44 * 16 ** (m - 1), 4 ** (2 * m - 1) * d_c),
                (m, "A constant law"))
        require(constants["E"][m] == Fraction(4 + 33 * 9 ** (m - 1), 3 ** (2 * m - 1) * d_c),
                (m, "E constant law"))

    actual = slices(tables[8][1], 8)["A"][5]
    correct = Fraction(4949, 3840)
    naive = -correct
    require(actual == correct and actual != naive, "j=8,m=2 sign hostile")
    return {
        "visible_odd_slot_checks": visible_counts,
        "visible_even_zero_checks": even_zero_checks,
        "C_bernoulli_constants": bernoulli_constants,
        "C_denominators": denominators,
        "j8_m2_correct_sign": str(correct),
        "naive_j_plus_k_sign": str(naive),
    }


def prediction_record():
    c5 = Fraction(46) * abs(BERNOULLI[10]) / factorial(10)
    c6 = Fraction(46) * abs(BERNOULLI[12]) / factorial(12)
    require(c5 == Fraction(23, 23950080), "P10 prediction changed")
    require(c6 == Fraction(15893, 653837184000), "P12 prediction changed")
    return {
        "status": "CONJECTURAL_BEYOND_J8",
        "first_new_nonterminal_test": "P_10,m=5,k=9",
        "P9_boundary": "k=9 is terminal and outside the slot law",
        "predicted_c5_C": str(c5),
        "first_Bernoulli_numerator_break": "P_12,m=6,B_12=-691/2730",
        "predicted_c6_C": str(c6),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    tables = {order: load_table(order) for order in TABLES}
    records = [
        {"table_shape": shape_record(tables)},
        {"C_faulhaber": faulhaber_record(tables)},
        {"slot_audit": slot_record(tables)},
        {"future_boundary": prediction_record()},
        {
            "evidence_scope": {
                "independent_table_level_exact": True,
                "interpolation_grid_rebuilt": False,
                "out_of_sample_grid_replayed": False,
            }
        },
    ]
    lines = [json.dumps(record, sort_keys=True, separators=(",", ":")) for record in records]
    record_digest = sha256("\n".join(lines).encode("ascii")).hexdigest()
    require(record_digest == EXPECTED_RECORD_DIGEST, "THM-3030 referee digest changed")
    transcript = "\n".join([
        "THM-3030 INDEPENDENT BERNOULLI-FAULHABER REFEREE",
        "method=direct frozen-table arithmetic;no primary verifier import",
        *lines,
        f"record_digest={record_digest}",
        "all_exact_table_checks=PASS",
    ]) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
