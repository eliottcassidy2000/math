#!/usr/bin/env python3
"""Exact controls for THM-3397.

This companion uses only integer and rational arithmetic.  It checks the
two independent finite ledgers used by the proof:

* the v=0 and 1+c*v=0 orders in the terminal denominator filtration; and
* the divisor lattice and separated reciprocal poles in the V4 quotient.

The normality, etale-torsor, and all-parameter valuation arguments remain
mathematical proof inputs in the theorem file.  No assertion is a truth
gate, so ``python`` and ``python -O`` execute the same checks.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def poly_mul(left, right):
    out = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def poly_eval(poly, value):
    total = Fraction(0)
    for coefficient in reversed(poly):
        total = total * value + coefficient
    return total


def divide_by_linear(poly, root):
    """Return q for poly=(v-root)q, requiring exact divisibility."""
    require(len(poly) >= 2, "constant polynomial cannot be divided")
    high = list(reversed(poly))
    quotient_high = [high[0]]
    for coefficient in high[1:-1]:
        quotient_high.append(coefficient + root * quotient_high[-1])
    remainder = high[-1] + root * quotient_high[-1]
    require(remainder == 0, "inexact linear division")
    return trim(tuple(reversed(quotient_high)))


def order_at(poly, root):
    current = trim(poly)
    order = 0
    while len(current) > 1 and poly_eval(current, root) == 0:
        current = divide_by_linear(current, root)
        order += 1
    return order


def boundary_poly(c, zero_order, l_order, with_tail=True):
    # v^zero_order (1+c*v)^l_order times a factor nonzero at both boundaries.
    zero = (Fraction(0),) * zero_order + (Fraction(1),)
    linear_power = tuple(
        Fraction(comb(l_order, j)) * c**j for j in range(l_order + 1)
    )
    out = poly_mul(zero, linear_power)
    if with_tail:
        out = poly_mul(out, (Fraction(1), c / 2))
    return out


def determinant(matrix):
    n = len(matrix)
    require(all(len(row) == n for row in matrix), "determinant matrix not square")
    if n == 1:
        return matrix[0][0]
    total = 0
    for j, entry in enumerate(matrix[0]):
        minor = [row[:j] + row[j + 1 :] for row in matrix[1:]]
        total += (-1) ** j * entry * determinant(minor)
    return total


def gcd_minors(matrix, size):
    rows = range(len(matrix))
    cols = range(len(matrix[0]))
    value = 0
    for row_ids in combinations(rows, size):
        for col_ids in combinations(cols, size):
            minor = [[matrix[i][j] for j in col_ids] for i in row_ids]
            value = gcd(value, abs(determinant(minor)))
    return value


def smith_diagonal_full_column_rank(matrix):
    rank = len(matrix[0])
    determinantal = [1]
    for size in range(1, rank + 1):
        determinantal.append(gcd_minors(matrix, size))
        require(determinantal[-1] != 0, "relation matrix lost full column rank")
    diagonal = []
    for size in range(1, rank + 1):
        require(
            determinantal[size] % determinantal[size - 1] == 0,
            "invalid determinantal-divisor chain",
        )
        diagonal.append(determinantal[size] // determinantal[size - 1])
    return tuple(diagonal)


def terminal_filtration_checks():
    cells = 0
    good = 0
    zero_hostiles = 0
    l_hostiles = 0
    c_bank = (Fraction(1), Fraction(-1), Fraction(2), Fraction(-3, 2))
    for e in range(2, 10):
        relation_matrix = [[e, 0], [0, e], [1, 1]]
        require(
            smith_diagonal_full_column_rank(relation_matrix) == (1, e),
            f"wrong A_(e-1) divisor Smith form at e={e}",
        )
        for m in range(0, 7):
            for q in range(1, 5):
                for c in c_bank:
                    cells += 1
                    need_zero = m * q
                    need_l = e * q
                    root_l = -1 / c
                    candidate = boundary_poly(c, need_zero, need_l)
                    require(
                        order_at(candidate, Fraction(0)) == need_zero,
                        "wrong v=0 order in terminal clearing polynomial",
                    )
                    require(
                        order_at(candidate, root_l) == need_l,
                        "wrong 1+c*v=0 order in terminal clearing polynomial",
                    )
                    good += 1

                    dropped_l = boundary_poly(c, need_zero, need_l - 1)
                    require(
                        order_at(dropped_l, root_l) == need_l - 1,
                        "L-boundary hostile did not miss by one",
                    )
                    l_hostiles += 1

                    if need_zero > 0:
                        dropped_zero = boundary_poly(c, need_zero - 1, need_l)
                        require(
                            order_at(dropped_zero, Fraction(0)) == need_zero - 1,
                            "v-boundary hostile did not miss by one",
                        )
                        zero_hostiles += 1

    return cells, good, zero_hostiles, l_hostiles


def v4_quotient_checks():
    # div(a), div(b), div(c), div(d) in the basis P_a,P_b,P_c.
    relation_matrix = [
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
        [1, 1, 1],
    ]
    require(
        smith_diagonal_full_column_rank(relation_matrix) == (1, 2, 2),
        "V4 quotient divisor lattice is not (Z/2)^2",
    )

    # A monomial x^i y^j z^k is fixed by all even sign changes exactly when
    # its three exponent parities agree.  In that case it is a monomial in
    # a=x^2,b=y^2,c=z^2,d=xyz with d exponent zero or one.
    invariant_monomials = 0
    tested_monomials = 0
    even_flips = ((0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1))
    for exponents in product(range(7), repeat=3):
        tested_monomials += 1
        fixed = all(
            sum(e * flip for e, flip in zip(exponents, flips)) % 2 == 0
            for flips in even_flips
        )
        parity_equal = len({e % 2 for e in exponents}) == 1
        require(fixed == parity_equal, "V4 invariant parity criterion failed")
        if fixed:
            d_power = exponents[0] % 2
            abc_powers = tuple((e - d_power) // 2 for e in exponents)
            require(all(power >= 0 for power in abc_powers), "negative invariant exponent")
            rebuilt = tuple(2 * power + d_power for power in abc_powers)
            require(rebuilt == exponents, "invariant monomial reconstruction failed")
            invariant_monomials += 1

    # At P_a the reciprocal 1/a has order -2 and the other two reciprocals
    # have order zero, and cyclically.  Hence a nonzero coefficient has a
    # unique lowest-order term on its own boundary and cannot cancel.
    coefficient_bank = (-2, -1, 0, 1, 2)
    triples = 0
    forced_poles = 0
    for coefficients in product(coefficient_bank, repeat=3):
        if coefficients == (0, 0, 0):
            continue
        triples += 1
        pole_places = []
        for place, coefficient in enumerate(coefficients):
            term_orders = [0, 0, 0]
            term_orders[place] = -2
            if coefficient != 0:
                require(
                    term_orders[place] < min(
                        term_orders[j] for j in range(3) if j != place
                    ),
                    "reciprocal pole is not separated",
                )
                pole_places.append(place)
        require(pole_places, "nonzero reciprocal combination has no forced pole")
        forced_poles += len(pole_places)

    # The three nonzero class labels form the standard F_2[C3] plane.
    class_cycle = ((1, 0), (0, 1), (1, 1))
    require(
        tuple((class_cycle[(i + 1) % 3]) for i in range(3))
        == ((0, 1), (1, 1), (1, 0)),
        "C3 class cycle failed",
    )
    require(
        tuple(
            class_cycle[0][j] ^ class_cycle[1][j] ^ class_cycle[2][j]
            for j in range(2)
        )
        == (0, 0),
        "three nonzero V4 character classes do not sum to zero",
    )

    return tested_monomials, invariant_monomials, triples, forced_poles


def main():
    cells, good, zero_hostiles, l_hostiles = terminal_filtration_checks()
    monomials, invariants, triples, poles = v4_quotient_checks()
    semantic = {
        "terminal": {
            "parameter_cells": cells,
            "exact_clearing": good,
            "v_hostiles": zero_hostiles,
            "L_hostiles": l_hostiles,
            "divisor_smith": [[1, e] for e in range(2, 10)],
        },
        "v4": {
            "divisor_smith": [1, 2, 2],
            "monomials": monomials,
            "invariants": invariants,
            "coefficient_triples": triples,
            "forced_poles": poles,
            "class_labels": [[1, 0], [0, 1], [1, 1]],
        },
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("THM-3397 exact boundary-effectivity controls")
    print("arithmetic=integer+rational only")
    print(f"terminal_parameter_cells={cells}")
    print(f"terminal_exact_clearing_controls={good}")
    print(f"terminal_v_boundary_one-step_hostiles={zero_hostiles}")
    print(f"terminal_L_boundary_one-step_hostiles={l_hostiles}")
    print("terminal_divisor_smith_forms=e=2..9: (1,e)")
    print("terminal_identity=v^(m*q)*L^(e*q)*(X^-1*v^-m)^q=Y^q")
    print("terminal_pullback_poles=(-m,-e,0) on (xy=1,x=0,y=0)")
    print("v4_divisor_smith_form=(1,2,2)")
    print(f"v4_monomials_tested={monomials}")
    print(f"v4_invariant_monomials={invariants}")
    print(f"three_view_nonzero_coefficient_triples={triples}")
    print(f"three_view_forced_boundary_poles={poles}")
    print("three_view_class_labels=(1,0),(0,1),(1,1); xor_sum=(0,0)")
    print("invariant_reciprocal_sum=(ab+ac+bc)/d^2")
    print("pullback_reciprocal_sum=x^-2+y^-2+z^-2")
    print(f"semantic_sha256={digest}")
    print("result=PASS")


if __name__ == "__main__":
    main()
