#!/usr/bin/env python3
"""Exact controls for THM-3260's ternary-Cantor residual theorem.

The proof in THM-3260 is symbolic.  This companion independently checks the
factorial-unit contraction, the two Cartier recurrences, the thickened
Cantor-window valuation lemma on a declared finite hostile universe, and the
full (not merely terminal-Bessel) moment polynomials at the first mixed
resonances.  It also records the minimal hostile to the broader digit rule.
"""

from fractions import Fraction
from itertools import product
from math import factorial


PRIME = 3
CARTIER_LENGTH = 8
DIRECT_GCD_LENGTH = 7
TERM_LENGTH = 6
EXACT_D_VALUES = (21, 57, 75, 165, 183, 219, 237)


def vp(integer: int, prime: int = PRIME) -> int:
    if integer == 0:
        raise ValueError("v_p(0) is not used in this script")
    answer = 0
    while integer % prime == 0:
        answer += 1
        integer //= prime
    return answer


def unit(integer: int, prime: int = PRIME) -> int:
    while integer % prime == 0:
        integer //= prime
    return integer % prime


def factorial_valuation(index: int) -> int:
    answer = 0
    while index:
        index //= PRIME
        answer += index
    return answer


def factorial_unit(index: int) -> int:
    """The unit part of index! modulo 3, computed without index!."""
    answer = 1
    for integer in range(1, index + 1):
        answer = answer * unit(integer) % PRIME
    return answer


def carry_potential(state) -> int:
    """Potential P(r,d,g) in Lemma 2.1's 12-state carry certificate."""
    remainder, double_carry, add_carry = state
    if remainder == 0:
        return 0
    return remainder - double_carry - add_carry


def carry_transition_certificate():
    """Check every local transition in the unbounded carry induction.

    A state (r,d,g) records the deficit r in A+B+C+r=M, the incoming
    carry d in B+B, and the incoming carry g in A+(B+B).  For a Cantor
    output digit e in {0,2}, the low input digits (a,b,c) determine the
    next state.  The checked inequality is precisely the induction step

        b-c-2(d'+g') + P(r',d',g') <= P(r,d,g).
    """
    maxima = []
    valid_transitions = 0
    passed = True
    for state in product(range(3), range(2), range(2)):
        remainder, double_carry, add_carry = state
        values = []
        for cantor_digit in (0, 2):
            for a_digit, b_digit, c_digit in product(range(3), repeat=3):
                numerator = (
                    a_digit + b_digit + c_digit + remainder - cantor_digit
                )
                if numerator < 0 or numerator % 3:
                    continue
                next_remainder = numerator // 3
                if next_remainder > 2:
                    continue
                doubled_digit = (2 * b_digit + double_carry) % 3
                next_double = (2 * b_digit + double_carry) // 3
                next_add = (a_digit + doubled_digit + add_carry) // 3
                next_state = (next_remainder, next_double, next_add)
                value = (
                    b_digit
                    - c_digit
                    - 2 * (next_double + next_add)
                    + carry_potential(next_state)
                )
                values.append(value)
                valid_transitions += 1
                passed &= value <= carry_potential(state)
        maxima.append((state, carry_potential(state), max(values), len(values)))
    return valid_transitions, passed, maxima


def ternary(index: int) -> str:
    if index == 0:
        return "0"
    digits = []
    while index:
        digits.append(str(index % PRIME))
        index //= PRIME
    return "".join(reversed(digits))


def is_cantor(index: int) -> bool:
    return set(ternary(index)) <= {"0", "2"}


def bessel_vu(n: int):
    """Valuation/unit arrays for C(n,j)=(n+j)!/(j!(n-j)!)."""
    valuations = [0]
    units = [1]
    valuation = 0
    residue = 1
    for j in range(1, n + 1):
        valuation += vp(n + j) + vp(n - j + 1) - vp(j)
        residue *= unit(n + j) * unit(n - j + 1) * unit(j)
        residue %= PRIME
        valuations.append(valuation)
        units.append(residue)
    return valuations, units


def trim(polynomial):
    answer = [coefficient % PRIME for coefficient in polynomial]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def add(left, right):
    size = max(len(left), len(right))
    return trim(
        [
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(size)
        ]
    )


def scale(polynomial, scalar: int):
    return trim([scalar * coefficient for coefficient in polynomial])


def shift(polynomial, amount: int = 1):
    return [0] * amount + polynomial


def frobenius_cube(polynomial):
    answer = [0] * (3 * (len(polynomial) - 1) + 1)
    for index, coefficient in enumerate(polynomial):
        answer[3 * index] = coefficient
    return trim(answer)


def polynomial_remainder(dividend, divisor):
    answer = trim(dividend)
    divisor = trim(divisor)
    while len(answer) >= len(divisor) and answer != [0]:
        offset = len(answer) - len(divisor)
        multiplier = answer[-1] * pow(divisor[-1], -1, PRIME) % PRIME
        for index, coefficient in enumerate(divisor):
            answer[index + offset] -= multiplier * coefficient
        answer = trim(answer)
    return answer


def polynomial_gcd(left, right):
    left = trim(left)
    right = trim(right)
    while right != [0]:
        left, right = right, polynomial_remainder(left, right)
    inverse = pow(left[-1], -1, PRIME)
    return scale(left, inverse)


def lower_hull_from_valuations(valuations):
    hull = []
    for point in enumerate(valuations):
        while len(hull) >= 2:
            old_slope = Fraction(
                hull[-1][1] - hull[-2][1], hull[-1][0] - hull[-2][0]
            )
            new_slope = Fraction(
                point[1] - hull[-1][1], point[0] - hull[-1][0]
            )
            if old_slope >= new_slope:
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def residual_edges_from_vu(valuations, units):
    hull = lower_hull_from_valuations(valuations)
    edges = {}
    for left, right in zip(hull, hull[1:]):
        i0, y0 = left
        i1, y1 = right
        slope = Fraction(y1 - y0, i1 - i0)
        step = slope.denominator
        residual = []
        for index in range(i0, i1 + 1, step):
            on_edge = valuations[index] == y0 + slope * (index - i0)
            residual.append(units[index] if on_edge else 0)
        edges[slope] = trim(residual)
    return hull, edges


def cantor_sectors(n: int):
    """Return P_n,Q_n in the normalization of THM-3260."""
    valuations, units = bessel_vu(n)
    p_polynomial = [
        units[2 * k] if valuations[2 * k] == k else 0
        for k in range(n // 2 + 1)
    ]
    if n == 0:
        return trim(p_polynomial), [0]
    valuations, units = bessel_vu(n - 1)
    q_polynomial = [
        units[2 * k + 1] if valuations[2 * k + 1] == k else 0
        for k in range(n // 2)
    ]
    return trim(p_polynomial), trim(q_polynomial)


def predicted_cartier(n: int):
    if n % 3 == 0:
        m = n // 3
        p_small, q_small = cantor_sectors(m)
        q_cube = frobenius_cube(q_small)
        return (
            add(frobenius_cube(p_small), scale(shift(q_cube), -1)),
            scale(shift(q_cube), -1),
        )
    m = (n - 2) // 3
    p_small, q_small = cantor_sectors(m)
    if m % 3 == 0:
        r_small = add(p_small, scale(q_small, -1))
    else:
        r_small = scale(add(p_small, q_small), -1)
    return (
        add(frobenius_cube(p_small), shift(frobenius_cube(r_small))),
        scale(frobenius_cube(p_small), -1),
    )


def bessel_coefficient(n: int, j: int) -> int:
    return factorial(n + j) // (factorial(j) * factorial(n - j))


def moment_coefficients(n: int, d: int):
    """Coefficients after v=Du, using the k-sum in THM-3260."""
    answer = []
    n_factorial = factorial(n)
    for j in range(n + 1):
        coefficient = 0
        for k in range(n - j + 1):
            coefficient += (
                (-1) ** (n - j - k)
                * n_factorial
                * d**k
                * bessel_coefficient(n - k, j)
                // factorial(k)
            )
        answer.append(coefficient)
    return answer


def integer_vu(coefficients):
    return [vp(value) for value in coefficients], [unit(value) for value in coefficients]


def termwise_cantor_control(max_length: int):
    checked = 0
    strict = True
    for n in range(2, 3**max_length, 2):
        if not is_cantor(n) or n % 3 != 2:
            continue
        d = n + 1
        exponent = vp(d)
        for delta in (1, 2):
            row = d - delta
            for j in range(row + 1):
                target = Fraction(j, 2) if delta == 1 else (
                    Fraction(0) if j == 0 else Fraction(j - 1, 2)
                )
                for k in range(row - j + 1):
                    value = (
                        exponent * k
                        - factorial_valuation(k)
                        + factorial_valuation(row - k + j)
                        - factorial_valuation(j)
                        - factorial_valuation(row - k - j)
                    )
                    strict &= value >= target
                    strict &= k == 0 or value > target
                    checked += 1
    return checked, strict


def exact_moment_controls():
    exact = True
    rows = []
    for d in EXACT_D_VALUES:
        row_data = []
        for n in (d - 2, d - 1):
            actual = moment_coefficients(n, d)
            actual_v, actual_u = integer_vu(actual)
            actual_hull, actual_edges = residual_edges_from_vu(actual_v, actual_u)

            terminal = [
                (-1) ** (n - j) * factorial(n) * bessel_coefficient(n, j)
                for j in range(n + 1)
            ]
            terminal_v, terminal_u = integer_vu(terminal)
            terminal_hull, terminal_edges = residual_edges_from_vu(
                terminal_v, terminal_u
            )
            exact &= actual_hull == terminal_hull
            exact &= actual_edges == terminal_edges
            row_data.append(actual_hull)
        rows.append((d, row_data))
    return exact, rows


def hostile_d201():
    d = 201
    polynomials = []
    hulls = []
    for n in (d - 2, d - 1):
        coefficients = moment_coefficients(n, d)
        valuations, units = integer_vu(coefficients)
        hull, edges = residual_edges_from_vu(valuations, units)
        hulls.append(hull)
        polynomials.append(edges[Fraction(1, 2)])
    left, right = polynomials
    return hulls, left, right, polynomial_gcd(left, right)


def sparse(polynomial):
    return [(index, coefficient) for index, coefficient in enumerate(polynomial) if coefficient]


def main():
    transitions, carry_passed, carry_maxima = carry_transition_certificate()
    print(
        f"carry_automaton_states=12 valid_transitions={transitions} "
        f"local_potential_certificate={carry_passed}"
    )
    print(f"carry_transition_maxima={carry_maxima}")

    contraction = True
    for m in range(512):
        for remainder in range(3):
            left = factorial_unit(3 * m + remainder)
            tail = -1 if remainder == 2 else 1
            right = (-1) ** m * factorial_unit(m) * tail % PRIME
            contraction &= left == right
    print(f"factorial_unit_contraction_n_lt_1536={contraction}")

    cantor_values = [
        n for n in range(3**CARTIER_LENGTH) if n % 2 == 0 and is_cantor(n)
    ]
    bessel_integral = True
    cartier = True
    for n in cantor_values:
        p_polynomial, q_polynomial = cantor_sectors(n)
        p_values, _ = bessel_vu(n)
        bessel_integral &= all(value >= Fraction(j, 2) for j, value in enumerate(p_values))
        if n:
            q_values, _ = bessel_vu(n - 1)
            bessel_integral &= all(
                value >= (Fraction(0) if j == 0 else Fraction(j - 1, 2))
                for j, value in enumerate(q_values)
            )
            cartier &= (p_polynomial, q_polynomial) == predicted_cartier(n)
    print(
        f"cantor_bessel_universe=all_{len(cantor_values)}_words_length_le_"
        f"{CARTIER_LENGTH} bessel_face_integrality={bessel_integral} "
        f"cartier_recursions={cartier}"
    )

    chamber_values = [
        n
        for n in cantor_values
        if n % 3 == 2 and len(ternary(n)) <= DIRECT_GCD_LENGTH
    ]
    direct_gcds = True
    simple_faces = True
    for n in chamber_values:
        p_polynomial, q_polynomial = cantor_sectors(n)
        direct_gcds &= polynomial_gcd(p_polynomial, q_polynomial) == [1]
        p_values, _ = bessel_vu(n)
        q_values, _ = bessel_vu(n - 1)
        simple_faces &= lower_hull_from_valuations(p_values) == [
            (0, 0),
            (n, n // 2),
        ]
        expected_q_hull = [(0, 0), (1, 0)]
        if n > 2:
            expected_q_hull.append((n - 1, (n - 2) // 2))
        simple_faces &= lower_hull_from_valuations(q_values) == expected_q_hull
    print(
        f"direct_chamber_words_length_le_{DIRECT_GCD_LENGTH}={len(chamber_values)} "
        f"simple_faces={simple_faces} residual_gcds_one={direct_gcds}"
    )

    term_count, termwise = termwise_cantor_control(TERM_LENGTH)
    print(
        f"thickened_cantor_terms_length_le_{TERM_LENGTH}={term_count} "
        f"all_on_or_above_and_k_positive_strict={termwise}"
    )

    exact, exact_rows = exact_moment_controls()
    print(
        f"full_moment_equals_terminal_faces_D={','.join(map(str, EXACT_D_VALUES))} "
        f"result={exact}"
    )
    print(f"full_moment_hulls={exact_rows}")

    hulls, left, right, common = hostile_d201()
    print(
        f"hostile_D=201 ternary={ternary(201)} unit_ternary={ternary(67)} "
        f"hulls={hulls} left_half={sparse(left)} right_half={sparse(right)} "
        f"gcd_degree={len(common)-1} gcd_sparse={sparse(common)}"
    )

    print(
        "broader_rule_hostiles="
        "D525:unit20111_3,D1731:unit210101_3 "
        "reason=more_than_one_digit_1"
    )


if __name__ == "__main__":
    main()
