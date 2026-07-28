#!/usr/bin/env python3
"""Independent exact audit for THM-2808.

This script deliberately does not use SymPy.  It has two independent parts:

1. a small exact ``Fraction`` implementation of Q[lambda], used to reduce
   D=x^a(x-1)^b(x-lambda)^c modulo the pole-stripped critical quadratic and
   check the collision quotient, boundary values, squarefreeness, and the
   double-critical-point division; and
2. an exhaustive permutation enumeration of products of two disjoint
   transpositions, modulo the full-cycle centralizer.  The three cycles are
   explicitly labelled, so repeated pole parts do not silently lose marked
   charts.  The exhaustive atlas is compared as a set with the noncrossing
   chord-gap parametrization.

The finite universe is all ordered positive triples through N=14.  The proof
in THM-2808 is all-degree; this is a hostile finite audit, not its replacement.
"""

import ast
from fractions import Fraction
from itertools import combinations, permutations
from math import comb
from pathlib import Path


ZERO = (Fraction(0),)
ONE = (Fraction(1),)


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    return any(
        isinstance(node, ast.Assert)
        for node in ast.walk(
            ast.parse(Path(path).read_text(encoding="utf-8"))
        )
    )


# ---------------------------------------------------------------------------
# Exact polynomials in lambda, stored coefficient-first.


def trim(poly):
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def poly_add(left, right):
    size = max(len(left), len(right))
    return trim(
        [
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(size)
        ]
    )


def poly_neg(poly):
    return tuple(-coefficient for coefficient in poly)


def poly_sub(left, right):
    return poly_add(left, poly_neg(right))


def poly_scale(poly, scalar):
    return trim([scalar * coefficient for coefficient in poly])


def poly_mul(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            result[left_index + right_index] += (
                left_coefficient * right_coefficient
            )
    return trim(result)


def monomial(coefficient, degree):
    return trim([Fraction(0)] * degree + [Fraction(coefficient)])


def poly_divmod(numerator, denominator):
    numerator_work = list(trim(numerator))
    denominator = trim(denominator)
    require(denominator != ZERO, "polynomial division by zero")
    quotient = [Fraction(0)] * max(
        1, len(numerator_work) - len(denominator) + 1
    )
    while (
        not (
            len(numerator_work) == 1
            and numerator_work[0] == 0
        )
        and len(numerator_work) >= len(denominator)
    ):
        shift = len(numerator_work) - len(denominator)
        coefficient = (
            numerator_work[-1] / denominator[-1]
        )
        quotient[shift] += coefficient
        for index, denominator_coefficient in enumerate(denominator):
            numerator_work[index + shift] -= (
                coefficient * denominator_coefficient
            )
        numerator_work = list(trim(numerator_work))
    return trim(quotient), trim(numerator_work)


def poly_derivative(poly):
    if len(poly) == 1:
        return ZERO
    return trim(
        [
            Fraction(index) * poly[index]
            for index in range(1, len(poly))
        ]
    )


def poly_gcd(left, right):
    left = trim(left)
    right = trim(right)
    while right != ZERO:
        _, remainder = poly_divmod(left, right)
        left, right = right, remainder
    if left == ZERO:
        return ZERO
    return poly_scale(left, Fraction(1, 1) / left[-1])


def poly_evaluate(poly, value):
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


# Polynomials in x whose coefficients are the lambda-polynomials above.


def xpoly_trim(poly):
    result = list(poly)
    while len(result) > 1 and result[-1] == ZERO:
        result.pop()
    return result


def xpoly_add(left, right):
    size = max(len(left), len(right))
    return xpoly_trim(
        [
            poly_add(
                left[index] if index < len(left) else ZERO,
                right[index] if index < len(right) else ZERO,
            )
            for index in range(size)
        ]
    )


def xpoly_neg(poly):
    return [poly_neg(coefficient) for coefficient in poly]


def xpoly_sub(left, right):
    return xpoly_add(left, xpoly_neg(right))


def xpoly_mul(left, right):
    result = [ZERO] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            result[left_index + right_index] = poly_add(
                result[left_index + right_index],
                poly_mul(left_coefficient, right_coefficient),
            )
    return xpoly_trim(result)


def xpoly_derivative(poly):
    if len(poly) == 1:
        return [ZERO]
    return xpoly_trim(
        [
            poly_scale(poly[index], Fraction(index))
            for index in range(1, len(poly))
        ]
    )


def xpoly_divmod_monic(numerator, denominator):
    numerator_work = xpoly_trim(numerator)
    denominator = xpoly_trim(denominator)
    require(denominator[-1] == ONE, "x-polynomial divisor is not monic")
    quotient = [ZERO] * max(
        1, len(numerator_work) - len(denominator) + 1
    )
    while len(numerator_work) >= len(denominator):
        shift = len(numerator_work) - len(denominator)
        coefficient = numerator_work[-1]
        quotient[shift] = poly_add(quotient[shift], coefficient)
        for index, denominator_coefficient in enumerate(denominator):
            numerator_work[index + shift] = poly_sub(
                numerator_work[index + shift],
                poly_mul(coefficient, denominator_coefficient),
            )
        numerator_work = xpoly_trim(numerator_work)
    return xpoly_trim(quotient), xpoly_trim(numerator_work)


def denominator_coefficients(a, b, c):
    """Return D=x^a(x-1)^b(x-lambda)^c coefficient-first in x."""
    degree = a + b + c
    coefficients = [ZERO] * (degree + 1)
    for first_index in range(b + 1):
        for second_index in range(c + 1):
            scalar = Fraction(
                comb(b, first_index)
                * comb(c, second_index)
                * (-1) ** (
                    b - first_index + c - second_index
                )
            )
            x_degree = a + first_index + second_index
            coefficients[x_degree] = poly_add(
                coefficients[x_degree],
                monomial(scalar, c - second_index),
            )
    return xpoly_trim(coefficients)


def critical_remainder(a, b, c):
    """Reduce D modulo E=x^2-sx+p using only its two-term recurrence."""
    degree = a + b + c
    sum_numerator_constant = a + c
    sum_numerator_linear = a + b
    critical_sum = (
        Fraction(sum_numerator_constant, degree),
        Fraction(sum_numerator_linear, degree),
    )
    critical_product = (
        Fraction(0),
        Fraction(a, degree),
    )

    # x^k = linear_coefficients[k] x + constant_coefficients[k] mod E.
    linear_coefficients = [ZERO, ONE]
    constant_coefficients = [ONE, ZERO]
    for _ in range(2, degree + 1):
        previous_linear = linear_coefficients[-1]
        previous_constant = constant_coefficients[-1]
        linear_coefficients.append(
            poly_add(
                poly_mul(previous_linear, critical_sum),
                previous_constant,
            )
        )
        constant_coefficients.append(
            poly_neg(
                poly_mul(previous_linear, critical_product)
            )
        )

    denominator = denominator_coefficients(a, b, c)
    secant_coefficient = ZERO
    critical_value = ZERO
    for x_degree, coefficient in enumerate(denominator):
        secant_coefficient = poly_add(
            secant_coefficient,
            poly_mul(
                coefficient, linear_coefficients[x_degree]
            ),
        )
        critical_value = poly_add(
            critical_value,
            poly_mul(
                coefficient, constant_coefficients[x_degree]
            ),
        )

    critical_sum_numerator = (
        Fraction(sum_numerator_constant),
        Fraction(sum_numerator_linear),
    )
    collision_discriminant = poly_sub(
        poly_mul(
            critical_sum_numerator, critical_sum_numerator
        ),
        monomial(Fraction(4 * degree * a), 1),
    )
    maxwell, collision_remainder = poly_divmod(
        secant_coefficient, collision_discriminant
    )

    return {
        "degree": degree,
        "critical_sum": critical_sum,
        "critical_product": critical_product,
        "collision": collision_discriminant,
        "denominator": denominator,
        "secant": secant_coefficient,
        "value": critical_value,
        "maxwell": maxwell,
        "collision_remainder": collision_remainder,
    }


def check_algebra(a, b, c):
    data = critical_remainder(a, b, c)
    degree = data["degree"]
    maxwell = data["maxwell"]
    collision = data["collision"]

    require(
        data["collision_remainder"] == ZERO,
        "collision factor does not divide the secant coefficient",
    )
    require(
        poly_mul(collision, maxwell) == data["secant"],
        "collision quotient does not reconstruct the secant coefficient",
    )
    require(
        len(data["secant"]) - 1 == degree - 1,
        "raw secant coefficient has the wrong degree",
    )
    require(
        len(maxwell) - 1 == degree - 3,
        "Maxwell quotient has the wrong degree",
    )
    require(
        poly_gcd(collision, poly_derivative(collision)) == ONE,
        "collision discriminant is not squarefree",
    )
    require(
        poly_gcd(maxwell, collision) == ONE,
        "Maxwell and collision factors meet",
    )
    require(
        poly_gcd(maxwell, poly_derivative(maxwell)) == ONE,
        "Maxwell quotient is not squarefree",
    )

    left_part = a + c
    right_part = b + c
    boundary_zero = (
        Fraction((-1) ** b)
        * Fraction(left_part) ** (left_part - 3)
        * Fraction(b) ** b
        / Fraction(degree) ** (degree - 1)
    )
    boundary_one = (
        Fraction((-1) ** (b + c - 1))
        * Fraction(a) ** a
        * Fraction(right_part) ** (right_part - 3)
        / Fraction(degree) ** (degree - 1)
    )
    leading_coefficient = (
        Fraction((-1) ** c)
        * Fraction(a + b) ** (a + b - 3)
        * Fraction(c) ** c
        / Fraction(degree) ** (degree - 1)
    )
    require(
        poly_evaluate(maxwell, Fraction(0)) == boundary_zero,
        "lambda=0 boundary formula",
    )
    require(
        poly_evaluate(maxwell, Fraction(1)) == boundary_one,
        "lambda=1 boundary formula",
    )
    require(
        maxwell[-1] == leading_coefficient,
        "Maxwell leading-coefficient formula",
    )
    require(
        poly_gcd(data["value"], maxwell) == ONE,
        "critical value vanishes at a Maxwell root",
    )

    # At a Maxwell root, D-value is constant modulo E and D-value has a
    # double zero at each E-root.  Verify E^2 division coefficientwise modulo
    # Q(lambda), using a monic x-division rather than a CAS quotient.
    critical_quadratic = [
        data["critical_product"],
        poly_neg(data["critical_sum"]),
        ONE,
    ]
    critical_quadratic_squared = xpoly_mul(
        critical_quadratic, critical_quadratic
    )
    shifted_denominator = list(data["denominator"])
    shifted_denominator[0] = poly_sub(
        shifted_denominator[0], data["value"]
    )
    simple_factor, square_remainder = xpoly_divmod_monic(
        shifted_denominator, critical_quadratic_squared
    )
    require(
        len(simple_factor) - 1 == degree - 4
        and simple_factor[-1] == ONE,
        "remaining simple factor has the wrong degree or leading term",
    )
    for coefficient in square_remainder:
        _, remainder = poly_divmod(coefficient, maxwell)
        require(
            remainder == ZERO,
            "critical quadratic square does not divide modulo Maxwell",
        )

    # Cleared derivative identity for F=(D-v)/D.
    pole_factor = denominator_coefficients(
        a - 1, b - 1, c - 1
    )
    left_side = xpoly_sub(
        xpoly_mul(
            xpoly_derivative(shifted_denominator),
            data["denominator"],
        ),
        xpoly_mul(
            shifted_denominator,
            xpoly_derivative(data["denominator"]),
        ),
    )
    right_side = xpoly_mul(
        [
            poly_scale(
                coefficient,
                Fraction(degree),
            )
            for coefficient in pole_factor
        ],
        critical_quadratic,
    )
    right_side = [
        poly_mul(data["value"], coefficient)
        for coefficient in right_side
    ]
    require(
        xpoly_sub(left_side, right_side) == [ZERO],
        "cleared response derivative identity",
    )

    return data


# ---------------------------------------------------------------------------
# Full-cycle centralizer and marked noncrossing chord audit.


def compose(left, right):
    return tuple(
        left[right[index]] for index in range(len(left))
    )


def involution(degree, chord_pair):
    permutation = list(range(degree))
    for first, second in chord_pair:
        permutation[first] = second
        permutation[second] = first
    return tuple(permutation)


def cycle_sets(permutation):
    unseen = set(range(len(permutation)))
    cycles = []
    while unseen:
        start = min(unseen)
        point = start
        cycle = []
        while point in unseen:
            unseen.remove(point)
            cycle.append(point)
            point = permutation[point]
        cycles.append(frozenset(cycle))
    return tuple(cycles)


def shifted_chord_pair(chord_pair, shift, degree):
    return tuple(
        sorted(
            tuple(
                sorted(
                    (
                        (first + shift) % degree,
                        (second + shift) % degree,
                    )
                )
            )
            for first, second in chord_pair
        )
    )


def canonical_unmarked_chords(chord_pair, degree):
    return min(
        shifted_chord_pair(chord_pair, shift, degree)
        for shift in range(degree)
    )


def canonical_marked_chords(chord_pair, labelled_cycles, degree):
    representatives = []
    for shift in range(degree):
        shifted_pair = shifted_chord_pair(
            chord_pair, shift, degree
        )
        shifted_cycles = tuple(
            tuple(
                sorted(
                    (point + shift) % degree
                    for point in cycle
                )
            )
            for cycle in labelled_cycles
        )
        representatives.append((shifted_pair, shifted_cycles))
    return min(representatives)


def exhaustive_chord_atlas(degree):
    full_cycle = tuple(
        (index + 1) % degree for index in range(degree)
    )
    marked_by_parts = {}
    unmarked_orbits = set()
    raw_noncrossing = set()
    for endpoints in combinations(range(degree), 4):
        first, second, third, fourth = endpoints
        pairings = (
            ((first, second), (third, fourth)),
            ((first, third), (second, fourth)),
            ((first, fourth), (second, third)),
        )
        for chord_pair in pairings:
            traced_cycles = cycle_sets(
                compose(
                    involution(degree, chord_pair), full_cycle
                )
            )
            if len(traced_cycles) != 3:
                continue
            normalized_pair = tuple(
                sorted(tuple(sorted(chord)) for chord in chord_pair)
            )
            raw_noncrossing.add(normalized_pair)
            unmarked_orbits.add(
                canonical_unmarked_chords(
                    normalized_pair, degree
                )
            )
            for cycle_order in permutations(range(3)):
                labelled_cycles = tuple(
                    traced_cycles[index] for index in cycle_order
                )
                parts = tuple(
                    len(cycle) for cycle in labelled_cycles
                )
                state = canonical_marked_chords(
                    normalized_pair, labelled_cycles, degree
                )
                marked_by_parts.setdefault(parts, set()).add(state)

    return {
        "marked": marked_by_parts,
        "unmarked": unmarked_orbits,
        "raw": raw_noncrossing,
    }


def gap_atlas(parts):
    degree = sum(parts)
    full_cycle = tuple(
        (index + 1) % degree for index in range(degree)
    )
    states = set()
    for central_index, central_part in enumerate(parts):
        left_index = (central_index + 1) % 3
        right_index = (central_index + 2) % 3
        left_part = parts[left_index]
        right_part = parts[right_index]
        for cut in range(1, central_part):
            endpoints = (
                0,
                left_part,
                left_part + cut,
                left_part + cut + right_part,
            )
            chord_pair = (
                (endpoints[0], endpoints[1]),
                (endpoints[2], endpoints[3]),
            )
            traced_cycles = cycle_sets(
                compose(
                    involution(degree, chord_pair), full_cycle
                )
            )
            labelled_cycles = [None, None, None]
            labelled_cycles[left_index] = frozenset(
                range(0, left_part)
            )
            labelled_cycles[right_index] = frozenset(
                range(
                    left_part + cut,
                    left_part + cut + right_part,
                )
            )
            labelled_cycles[central_index] = frozenset(
                set(range(degree))
                - labelled_cycles[left_index]
                - labelled_cycles[right_index]
            )
            require(
                set(traced_cycles) == set(labelled_cycles),
                "gap trace does not give the claimed labelled cycles",
            )
            states.add(
                canonical_marked_chords(
                    chord_pair,
                    tuple(labelled_cycles),
                    degree,
                )
            )
    return states


def half_turn_fixed_count(raw_chords, degree):
    if degree % 2:
        return 0
    half_turn = degree // 2
    return sum(
        shifted_chord_pair(
            chord_pair, half_turn, degree
        )
        == chord_pair
        for chord_pair in raw_chords
    )


def check_chords(degree):
    exhaustive = exhaustive_chord_atlas(degree)
    expected_partitions = comb(degree - 1, 2)
    require(
        len(exhaustive["marked"]) == expected_partitions,
        "not every ordered positive pole partition occurs",
    )
    require(
        len(exhaustive["raw"]) == 2 * comb(degree, 4),
        "wrong number of raw noncrossing chord pairs",
    )

    marked_total = 0
    for parts, exhaustive_states in exhaustive["marked"].items():
        explicit_states = gap_atlas(parts)
        require(
            explicit_states == exhaustive_states,
            "gap atlas differs from the exhaustive centralizer quotient",
        )
        require(
            len(explicit_states) == degree - 3,
            "marked Nielsen count is not N-3",
        )
        marked_total += len(explicit_states)

    fixed_counts = []
    for shift in range(degree):
        fixed_counts.append(
            sum(
                shifted_chord_pair(
                    chord_pair, shift, degree
                )
                == chord_pair
                for chord_pair in exhaustive["raw"]
            )
        )
    expected_half_turn = (
        degree * (degree - 2) // 4
        if degree % 2 == 0
        else 0
    )
    require(
        half_turn_fixed_count(
            exhaustive["raw"], degree
        )
        == expected_half_turn,
        "half-turn fixed-pair count",
    )
    for shift, fixed_count in enumerate(fixed_counts):
        expected = 0
        if shift == 0:
            expected = 2 * comb(degree, 4)
        elif degree % 2 == 0 and shift == degree // 2:
            expected = expected_half_turn
        require(
            fixed_count == expected,
            "unexpected nontrivial rotational stabilizer",
        )

    burnside_count = Fraction(sum(fixed_counts), degree)
    formula_count = Fraction(
        (degree - 1) * (degree - 2) * (degree - 3),
        12,
    )
    if degree % 2 == 0:
        formula_count += Fraction(degree - 2, 4)
    require(
        burnside_count.denominator == 1,
        "Burnside count is not integral",
    )
    require(
        len(exhaustive["unmarked"]) == burnside_count,
        "explicit unmarked orbit count differs from Burnside",
    )
    require(
        burnside_count == formula_count,
        "unmarked H3 formula",
    )

    return {
        "partitions": expected_partitions,
        "marked": marked_total,
        "unmarked": int(burnside_count),
        "half_turn": expected_half_turn,
    }


def main():
    path = Path(__file__)
    require(not has_asserts(path), "truth-bearing Python assert found")

    boundary = critical_remainder(1, 1, 1)
    require(
        boundary["maxwell"] == (Fraction(-1, 18),),
        "N=3 empty-boundary Maxwell constant",
    )

    print("THM-2808 INDEPENDENT EXACT AUDIT")
    print("algebra: custom Fraction Q[lambda] engine; no CAS")
    print("Nielsen: exhaustive labelled cycles / full-cycle centralizer")
    print("universe: N=3 boundary; all ordered positive triples, 4<=N<=14")
    print("N | ordered triples | marked charts | unmarked H3 | half-turn fixed")

    algebra_cases = 0
    algebra_roots = 0
    total_marked = 0
    for degree in range(4, 15):
        for a in range(1, degree - 1):
            for b in range(1, degree - a):
                c = degree - a - b
                data = check_algebra(a, b, c)
                algebra_cases += 1
                algebra_roots += len(data["maxwell"]) - 1

        chord_counts = check_chords(degree)
        total_marked += chord_counts["marked"]
        require(
            chord_counts["marked"]
            == comb(degree - 1, 2) * (degree - 3),
            "aggregate marked count",
        )
        print(
            f"{degree:2d} | {chord_counts['partitions']:15d}"
            f" | {chord_counts['marked']:13d}"
            f" | {chord_counts['unmarked']:11d}"
            f" | {chord_counts['half_turn']:15d}"
        )

    require(
        algebra_roots == total_marked,
        "aggregate Maxwell degree and marked atlas disagree",
    )
    print("N=3: Q=-1/18, no Maxwell root, no marked class")
    print(f"ordered algebra cases checked: {algebra_cases}")
    print(f"simple Maxwell roots / marked charts checked: {algebra_roots}")
    print("collision, boundary, E^2, response, and chord-orbit gates: PASS")
    print("assert_nodes=0")


if __name__ == "__main__":
    main()
