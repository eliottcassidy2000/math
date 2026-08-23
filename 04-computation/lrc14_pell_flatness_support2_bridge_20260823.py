#!/usr/bin/env python3
"""Exact finite audit of the support-two flatness / Pell selector bridge.

This is a scratch theorem-control surface.  It does not use reserved
THM-3742--3745 as proved dependencies and makes no LRC(14) claim.
"""

from collections import Counter
from fractions import Fraction
from math import gcd, isqrt


CAP = 356
PRIME = 13


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def pell_numbers(last_index: int) -> list[int]:
    values = [0, 1]
    while len(values) <= last_index:
        values.append(2 * values[-1] + values[-2])
    return values


def mat_mul(left, right, modulus=None):
    answer = tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )
    if modulus is not None:
        answer = tuple(tuple(entry % modulus for entry in row) for row in answer)
    return answer


def mat_vec(matrix, vector, modulus=None):
    answer = tuple(sum(matrix[i][j] * vector[j] for j in range(2)) for i in range(2))
    if modulus is not None:
        answer = tuple(entry % modulus for entry in answer)
    return answer


def projective_ratio(first: int, second: int):
    first %= PRIME
    second %= PRIME
    require((first, second) != (0, 0), "primitive ratio cannot vanish mod 13")
    if second == 0:
        return "inf"
    return first * pow(second, -1, PRIME) % PRIME


def distance(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def margin(speeds: tuple[int, ...], time: Fraction) -> Fraction:
    return min(distance(time * speed) for speed in speeds)


def exact_loneliness(speeds: tuple[int, ...]) -> tuple[Fraction, tuple[Fraction, ...]]:
    """Exact lower-envelope maximum on [0,1/2]."""
    denominators = {2 * speed for speed in speeds}
    for index, first in enumerate(speeds):
        for second in speeds[index + 1 :]:
            denominators.add(first + second)
            denominators.add(abs(first - second))
    best_numerator = -1
    best_denominator = 1
    arguments: list[Fraction] = []
    for denominator in denominators:
        if denominator == 0:
            continue
        for numerator in range(denominator // 2 + 1):
            margin_numerator = denominator
            for speed in speeds:
                residue = (numerator * speed) % denominator
                margin_numerator = min(margin_numerator, residue, denominator - residue)
            comparison_left = margin_numerator * best_denominator
            comparison_right = best_numerator * denominator
            if comparison_left > comparison_right:
                best_numerator = margin_numerator
                best_denominator = denominator
                arguments = [Fraction(numerator, denominator)]
            elif comparison_left == comparison_right:
                arguments.append(Fraction(numerator, denominator))
    return Fraction(best_numerator, best_denominator), tuple(sorted(set(arguments)))


def pell_prefix_formula(length: int, pell: list[int], states: list[tuple[int, int]]) -> Fraction:
    """Candidate THM-3744 value, checked here only in a finite universe."""
    if length % 2 == 0:
        return Fraction(
            pell[length] - pell[length - 1] + 1,
            2 * (pell[length] + 1),
        )
    if length % 4 == 1:
        depth = (length - 1) // 4
        x_coordinate = states[depth][0]
        return Fraction(x_coordinate, 2 * pell[2 * depth + 1])
    depth = (length - 3) // 4
    x_coordinate = states[depth + 1][0]
    return Fraction(pell[2 * depth + 1], x_coordinate)


def representable_by_adjacent_generators(value: int, degree: int) -> bool:
    for first_coefficient in range(value // degree + 1):
        remainder = value - first_coefficient * degree
        if remainder % (degree + 1) == 0:
            return True
    return False


def main() -> None:
    pell = pell_numbers(40)

    # Square-triangular state and half-Hadamard compiler.
    update = ((3, 4), (2, 3))  # (x,s) -> (3x+4s,2x+3s)
    state = (1, 0)
    states = [state]
    compiler_rows = []
    for depth in range(1, 16):
        state = mat_vec(update, state)
        states.append(state)
        x_coordinate, even_coordinate = state
        require(even_coordinate % 2 == 0, f"even Pell coordinate depth {depth}")
        square_root = even_coordinate // 2
        degree = (x_coordinate + 1) // 2
        lower = x_coordinate - even_coordinate
        upper = x_coordinate + even_coordinate
        require((lower, upper) == (pell[2 * depth - 1], pell[2 * depth + 1]),
                f"half-angle Pell pair depth {depth}")
        require(degree * (degree - 1) // 2 == square_root * square_root,
                f"square triangular conductor delta depth {depth}")
        require(lower + upper == 4 * degree - 2, f"relation mass depth {depth}")
        require(gcd(lower, upper) == 1 and lower % 2 == upper % 2 == 1,
                f"primitive odd pair depth {depth}")
        compiler_rows.append((depth, degree, square_root, lower, upper, lower + upper))

    compiler_in_deck = [
        row for row in compiler_rows if row[3] < row[4] and row[5] <= CAP
    ]
    require(
        compiler_in_deck
        == [(1, 2, 1, 1, 5, 6), (2, 9, 6, 5, 29, 34), (3, 50, 35, 29, 169, 198)],
        "exact Pell compiler / flatness deck intersection",
    )

    # Central sign and projective half-angle orbit modulo 13.
    identity = ((1, 0), (0, 1))
    negative_identity = ((PRIME - 1, 0), (0, PRIME - 1))
    power = identity
    first_identity = None
    first_negative_identity = None
    for exponent in range(1, 50):
        power = mat_mul(power, update, PRIME)
        if power == negative_identity and first_negative_identity is None:
            first_negative_identity = exponent
        if power == identity:
            first_identity = exponent
            break
    require((first_negative_identity, first_identity) == (7, 14), "mod-13 central order")

    state_mod = (1, 0)
    projective_cycle = []
    degree_residues = []
    inverse_two = pow(2, -1, PRIME)
    for depth in range(14):
        x_coordinate, even_coordinate = state_mod
        lower = x_coordinate - even_coordinate
        upper = x_coordinate + even_coordinate
        projective_cycle.append(projective_ratio(lower, upper))
        degree_residues.append((x_coordinate + 1) * inverse_two % PRIME)
        state_mod = mat_vec(update, state_mod, PRIME)
    require(projective_cycle[:7] == [1, 8, 6, "inf", 0, 11, 5], "projective cycle")
    require(projective_cycle[7:] == projective_cycle[:7], "projective period seven")
    require(
        all(degree_residues[depth + 7] == (1 - degree_residues[depth]) % PRIME
            for depth in range(7)),
        "central sign complements conductor degree residue",
    )
    pell_projective_orbit = set(projective_cycle[:7])

    # Complete support-two ratio deck from the l1<=356 flatness branch.
    ratio_deck = [
        (first, second)
        for first in range(1, CAP)
        for second in range(first + 1, CAP + 1 - first)
        if gcd(first, second) == 1
    ]
    require(len(ratio_deck) == 19314, "flatness ratio deck cardinality")
    distribution = Counter(projective_ratio(*pair) for pair in ratio_deck)
    projective_selected = sum(
        count for residue, count in distribution.items() if residue in pell_projective_orbit
    )
    require((projective_selected, len(ratio_deck) - projective_selected) == (9651, 9663),
            "projective selector split")

    # The mod-13 quotient loses midpoint parity; the primitive ratio quotient
    # loses common scale and hence tie count/locations.
    require(projective_ratio(1, 5) == projective_ratio(1, 18) == 8,
            "same projective Pell observer")
    require((1 % 2 == 1 and 5 % 2 == 1) and not (1 % 2 == 1 and 18 % 2 == 1),
            "projective parity hostile")
    require(gcd(1, 5) == 1 and gcd(2, 10) == 2, "scale hostile gcd")
    ties_primitive = (Fraction(1, 2),)
    ties_scaled = (Fraction(1, 4), Fraction(3, 4))

    # Monomial conductor boundary A=k[b^d,b^(d+1)].  This proves only the
    # semigroup boundary needed as a control, not reserved THM-3745's general F.
    for degree in range(2, 31):
        conductor = degree * (degree - 1)
        gaps = [
            value for value in range(conductor)
            if not representable_by_adjacent_generators(value, degree)
        ]
        require(len(gaps) == degree * (degree - 1) // 2, f"triangular delta d={degree}")
        require(not representable_by_adjacent_generators(conductor - 1, degree),
                f"sharp conductor predecessor d={degree}")
        require(all(representable_by_adjacent_generators(value, degree)
                    for value in range(conductor, conductor + 2 * degree + 2)),
                f"conductor tail d={degree}")
    square_delta_degrees = [
        degree
        for degree in range(2, CAP + 1)
        if isqrt(degree * (degree - 1) // 2) ** 2 == degree * (degree - 1) // 2
    ]
    require(square_delta_degrees == [2, 9, 50, 289], "square triangular degrees <=356")

    # Equating relation l1 mass with conductor degree is only a scalar
    # projection.  Its fibres are already nontrivial and grow rapidly.
    scalar_shell_counts = {
        degree: sum(first + second == degree for first, second in ratio_deck)
        for degree in square_delta_degrees
    }
    require(scalar_shell_counts == {2: 0, 9: 3, 50: 10, 289: 136},
            "untyped scalar-shell hostile")

    # Reserved THM-3744 has no live proved statement.  Independently replay
    # its candidate prefix values through the full thirteen-speed prefix.
    prefix_rows = []
    for length in range(1, 14):
        speeds = tuple(pell[1 : length + 1])
        value, arguments = exact_loneliness(speeds)
        predicted = pell_prefix_formula(length, pell, states)
        require(value == predicted, f"Pell-prefix candidate formula length {length}")
        require(arguments == (value,), f"unique prefix endpoint length {length}")
        require(value > Fraction(1, 14), f"finite Pell prefix is safely above LRC14 floor {length}")
        prefix_rows.append((length, pell[length], value))

    # In THM-2052's notation, every support-two row here is automatically in
    # W_(91^6,3); it cannot be the outside-W rank-twelve increment.
    require(CAP < 91**6, "support-two flatness row lies in rank code")

    print("PELL_COMPILER_ROWS_IN_FLATNESS_DECK", compiler_in_deck)
    print("MOD13_MATRIX_CENTRAL_ORDER", first_negative_identity, first_identity)
    print("MOD13_PROJECTIVE_CYCLE", projective_cycle[:7])
    print("MOD13_CONDUCTOR_DEGREE_RESIDUES_14", degree_residues)
    ordered_distribution = [distribution[index] for index in range(PRIME)] + [distribution["inf"]]
    print("RATIO_DECK_SIZE_AND_MOD13_DISTRIBUTION", len(ratio_deck), ordered_distribution)
    print("PELL_PROJECTIVE_SELECTED_VS_COMPLEMENT", projective_selected,
          len(ratio_deck) - projective_selected)
    print("PROJECTIVE_PARITY_HOSTILE", (1, 5), (1, 18), projective_ratio(1, 5))
    print("RATIO_SCALE_TIE_HOSTILE", (1, 5), ties_primitive, (2, 10), ties_scaled)
    print("MONOMIAL_SQUARE_DELTA_DEGREES_LE_356", square_delta_degrees)
    print("UNTYPED_L1_SHELL_FIBRES", scalar_shell_counts)
    print("PELL_PREFIX_FINITE_EXACT_ROWS", prefix_rows)
    print("THM2052_SUPPORT2_BRANCH", "inside_W_only")
    print("RESULT: PASS; NO_LRC14_CLAIM")


if __name__ == "__main__":
    main()

