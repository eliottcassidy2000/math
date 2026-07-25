#!/usr/bin/env python3
"""Optimization-safe exact controls for THM-2346.

This script verifies the finite uniform ANOVA formulas, the blockwise
Mobius-to-global bridge, the exact score-table gauge hostile, tensor-power
rank cases, equal-prime orbit counts, the sharply typed tournament boundary,
and the THM-2339 two-token word-metric matrices.  It is not a knot-distance
computation.  Every load-bearing check uses ``require`` so normal and
optimized Python execute the same verification.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb


Q = Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def powerset_masks(n: int):
    return range(1 << n)


def mask_coordinates(mask: int, n: int) -> tuple[int, ...]:
    return tuple(i for i in range(n) if mask & (1 << i))


def is_subset(left: int, right: int) -> bool:
    return left & ~right == 0


def product_q(values) -> Q:
    answer = Q(1)
    for value in values:
        answer *= value
    return answer


def all_assignments(m: int, n: int):
    return product(range(m), repeat=n)


def conditional_mean(
    values: dict[tuple[int, ...], Q],
    active: tuple[int, ...],
    active_values: tuple[int, ...],
) -> Q:
    fixed = dict(zip(active, active_values))
    compatible = [
        value
        for assignment, value in values.items()
        if all(assignment[s] == colour for s, colour in fixed.items())
    ]
    return sum(compatible, Q(0)) / len(compatible)


def anova_components(
    values: dict[tuple[int, ...], Q],
    m: int,
    n: int,
) -> dict[int, dict[tuple[int, ...], Q]]:
    components: dict[int, dict[tuple[int, ...], Q]] = {}
    for mask in powerset_masks(n):
        active = mask_coordinates(mask, n)
        table: dict[tuple[int, ...], Q] = {}
        for active_values in all_assignments(m, len(active)):
            total = Q(0)
            for submask in powerset_masks(len(active)):
                subset_positions = mask_coordinates(submask, len(active))
                subset = tuple(active[i] for i in subset_positions)
                subset_values = tuple(
                    active_values[i] for i in subset_positions
                )
                total += (
                    (-1) ** (len(active) - len(subset))
                    * conditional_mean(values, subset, subset_values)
                )
            table[active_values] = total
        components[mask] = table
    return components


def component_value(
    table: dict[tuple[int, ...], Q],
    mask: int,
    assignment: tuple[int, ...],
    n: int,
) -> Q:
    coordinates = mask_coordinates(mask, n)
    return table[tuple(assignment[s] for s in coordinates)]


def check_zero_marginals(
    components: dict[int, dict[tuple[int, ...], Q]],
    m: int,
    n: int,
) -> None:
    for mask, table in components.items():
        active = mask_coordinates(mask, n)
        for position in range(len(active)):
            other_positions = tuple(
                i for i in range(len(active)) if i != position
            )
            for other_values in all_assignments(m, len(other_positions)):
                fixed = dict(zip(other_positions, other_values))
                marginal = Q(0)
                for colour in range(m):
                    point = tuple(
                        colour if i == position else fixed[i]
                        for i in range(len(active))
                    )
                    marginal += table[point]
                require(
                    marginal == 0,
                    f"nonzero marginal at mask={mask}, position={position}",
                )


def inner_product(
    first: dict[tuple[int, ...], Q],
    second: dict[tuple[int, ...], Q],
) -> Q:
    require(first.keys() == second.keys(), "inner-product domain mismatch")
    return sum((first[x] * second[x] for x in first), Q(0)) / len(first)


def lifted_component(
    table: dict[tuple[int, ...], Q],
    mask: int,
    m: int,
    n: int,
) -> dict[tuple[int, ...], Q]:
    return {
        assignment: component_value(table, mask, assignment, n)
        for assignment in all_assignments(m, n)
    }


def subset_mobius(
    subset_scores: dict[tuple[int, int], Q],
    m: int,
    n: int,
) -> dict[tuple[int, int], Q]:
    coefficients: dict[tuple[int, int], Q] = {}
    for colour in range(m):
        for mask in powerset_masks(n):
            if mask == 0:
                continue
            total = Q(0)
            submask = mask
            while True:
                total += (
                    (-1) ** (mask.bit_count() - submask.bit_count())
                    * subset_scores[colour, submask]
                )
                if submask == 0:
                    break
                submask = (submask - 1) & mask
            coefficients[colour, mask] = total
    return coefficients


def allocation_values_from_mobius(
    coefficients: dict[tuple[int, int], Q],
    m: int,
    n: int,
) -> dict[tuple[int, ...], Q]:
    values: dict[tuple[int, ...], Q] = {}
    for assignment in all_assignments(m, n):
        colour_masks = [0] * m
        for token, colour in enumerate(assignment):
            colour_masks[colour] |= 1 << token
        total = Q(0)
        for colour in range(m):
            for mask in powerset_masks(n):
                if mask and is_subset(mask, colour_masks[colour]):
                    total += coefficients[colour, mask]
        values[assignment] = total
    return values


def lambda_coefficient(
    coefficients: dict[tuple[int, int], Q],
    colour: int,
    mask: int,
    m: int,
    n: int,
) -> Q:
    return sum(
        (
            Q(1, m) ** (larger.bit_count() - mask.bit_count())
            * coefficients[colour, larger]
        )
        for larger in powerset_masks(n)
        if larger and is_subset(mask, larger)
    )


def centred_indicator(colour: int, value: int, m: int) -> Q:
    return Q(int(colour == value)) - Q(1, m)


def bridge_component(
    coefficients: dict[tuple[int, int], Q],
    mask: int,
    point: tuple[int, ...],
    m: int,
    n: int,
) -> Q:
    active = mask_coordinates(mask, n)
    if not active:
        return sum(
            (
                Q(1, m) ** larger.bit_count()
                * coefficients[colour, larger]
            )
            for colour in range(m)
            for larger in powerset_masks(n)
            if larger
        )
    return sum(
        (
            lambda_coefficient(
                coefficients, colour, mask, m, n
            )
            * product_q(
                centred_indicator(colour, value, m) for value in point
            )
        )
        for colour in range(m)
    )


def deterministic_scores(m: int, n: int) -> dict[tuple[int, int], Q]:
    scores: dict[tuple[int, int], Q] = {}
    for colour in range(m):
        scores[colour, 0] = Q(0)
        for mask in powerset_masks(n):
            if mask:
                numerator = (
                    (colour + 2) * (mask + 3)
                    + (colour + 1) * mask.bit_count() ** 2
                    - 2 * mask
                )
                scores[colour, mask] = Q(numerator, colour + 3)
    return scores


def check_projection_and_bridge() -> None:
    cases = 0
    for m, n in ((2, 3), (3, 3), (4, 3), (3, 4)):
        scores = deterministic_scores(m, n)
        coefficients = subset_mobius(scores, m, n)
        values = allocation_values_from_mobius(coefficients, m, n)
        components = anova_components(values, m, n)
        check_zero_marginals(components, m, n)

        for assignment, value in values.items():
            reconstructed = sum(
                (
                    component_value(table, mask, assignment, n)
                    for mask, table in components.items()
                ),
                Q(0),
            )
            require(
                reconstructed == value,
                f"ANOVA reconstruction failed at m={m}, n={n}",
            )

        lifted = {
            mask: lifted_component(table, mask, m, n)
            for mask, table in components.items()
        }
        masks = list(powerset_masks(n))
        for i, left in enumerate(masks):
            for right in masks[i + 1 :]:
                require(
                    inner_product(lifted[left], lifted[right]) == 0,
                    f"ANOVA orthogonality failed at m={m}, n={n}",
                )

        for mask, table in components.items():
            active = mask_coordinates(mask, n)
            for point in all_assignments(m, len(active)):
                require(
                    table[point]
                    == bridge_component(
                        coefficients, mask, point, m, n
                    ),
                    f"Mobius bridge failed at m={m}, n={n}, mask={mask}",
                )

        for mask, table in components.items():
            if mask.bit_count() == 2:
                for point in all_assignments(m, 2):
                    require(
                        table[point] == table[point[::-1]],
                        "allocation pair tensor is not symmetric",
                    )
        cases += 1

    print("projection and Mobius-to-ANOVA bridge")
    print(f"exact_cases={cases}")
    print("reconstruction=PASS zero_marginals=PASS orthogonality=PASS")
    print("allocation_pair_symmetry=PASS")


def matrix_rank(matrix: list[list[Q]]) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (
                row
                for row in range(pivot_row, rows)
                if work[row][column] != 0
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [
            entry / pivot_value for entry in work[pivot_row]
        ]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    x - factor * y
                    for x, y in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def tensor_power(vector: tuple[Q, ...], k: int) -> tuple[Q, ...]:
    entries = (Q(1),)
    for _ in range(k):
        entries = tuple(x * y for x in entries for y in vector)
    return entries


def check_gram_ranks() -> None:
    printed = []
    for m in range(2, 7):
        contrasts = [
            tuple(Q(int(i == colour)) - Q(1, m) for i in range(m))
            for colour in range(m)
        ]
        for k in range(1, 6):
            columns = [tensor_power(vector, k) for vector in contrasts]
            matrix = [
                [columns[column][row] for column in range(m)]
                for row in range(m**k)
            ]
            observed_rank = matrix_rank(matrix)
            expected_rank = m - 1 if k == 1 else (1 if m == 2 else m)
            require(
                observed_rank == expected_rank,
                f"tensor rank mismatch at m={m}, k={k}",
            )

            a = Q(m - 1, m)
            b = Q(-1, m)
            eta_perp = a**k - b**k
            eta_one = a**k + (m - 1) * b**k
            require(
                (eta_perp == 0) + (eta_one == 0)
                == m - observed_rank,
                f"Gram eigenvalue multiplicity mismatch at m={m}, k={k}",
            )
        printed.append(
            f"m={m}: ranks k=1..5 "
            + ",".join(
                str(
                    m - 1 if k == 1 else (1 if m == 2 else m)
                )
                for k in range(1, 6)
            )
        )

    print("diagonal Veronese rank bank")
    for line in printed:
        print(line)


def check_two_colour_gauge_hostile() -> None:
    m = 2
    n = 2
    coefficients = {
        (colour, mask): Q(0)
        for colour in range(m)
        for mask in (1, 2, 3)
    }
    coefficients[0, 3] = Q(1)
    coefficients[1, 3] = Q(-1)
    values = allocation_values_from_mobius(coefficients, m, n)
    expected = {
        (0, 0): Q(1),
        (0, 1): Q(0),
        (1, 0): Q(0),
        (1, 1): Q(-1),
    }
    require(values == expected, "two-colour gauge hostile table mismatch")
    components = anova_components(values, m, n)
    require(
        all(value == 0 for value in components[3].values()),
        "blockwise pair hostile retained a global pair interaction",
    )
    expected_unary = {(0,): Q(1, 2), (1,): Q(-1, 2)}
    require(
        components[1] == expected_unary
        and components[2] == expected_unary,
        "two-colour gauge hostile unary fields mismatch",
    )
    require(
        lambda_coefficient(coefficients, 0, 3, m, n)
        + lambda_coefficient(coefficients, 1, 3, m, n)
        == 0,
        "two-colour parity gauge criterion failed",
    )

    print("two-colour global-gauge hostile")
    print("block_pair_coefficients=(1,-1)")
    print("global_table=((1,0),(0,-1))")
    print("global_pair=0 unary_fields=((1/2,-1/2),(1/2,-1/2))")


def histogram(
    assignment: tuple[int, ...],
    token_types: tuple[int, ...],
    type_count: int,
    m: int,
) -> tuple[tuple[int, ...], ...]:
    rows = [[0] * m for _ in range(type_count)]
    for token_type, colour in zip(token_types, assignment):
        rows[token_type][colour] += 1
    return tuple(tuple(row) for row in rows)


def check_equal_prime_quotient() -> None:
    multiplicities = (2, 3)
    token_types = tuple(
        token_type
        for token_type, multiplicity in enumerate(multiplicities)
        for _ in range(multiplicity)
    )
    m = 3
    n = len(token_types)
    histograms = {
        histogram(assignment, token_types, len(multiplicities), m)
        for assignment in all_assignments(m, n)
    }
    expected_count = product_q(
        Q(comb(multiplicity + m - 1, m - 1))
        for multiplicity in multiplicities
    )
    require(
        len(histograms) == expected_count,
        "equal-prime contingency-table orbit count mismatch",
    )

    def energy_from_counts(table: tuple[tuple[int, ...], ...]) -> Q:
        return sum(
            Q((token_type + 1) * (colour + 2) * count * count)
            for token_type, row in enumerate(table)
            for colour, count in enumerate(row)
        )

    by_histogram: dict[tuple[tuple[int, ...], ...], Q] = {}
    for assignment in all_assignments(m, n):
        table = histogram(
            assignment, token_types, len(multiplicities), m
        )
        value = energy_from_counts(table)
        if table in by_histogram:
            require(
                by_histogram[table] == value,
                "energy did not descend to equal-prime quotient",
            )
        by_histogram[table] = value

    print("equal-prime quotient")
    print(
        f"multiplicities={multiplicities} colours={m} "
        f"orbit_count={len(histograms)} expected={expected_count}"
    )
    print("contingency_energy_descent=PASS")


def transpose(matrix: tuple[tuple[Q, ...], ...]):
    return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix)))


def check_tournament_boundary() -> None:
    cyclic = (
        (Q(0), Q(1), Q(-1)),
        (Q(-1), Q(0), Q(1)),
        (Q(1), Q(-1), Q(0)),
    )
    require(
        all(sum(row, Q(0)) == 0 for row in cyclic),
        "cyclic tournament control lacks zero rows",
    )
    require(
        all(
            cyclic[i][j] == -cyclic[j][i]
            for i in range(3)
            for j in range(3)
        ),
        "cyclic tournament control is not skew",
    )
    cyclic_values = {
        (i, j): cyclic[i][j] for i, j in all_assignments(3, 2)
    }
    cyclic_anova = anova_components(cyclic_values, 3, 2)
    require(
        all(value == 0 for value in cyclic_anova[0].values())
        and all(value == 0 for value in cyclic_anova[1].values())
        and all(value == 0 for value in cyclic_anova[2].values()),
        "skew tournament control leaked to lower ANOVA orders",
    )
    require(
        cyclic_anova[3] == cyclic_values,
        "skew tournament pair component mismatch",
    )

    print("tournament typing controls")
    print("vertices=block_colours pair_observable=ordered_token_pair_skew")
    print("three_colour_cyclic_skew_positive_control=PASS")
    print("allocation_pair_alternating_sector=0")


def matrix_anova(matrix: tuple[tuple[Q, ...], ...]):
    m = len(matrix)
    values = {
        (row, column): matrix[row][column]
        for row, column in all_assignments(m, 2)
    }
    components = anova_components(values, m, 2)
    mean = components[0][()]
    unary_row = tuple(components[1][(i,)] for i in range(m))
    unary_column = tuple(components[2][(i,)] for i in range(m))
    pair = tuple(
        tuple(components[3][(i, j)] for j in range(m))
        for i in range(m)
    )
    return mean, unary_row, unary_column, pair


def squared_uniform_norm(matrix: tuple[tuple[Q, ...], ...]) -> Q:
    m = len(matrix)
    return sum((entry * entry for row in matrix for entry in row), Q(0)) / (
        m * m
    )


def check_thm2339_hostile() -> None:
    e0 = ((Q(1), Q(2)), (Q(2), Q(1)))
    e1 = ((Q(0), Q(2)), (Q(2), Q(1)))
    mean0, row0, column0, pair0 = matrix_anova(e0)
    mean1, row1, column1, pair1 = matrix_anova(e1)

    expected_pair0 = (
        (Q(-1, 2), Q(1, 2)),
        (Q(1, 2), Q(-1, 2)),
    )
    expected_pair1 = (
        (Q(-3, 4), Q(3, 4)),
        (Q(3, 4), Q(-3, 4)),
    )
    require(
        (mean0, row0, column0, pair0)
        == (
            Q(3, 2),
            (Q(0), Q(0)),
            (Q(0), Q(0)),
            expected_pair0,
        ),
        "THM-2339 d0 ANOVA mismatch",
    )
    require(
        (mean1, row1, column1, pair1)
        == (
            Q(5, 4),
            (Q(-1, 4), Q(1, 4)),
            (Q(-1, 4), Q(1, 4)),
            expected_pair1,
        ),
        "THM-2339 d1 ANOVA mismatch",
    )
    difference_pair = tuple(
        tuple(pair1[i][j] - pair0[i][j] for j in range(2))
        for i in range(2)
    )
    require(
        difference_pair
        == (
            (Q(-1, 4), Q(1, 4)),
            (Q(1, 4), Q(-1, 4)),
        ),
        "THM-2339 pair difference mismatch",
    )
    require(
        squared_uniform_norm(pair0) == Q(1, 4)
        and squared_uniform_norm(pair1) == Q(9, 16),
        "THM-2339 interaction norm mismatch",
    )
    require(
        pair0 == transpose(pair0) and pair1 == transpose(pair1),
        "THM-2339 hostile pair is not symmetric",
    )
    require(
        min(
            assignment for assignment, value in {
                (i, j): e0[i][j] for i, j in all_assignments(2, 2)
            }.items()
            if value == min(entry for row in e0 for entry in row)
        )
        == (0, 0),
        "THM-2339 d0 minimum positive control failed",
    )
    opt0 = {
        (i, j)
        for i, j in all_assignments(2, 2)
        if e0[i][j] == min(entry for row in e0 for entry in row)
    }
    opt1 = {
        (i, j)
        for i, j in all_assignments(2, 2)
        if e1[i][j] == min(entry for row in e1 for entry in row)
    }
    require(opt0 == {(0, 0), (1, 1)}, "THM-2339 d0 optima mismatch")
    require(opt1 == {(0, 0)}, "THM-2339 d1 optimum mismatch")

    print("THM-2339 two-token hostile")
    print("E0=((1,2),(2,1))")
    print("ANOVA0=3/2-(1/2)h_p*h_q pair_L2_sq=1/4")
    print("E1=((0,2),(2,1))")
    print(
        "ANOVA1=5/4-(1/4)h_p-(1/4)h_q-(3/4)h_p*h_q "
        "pair_L2_sq=9/16"
    )
    print("delta_pair=-(1/4)h_p*h_q")
    print("opt0={(a,a),(b,b)} opt1={(a,a)}")


def main() -> None:
    print("THM-2346 exact global allocation ANOVA controls")
    check_projection_and_bridge()
    check_gram_ranks()
    check_two_colour_gauge_hostile()
    check_equal_prime_quotient()
    check_tournament_boundary()
    check_thm2339_hostile()
    print("all exact controls passed")


if __name__ == "__main__":
    main()
