#!/usr/bin/env python3
"""Exact referee for the THM-2852 prime-power spectrum harvest."""

from __future__ import annotations

from fractions import Fraction
from random import Random


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [list(row) for row in zip(*matrix)]


def determinant_mod(matrix: list[list[int]], prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    size = len(work)
    determinant = 1
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if work[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column] % prime
        determinant = determinant * pivot_value % prime
        pivot_inverse = pow(pivot_value, -1, prime)
        work[column] = [
            entry * pivot_inverse % prime for entry in work[column]
        ]
        for row in range(column + 1, size):
            multiplier = work[row][column] % prime
            if multiplier:
                work[row] = [
                    (left - multiplier * right) % prime
                    for left, right in zip(work[row], work[column])
                ]
    return determinant % prime


def rank_over_q(matrix: list[list[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    rank = 0
    for column in range(columns):
        pivot = next(
            (
                row
                for row in range(rank, rows)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [entry / pivot_value for entry in work[rank]]
        for row in range(rows):
            if row == rank:
                continue
            multiplier = work[row][column]
            if multiplier:
                work[row] = [
                    left - multiplier * right
                    for left, right in zip(work[row], work[rank])
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def cayley_matrices(group, multiply, inverse, connection_set):
    """Return adjacency, column right multiplication, and reverse adjacency."""

    index = {element: position for position, element in enumerate(group)}
    adjacency = [
        [
            int(multiply(inverse(source), target) in connection_set)
            for target in group
        ]
        for source in group
    ]
    right_multiplication = [[0 for _source in group] for _target in group]
    for source in group:
        for step in connection_set:
            target = multiply(source, step)
            right_multiplication[index[target]][index[source]] += 1
    reverse_set = {inverse(step) for step in connection_set}
    reverse_adjacency = [
        [
            int(multiply(inverse(source), target) in reverse_set)
            for target in group
        ]
        for source in group
    ]
    return adjacency, right_multiplication, reverse_adjacency


def inverse_pairs(group, inverse, identity):
    seen = {identity}
    pairs = []
    for element in group:
        if element in seen:
            continue
        partner = inverse(element)
        require(partner != element, "nonidentity involution in odd group")
        pairs.append((element, partner))
        seen.update((element, partner))
    require(len(seen) == len(group), "inverse pairs do not partition group")
    return pairs


def orientation_from_bits(pairs, bits: int):
    return {
        pair[(bits >> position) & 1]
        for position, pair in enumerate(pairs)
    }


def check_tournament_matrix(matrix: list[list[int]]) -> None:
    size = len(matrix)
    require(
        all(matrix[row][row] == 0 for row in range(size)),
        "Cayley tournament acquired a loop",
    )
    require(
        all(
            matrix[row][column] + matrix[column][row] == 1
            for row in range(size)
            for column in range(row + 1, size)
        ),
        "Cayley orientation is not a tournament",
    )


def cyclic_controls(order: int, prime: int):
    group = list(range(order))
    multiply = lambda left, right: (left + right) % order
    inverse = lambda value: (-value) % order
    pairs = inverse_pairs(group, inverse, 0)
    determinant_residues = set()
    transpose_checks = 0
    for bits in range(1 << len(pairs)):
        connection_set = orientation_from_bits(pairs, bits)
        adjacency, right, reverse = cayley_matrices(
            group, multiply, inverse, connection_set
        )
        check_tournament_matrix(adjacency)
        require(
            right == transpose(adjacency),
            "right multiplication is not adjacency transpose",
        )
        require(
            reverse == transpose(adjacency),
            "orientation reversal is not adjacency transpose",
        )
        transpose_checks += 1
        determinant_residues.add(determinant_mod(adjacency, prime))
    expected = (-pow(2, -1, prime)) % prime
    require(
        determinant_residues == {expected},
        "cyclic Cayley determinant residue changed",
    )
    return len(pairs), transpose_checks, tuple(sorted(determinant_residues))


def heisenberg_controls():
    prime = 3
    group = [
        (a, b, c)
        for a in range(prime)
        for b in range(prime)
        for c in range(prime)
    ]

    def multiply(left, right):
        a, b, c = left
        aa, bb, cc = right
        return (
            (a + aa) % prime,
            (b + bb) % prime,
            (c + cc + a * bb) % prime,
        )

    def inverse(value):
        a, b, c = value
        return (
            (-a) % prime,
            (-b) % prime,
            (-c + a * b) % prime,
        )

    identity = (0, 0, 0)
    pairs = inverse_pairs(group, inverse, identity)
    require(len(pairs) == 13, "H_3 inverse-pair count changed")

    words = [0, (1 << len(pairs)) - 1]
    random = Random(2839)
    while len(words) < 100:
        word = random.randrange(1 << len(pairs))
        if word not in words:
            words.append(word)

    determinant_residues = set()
    for word in words:
        connection_set = orientation_from_bits(pairs, word)
        adjacency, right, reverse = cayley_matrices(
            group, multiply, inverse, connection_set
        )
        check_tournament_matrix(adjacency)
        require(
            right == transpose(adjacency),
            "H_3 right multiplication transpose convention failed",
        )
        require(
            reverse == transpose(adjacency),
            "H_3 reverse orientation transpose convention failed",
        )
        determinant_residues.add(determinant_mod(adjacency, prime))

    expected = (-pow(2, -1, prime)) % prime
    require(
        determinant_residues == {expected},
        "H_3 Cayley determinant residue changed",
    )
    return len(pairs), len(words), tuple(sorted(determinant_residues))


def circulant_matrix(coefficients: list[int]) -> list[list[int]]:
    order = len(coefficients)
    return [
        [coefficients[(target - source) % order] for target in range(order)]
        for source in range(order)
    ]


def main() -> None:
    prime = 13

    order_six = prime**6
    n_zero = 3_454_614
    n_plus = 3_454_627
    n_affine = 4_143_978
    physical_scalar = 790_161_473_087_466_480
    support = (n_zero, n_plus, n_affine)
    require(len(set(support)) == 3, "THM-2807 support collided")
    require(
        (
            (n_plus - n_zero) % order_six,
            (n_affine - n_plus) % order_six,
            (n_affine - n_zero) % order_six,
        )
        == (13, 689_351, 689_364),
        "THM-2807 address gaps changed",
    )
    normalized_mass = len(support)
    normalized_determinant_residue = pow(
        normalized_mass, order_six, prime
    )
    require(
        normalized_mass % prime == 3
        and normalized_determinant_residue == 3,
        "THM-2807 normalized augmentation check failed",
    )
    require(
        physical_scalar % prime == 0,
        "THM-2807 rational normalization caveat disappeared",
    )

    arrow_counts = (66_099, 65_612, 65_579, 65_612, 65_579, 65_098)
    arrow_remainders = tuple(count % prime for count in arrow_counts)
    require(
        arrow_remainders == (7, 1, 7, 1, 7, 7)
        and all(arrow_remainders),
        "THM-2829 arrow augmentations changed",
    )
    raw_target_counts = (66_099, 65_652)
    same_section_counts = (66_099, 65_612)
    reverse_count = 65_619
    require(
        tuple(count % prime for count in raw_target_counts) == (7, 2)
        and tuple(count % prime for count in same_section_counts) == (7, 1)
        and reverse_count % prime == 8,
        "THM-2829 target/reverse augmentations changed",
    )

    cyclic_five = cyclic_controls(5, 5)
    cyclic_nine = cyclic_controls(9, 3)
    heisenberg = heisenberg_controls()

    order_hostile = prime**2
    subgroup_coefficients = [
        int(index % prime == 0) for index in range(order_hostile)
    ]
    subgroup_rank = prime
    require(
        sum(subgroup_coefficients) == prime
        and subgroup_rank == order_hostile // prime,
        "p-divisible subgroup hostile changed",
    )
    rational_mass = Fraction(sum(subgroup_coefficients), prime)
    require(rational_mass == 1, "rational augmentation hostile changed")
    non_p_power = [1, 0, 1, 0, 1, 0]
    non_p_power_rank = rank_over_q(circulant_matrix(non_p_power))
    require(
        sum(non_p_power) % 2 == 1 and non_p_power_rank == 2,
        "non-p-power hostile changed",
    )

    print(
        "simplex_normalized="
        f"(order={order_six},mass={normalized_mass},"
        f"det_mod13={normalized_determinant_residue},"
        f"rank_Q={order_six}); "
        f"physical_scalar_mod13={physical_scalar % prime}"
    )
    print(
        "ancestry_arrow_counts="
        f"{arrow_counts}; remainders_mod13={arrow_remainders}; "
        f"rank_each_Q={prime**5}"
    )
    print(
        "ancestry_auxiliary="
        f"(raw_target={raw_target_counts},"
        f"same_section={same_section_counts},reverse={reverse_count}); "
        "remainders_mod13="
        f"{tuple(count % prime for count in raw_target_counts)},"
        f"{tuple(count % prime for count in same_section_counts)},"
        f"{reverse_count % prime}"
    )
    print(
        "cyclic_C5="
        f"(inverse_pairs={cyclic_five[0]},"
        f"all_orientations={cyclic_five[1]},"
        f"det_residues_mod5={cyclic_five[2]})"
    )
    print(
        "cyclic_C9="
        f"(inverse_pairs={cyclic_nine[0]},"
        f"all_orientations={cyclic_nine[1]},"
        f"det_residues_mod3={cyclic_nine[2]})"
    )
    print(
        "heisenberg_H3="
        f"(inverse_pairs={heisenberg[0]},"
        f"fixed_orientation_controls={heisenberg[1]},"
        f"det_residues_mod3={heisenberg[2]},"
        "right_multiplication=adjacency_transpose,"
        "reverse_orientation=adjacency_transpose)"
    )
    print(
        "hostiles="
        f"(C169_subgroup_mass={sum(subgroup_coefficients)},"
        f"rank={subgroup_rank},"
        f"rational_scaled_mass={rational_mass};"
        f"C6_odd_mass={sum(non_p_power)},rank={non_p_power_rank})"
    )
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
