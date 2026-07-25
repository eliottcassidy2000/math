#!/usr/bin/env python3
"""Exact companion for THM-2249.

Universe:
  * every 3x3 nonnegative integer transport table with all margins N,
    for 1 <= N <= 12;
  * every bijection between two labelled six-vertex, three-block cores;
  * all 2^3 source and 2^3 target internal orientations for those cores.

The proof of THM-2249 is uniform in N.  The bounded census is an independent
hostile audit of the universal zero-layer lemma through quotient order five,
the directed-triangle formula, its sharp floor, and the claim that the forced
quotient pairs are a literal subset of the reversal cost.
"""

from itertools import permutations, product


Z3 = range(3)
PERMS = tuple(permutations(Z3))


def beats(a: int, b: int) -> bool:
    """Arc relation a -> b in the directed triangle."""
    return (b - a) % 3 == 1


def rotation(a: int) -> tuple[int, int, int]:
    return tuple((i + a) % 3 for i in Z3)


def reflection(a: int) -> tuple[int, int, int]:
    return tuple((a - i) % 3 for i in Z3)


ROTATIONS = tuple(rotation(a) for a in Z3)
REFLECTIONS = tuple(reflection(a) for a in Z3)


def pair_frustration(
    sigma: tuple[int, int, int], tau: tuple[int, int, int]
) -> int:
    return sum(
        tau[(i + 1) % 3] == (sigma[i] - 1) % 3
        for i in Z3
    )


def frustration(x: tuple[tuple[int, ...], ...]) -> int:
    return sum(
        x[i][j] * x[(i + 1) % 3][(j - 1) % 3]
        for i in Z3
        for j in Z3
    )


def transport_tables(n: int):
    """Enumerate every 3x3 integer table with all margins n."""
    for a00 in range(n + 1):
        for a01 in range(n - a00 + 1):
            a02 = n - a00 - a01
            for a10 in range(n + 1):
                for a11 in range(n - a10 + 1):
                    a12 = n - a10 - a11
                    a20 = n - a00 - a10
                    a21 = n - a01 - a11
                    a22 = n - a02 - a12
                    if min(a20, a21, a22) < 0:
                        continue
                    if a20 + a21 + a22 != n:
                        continue
                    yield (
                        (a00, a01, a02),
                        (a10, a11, a12),
                        (a20, a21, a22),
                    )


def decompose_into_permutations(
    x: tuple[tuple[int, ...], ...], n: int
) -> dict[tuple[int, int, int], int]:
    """Lexicographic Hall decomposition of an integer transport."""
    work = [list(row) for row in x]
    counts = {sigma: 0 for sigma in PERMS}
    for _ in range(n):
        sigma = next(
            (
                candidate
                for candidate in PERMS
                if all(work[i][candidate[i]] > 0 for i in Z3)
            ),
            None,
        )
        assert sigma is not None
        counts[sigma] += 1
        for i in Z3:
            work[i][sigma[i]] -= 1
    assert all(value == 0 for row in work for value in row)
    return counts


def decomposition_formula(
    counts: dict[tuple[int, int, int], int]
) -> int:
    r = [counts[sigma] for sigma in ROTATIONS]
    s = [counts[sigma] for sigma in REFLECTIONS]
    r_total = sum(r)
    s_total = sum(s)
    return (
        3 * (r[0] * r[1] + r[1] * r[2] + r[2] * r[0])
        + 2 * r_total * s_total
        + 3 * sum(value * value for value in s)
    )


def sharp_floor(n: int) -> int:
    if n == 1:
        return 3
    return min(3 * (n - 1), 2 * n + 1)


def matrix_of_permutation(
    sigma: tuple[int, int, int], multiplicity: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(multiplicity if sigma[i] == j else 0 for j in Z3)
        for i in Z3
    )


def add_matrices(*matrices):
    return tuple(
        tuple(sum(matrix[i][j] for matrix in matrices) for j in Z3)
        for i in Z3
    )


def tournament_from_mask(order: int, mask: int):
    adjacency = [
        [False for _ in range(order)] for _ in range(order)
    ]
    bit = 0
    for i in range(order):
        for j in range(i + 1, order):
            if (mask >> bit) & 1:
                adjacency[i][j] = True
            else:
                adjacency[j][i] = True
            bit += 1
    return adjacency


def quotient_pair_frustration(adjacency, sigma, tau) -> int:
    order = len(adjacency)
    return sum(
        adjacency[i][k] and adjacency[tau[k]][sigma[i]]
        for i in range(order)
        for k in range(order)
    )


def is_automorphism(adjacency, sigma) -> bool:
    order = len(adjacency)
    return all(
        adjacency[i][k] == adjacency[sigma[i]][sigma[k]]
        for i in range(order)
        for k in range(i + 1, order)
    )


def audit_universal_zero_layer() -> tuple[int, int, int]:
    tournament_count = 0
    permutation_checks = 0
    automorphism_pair_checks = 0
    for order in range(2, 6):
        sym = tuple(permutations(range(order)))
        for mask in range(1 << (order * (order - 1) // 2)):
            tournament_count += 1
            adjacency = tournament_from_mask(order, mask)
            automorphisms = []
            for sigma in sym:
                diagonal_zero = (
                    quotient_pair_frustration(adjacency, sigma, sigma) == 0
                )
                assert diagonal_zero == is_automorphism(adjacency, sigma)
                permutation_checks += 1
                if diagonal_zero:
                    automorphisms.append(sigma)
            for index, sigma in enumerate(automorphisms):
                for tau in automorphisms[index + 1 :]:
                    assert (
                        quotient_pair_frustration(adjacency, sigma, tau)
                        + quotient_pair_frustration(adjacency, tau, sigma)
                        > 0
                    )
                    automorphism_pair_checks += 1
    return tournament_count, permutation_checks, automorphism_pair_checks


def core_arc(
    vertex_a: tuple[int, int],
    vertex_b: tuple[int, int],
    internal_bits: tuple[int, int, int],
) -> bool:
    block_a, local_a = vertex_a
    block_b, local_b = vertex_b
    if block_a != block_b:
        return beats(block_a, block_b)
    assert local_a != local_b
    forward = bool(internal_bits[block_a])
    return forward if (local_a, local_b) == (0, 1) else not forward


def audit_six_vertex_bijections() -> tuple[int, int]:
    vertices = tuple((block, local) for block in Z3 for local in range(2))
    source_pairs = tuple(
        (u, v)
        for index, u in enumerate(vertices)
        for v in vertices[index + 1 :]
    )
    bijection_count = 0
    orientation_checks = 0

    for target_order in permutations(range(6)):
        bijection_count += 1
        image = {
            vertices[source_index]: vertices[target_index]
            for source_index, target_index in enumerate(target_order)
        }
        table = tuple(
            tuple(
                sum(
                    image[(source_block, local)][0] == target_block
                    for local in range(2)
                )
                for target_block in Z3
            )
            for source_block in Z3
        )
        forced = 0
        for u, v in source_pairs:
            image_u, image_v = image[u], image[v]
            if u[0] == v[0] or image_u[0] == image_v[0]:
                continue
            if beats(u[0], v[0]) != beats(image_u[0], image_v[0]):
                forced += 1
        assert forced == frustration(table)

        for source_bits in product((0, 1), repeat=3):
            for target_bits in product((0, 1), repeat=3):
                full_cost = sum(
                    core_arc(u, v, source_bits)
                    != core_arc(image[u], image[v], target_bits)
                    for u, v in source_pairs
                )
                assert full_cost >= forced
                orientation_checks += 1

    return bijection_count, orientation_checks


def main() -> None:
    print("THM-2249 DIRECTED-TRIANGLE FRUSTRATION AUDIT")
    print("pair table rows/cols R0 R1 R2 F0 F1 F2")
    ordered = ROTATIONS + REFLECTIONS
    table = tuple(
        tuple(pair_frustration(sigma, tau) for tau in ordered)
        for sigma in ordered
    )
    expected = (
        (0, 3, 0, 1, 1, 1),
        (0, 0, 3, 1, 1, 1),
        (3, 0, 0, 1, 1, 1),
        (1, 1, 1, 3, 0, 0),
        (1, 1, 1, 0, 3, 0),
        (1, 1, 1, 0, 0, 3),
    )
    assert table == expected
    for row in table:
        print(" ".join(map(str, row)))

    tournaments, permutation_checks, automorphism_pairs = (
        audit_universal_zero_layer()
    )
    print(
        "universal_zero_layer "
        f"orders=2..5 tournaments={tournaments} "
        f"permutation_checks={permutation_checks} "
        f"automorphism_pair_checks={automorphism_pairs}"
    )

    for n in range(1, 13):
        tables = tuple(transport_tables(n))
        rotation_matrices = {
            matrix_of_permutation(sigma, n) for sigma in ROTATIONS
        }
        zero = []
        nonzero_values = []
        for x in tables:
            counts = decompose_into_permutations(x, n)
            value = frustration(x)
            assert value == decomposition_formula(counts)
            if value == 0:
                zero.append(x)
            else:
                nonzero_values.append(value)
        assert set(zero) == rotation_matrices
        observed_floor = min(nonzero_values)
        assert observed_floor == sharp_floor(n)

        reflection_witness = add_matrices(
            matrix_of_permutation(ROTATIONS[0], n - 1),
            matrix_of_permutation(REFLECTIONS[0], 1),
        )
        assert frustration(reflection_witness) == 2 * n + 1
        if n >= 2:
            rotation_witness = add_matrices(
                matrix_of_permutation(ROTATIONS[0], n - 1),
                matrix_of_permutation(ROTATIONS[1], 1),
            )
            assert frustration(rotation_witness) == 3 * (n - 1)

        print(
            f"N={n:2d} tables={len(tables):4d} zero={len(zero)} "
            f"min_positive={observed_floor:2d} sharp_floor={sharp_floor(n):2d}"
        )

    bijections, orientation_checks = audit_six_vertex_bijections()
    print(
        f"N=2 bijections={bijections} "
        f"internal_orientation_checks={orientation_checks}"
    )
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
