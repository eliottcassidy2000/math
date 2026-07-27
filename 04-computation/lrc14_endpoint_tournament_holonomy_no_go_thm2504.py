#!/usr/bin/env python3
"""Exact companion for THM-2504.

Audits which endpoint-mask and root-label relations can define an intrinsic
tournament.  All graph arithmetic is exact over GF(2); no third-party package
is used.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb


MASK_BITS = 7
MASK_COUNT = 1 << MASK_BITS
ROOT_COUNT = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def word(mask: int) -> tuple[int, ...]:
    """Seven-bit word in the THM-2502 role order, most significant first."""

    return tuple((mask >> (MASK_BITS - 1 - index)) & 1 for index in range(MASK_BITS))


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    """Move coordinate i to coordinate permutation[i]."""

    result = 0
    for index in range(MASK_BITS):
        if mask & (1 << index):
            result |= 1 << permutation[index]
    return result


def swapping_permutation(left: int, right: int) -> tuple[int, ...]:
    """An involutive coordinate permutation swapping equal-weight masks."""

    left_only = [
        index
        for index in range(MASK_BITS)
        if (left & (1 << index)) and not (right & (1 << index))
    ]
    right_only = [
        index
        for index in range(MASK_BITS)
        if (right & (1 << index)) and not (left & (1 << index))
    ]
    require(len(left_only) == len(right_only), "equal-weight swap cardinality")
    permutation = list(range(MASK_BITS))
    for first, second in zip(left_only, right_only):
        permutation[first] = second
        permutation[second] = first
    return tuple(permutation)


def gf2_rank(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    rank = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(rows):
            if row != rank and work[row][column]:
                work[row] = [
                    entry ^ pivot_entry
                    for entry, pivot_entry in zip(work[row], work[rank])
                ]
        rank += 1
    return rank


def gf2_solve(matrix: list[list[int]], target: list[int]) -> list[int]:
    """Return the zero-free-variable solution of a consistent GF(2) system."""

    work = [row[:] + [entry] for row, entry in zip(matrix, target)]
    rows = len(work)
    columns = len(matrix[0]) if rows else 0
    rank = 0
    pivots: list[int] = []
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(rows):
            if row != rank and work[row][column]:
                work[row] = [
                    entry ^ pivot_entry
                    for entry, pivot_entry in zip(work[row], work[rank])
                ]
        pivots.append(column)
        rank += 1
    for row in range(rank, rows):
        require(not work[row][-1], "inconsistent GF(2) system")
    solution = [0] * columns
    for row, column in enumerate(pivots):
        solution[column] = work[row][-1]
    return solution


def tile_edges(vertex_count: int) -> list[tuple[int, int]]:
    """Edges of K_n with the numeric successor path removed."""

    return [
        (left, right)
        for left in range(vertex_count)
        for right in range(left + 2, vertex_count)
    ]


def laplacian(
    vertex_count: int,
    edges: list[tuple[int, int]],
) -> list[list[int]]:
    matrix = [[0] * vertex_count for _ in range(vertex_count)]
    for left, right in edges:
        matrix[left][left] ^= 1
        matrix[right][right] ^= 1
        matrix[left][right] ^= 1
        matrix[right][left] ^= 1
    return matrix


def boundary(
    vertex_count: int,
    edges: list[tuple[int, int]],
    cochain: list[int],
) -> list[int]:
    result = [0] * vertex_count
    for value, (left, right) in zip(cochain, edges):
        if value:
            result[left] ^= 1
            result[right] ^= 1
    return result


def cut_cycle_decomposition(
    vertex_count: int,
    cochain: list[int],
) -> tuple[list[int], list[int], list[int], int]:
    """Decompose a tile cochain as cut + cycle when the bicycle is zero."""

    edges = tile_edges(vertex_count)
    require(len(cochain) == len(edges), "cochain length")
    graph_laplacian = laplacian(vertex_count, edges)
    rank = gf2_rank(graph_laplacian)
    require(rank == vertex_count - 1, "tile graph Laplacian corank")

    cochain_boundary = boundary(vertex_count, edges, cochain)
    # Fix the constant-potential gauge by p_0=0.  The remaining principal
    # minor is invertible exactly because the bicycle space is zero here.
    reduced = [row[1:] for row in graph_laplacian[1:]]
    potential = [0] + gf2_solve(reduced, cochain_boundary[1:])
    cut = [potential[left] ^ potential[right] for left, right in edges]
    cycle = [entry ^ cut_entry for entry, cut_entry in zip(cochain, cut)]
    require(not any(boundary(vertex_count, edges, cycle)), "cycle boundary")
    require(
        all(entry == (cut_entry ^ cycle_entry) for entry, cut_entry, cycle_entry in zip(cochain, cut, cycle)),
        "cut-cycle reconstruction",
    )
    require(gf2_rank(reduced) == vertex_count - 1, "unique cut-cycle split")
    return potential, cut, cycle, rank


def first_failure_data() -> tuple[Counter, Counter]:
    cells: Counter = Counter()
    layers: Counter = Counter()
    for mask in range(MASK_COUNT):
        bits = word(mask)
        first_failure = next(
            (index + 1 for index, bit in enumerate(bits[:5]) if bit == 0),
            0,
        )
        blocker_word = bits[5:]
        cells[(first_failure, blocker_word)] += 1
        layers[first_failure] += 1
    return cells, layers


def cyclic_autocorrelation(support: set[int], modulus: int) -> tuple[int, ...]:
    return tuple(
        sum(int((point + shift) % modulus in support) for point in support)
        for shift in range(modulus)
    )


def main() -> None:
    masks = range(MASK_COUNT)
    mask_pairs = list(combinations(masks, 2))
    require(len(mask_pairs) == 8128, "mask pair count")

    # Every pair is swapped by a translation of the abstract untyped cube.
    for left, right in mask_pairs:
        difference = left ^ right
        require((left ^ difference, right ^ difference) == (right, left), "translation swap")

    # With polarity fixed but coordinate names forgotten, every equal-weight
    # pair is swapped by an explicit coordinate involution.
    equal_weight_pairs = 0
    for left, right in mask_pairs:
        if left.bit_count() != right.bit_count():
            continue
        equal_weight_pairs += 1
        permutation = swapping_permutation(left, right)
        require(permute_mask(left, permutation) == right, "left coordinate swap")
        require(permute_mask(right, permutation) == left, "right coordinate swap")
    equal_weight_formula = sum(comb(comb(MASK_BITS, weight), 2) for weight in range(MASK_BITS + 1))
    require(equal_weight_pairs == equal_weight_formula == 1652, "equal-weight ties")

    parity_oriented = sum(
        int((left.bit_count() - right.bit_count()) & 1)
        for left, right in mask_pairs
    )
    parity_ties = len(mask_pairs) - parity_oriented
    require(parity_oriented == 4096 and parity_ties == 4032, "Hamming-parity census")

    cells, layers = first_failure_data()
    expected_cell_sizes = (1, 16, 8, 4, 2, 1)
    expected_layer_sizes = (4, 64, 32, 16, 8, 4)
    for blocker_word in ((0, 0), (0, 1), (1, 0), (1, 1)):
        require(
            tuple(cells[(label, blocker_word)] for label in range(6)) == expected_cell_sizes,
            "24-cell sizes",
        )
    require(tuple(layers[label] for label in range(6)) == expected_layer_sizes, "layer sizes")
    cell_ties = sum(comb(size, 2) for size in cells.values())
    layer_ties = sum(comb(size, 2) for size in layers.values())
    require(cell_ties == 620 and layer_ties == 2672, "cell/layer ties")

    # The THM-2502 lex tournament is exactly numeric order in this bit gauge.
    lex_reversal = []
    for left, right in mask_pairs:
        first_difference = next(
            index
            for index, (left_bit, right_bit) in enumerate(zip(word(left), word(right)))
            if left_bit != right_bit
        )
        lex_points_forward = word(left)[first_difference] < word(right)[first_difference]
        lex_reversal.append(int(not lex_points_forward))
    require(not any(lex_reversal), "lex tournament equals numeric order")

    edges_128 = tile_edges(MASK_COUNT)
    parity_cut_128 = [
        (left.bit_count() ^ right.bit_count()) & 1
        for left, right in edges_128
    ]
    potential_128, cut_128, cycle_128, rank_128 = cut_cycle_decomposition(
        MASK_COUNT,
        parity_cut_128,
    )
    require(potential_128 == [mask.bit_count() & 1 for mask in masks], "parity potential")
    require(cut_128 == parity_cut_128 and not any(cycle_128), "parity is an exact cut")
    require(len(edges_128) == 8001 and sum(parity_cut_128) == 4011, "128 tile census")

    # The uniform common-chart packet model explicitly constructed in
    # THM-2457 is bidirected on every pair, so comparison of opposite entries
    # produces no tournament edge.  This is not a full scalar LRC row.
    uniform_entry = Fraction(1, MASK_COUNT**2)
    uniform_matrix = [[uniform_entry] * MASK_COUNT for _ in masks]
    uniform_asymmetry_ties = sum(
        int(uniform_matrix[left][right] == uniform_matrix[right][left])
        for left, right in mask_pairs
    )
    require(uniform_asymmetry_ties == 8128, "uniform co-support ties")
    sparse_loop_model = [[0] * MASK_COUNT for _ in masks]
    sparse_loop_model[5][5] = 1
    require(
        all(
            sparse_loop_model[left][right] == sparse_loop_model[right][left] == 0
            for left, right in mask_pairs
        ),
        "sparse loop model has no off-diagonal relation",
    )

    edges_13 = tile_edges(ROOT_COUNT)
    cyclic_reversal_13 = [
        int(right - left >= 7)
        for left, right in edges_13
    ]
    _, cyclic_cut_13, cyclic_cycle_13, rank_13 = cut_cycle_decomposition(
        ROOT_COUNT,
        cyclic_reversal_13,
    )
    require(len(edges_13) == 66 and sum(cyclic_reversal_13) == 21, "13-root cyclic support")
    require(sum(cyclic_cut_13) == 38 and sum(cyclic_cycle_13) == 35, "13-root split")

    triangle = ((0, 2), (2, 7), (0, 7))
    edge_value = dict(zip(edges_13, cyclic_reversal_13))
    triangle_holonomy = sum(edge_value[edge] for edge in triangle) & 1
    require(triangle_holonomy == 1, "nonzero cyclic triangle holonomy")

    global_converse_fixed_path_13 = [entry ^ 1 for entry in cyclic_reversal_13]
    _, converse_cut_13, converse_cycle_13, _ = cut_cycle_decomposition(
        ROOT_COUNT,
        global_converse_fixed_path_13,
    )
    require(sum(converse_cut_13) == 30 and sum(converse_cycle_13) == 39, "fixed-path converse split")

    # A physical relabeling u -> -u also transports the numeric Hamiltonian
    # path.  Transport the whole decomposition edge by edge; its support is
    # isomorphic to the original 38/35 split.  Thus 30/39 is deliberately a
    # fixed-path gauge calculation, not a gauge-invariant support count.
    inversion = lambda vertex: -vertex % ROOT_COUNT
    transported_edges_13 = [
        tuple(sorted((inversion(left), inversion(right))))
        for left, right in edges_13
    ]
    require(len(set(transported_edges_13)) == len(edges_13), "transported tile graph")
    require(
        gf2_rank(laplacian(ROOT_COUNT, transported_edges_13)) == ROOT_COUNT - 1,
        "transported Laplacian rank",
    )
    require(
        not any(boundary(ROOT_COUNT, transported_edges_13, cyclic_cycle_13)),
        "transported cycle boundary",
    )

    quadratic_residues = {pow(value, 2, ROOT_COUNT) for value in range(1, ROOT_COUNT)}
    require(quadratic_residues == {1, 3, 4, 9, 10, 12}, "quadratic residue set")
    require(quadratic_residues == {-value % ROOT_COUNT for value in quadratic_residues}, "QR negation symmetry")

    # THM-2458 four-root guard AP with start 0 and step 5.
    guard_support = {0, 5, 10, 2}
    guard_autocorrelation = cyclic_autocorrelation(guard_support, ROOT_COUNT)
    expected_autocorrelation = (4, 0, 1, 2, 0, 3, 0, 0, 3, 0, 2, 1, 0)
    require(guard_autocorrelation == expected_autocorrelation, "guard autocorrelation")
    require(
        all(guard_autocorrelation[shift] == guard_autocorrelation[-shift % ROOT_COUNT] for shift in range(ROOT_COUNT)),
        "guard overlap negation symmetry",
    )

    print("THM-2504 ENDPOINT TOURNAMENT-HOLONOMY NO-GO AUDIT")
    print(
        f"mask_pairs={len(mask_pairs)};translation_swaps={len(mask_pairs)};"
        f"equal_weight_ties={equal_weight_pairs}"
    )
    print(f"parity_oriented={parity_oriented};parity_ties={parity_ties}")
    print(
        "first_failure_layers="
        + ",".join(map(str, expected_layer_sizes))
        + f";layer_ties={layer_ties};cell_ties={cell_ties}"
    )
    print("lex_reversal_support=0;lex_cycle_support=0")
    print(
        f"tile128=edges:{len(edges_128)},laplacian_rank:{rank_128},bicycle_dim:0,"
        f"parity_cut_support:{sum(parity_cut_128)},cycle_support:{sum(cycle_128)}"
    )
    print(f"cosupport_uniform_asymmetry_ties={uniform_asymmetry_ties};sparse_loop_offdiagonal=0")
    print(
        f"root13=edges:{len(edges_13)},laplacian_rank:{rank_13},bicycle_dim:0,"
        f"cyclic_support:{sum(cyclic_reversal_13)},cut_support:{sum(cyclic_cut_13)},"
        f"cycle_support:{sum(cyclic_cycle_13)}"
    )
    print(
        f"root13_global_converse_fixed_path=cut_support:{sum(converse_cut_13)},"
        f"cycle_support:{sum(converse_cycle_13)};triangle_0_2_7={triangle_holonomy}"
    )
    print(
        f"root13_inversion_transported_path=cut_support:{sum(cyclic_cut_13)},"
        f"cycle_support:{sum(cyclic_cycle_13)}"
    )
    print("quadratic_residue_negation_symmetry=PASS")
    print("guard_autocorrelation=" + ",".join(map(str, guard_autocorrelation)) + ";symmetric=PASS")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
