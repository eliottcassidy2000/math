#!/usr/bin/env python3
"""Exact finite companion for THM-2703.

The theorem is a general representation/lattice argument.  This script
exhausts small weighted fixed-subtree/three-arm templates, checks the mod-two
standard-block formula directly, audits continuant parity, and verifies the
sharp D4 and odd-arm controls.  No floating point arithmetic is used.
"""

from __future__ import annotations

from itertools import product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def gf2_rank(rows: list[int], columns: int) -> int:
    work = rows[:]
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(work))
                      if (work[i] >> column) & 1), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for i in range(len(work)):
            if i != rank and ((work[i] >> column) & 1):
                work[i] ^= work[rank]
        rank += 1
    return rank


def gf2_rows(matrix: list[list[int]]) -> list[int]:
    return [sum((entry & 1) << j for j, entry in enumerate(row))
            for row in matrix]


def gf2_nullity(matrix: list[list[int]]) -> int:
    return len(matrix[0]) - gf2_rank(gf2_rows(matrix), len(matrix[0]))


def chain_matrix(weights: tuple[int, ...]) -> list[list[int]]:
    n = len(weights)
    matrix = [[0] * n for _ in range(n)]
    for i, weight in enumerate(weights):
        matrix[i][i] = -weight
        if i:
            matrix[i - 1][i] = 1
            matrix[i][i - 1] = 1
    return matrix


def build_tree(
    fixed_weights: tuple[int, ...],
    arm_specs: tuple[tuple[int, tuple[int, ...]], ...],
) -> tuple[list[list[int]], list[int], tuple[tuple[int, ...], ...]]:
    """Build a fixed path with three cyclic copies of every rooted chain."""
    size = len(fixed_weights) + 3 * sum(len(weights) for _, weights in arm_specs)
    matrix = [[0] * size for _ in range(size)]
    permutation = list(range(size))

    for i, weight in enumerate(fixed_weights):
        matrix[i][i] = -weight
        if i:
            matrix[i - 1][i] = matrix[i][i - 1] = 1

    cursor = len(fixed_weights)
    representative_arms: list[tuple[int, ...]] = []
    for attachment, weights in arm_specs:
        representative_arms.append(weights)
        copies: list[list[int]] = []
        for _ in range(3):
            indices = list(range(cursor, cursor + len(weights)))
            cursor += len(weights)
            copies.append(indices)
            for j, (index, weight) in enumerate(zip(indices, weights)):
                matrix[index][index] = -weight
                if j:
                    matrix[indices[j - 1]][index] = 1
                    matrix[index][indices[j - 1]] = 1
            matrix[attachment][indices[0]] = 1
            matrix[indices[0]][attachment] = 1
        for copy in range(3):
            for position in range(len(weights)):
                permutation[copies[copy][position]] = copies[(copy + 1) % 3][position]

    require(cursor == size, "tree cursor")
    require(sorted(permutation) == list(range(size)), "tree permutation")
    return matrix, permutation, tuple(representative_arms)


def standard_kernel_dimension(matrix: list[list[int]], permutation: list[int]) -> int:
    """Dimension of ker(M) intersect ker(1+sigma+sigma^2) over F2."""
    n = len(matrix)
    inverse = [0] * n
    for i, image in enumerate(permutation):
        inverse[image] = i
    inverse2 = [inverse[inverse[i]] for i in range(n)]
    norm_rows = []
    for row in range(n):
        mask = (1 << row) ^ (1 << inverse[row]) ^ (1 << inverse2[row])
        norm_rows.append(mask)
    stacked = gf2_rows(matrix) + norm_rows
    return n - gf2_rank(stacked, n)


def verify_template(
    fixed_weights: tuple[int, ...],
    arm_specs: tuple[tuple[int, tuple[int, ...]], ...],
) -> int:
    matrix, permutation, arms = build_tree(fixed_weights, arm_specs)
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            require(matrix[permutation[i]][permutation[j]] == matrix[i][j],
                    "C3 matrix invariance")
    observed = standard_kernel_dimension(matrix, permutation)
    predicted = 2 * sum(gf2_nullity(chain_matrix(weights)) for weights in arms)
    require(observed == predicted, "standard arm-block multiplicity")
    return observed


def continuant(weights: tuple[int, ...]) -> int:
    d_previous_previous = 1
    d_previous = weights[0]
    for weight in weights[1:]:
        d_previous_previous, d_previous = (
            d_previous,
            weight * d_previous - d_previous_previous,
        )
    return d_previous


def mat2_mul(a: tuple[int, int, int, int],
             b: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return (
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    )


def main() -> None:
    template_count = 0
    positive_standard_count = 0

    # One fixed vertex and one freely moving arm orbit.
    for center_weight in range(1, 5):
        for length in range(1, 6):
            for weights in product(range(1, 5), repeat=length):
                observed = verify_template((center_weight,), ((0, weights),))
                template_count += 1
                positive_standard_count += observed > 0

    # A nontrivial fixed subtree and two independently moving arm orbits.
    arm_bank = tuple(
        weights
        for length in range(1, 3)
        for weights in product(range(1, 4), repeat=length)
    )
    for fixed_weights in product(range(1, 4), repeat=2):
        for left_arm in arm_bank:
            for right_arm in arm_bank:
                observed = verify_template(
                    fixed_weights,
                    ((0, left_arm), (1, right_arm)),
                )
                template_count += 1
                positive_standard_count += observed > 0

    # Negative continued-fraction numerators are the chain determinants.
    continuant_count = 0
    even_continuant_count = 0
    for length in range(1, 7):
        for weights in product(range(2, 7), repeat=length):
            determinant = continuant(weights)
            transfer = (1, 0, 0, 1)
            for weight in weights:
                transfer = mat2_mul(transfer, (weight, -1, 1, 0))
            require(transfer[0] == determinant, "continuant transfer numerator")
            direct = int(sp.Matrix(chain_matrix(weights)).det())
            require(abs(direct) == determinant, "chain determinant")
            require((determinant % 2 == 0) ==
                    (gf2_nullity(chain_matrix(weights)) > 0),
                    "even continuant criterion")
            continuant_count += 1
            even_continuant_count += determinant % 2 == 0

    for length in range(1, 13):
        require(continuant((2,) * length) == length + 1, "A-chain determinant")
        require((continuant((2,) * length) % 2 == 0) == (length % 2 == 1),
                "A-chain parity")

    # D4: one fixed central node and a triple of moving -2 arms.
    d4_matrix, d4_permutation, _ = build_tree((2,), ((0, (2,)),))
    d4 = sp.Matrix(d4_matrix)
    d4_smith = smith_normal_form(d4, domain=sp.ZZ)
    d4_diagonal = tuple(abs(int(d4_smith[i, i])) for i in range(4))
    require(int(d4.det()) == 4, "D4 determinant")
    require(d4_diagonal == (1, 1, 2, 2), "D4 Smith form")
    require(standard_kernel_dimension(d4_matrix, d4_permutation) == 2,
            "D4 one standard plane")

    # Same three-arm shape with odd arms: no standard plane.
    odd_matrix, odd_permutation, _ = build_tree((2,), ((0, (3,)),))
    odd = sp.Matrix(odd_matrix)
    require(abs(int(odd.det())) == 27, "odd-arm determinant control")
    require(standard_kernel_dimension(odd_matrix, odd_permutation) == 0,
            "odd-arm triality has no standard plane")

    print("THM-2703 C3 BOUNDARY-TREE ARM AUDIT")
    print(f"weighted_tree_templates_checked={template_count}")
    print(f"templates_with_standard_kernel={positive_standard_count}")
    print(f"continuant_chains_checked={continuant_count}")
    print(f"even_continuant_chains={even_continuant_count}")
    print("standard_multiplicity=sum_arm_nullities")
    print(f"d4_det={int(d4.det())} d4_smith={d4_diagonal} standard_dim=2")
    print(f"odd_arm_det={abs(int(odd.det()))} standard_dim=0")
    print("all_minus_two_chain_det=length+1 parity_even_iff_length_odd")
    print("surface_corollary_requires_full_rank_boundary_lattice")
    print("unit_squareclass_branch_not_tested")
    print("PASS")


if __name__ == "__main__":
    main()
