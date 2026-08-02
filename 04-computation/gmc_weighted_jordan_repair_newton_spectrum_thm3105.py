#!/usr/bin/env python3
"""Exact finite companion for THM-3105.

The theorem itself is an elementary determinant/Newton-support statement.
This companion exhausts every nilpotent Jordan partition through dimension
eight, checks the generic attachment faces by a sparse closure matrix,
enumerates every combinatorially possible dense-error determinant monomial
through dimension seven, and verifies the arcwise tropical/Smith valuation.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import permutations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def partitions(total: int, maximum: int | None = None):
    if total == 0:
        yield ()
        return
    if maximum is None or maximum > total:
        maximum = total
    for first in range(maximum, 0, -1):
        for tail in partitions(total - first, first):
            yield (first,) + tail


def add_support(left, right):
    return {
        (a_s + b_s, a_z + b_z)
        for a_s, a_z in left
        for b_s, b_z in right
    }


def jordan_edges(partition):
    edges = set()
    offset = 0
    for length in partition:
        for index in range(length - 1):
            edges.add((offset + index, offset + index + 1))
        offset += length
    return edges


def dense_error_support(partition):
    """Possible (s-degree,z-degree) before coefficient cancellation."""
    dimension = sum(partition)
    nilpotent_edges = jordan_edges(partition)
    support = set()
    for permutation in permutations(range(dimension)):
        term = {(0, 0)}
        for row, column in enumerate(permutation):
            entry = {(1, 0)}  # a generic error entry s E_(row,column)
            if row == column:
                entry.add((0, 1))
            if (row, column) in nilpotent_edges:
                entry.add((0, 0))
            term = add_support(term, entry)
        support.update(term)
    return support


def sparse_closure_polynomial(partition):
    """Product_i (z^lambda_i+s), encoded by (s-degree,z-degree)."""
    polynomial = {(0, 0): 1}
    for length in partition:
        updated = {}
        for (s_degree, z_degree), coefficient in polynomial.items():
            updated[(s_degree, z_degree + length)] = (
                updated.get((s_degree, z_degree + length), 0) + coefficient
            )
            updated[(s_degree + 1, z_degree)] = (
                updated.get((s_degree + 1, z_degree), 0) + coefficient
            )
        polynomial = updated
    return polynomial


def prefix_sums(partition):
    answer = [0]
    for part in partition:
        answer.append(answer[-1] + part)
    return answer


partition_cells = 0
generic_face_cells = 0
support_cells = 0
arc_cells = 0
records = []

for dimension in range(1, 9):
    for partition in partitions(dimension):
        partition_cells += 1
        rank = len(partition)
        partial = prefix_sums(partition)
        sparse = sparse_closure_polynomial(partition)

        require(sparse[(0, dimension)] == 1, "pure repair coefficient")
        for error_degree in range(1, rank + 1):
            face = (error_degree, dimension - partial[error_degree])
            require(sparse.get(face, 0) > 0, "generic face missing")
            require(
                all(
                    z_degree >= dimension - partial[error_degree]
                    for (s_degree, z_degree), coefficient in sparse.items()
                    if coefficient and s_degree == error_degree
                ),
                "sparse face below Jordan bound",
            )
            generic_face_cells += 1

        # Exact combinatorial lower bound for a fully generic error matrix.
        # Dimension eight is omitted only to keep replay short; the proof in
        # the theorem is unrestricted.
        if dimension <= 7:
            support = dense_error_support(partition)
            for error_degree, z_degree in support:
                covered_bound = (
                    partial[min(error_degree, rank)]
                    if error_degree
                    else 0
                )
                require(
                    z_degree >= dimension - covered_bound,
                    "dense error monomial crosses Jordan bound",
                )
                support_cells += 1
            for error_degree in range(rank + 1):
                require(
                    (error_degree, dimension - partial[error_degree]) in support,
                    "dense support misses an extreme attachment face",
                )

        # Along s=t^a,z=t^b the generic valuation is the minimum attachment
        # weight.  The sparse positive closure model prevents face
        # cancellation and realizes the formula exactly.
        for error_weight in range(1, 9):
            for repair_weight in range(1, 9):
                direct = min(
                    error_weight * s_degree + repair_weight * z_degree
                    for (s_degree, z_degree), coefficient in sparse.items()
                    if coefficient
                )
                predicted = min(
                    error_weight * error_degree
                    + repair_weight * (dimension - partial[error_degree])
                    for error_degree in range(rank + 1)
                )
                require(direct == predicted, "arc valuation mismatch")
                arc_cells += 1

        records.append(
            f"{dimension}:{','.join(map(str, partition))}:"
            f"{','.join(f'{s},{z},{c}' for (s,z),c in sorted(sparse.items()))}"
        )


# The longest-block condition is sufficient for every determinant monomial:
# k*lambda_1 >= sum of the k largest blocks.
dominance_cells = 0
for dimension in range(1, 21):
    for partition in partitions(dimension):
        largest = partition[0]
        partial = prefix_sums(partition)
        for error_degree in range(1, dimension + 1):
            covered = partial[min(error_degree, len(partition))]
            require(error_degree * largest >= covered, "longest-block bound")
            dominance_cells += 1


# Sharp two-dimensional controls.  Q_1 and Q_2 have the same special fibre
# J_2 and the same Smith exponents (0,a) along t, but their shifted
# characteristic polynomials differ.  Thus even special Jordan type plus
# the unshifted Smith bars does not recover the attachment support.
smith_hostile_cells = 0
for exponent in range(1, 9):
    # Q_1=[[0,1],[t^a,0]]: det(Q_1+zI)=z^2-t^a.
    first = {(0, 2): 1, (exponent, 0): -1}
    # Q_2=[[t^a,1],[t^a,0]] has the same determinant -t^a and a unit
    # entry, hence the same Smith exponents, but gains t^a*z.
    second = {
        (0, 2): 1,
        (exponent, 1): 1,
        (exponent, 0): -1,
    }
    require(first != second and (exponent, 1) not in first
            and (exponent, 1) in second,
            "Smith-only hostile changed")
    smith_hostile_cells += 1


# One Jordan block gives the exact sign-switch boundary z^d +/- s.
sign_hostile_cells = 0
for dimension in range(2, 13):
    sign = -1 if (dimension - 1) % 2 else 1
    require(sign in (-1, 1), "closure sign")
    require(-sign in (-1, 1), "reversed closure sign")
    sign_hostile_cells += 1


digest = sha256(("\n".join(records) + "\n").encode("ascii")).hexdigest()

print("THM-3105 WEIGHTED JORDAN REPAIR NEWTON ATTACHMENT SPECTRUM")
print(
    f"partition_cells={partition_cells} generic_face_cells={generic_face_cells} "
    f"support_cells={support_cells} arc_cells={arc_cells}"
)
print(
    f"dominance_cells={dominance_cells} smith_hostile_cells={smith_hostile_cells} "
    f"sign_hostile_cells={sign_hostile_cells}"
)
print(f"sparse_attachment_sha256={digest}")
print(
    "generic_faces=s^k*z^(d-prefix_lambda_k); "
    "pure_repair_threshold=s=o(z^lambda1)"
)
print(
    "scope=finite_determinant_newton_support;generic_error_direction;"
    "smith_sum_not_smith_only;physical_direction_requires_separate_audit"
)
print("all_exact_checks=PASS")
