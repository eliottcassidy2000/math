#!/usr/bin/env python3
"""Finite-exact scout for HYP-9080.

For a tournament C and x in C, this checks the two independent forms

  chi_x = 2[(H-disc)(C) - (H-disc)(C-x)]
        = 2 Delta H_x - E_odd(C-x;u_x) + disc(C-x).

The local-unavoidability and averaged inequalities remain conjectural.  Every
gate uses ``require`` and therefore remains active under ``python -O``.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import permutations
from pathlib import Path


N7_REPS = Path("05-knowledge/results/n7_class_reps_boxeph_S152.txt")
N7_REPS_NORMALIZED_SHA256 = "151070ab6234a9fd9f449f05f93feab64a04b17ad8dae71ba7b7a19679ccd262"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def pair_bit(n: int, i: int, j: int) -> int:
    require(0 <= i < j < n, "pair-bit domain")
    return i * (2 * n - i - 1) // 2 + j - i - 1


def tournament_matrix(n: int, mask: int) -> list[list[int]]:
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            value = 1 if (mask >> pair_bit(n, i, j)) & 1 else -1
            matrix[i][j] = value
            matrix[j][i] = -value
    return matrix


def determinant(matrix: list[list[int]]) -> int:
    """Fraction-free Bareiss determinant."""
    n = len(matrix)
    if n == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(n - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, n) if work[row][pivot_index]),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, n):
            for column in range(pivot_index + 1, n):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                require(numerator % previous == 0, "Bareiss exact division")
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def principal(matrix: list[list[int]], vertices: tuple[int, ...]) -> list[list[int]]:
    return [[matrix[i][j] for j in vertices] for i in vertices]


def identity_plus(skew: list[list[int]]) -> list[list[int]]:
    return [
        [skew[i][j] + int(i == j) for j in range(len(skew))]
        for i in range(len(skew))
    ]


def discriminant(matrix: list[list[int]], vertices: tuple[int, ...]) -> int:
    size = len(vertices)
    require(size >= 1, "empty discriminant convention not used")
    det = determinant(identity_plus(principal(matrix, vertices)))
    denominator = 1 << (size - 1)
    require(det % denominator == 0, "discriminant normalization")
    return det // denominator


def hamiltonian_subset_counts(matrix: list[list[int]]) -> list[int]:
    """H for every nonempty induced labelled subset in one DP."""
    n = len(matrix)
    endpoint = [[0] * n for _ in range(1 << n)]
    totals = [0] * (1 << n)
    for vertex in range(n):
        endpoint[1 << vertex][vertex] = 1
        totals[1 << vertex] = 1
    for subset in range(1, 1 << n):
        for last in range(n):
            count = endpoint[subset][last]
            if not count:
                continue
            for nxt in range(n):
                if not (subset >> nxt) & 1 and matrix[last][nxt] == 1:
                    endpoint[subset | (1 << nxt)][nxt] += count
        totals[subset] = sum(endpoint[subset])
    return totals


def rooted_odd_energy(
    matrix: list[list[int]], deleted: int, base_vertices: tuple[int, ...]
) -> int:
    """Determinant-response evaluation independent of the deletion formula."""
    base = principal(matrix, base_vertices)
    M = identity_plus(base)
    root = tuple(matrix[vertex][deleted] for vertex in base_vertices)
    rank_one = [
        [M[i][j] + root[i] * root[j] for j in range(len(root))]
        for i in range(len(root))
    ]
    numerator = determinant(rank_one) - determinant(M)
    denominator = 1 << (len(root) - 1)
    require(numerator % denominator == 0, "rooted energy normalization")
    return numerator // denominator


def hamiltonian_paths(matrix: list[list[int]], vertices: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        path
        for path in permutations(vertices)
        if all(matrix[path[i]][path[i + 1]] == 1 for i in range(len(path) - 1))
    )


def audit_deletion_incidence(
    matrix: list[list[int]], deleted: int, delta_h: int
) -> None:
    """Directly replay the THM-4094 fibers for the small labelled universe."""
    n = len(matrix)
    base_vertices = tuple(vertex for vertex in range(n) if vertex != deleted)
    old_paths = hamiltonian_paths(matrix, base_vertices)
    new_paths = hamiltonian_paths(matrix, tuple(range(n)))
    old_set = set(old_paths)
    fiber_sizes = {path: 0 for path in old_paths}
    orphans = 0
    for path in new_paths:
        reduced = tuple(vertex for vertex in path if vertex != deleted)
        if reduced in old_set:
            fiber_sizes[reduced] += 1
        else:
            orphans += 1
    require(all(size >= 1 for size in fiber_sizes.values()), "left-total deletion incidence")
    excess = sum(size - 1 for size in fiber_sizes.values())
    require(delta_h == excess + orphans, "Hamiltonian deletion-deficit identity")


def is_strong(matrix: list[list[int]]) -> bool:
    n = len(matrix)
    for start in range(n):
        seen = {start}
        stack = [start]
        while stack:
            current = stack.pop()
            for nxt in range(n):
                if nxt not in seen and matrix[current][nxt] == 1:
                    seen.add(nxt)
                    stack.append(nxt)
        if len(seen) != n:
            return False
    return True


def evaluate_tournament(n: int, mask: int, incidence_audit: bool) -> dict[str, object]:
    matrix = tournament_matrix(n, mask)
    vertices = tuple(range(n))
    full_subset = (1 << n) - 1
    H_by_subset = hamiltonian_subset_counts(matrix)
    H_full = H_by_subset[full_subset]
    disc_full = discriminant(matrix, vertices)
    slack_full = H_full - disc_full
    chis: list[int] = []
    deletion_slacks: list[int] = []
    root_gates = 0
    incidence_gates = 0
    for deleted in range(n):
        base_vertices = tuple(vertex for vertex in vertices if vertex != deleted)
        H_deleted = H_by_subset[full_subset ^ (1 << deleted)]
        disc_deleted = discriminant(matrix, base_vertices)
        slack_deleted = H_deleted - disc_deleted
        delta_h = H_full - H_deleted
        chi = 2 * (slack_full - slack_deleted)
        require(chi % 2 == 0, "chi parity")

        odd_energy = rooted_odd_energy(matrix, deleted, base_vertices)
        require(odd_energy == 2 * disc_full - disc_deleted, "rooted deletion response")
        boundary_chi = 2 * delta_h - odd_energy + disc_deleted
        require(boundary_chi == chi, "two chi formulas")
        require(delta_h >= 0, "Hamiltonian deletion monotonicity")
        root_gates += 2

        if incidence_audit:
            audit_deletion_incidence(matrix, deleted, delta_h)
            incidence_gates += 1
        chis.append(chi)
        deletion_slacks.append(slack_deleted)

    average_halfcharge = sum(chis) // 2
    require(2 * average_halfcharge == sum(chis), "average charge parity")
    require(
        average_halfcharge == n * slack_full - sum(deletion_slacks),
        "global deletion-charge identity",
    )
    return {
        "H": H_full,
        "disc": disc_full,
        "slack": slack_full,
        "chis": tuple(chis),
        "average_halfcharge": average_halfcharge,
        "strong": is_strong(matrix),
        "root_gates": root_gates,
        "incidence_gates": incidence_gates,
    }


def census(n: int, masks: tuple[int, ...], incidence_audit: bool) -> dict[str, object]:
    all_negative = 0
    negative_individuals = 0
    first_negative: tuple[int, int] | None = None
    min_chi: int | None = None
    min_chi_witness: tuple[int, int] | None = None
    min_average: int | None = None
    min_best: int | None = None
    min_best_witness: int | None = None
    strong_count = 0
    strong_min_best: int | None = None
    strong_min_best_witness: tuple[int, int, int, tuple[int, ...]] | None = None
    root_gates = 0
    incidence_gates = 0

    for mask in masks:
        row = evaluate_tournament(n, mask, incidence_audit)
        chis = row["chis"]
        require(isinstance(chis, tuple), "chi tuple")
        best = max(chis)
        worst = min(chis)
        all_negative += best < 0
        negative_individuals += sum(chi < 0 for chi in chis)
        if worst < 0 and first_negative is None:
            first_negative = (mask, chis.index(worst))
        if min_chi is None or worst < min_chi:
            min_chi = worst
            min_chi_witness = (mask, chis.index(worst))
        average = row["average_halfcharge"]
        require(isinstance(average, int), "average integer")
        if min_average is None or average < min_average:
            min_average = average
        if min_best is None or best < min_best:
            min_best = best
            min_best_witness = mask
        if row["strong"]:
            strong_count += 1
            if strong_min_best is None or best < strong_min_best:
                strong_min_best = best
                strong_min_best_witness = (
                    mask,
                    row["H"],
                    row["disc"],
                    chis,
                )
        root_gates += int(row["root_gates"])
        incidence_gates += int(row["incidence_gates"])

    require(min_chi is not None and min_average is not None and min_best is not None, "nonempty census")
    return {
        "count": len(masks),
        "all_negative": all_negative,
        "negative_individuals": negative_individuals,
        "first_negative": first_negative,
        "min_chi": min_chi,
        "min_chi_witness": min_chi_witness,
        "min_average": min_average,
        "min_best": min_best,
        "min_best_witness": min_best_witness,
        "strong_count": strong_count,
        "strong_min_best": strong_min_best,
        "strong_min_best_witness": strong_min_best_witness,
        "root_gates": root_gates,
        "incidence_gates": incidence_gates,
    }


def normalized_n7_reps() -> tuple[int, ...]:
    raw = N7_REPS.read_bytes()
    normalized = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    require(sha256(normalized).hexdigest() == N7_REPS_NORMALIZED_SHA256, "n7 rep hash")
    reps = tuple(int(line) for line in normalized.decode("ascii").splitlines() if line.strip())
    require(len(reps) == 456 and len(set(reps)) == 456, "n7 representative universe")
    return reps


def main() -> None:
    print("HYP-9080 tournament deletion-slack local-unavoidability scout")
    total_root_gates = 0
    total_incidence_gates = 0
    labelled_rows: dict[int, dict[str, object]] = {}
    for n in range(2, 7):
        masks = tuple(range(1 << (n * (n - 1) // 2)))
        row = census(n, masks, incidence_audit=n <= 5)
        labelled_rows[n] = row
        total_root_gates += int(row["root_gates"])
        total_incidence_gates += int(row["incidence_gates"])
        print(
            f"labelled n={n} count={row['count']} all_negative={row['all_negative']} "
            f"negative_individuals={row['negative_individuals']} min_chi={row['min_chi']} "
            f"min_chi_witness={row['min_chi_witness']} "
            f"min_average_halfcharge={row['min_average']} min_best_chi={row['min_best']}"
        )

    require(labelled_rows[4]["first_negative"] == (2, 3), "first negative chi witness")
    witness = evaluate_tournament(4, 2, incidence_audit=True)
    require(
        (witness["H"], witness["disc"], witness["slack"], witness["chis"])
        == (3, 2, 1, (2, 2, 2, -2)),
        "source-over-C3 hostile",
    )
    require(labelled_rows[6]["min_chi_witness"] == (72, 5), "order-six minimum witness")

    reps = normalized_n7_reps()
    n7 = census(7, reps, incidence_audit=False)
    total_root_gates += int(n7["root_gates"])
    total_incidence_gates += int(n7["incidence_gates"])
    print(
        f"n7_isomorphism_reps count={n7['count']} strong={n7['strong_count']} "
        f"all_negative={n7['all_negative']} negative_individuals={n7['negative_individuals']} "
        f"min_chi={n7['min_chi']} min_chi_witness={n7['min_chi_witness']} "
        f"min_average_halfcharge={n7['min_average']} min_best_chi={n7['min_best']}"
    )
    require(n7["strong_count"] == 353, "order-seven strong class count")
    require(n7["min_chi_witness"] == (664, 6), "order-seven minimum witness")
    require(
        n7["strong_min_best_witness"] == (171, 27, 4, (44, 34, 20, 44, 20, 20, 34)),
        "strong order-seven minimum best certificate",
    )
    print(
        "n7_strong_min_best_chi=44 witness_mask=171 H=27 disc=4 slack=23 "
        "chis=44,34,20,44,20,20,34"
    )

    require(all(row["all_negative"] == 0 for row in labelled_rows.values()), "labelled local unavoidability")
    require(all(row["min_average"] >= 0 for row in labelled_rows.values()), "labelled average charge")
    require(n7["all_negative"] == 0 and n7["min_average"] >= 0, "n7 finite frontier")
    print("first_negative_individual=n4 mask=2 x=3 H_disc_slack=3,2,1 deletion_slack=2 chi=-2")
    print(f"rooted_response_gates={total_root_gates}")
    print(f"full_deletion_incidence_gates={total_incidence_gates} through_n=5")
    print("result=FINITE-EXACT_PASS")
    print("local_unavoidability=OPEN average_form=OPEN H_ge_disc=OPEN")


if __name__ == "__main__":
    main()
