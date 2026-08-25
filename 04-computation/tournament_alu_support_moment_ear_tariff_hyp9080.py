#!/usr/bin/env python3
"""Exact support-moment and mixed-ear scout for HYP-9080.

This companion has two independent finite roles.

1.  It reconstructs the OCF support atoms by Boolean Mobius inversion of
    induced-subtournament Hamiltonian counts and the Pfaffian-square support
    atoms by exact principal determinants.  It checks the local and averaged
    deletion-charge moment identities.
2.  It builds the full Hamiltonian ear cut from Start/End/Q data and the
    discriminant ear cut from the exact inverse response R=(I-K^2)^(-1).
    It exhausts every nonconstant ear over every labelled order-six base.

The local-unavoidability, averaged-unavoidability, and combined cut-tariff
inequalities remain open.  All gates use ``require`` and survive ``-O``.
"""

from __future__ import annotations

from fractions import Fraction
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
    """H on every induced subset, with the OCF empty convention H(empty)=1."""
    n = len(matrix)
    endpoint = [[0] * n for _ in range(1 << n)]
    totals = [0] * (1 << n)
    totals[0] = 1
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


def mobius_atoms(zeta_values: list[int], n: int) -> list[int]:
    atoms = zeta_values[:]
    for bit in range(n):
        for mask in range(1 << n):
            if (mask >> bit) & 1:
                atoms[mask] -= atoms[mask ^ (1 << bit)]
    return atoms


def pfaffian_square_atoms(matrix: list[list[int]]) -> list[int]:
    n = len(matrix)
    atoms = [0] * (1 << n)
    for mask in range(1 << n):
        if mask.bit_count() % 2:
            continue
        vertices = tuple(i for i in range(n) if (mask >> i) & 1)
        value = determinant(principal(matrix, vertices))
        require(value >= 0, "skew principal determinant is a Pfaffian square")
        atoms[mask] = value
    return atoms


def support_moment_audit(n: int, mask: int) -> int:
    matrix = tournament_matrix(n, mask)
    all_mask = (1 << n) - 1
    vertices = tuple(range(n))
    h = hamiltonian_subset_counts(matrix)
    c = mobius_atoms(h, n)
    require(c[0] == 1 and all(value >= 0 for value in c), "OCF Mobius support atoms")

    p = pfaffian_square_atoms(matrix)
    pf_sum = sum(p)
    require(pf_sum % (1 << (n - 1)) == 0, "full Pfaffian normalization")
    disc_full = pf_sum // (1 << (n - 1))
    require(disc_full == discriminant(matrix, vertices), "principal-minor discriminant")

    h_deletions = [h[all_mask ^ (1 << x)] for x in range(n)]
    disc_deletions = [
        discriminant(matrix, tuple(v for v in vertices if v != x))
        for x in range(n)
    ]
    slack_full = h[all_mask] - disc_full
    deletion_slacks = [
        h_deletions[x] - disc_deletions[x]
        for x in range(n)
    ]

    ocf_moment = sum(state.bit_count() * c[state] for state in range(1 << n))
    pf_numerator = sum(
        (2 * state.bit_count() - n) * p[state]
        for state in range(1 << n)
    )
    pf_moment = Fraction(pf_numerator, 1 << (n - 1))
    require(pf_moment.denominator == 1, "Pfaffian support moment integrality")
    average_charge = ocf_moment - pf_moment
    require(
        average_charge == n * slack_full - sum(deletion_slacks),
        "averaged OCF/Pfaffian support-moment identity",
    )

    for x in range(n):
        ocf_incident = sum(c[state] for state in range(1 << n) if (state >> x) & 1)
        pf_incident = sum(p[state] for state in range(1 << n) if (state >> x) & 1)
        pf_avoiding = pf_sum - pf_incident
        local_charge = ocf_incident - Fraction(
            pf_incident - pf_avoiding,
            1 << (n - 1),
        )
        require(local_charge.denominator == 1, "local support charge integrality")
        require(
            local_charge == slack_full - deletion_slacks[x],
            "local OCF/Pfaffian support-moment identity",
        )
    return n + 2


def normalized_n7_reps() -> tuple[int, ...]:
    raw = N7_REPS.read_bytes()
    normalized = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    require(sha256(normalized).hexdigest() == N7_REPS_NORMALIZED_SHA256, "n7 rep hash")
    reps = tuple(int(line) for line in normalized.decode("ascii").splitlines() if line.strip())
    require(len(reps) == 456 and len(set(reps)) == 456, "n7 representative universe")
    return reps


def rational_inverse(matrix: list[list[int]]) -> list[list[Fraction]]:
    n = len(matrix)
    work = [
        [Fraction(matrix[i][j]) for j in range(n)]
        + [Fraction(int(i == j)) for j in range(n)]
        for i in range(n)
    ]
    for column in range(n):
        pivot = next((row for row in range(column, n) if work[row][column]), None)
        require(pivot is not None, "invertible response matrix")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(n):
            if row == column:
                continue
            scale = work[row][column]
            if scale:
                work[row] = [
                    work[row][j] - scale * work[column][j]
                    for j in range(2 * n)
                ]
    return [row[n:] for row in work]


def discriminant_response(matrix: list[list[int]]) -> list[list[Fraction]]:
    n = len(matrix)
    response = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            response[i][j] = int(i == j) - sum(
                matrix[i][k] * matrix[k][j]
                for k in range(n)
            )
    inverse = rational_inverse(response)
    require(
        all(inverse[i][j] == inverse[j][i] for i in range(n) for j in range(n)),
        "symmetric discriminant response",
    )
    return inverse


def ear_hamiltonian_response(
    matrix: list[list[int]],
) -> tuple[int, list[list[Fraction]], list[Fraction]]:
    n = len(matrix)
    start = [0] * n
    end = [0] * n
    q = [[0] * n for _ in range(n)]
    for word in permutations(range(n)):
        bad = [
            i
            for i in range(n - 1)
            if matrix[word[i]][word[i + 1]] != 1
        ]
        if not bad:
            start[word[0]] += 1
            end[word[-1]] += 1
            for i in range(n - 1):
                q[word[i]][word[i + 1]] += 1
        elif len(bad) == 1:
            i = bad[0]
            q[word[i]][word[i + 1]] += 1
    h_base = sum(start)
    require(h_base == sum(end), "Hamiltonian endpoint totals")
    w = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    field = [Fraction(0) for _ in range(n)]
    for i in range(n):
        for j in range(n):
            w[i][j] = Fraction(q[i][j] + q[j][i], 2)
        column = sum(q[j][i] for j in range(n))
        row = sum(q[i][j] for j in range(n))
        field[i] = Fraction(start[i] - end[i]) + Fraction(column - row, 2)
    require(sum(field) == 0, "ear orientation field has zero sum")
    return h_base, w, field


def extend_by_ear(matrix: list[list[int]], subset: int) -> list[list[int]]:
    n = len(matrix)
    extended = [row[:] + [0] for row in matrix] + [[0] * (n + 1)]
    for vertex in range(n):
        root_sign = 1 - 2 * ((subset >> vertex) & 1)
        extended[vertex][n] = root_sign
        extended[n][vertex] = -root_sign
    return extended


def ear_packet(
    matrix: list[list[int]],
) -> tuple[int, list[list[Fraction]], list[Fraction], Fraction]:
    n = len(matrix)
    vertices = tuple(range(n))
    h_base, w, field = ear_hamiltonian_response(matrix)
    d = discriminant(matrix, vertices)
    response = discriminant_response(matrix)
    one_response = sum(response[i][j] for i in range(n) for j in range(n))
    g_empty = Fraction(d, 2) * (1 + one_response)
    require(g_empty.denominator == 1, "constant-ear discriminant integrality")
    a0 = g_empty - d
    combined = [
        [w[i][j] + 2 * d * response[i][j] for j in range(n)]
        for i in range(n)
    ]
    require(
        all(combined[i][j] == combined[j][i] for i in range(n) for j in range(n)),
        "combined ear weights symmetric",
    )
    return h_base, combined, field, a0


def cut_value(weights: list[list[Fraction]], subset: int) -> Fraction:
    n = len(weights)
    return sum(
        weights[i][j]
        for i in range(n)
        for j in range(i + 1, n)
        if ((subset >> i) & 1) != ((subset >> j) & 1)
    )


def ear_halfcharge(
    combined: list[list[Fraction]],
    field: list[Fraction],
    a0: Fraction,
    subset: int,
) -> Fraction:
    return cut_value(combined, subset) + sum(
        field[i] for i in range(len(field)) if (subset >> i) & 1
    ) - a0


def direct_ear_audit(base_order: int) -> int:
    gates = 0
    for mask in range(1 << (base_order * (base_order - 1) // 2)):
        matrix = tournament_matrix(base_order, mask)
        vertices = tuple(range(base_order))
        h_base, combined, field, a0 = ear_packet(matrix)
        d_base = discriminant(matrix, vertices)
        for subset in range(1 << base_order):
            extended = extend_by_ear(matrix, subset)
            full_vertices = tuple(range(base_order + 1))
            h_full = hamiltonian_subset_counts(extended)[(1 << (base_order + 1)) - 1]
            disc_full = discriminant(extended, full_vertices)
            direct = (h_full - disc_full) - (h_base - d_base)
            formula = ear_halfcharge(combined, field, a0, subset)
            require(formula.denominator == 1 and formula == direct, "direct ear cut formula")
            gates += 1
    return gates


def exhaustive_order_six_ears() -> dict[str, object]:
    n = 6
    all_nonconstant_ears = 0
    min_combined_edge: Fraction | None = None
    edge_witness: tuple[int, int, int] | None = None
    negative_edge_bases = 0
    min_halfcharge: Fraction | None = None
    charge_witness: tuple[int, int] | None = None
    negative_charge_bases = 0
    zero_min_bases = 0
    min_tariff_gap: Fraction | None = None
    tariff_witness: tuple[int, int] | None = None

    for mask in range(1 << (n * (n - 1) // 2)):
        matrix = tournament_matrix(n, mask)
        _, combined, field, a0 = ear_packet(matrix)
        base_min_edge = min(combined[i][j] for i in range(n) for j in range(i + 1, n))
        if min_combined_edge is None or base_min_edge < min_combined_edge:
            min_combined_edge = base_min_edge
            edge = next(
                (i, j)
                for i in range(n)
                for j in range(i + 1, n)
                if combined[i][j] == base_min_edge
            )
            edge_witness = (mask, edge[0], edge[1])
        negative_edge_bases += base_min_edge < 0

        base_min_charge: Fraction | None = None
        base_min_gap: Fraction | None = None
        for subset in range(1, (1 << n) - 1):
            halfcharge = ear_halfcharge(combined, field, a0, subset)
            require(halfcharge.denominator == 1, "mixed-ear halfcharge integrality")
            field_sum = sum(field[i] for i in range(n) if (subset >> i) & 1)
            gap = cut_value(combined, subset) - a0 - abs(field_sum)
            require(gap.denominator == 1, "combined tariff gap integrality")
            if base_min_charge is None or halfcharge < base_min_charge:
                base_min_charge = halfcharge
            if base_min_gap is None or gap < base_min_gap:
                base_min_gap = gap
            if min_halfcharge is None or halfcharge < min_halfcharge:
                min_halfcharge = halfcharge
                charge_witness = (mask, subset)
            if min_tariff_gap is None or gap < min_tariff_gap:
                min_tariff_gap = gap
                tariff_witness = (mask, subset)
            all_nonconstant_ears += 1
        require(base_min_charge is not None and base_min_gap is not None, "nonempty mixed-ear fibre")
        negative_charge_bases += base_min_charge < 0
        zero_min_bases += base_min_charge == 0
        require(base_min_charge == base_min_gap, "complement-paired tariff equivalence")

    require(
        min_combined_edge is not None
        and edge_witness is not None
        and min_halfcharge is not None
        and charge_witness is not None
        and min_tariff_gap is not None
        and tariff_witness is not None,
        "nonempty order-six census",
    )
    return {
        "bases": 1 << 15,
        "all_nonconstant_ears": all_nonconstant_ears,
        "min_combined_edge": min_combined_edge,
        "edge_witness": edge_witness,
        "negative_edge_bases": negative_edge_bases,
        "min_halfcharge": min_halfcharge,
        "charge_witness": charge_witness,
        "negative_charge_bases": negative_charge_bases,
        "zero_min_bases": zero_min_bases,
        "min_tariff_gap": min_tariff_gap,
        "tariff_witness": tariff_witness,
    }


def c3_firewall() -> tuple[int, tuple[int, ...]]:
    matrix = tournament_matrix(3, 2)
    _, combined, field, a0 = ear_packet(matrix)
    charges = tuple(
        int(2 * ear_halfcharge(combined, field, a0, subset))
        for subset in range(1 << 3)
    )
    require(charges == (-2, 4, 4, 4, 4, 4, 4, -2), "C3 ear firewall")
    return 2, charges


def main() -> None:
    print("HYP-9080 OCF/Pfaffian support moments and mixed-ear tariff scout")
    moment_gates = 0
    moment_objects = 0
    for n in range(2, 7):
        count = 1 << (n * (n - 1) // 2)
        for mask in range(count):
            moment_gates += support_moment_audit(n, mask)
        moment_objects += count
        print(f"support_moment_labelled_n={n} count={count} status=PASS")

    reps = normalized_n7_reps()
    for mask in reps:
        moment_gates += support_moment_audit(7, mask)
    moment_objects += len(reps)
    print(f"support_moment_n7_reps={len(reps)} status=PASS")
    print(f"support_moment_objects={moment_objects} gates={moment_gates}")

    direct_gates = sum(direct_ear_audit(n) for n in range(2, 5))
    print(f"direct_ear_cut_gates={direct_gates} through_base_order=4")

    c3_mask, c3_charges = c3_firewall()
    print(
        f"c3_base_mask={c3_mask} ear_chis="
        + ",".join(str(value) for value in c3_charges)
    )

    row = exhaustive_order_six_ears()
    print(
        f"base_order=6 labelled_bases={row['bases']} "
        f"all_nonconstant_ears={row['all_nonconstant_ears']} "
        f"min_combined_edge={row['min_combined_edge']} "
        f"edge_witness={row['edge_witness']} "
        f"negative_edge_bases={row['negative_edge_bases']}"
    )
    print(
        f"min_nonconstant_chi={2 * row['min_halfcharge']} "
        f"chi_witness={row['charge_witness']} "
        f"negative_chi_bases={row['negative_charge_bases']} "
        f"zero_min_bases={row['zero_min_bases']}"
    )
    print(
        f"min_combined_tariff_gap={row['min_tariff_gap']} "
        f"tariff_witness={row['tariff_witness']}"
    )
    require(row["all_nonconstant_ears"] == 2_031_616, "order-six ear universe")
    require(row["min_combined_edge"] == Fraction(1, 2), "order-six edge minimum")
    require(row["negative_edge_bases"] == 0, "order-six combined edge sign")
    require(row["min_halfcharge"] == 0, "order-six mixed-ear charge minimum")
    require(row["negative_charge_bases"] == 0, "order-six mixed-ear nonnegativity")
    require(row["zero_min_bases"] == 10_368, "order-six zero-minimum base count")
    require(row["min_tariff_gap"] == 0, "order-six tariff gap minimum")
    print("result=FINITE-EXACT_PASS")
    print("ALU=OPEN combined_ear_tariff=OPEN LU=OPEN H_ge_disc=OPEN")


if __name__ == "__main__":
    main()
