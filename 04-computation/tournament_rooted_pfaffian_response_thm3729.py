#!/usr/bin/env python3
"""Exact probe for the rooted-Pfaffian form of THM-1950's response s.

For a skew tournament matrix K and root vector u this verifies

  det(I+K) u^T(I+K)^(-1)u
    = sum_{S odd} Pf([[K[S],u[S]],[-u[S]^T,0]])^2.

The left side is evaluated independently by a rank-one determinant lemma.
"""

from __future__ import annotations

from itertools import combinations, product
from hashlib import sha256


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def determinant(matrix: list[list[int]]) -> int:
    """Fraction-free Bareiss determinant."""
    n = len(matrix)
    if n == 0:
        return 1
    a = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot_value - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss exact division")
                a[i][j] = numerator // previous
        previous = pivot_value
    return sign * a[n - 1][n - 1]


def pfaffian(matrix: list[list[int]]) -> int:
    n = len(matrix)
    require(n % 2 == 0, "Pfaffian dimension")
    if n == 0:
        return 1
    total = 0
    for j in range(1, n):
        minor = [
            [matrix[r][c] for c in range(1, n) if c != j]
            for r in range(1, n)
            if r != j
        ]
        total += (-1 if j % 2 == 0 else 1) * matrix[0][j] * pfaffian(minor)
    return total


def subsets(n: int, parity: int):
    for size in range(parity, n + 1, 2):
        yield from combinations(range(n), size)


def principal(matrix: list[list[int]], indices: tuple[int, ...]) -> list[list[int]]:
    return [[matrix[i][j] for j in indices] for i in indices]


def rooted_matrix(K: list[list[int]], indices: tuple[int, ...], u: tuple[int, ...]):
    size = len(indices)
    out = [[0] * (size + 1) for _ in range(size + 1)]
    for i, vi in enumerate(indices):
        for j, vj in enumerate(indices):
            out[i][j] = K[vi][vj]
        out[i][size] = u[vi]
        out[size][i] = -u[vi]
    return out


def energies(K: list[list[int]], u: tuple[int, ...]) -> tuple[int, int]:
    n = len(K)
    even = sum(pfaffian(principal(K, S)) ** 2 for S in subsets(n, 0))
    odd = sum(pfaffian(rooted_matrix(K, S, u)) ** 2 for S in subsets(n, 1))
    return even, odd


def tournament_matrix(n: int, mask: int) -> list[list[int]]:
    K = [[0] * n for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            value = 1 if (mask >> bit) & 1 else -1
            K[i][j] = value
            K[j][i] = -value
            bit += 1
    return K


def identity(n: int) -> list[list[int]]:
    return [[int(i == j) for j in range(n)] for i in range(n)]


def add(A: list[list[int]], B: list[list[int]]) -> list[list[int]]:
    return [[a + b for a, b in zip(arow, brow)] for arow, brow in zip(A, B)]


def outer(u: tuple[int, ...]) -> list[list[int]]:
    return [[x * y for y in u] for x in u]


def det_response(K: list[list[int]], u: tuple[int, ...]) -> tuple[int, int]:
    M = add(identity(len(K)), K)
    det_m = determinant(M)
    response_numerator = determinant(add(M, outer(u))) - det_m
    return det_m, response_numerator


def switch(K: list[list[int]], signs: tuple[int, ...]) -> list[list[int]]:
    return [
        [signs[i] * K[i][j] * signs[j] for j in range(len(K))]
        for i in range(len(K))
    ]


def hamiltonian_paths(K: list[list[int]]) -> int:
    n = len(K)
    if n == 0:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and K[last][nxt] == 1:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def induced(K: list[list[int]], indices: tuple[int, ...]) -> list[list[int]]:
    return [[K[i][j] for j in indices] for i in indices]


def sub_hamiltonian_square_energies(K: list[list[int]]) -> tuple[int, int]:
    n = len(K)
    even = sum(hamiltonian_paths(induced(K, S)) ** 2 for S in subsets(n, 0))
    odd = sum(hamiltonian_paths(induced(K, S)) ** 2 for S in subsets(n, 1))
    return even, odd


def main() -> None:
    print("== rooted Pfaffian response / parity-energy probe ==")
    total_tournaments = 0
    min_slack_by_n: list[tuple[int, int, int]] = []
    min_subpath_slack_by_n: list[tuple[int, int, int]] = []
    first_subpath_hostile = None
    semantic_rows: list[tuple[int, int, int, int, int]] = []

    for n in range(1, 7):
        masks = 1 << (n * (n - 1) // 2)
        min_slack = None
        min_subpath_slack = None
        equality_count = 0
        for mask in range(masks):
            K = tournament_matrix(n, mask)
            u = (1,) * n
            even, odd = energies(K, u)
            det_m, response = det_response(K, u)
            require(even == det_m, f"even principal identity n={n} mask={mask}")
            require(odd == response, f"odd rooted identity n={n} mask={mask}")
            require(even > 0 and odd > 0, "positive energies")
            H = hamiltonian_paths(K)
            slack = (1 << (n - 1)) * H - max(even, odd)
            require(slack >= 0, f"finite H/P control n={n} mask={mask}")
            sub_even, sub_odd = sub_hamiltonian_square_energies(K)
            subpath_slack = (1 << (n - 1)) * H - max(sub_even, sub_odd)
            if subpath_slack < 0 and first_subpath_hostile is None:
                first_subpath_hostile = (
                    n,
                    mask,
                    H,
                    sub_even,
                    sub_odd,
                    even,
                    odd,
                )
            if slack == 0:
                equality_count += 1
            if min_slack is None or slack < min_slack:
                min_slack = slack
            if min_subpath_slack is None or subpath_slack < min_subpath_slack:
                min_subpath_slack = subpath_slack
            total_tournaments += 1
        require(min_slack is not None, "nonempty census")
        require(min_subpath_slack is not None, "nonempty subpath census")
        min_slack_by_n.append((n, min_slack, equality_count))
        min_subpath_slack_by_n.append((n, min_subpath_slack, masks))
        semantic_rows.append(
            (n, masks, min_slack, min_subpath_slack, equality_count, total_tournaments)
        )
        print(
            f"n={n} tournaments={masks} min_scaled_H_minus_parity_energy={min_slack} "
            f"min_scaled_H_minus_subpath_square_energy={min_subpath_slack} "
            f"equalities={equality_count}"
        )

    # A non-Boolean root vector checks the general identity independently.
    K5 = tournament_matrix(5, 0b1011010011)
    u5 = (2, -1, 3, 0, -2)
    even5, odd5 = energies(K5, u5)
    det5, response5 = det_response(K5, u5)
    require(even5 == det5 and odd5 == response5, "integer root-vector control")

    # Switching covariance needs the root vector to move with the gauge.
    signs = (1, -1, -1, 1, -1)
    K5s = switch(K5, signs)
    u5s = tuple(signs[i] * u5[i] for i in range(5))
    require(energies(K5s, u5s) == (even5, odd5), "switching covariance")
    reset_odd = energies(K5s, u5)[1]
    require(reset_odd != odd5, "fixed-root hostile must change rooted energy")
    ones5 = (1,) * 5
    sign_root_base = energies(K5, ones5)
    sign_root_covariant = energies(K5s, signs)
    sign_root_reset = energies(K5s, ones5)
    require(sign_root_covariant == sign_root_base, "sign-root switching covariance")
    require(sign_root_reset != sign_root_base, "sign-root reset hostile")

    # Averaging every incident sign word kills the cross terms in
    # u^T adj(I+K)u.  The surviving diagonal cofactors are exactly the
    # determinants of the vertex deletions.
    sign_root_odd_numerators = [
        energies(K5, tuple(root))[1]
        for root in product((-1, 1), repeat=5)
    ]
    deletion_determinants = [
        determinant(
            add(
                identity(4),
                induced(K5, tuple(j for j in range(5) if j != vertex)),
            )
        )
        for vertex in range(5)
    ]
    require(
        sum(sign_root_odd_numerators)
        == (1 << 5) * sum(deletion_determinants),
        "sign-root average equals deletion-cofactor sum",
    )
    sign_root_average_numerator = sum(sign_root_odd_numerators) // (1 << 5)
    require(
        sign_root_average_numerator == sum(deletion_determinants),
        "integral sign-root average",
    )

    # Canonical controls: transitive T5 and directed C3.
    transitive5 = tournament_matrix(5, (1 << 10) - 1)
    cycle3 = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]
    transitive_energy = energies(transitive5, (1,) * 5)
    cycle_energy = energies(cycle3, (1,) * 3)
    require(transitive_energy == (16, 16), "transitive parity means")
    require(cycle_energy == (4, 12), "C3 parity means")

    print(f"general_root_u={u5} even={even5} odd={odd5}")
    print(f"switching_covariant_odd={odd5} fixed_root_hostile_odd={reset_odd}")
    print(
        f"sign_root_base={sign_root_base} covariant={sign_root_covariant} "
        f"reset={sign_root_reset}"
    )
    print(
        f"sign_root_average_odd_numerator={sign_root_average_numerator} "
        f"deletion_determinants={deletion_determinants}"
    )
    print(f"transitive_T5_energies={transitive_energy} C3_energies={cycle_energy}")
    require(first_subpath_hostile is not None, "expected subpath-square hostile")
    print(f"first_subpath_square_hostile={first_subpath_hostile}")
    print(f"total_exhaustive_tournaments={total_tournaments}")
    digest = sha256(
        repr(
            (
                semantic_rows,
                even5,
                odd5,
                reset_odd,
                sign_root_base,
                sign_root_reset,
                sign_root_average_numerator,
                deletion_determinants,
                first_subpath_hostile,
            )
        ).encode()
    ).hexdigest()
    print(f"semantic_sha256={digest}")
    print("ALL GATES PASS")


if __name__ == "__main__":
    main()
