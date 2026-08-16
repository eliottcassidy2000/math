#!/usr/bin/env python3
"""Exact finite companion for THM-3503.

The universal ultrametric and dimension statements are proved in the theorem.
This companion audits their finite Rule-30 realization, the period/innovation
indexing, and the sharp scalar boundary models.  It intentionally contains no
``assert`` statements, so normal and optimized execution use identical gates.
"""

from __future__ import annotations

import hashlib


MAX_WIDTH = 24
EXPECTED_V = (1, 3, 4, 6, 7, 9, 15, 16, 24)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def phi(x: int, mask: int) -> int:
    return (x ^ ((x << 1) | (x << 2))) & mask


def valuation2(x: int) -> int:
    require(x != 0, "valuation requires a nonzero integer")
    return (x & -x).bit_length() - 1


def seed_period(width: int) -> tuple[int, list[int]]:
    require(width >= 1, "width must be positive")
    mask = (1 << width) - 1
    row = 1
    orbit = [row]
    for period in range(1, 1 << width):
        row = phi(row, mask)
        if row == 1:
            return period, orbit
        orbit.append(row)
    raise RuntimeError("seed did not return inside the finite quotient")


def audit_period_tower() -> tuple[list[int], list[int], list[int]]:
    periods: list[int] = []
    eps: list[int] = []
    distinct_checks = 0
    for width in range(1, MAX_WIDTH + 1):
        period, orbit = seed_period(width)
        periods.append(period)
        require(period & (period - 1) == 0, "seed period is not a power of two")
        require(len(set(orbit)) == period, "seed orbit is not embedded")
        distinct_checks += len(orbit)

        mask_next = (1 << (width + 1)) - 1
        row = 1
        for _ in range(period):
            row = phi(row, mask_next)
        epsilon = (row >> width) & 1
        eps.append(epsilon)
        if width < MAX_WIDTH:
            # Checked after the next period has been constructed below.
            pass

    for width in range(1, MAX_WIDTH):
        require(
            periods[width] == (periods[width - 1] << eps[width - 1]),
            "period-lift recursion failed",
        )

    innovations = [k for k, bit in enumerate(eps, start=1) if bit]
    v = innovations[: len(EXPECTED_V)]
    require(tuple(v) == EXPECTED_V, "innovation prefix mismatch")
    print(f"period_tower widths=1..{MAX_WIDTH} orbit_cells={distinct_checks}")
    print("innovation_prefix=" + ",".join(map(str, v)))
    return periods, eps, innovations


def audit_exact_metric(periods: list[int], innovations: list[int]) -> int:
    checks = 0
    # v_m=kappa_(m+1).  Keep one extra bit so the first differing coordinate
    # is visible, and exhaust every phase in the finite seed cycle.
    for m in range(len(innovations)):
        v_m = innovations[m]
        width = v_m + 1
        if width > MAX_WIDTH:
            break
        period = periods[width - 1]
        shift = 1 << m
        require(period >= 2 * shift, "innovation level has wrong phase period")
        mask = (1 << width) - 1
        row = 1
        rows = [row]
        for _ in range(period + shift):
            row = phi(row, mask)
            rows.append(row)
        for t in range(period):
            difference = rows[t + shift] ^ rows[t]
            require(difference != 0, "dyadic phase shift collapsed")
            require(valuation2(difference) == v_m, "exact metric regrading failed")
            checks += 1

        # Every odd multiple of 2^m has the same first differing state bit.
        for odd in range(1, period // shift, 2):
            for t in range(shift):
                difference = rows[t + odd * shift] ^ rows[t]
                require(valuation2(difference) == v_m, "odd-multiple metric failed")
                checks += 1

    print(f"exact_metric_checks={checks}")
    return checks


def audit_covering_counts(periods: list[int]) -> int:
    checks = 0
    max_period = periods[-1]
    mask = (1 << MAX_WIDTH) - 1
    rows = []
    row = 1
    for _ in range(max_period):
        rows.append(row)
        row = phi(row, mask)
    for width, expected in enumerate(periods, start=1):
        residues = {x & ((1 << width) - 1) for x in rows}
        require(len(residues) == expected, "covering number is not the seed period")
        checks += len(residues)
    print(f"covering_count_cells={checks}")
    return checks


def audit_scalar_boundaries() -> None:
    # The period-floor/no-111 constraints permit both extremal asymptotic
    # densities: isolated Mersenne innovations have dimension zero, while the
    # periodic word 110 has dimension 2/3.  These are abstract scalar towers,
    # not Rule-30 realizations.
    horizon = 4096
    sparse = [0] * horizon
    k = 1
    while k < horizon:
        sparse[k] = 1
        k = 2 * k + 1
    dense = [1 if i % 3 in (1, 2) else 0 for i in range(horizon)]

    for word in (sparse, dense):
        require(all(not (word[i] and word[i + 1] and word[i + 2]) for i in range(horizon - 2)),
                "synthetic tower violates no-111")
        require(sum(word) > 0, "synthetic tower has no innovations")

    sparse_count = sum(sparse)
    dense_count = sum(dense)
    require(sparse_count <= 13, "sparse tower count mismatch")
    require(abs(dense_count / horizon - 2 / 3) < 1 / horizon, "dense tower density mismatch")
    digest = hashlib.sha256(bytes(sparse) + bytes(dense)).hexdigest()
    print(f"scalar_hostiles horizon={horizon} sparse={sparse_count} dense={dense_count}")
    print(f"scalar_hostile_sha256={digest}")


def main() -> None:
    periods, _eps, innovations = audit_period_tower()
    metric_checks = audit_exact_metric(periods, innovations)
    cover_cells = audit_covering_counts(periods)
    audit_scalar_boundaries()
    require(metric_checks > 0 and cover_cells > 0, "audit universes were empty")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
