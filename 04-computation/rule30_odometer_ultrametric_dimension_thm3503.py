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
    # Exhaust every unordered pair in the largest audited phase quotient.
    # This simultaneously checks every available dyadic valuation and every
    # nontrivial odd multiplier, rather than duplicating antipodal pairs.
    period = periods[-1]
    require(period & (period - 1) == 0, "largest phase period is not dyadic")
    mask = (1 << MAX_WIDTH) - 1
    rows: list[int] = []
    row = 1
    for _ in range(period):
        rows.append(row)
        row = phi(row, mask)
    require(row == 1, "largest phase orbit did not close")

    checks = 0
    for s in range(period):
        for t in range(s + 1, period):
            m = valuation2(t - s)
            require(m < len(innovations), "missing innovation depth")
            difference = rows[t] ^ rows[s]
            require(difference != 0, "distinct phases collapsed")
            require(
                valuation2(difference) == innovations[m],
                "exact all-pairs metric regrading failed",
            )
            checks += 1
    require(checks == period * (period - 1) // 2, "pair universe is incomplete")

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

    sparse_positions = [i for i, bit in enumerate(sparse) if bit]
    dense_positions = [i for i, bit in enumerate(dense) if bit]
    for positions in (sparse_positions, dense_positions):
        # The r-th innovation obeys the exact dyadic strict ceiling inherited
        # from the period floor: kappa_r < 2^r.
        for r, position in enumerate(positions, start=1):
            require(position < (1 << r), "synthetic tower violates period floor")

    sparse_count = len(sparse_positions)
    dense_count = sum(dense)
    require(sparse_count == 12, "sparse tower count mismatch")
    require(
        abs(3 * dense_count - 2 * horizon) < 3,
        "dense tower density mismatch",
    )
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
