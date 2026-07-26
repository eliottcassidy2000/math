#!/usr/bin/env python3
"""Exact algebra and finite-clock audit for THM-2449.

The Fourier checks use the tensor cyclotomic basis
  1,zeta_7,...,zeta_7^5  times  1,zeta_13,...,zeta_13^11.
No floating-point approximation or external package is used.
"""

from fractions import Fraction as F
from math import gcd
from random import Random


P = 7
Q = 13
Interval = tuple[F, F]
IntervalSet = tuple[Interval, ...]
Matrix = tuple[tuple[F, ...], ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def anchored_matrix(a: F, w: tuple[F, ...]) -> Matrix:
    require(len(w) == Q and w[0] == 0, "anchored profile")
    return tuple(
        tuple((a + w[s]) if ell == 0 else w[s] for s in range(Q))
        for ell in range(P)
    )


def rectangle(
    matrix: Matrix,
    ell: int,
    s: int,
) -> F:
    return (
        matrix[0][s]
        - matrix[ell][s]
        - matrix[0][0]
        + matrix[ell][0]
    )


def additive(matrix: Matrix) -> bool:
    return all(
        rectangle(matrix, ell, s) == 0
        for ell in range(1, P)
        for s in range(1, Q)
    )


def mixed_remainder(
    matrix: Matrix,
    kappa: int,
    b: int,
) -> tuple[tuple[F, ...], ...]:
    """Reduce the mixed Fourier numerator in the 6 by 12 tensor basis."""

    require(1 <= kappa < P and 1 <= b < Q, "charged colours")
    coefficients = [[F(0) for _ in range(Q)] for _ in range(P)]
    for ell in range(P):
        for s in range(Q):
            coefficients[(kappa * ell) % P][(b * s) % Q] += matrix[ell][s]

    # zeta_p^(p-1)=-(1+...+zeta_p^(p-2)) in each tensor factor.
    return tuple(
        tuple(
            coefficients[i][j]
            - coefficients[P - 1][j]
            - coefficients[i][Q - 1]
            + coefficients[P - 1][Q - 1]
            for j in range(Q - 1)
        )
        for i in range(P - 1)
    )


def mixed_zero(matrix: Matrix, kappa: int, b: int) -> bool:
    return all(
        value == 0
        for row in mixed_remainder(matrix, kappa, b)
        for value in row
    )


def add_cell(matrix: Matrix, ell: int, s: int, value: F) -> Matrix:
    return tuple(
        tuple(
            entry + (value if (row, column) == (ell, s) else 0)
            for column, entry in enumerate(entries)
        )
        for row, entries in enumerate(matrix)
    )


def scale_add(first: Matrix, second: Matrix, scalar: F) -> Matrix:
    return tuple(
        tuple(
            first[ell][s] + scalar * second[ell][s]
            for s in range(Q)
        )
        for ell in range(P)
    )


def grid_mask_source(denominator: int, mask: int) -> IntervalSet:
    intervals: list[Interval] = []
    index = 0
    while index < denominator:
        if not (mask >> index) & 1:
            index += 1
            continue
        start = index
        while index < denominator and (mask >> index) & 1:
            index += 1
        intervals.append((F(start, denominator), F(index, denominator)))
    return tuple(intervals)


def measure(intervals: IntervalSet) -> F:
    return sum((right - left for left, right in intervals), F(0))


def overlap_length(
    first: IntervalSet,
    second: IntervalSet,
    clock: int,
) -> F:
    """Integral 1_first(x) 1_second({clock*x}) dx by literal branches."""

    total = F(0)
    for left, right in first:
        for prefix in range(clock):
            for target_left, target_right in second:
                pull_left = (prefix + target_left) / clock
                pull_right = (prefix + target_right) / clock
                total += max(
                    F(0),
                    min(right, pull_right) - max(left, pull_left),
                )
    return total


class Counts:
    hostile_profiles = 0
    hostile_representative_modes = 0
    hostile_all_mode_checks = 0
    one_defect_controls = 0
    one_defect_modes = 0
    random_matrices = 0
    mixing_pairs = 0


def audit_matrix_dichotomy() -> None:
    # CRT identifies the Galois group with all independently charged pairs.
    orbit = {
        (unit % P, unit % Q)
        for unit in range(1, P * Q)
        if gcd(unit, P * Q) == 1
    }
    require(
        orbit
        == {
            (kappa, b)
            for kappa in range(1, P)
            for b in range(1, Q)
        },
        "CRT Galois orbit",
    )

    # Exhaust all binary replica profiles.  These include the flat w=0
    # boundary and every nonflat binary target row.
    for bits in range(1 << (Q - 1)):
        w = (F(0),) + tuple(
            F((bits >> index) & 1) for index in range(Q - 1)
        )
        matrix = anchored_matrix(F(1 + bits % 3), w)
        require(additive(matrix), "replica interaction")
        require(
            matrix[0][0] > 0
            and all(matrix[ell][0] == 0 for ell in range(1, P)),
            "anchored replica column",
        )
        require(mixed_zero(matrix, 1, 1), "replica mixed mode")
        Counts.hostile_representative_modes += 1
        # Sixteen evenly spaced profiles replay the complete Galois orbit;
        # the one-mode exhaustive loop is the independent large hostile core.
        if bits % 257 == 0:
            for kappa in range(1, P):
                for b in range(1, Q):
                    require(mixed_zero(matrix, kappa, b), "replica orbit mode")
                    Counts.hostile_all_mode_checks += 1
        Counts.hostile_profiles += 1

    # A defect in any one of the 6*12 decisive rectangles makes every one
    # of the 72 mixed modes nonzero.
    base_w = tuple(F(0) if s == 0 else F(s % 3) for s in range(Q))
    base = anchored_matrix(F(2), base_w)
    for ell in range(1, P):
        for s in range(1, Q):
            matrix = add_cell(base, ell, s, F(1))
            require(rectangle(matrix, ell, s) == -1, "one-cell rectangle")
            require(not additive(matrix), "one-cell interaction")
            for kappa in range(1, P):
                for b in range(1, Q):
                    require(
                        not mixed_zero(matrix, kappa, b),
                        "one defect lost a charged colour",
                    )
                    Counts.one_defect_modes += 1
            Counts.one_defect_controls += 1

    # Secondary signed audit: additive iff one/every mixed mode vanishes.
    rng = Random(2_446)
    for _ in range(256):
        matrix = tuple(
            tuple(F(rng.randrange(-2, 4)) for _ in range(Q))
            for _ in range(P)
        )
        flags = [
            mixed_zero(matrix, kappa, b)
            for kappa in range(1, P)
            for b in range(1, Q)
        ]
        require(all(flags) == additive(matrix), "ANOVA/Fourier equivalence")
        require(not any(flags) or all(flags), "Galois all-or-all")
        Counts.random_matrices += 1


def audit_scalar_mixing() -> None:
    target = (
        (F(1, 11), F(4, 11)),
        (F(7, 11), F(10, 11)),
    )
    target_mass = measure(target)

    # Exhaust every union of cells through denominator seven.
    for D0 in range(1, 8):
        for mask in range(1 << D0):
            source = grid_mask_source(D0, mask)
            source_mass = measure(source)
            for N in (1, 2, 3):
                for t in (1, 2):
                    next_N = N + D0 * t
                    first = overlap_length(source, target, N)
                    second = overlap_length(source, target, next_N)
                    require(
                        N * (first - source_mass * target_mass)
                        == next_N * (second - source_mass * target_mass),
                        "grid scalar covariance",
                    )
                    Counts.mixing_pairs += 1

    # Exhaust single intervals on the 13 and 26 grids.  Here
    # R=13*N and the reduced clocks still differ by D0.
    for D0 in (1, 2):
        D = 13 * D0
        for left in range(D):
            for right in range(left + 1, D + 1):
                source = ((F(left, D), F(right, D)),)
                source_mass = measure(source)
                for N in (1, 2):
                    next_N = N + D0
                    clock = 13 * N
                    next_clock = 13 * next_N
                    first = overlap_length(source, target, clock)
                    second = overlap_length(source, target, next_clock)
                    require(
                        clock * (first - source_mass * target_mass)
                        == next_clock * (second - source_mass * target_mass),
                        "13-grid scalar covariance",
                    )
                    Counts.mixing_pairs += 1


def audit_one_exception_boundary() -> tuple[bool, bool, bool]:
    # Delta_R=1-13/R vanishes at exactly R=13.  The surrounding matrices
    # remain nonnegative for every displayed clock.
    w = (F(0),) + (F(2),) * (Q - 1)
    hostile = anchored_matrix(F(3), w)
    mean = add_cell(hostile, 1, 1, F(1))
    discrepancy = add_cell(
        tuple(tuple(F(0) for _ in range(Q)) for _ in range(P)),
        1,
        1,
        F(-13),
    )
    statuses: list[bool] = []
    for clock in (13, 13**2, 13**3):
        matrix = scale_add(mean, discrepancy, F(1, clock))
        require(all(entry >= 0 for row in matrix for entry in row), "positive tail")
        zero = mixed_zero(matrix, 1, 1)
        require(zero == additive(matrix), "exception ANOVA")
        statuses.append(zero)
    require(statuses == [True, False, False], "one exception boundary")
    return tuple(statuses)  # type: ignore[return-value]


def main() -> None:
    audit_matrix_dichotomy()
    audit_scalar_mixing()
    exception = audit_one_exception_boundary()

    print("THM-2449 exact owner delta-replica companion")
    print("galois_mixed_orbit_size=72")
    print(
        "anchored_replica_hostiles:",
        f"profiles={Counts.hostile_profiles}",
        f"representative_zero_checks={Counts.hostile_representative_modes}",
        f"all_mode_checks={Counts.hostile_all_mode_checks}",
    )
    print(
        "one_rectangle_controls:",
        f"rectangles={Counts.one_defect_controls}",
        f"mixed_nonzero_checks={Counts.one_defect_modes}",
    )
    print("secondary_random_signed_matrices:", Counts.random_matrices)
    print("exact_scalar_mixing_pairs:", Counts.mixing_pairs)
    print("single_exception_clocks_13_169_2197:", exception)
    print(
        "VERIFIED: one owner rectangle defect forces all 72 charged "
        "source-target colours; persistent failure is delta plus six replicas."
    )


if __name__ == "__main__":
    main()
