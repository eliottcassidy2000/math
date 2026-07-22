#!/usr/bin/env python3
"""Exact referee for the mod-7 fourfold-moment obstruction.

Universe and hypotheses.  There are seven danger combs of radius ``1/14``.
Exactly four speeds are divisible by 7 and three are not.  The guard comb
``E_h`` has radius ``1/7``, with ``h`` odd and ``7`` not dividing ``h``.
Write ``H=1_Eh`` and let ``U,V`` count the divisible and nondivisible danger
combs.  Under containment, ``H=0`` implies ``U+V>=1``.

The two imported analytic facts are named explicitly: THM-1234 gives the
``44/273`` five-comb pair floor, and THM-2080 gives the reduced mixed-overlap
charge formula.  Everything specific to the new obstruction is checked here
using only Python's standard library and exact rational arithmetic:

* the seven-shift orbit identity, including its excluded null boundaries;
* the sharp top-three negative-charge census and its analytic tail;
* the twelve-subset averaging identity defining ``P``;
* all 39 admissible Boolean/count states of the pointwise majorant; and
* the integrated upper bound and its exact contradiction gap.

Reproduce with either command (the byte streams must agree):

    python 04-computation/lrc14_mod7_fourfold_moment_obstruction_codex_20260722.py
    python -O 04-computation/lrc14_mod7_fourfold_moment_obstruction_codex_20260722.py
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, product
from math import comb, gcd


SEVEN = 7
GUARD_RADIUS = F(1, 7)
DANGER_MEASURE = F(1, 7)
GUARD_MEASURE = F(2, 7)
PAIR_FLOOR = F(44, 273)
ALPHA = (
    -F(91, 12),
    -F(65, 12),
    -F(121, 20),
    -F(443, 60),
    -F(113, 12),
)
EXPECTED_SLACK_H0 = (
    (None, 0, 770, 1400),
    (0, 700, 1260, 1680),
    (364, 854, 1204, 1414),
    (434, 714, 854, 854),
    (210, 280, 210, 0),
)
EXPECTED_SLACK_H1 = (
    (0, 350, 560, 630),
    (0, 140, 140, 0),
    (630, 560, 350, 0),
    (1260, 980, 560, 0),
    (1890, 1400, 770, 0),
)


def require(condition: bool, message: str) -> None:
    """An optimization-proof replacement for ``assert``."""
    if not condition:
        raise RuntimeError(message)


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def guard_at_phase(x: F) -> int:
    return int(circle_distance(x) < GUARD_RADIUS)


def audit_seven_shift_independence() -> tuple[int, int, int, int]:
    """Audit the orbit proof of ``E[H phi(U)]=(2/7)E[phi(U)]``.

    If ``x={h t}``, the seven shifts replace it by ``x+h*k/7``.  Since
    ``h`` is invertible mod 7, exactly two of these seven phases lie in the
    open guard interval on every component of the complement of the boundary
    grid.  A 7-divisible danger comb, hence ``U``, is unchanged by the same
    shifts.  Testing the five coordinate functions on ``U in {0,...,4}``
    proves the identity for every function of U.  Boundary phases form a
    finite null set and are displayed separately rather than silently used.
    """

    orbit_checks = 0
    basis_coefficient_checks = 0
    boundary_checks = 0
    patterns: set[tuple[int, ...]] = set()

    for h_residue in range(1, SEVEN):
        require(gcd(h_residue, SEVEN) == 1, "nonzero mod-7 residue not invertible")
        for cell in range(SEVEN):
            # Midpoint of one component between consecutive guard boundaries.
            phase = F(2 * cell + 1, 14)
            flags = tuple(
                guard_at_phase(phase + F(h_residue * shift, SEVEN))
                for shift in range(SEVEN)
            )
            require(sum(flags) == 2, "seven-shift guard count changed")
            patterns.add(flags)
            orbit_checks += 1

            # U is constant on the orbit.  These formal coefficient vectors
            # test the basis 1_{U=j}; linearity then handles arbitrary phi.
            for u in range(5):
                guarded = tuple(sum(flags) if j == u else 0 for j in range(5))
                unguarded = tuple(SEVEN if j == u else 0 for j in range(5))
                require(
                    guarded
                    == tuple(GUARD_MEASURE * coefficient for coefficient in unguarded),
                    "periodic-function independence coefficient failed",
                )
                basis_coefficient_checks += 5

            boundary_phase = F(cell, SEVEN)
            boundary_count = sum(
                guard_at_phase(boundary_phase + F(h_residue * shift, SEVEN))
                for shift in range(SEVEN)
            )
            require(boundary_count == 1, "open-boundary control changed")
            boundary_checks += 1

    # Necessity control: if 7 divides h, the orbit does not move the guard.
    excluded_flags = tuple(guard_at_phase(F(1, 14)) for _ in range(SEVEN))
    require(sum(excluded_flags) == 7, "7|h necessity control unexpectedly passed")

    # The other half of the proof is exact integer periodicity:
    # q=7r gives q(t+k/7)=qt+rk, with rk integral.
    period_checks = 0
    for r in range(1, 29):
        for shift in range(SEVEN):
            increment = F(SEVEN * r * shift, SEVEN)
            require(increment.denominator == 1, "divisible comb lost 1/7 periodicity")
            period_checks += 1

    return orbit_checks, basis_coefficient_checks, boundary_checks, period_checks


def fold_correction(a: int, b: int) -> F:
    x = F(a % 14, 14)
    y = F(b % 7, 7)
    return min(x, y) + max(F(0), x + y - 1) - 2 * x * y


def reduced_charge(a: int, b: int) -> F:
    require(a > 0 and b > 0 and gcd(a, b) == 1, "charge pair is not reduced")
    return F(2, a * b) * fold_correction(a, b)


def audit_three_charge_floor() -> tuple[int, int, tuple[tuple[int, int, F], ...]]:
    """Certify the three worst distinct charges for the allowed parity class."""

    residue_corrections = tuple(
        fold_correction(a_residue, b_residue)
        for a_residue in range(14)
        for b_residue in range(7)
    )
    require(min(residue_corrections) == -F(6, 49), "fold correction minimum changed")
    require(min(residue_corrections) >= -F(1, 8), "analytic correction floor failed")

    # For h odd and 7 not dividing h, reduction (h,q)=g(a,b) makes a odd.
    # For a nondivisible q it also makes 7 divide neither a nor b.  Distinct q
    # give distinct reduced pairs.  Enumerate the complete ab<=20 head.
    rows: list[tuple[int, int, F]] = []
    reduced_count = 0
    for a in range(1, 21):
        for b in range(1, 20 // a + 1):
            if a % 2 == 0 or a % 7 == 0 or b % 7 == 0 or gcd(a, b) != 1:
                continue
            reduced_count += 1
            charge = reduced_charge(a, b)
            if charge < 0:
                rows.append((a, b, -charge))
    rows.sort(key=lambda row: (-row[2], row[0], row[1]))

    expected_top = (
        (1, 6, F(5, 294)),
        (11, 1, F(8, 539)),
        (1, 5, F(3, 245)),
    )
    require(tuple(rows[:3]) == expected_top, "top-three charge rows changed")
    require(rows[3][2] == expected_top[2][2], "third-place tie disappeared")

    # For ab>=21, -charge <= 1/(4ab) <= 1/84, below third place.
    require(F(1, 84) < expected_top[2][2], "analytic tail can enter the top three")
    return reduced_count, len(rows), tuple(rows)


def pair_average_polynomial(u: int, v: int) -> F:
    return F(comb(u, 2), 2) + F(u * v, 2) + F(comb(v, 2), 3)


def audit_five_subset_average() -> tuple[int, int]:
    """Audit P as the average pair count over all 3+2 five-subsets."""

    divisible = tuple(range(4))
    nondivisible = tuple(range(4, 7))
    subsets = tuple(
        first + second
        for first in combinations(divisible, 3)
        for second in combinations(nondivisible, 2)
    )
    require(len(subsets) == 12, "3+2 subset count changed")

    state_checks = 0
    for flags in product((0, 1), repeat=7):
        u = sum(flags[index] for index in divisible)
        v = sum(flags[index] for index in nondivisible)
        total = sum(
            sum(flags[i] * flags[j] for i, j in combinations(subset, 2))
            for subset in subsets
        )
        average = F(total, len(subsets))
        require(average == pair_average_polynomial(u, v), "five-subset P identity failed")
        state_checks += 1

    # Direct coefficient census: DD, DV, VV pairs occur in 6, 6, 4 subsets.
    type_counts: dict[str, set[int]] = {"DD": set(), "DV": set(), "VV": set()}
    for i, j in combinations(range(7), 2):
        occurrence = sum(i in subset and j in subset for subset in subsets)
        kind = "DD" if j < 4 else ("VV" if i >= 4 else "DV")
        type_counts[kind].add(occurrence)
    require(type_counts == {"DD": {6}, "DV": {6}, "VV": {4}}, "subset coefficients")
    return len(subsets), state_checks


def majorant_left(h_value: int, u: int, v: int) -> F:
    return (
        pair_average_polynomial(u, v)
        + F(65, 12) * h_value
        - F(65, 42) * u
        - F(13, 6) * v
        + ALPHA[u] * (F(h_value) - GUARD_MEASURE)
        + F(4, 3) * h_value * v
    )


def audit_pointwise_majorant() -> tuple[int, tuple[tuple[object, ...], ...], tuple[tuple[object, ...], ...]]:
    tables: list[tuple[tuple[object, ...], ...]] = []
    admissible_checks = 0
    for h_value in (0, 1):
        rows: list[tuple[object, ...]] = []
        for u in range(5):
            row: list[object] = []
            for v in range(4):
                slack = -420 * majorant_left(h_value, u, v)
                require(slack.denominator == 1, "scaled majorant slack is not integral")
                if h_value == 0 and u == 0 and v == 0:
                    require(slack == -910, "uncovered hostile control changed")
                    row.append(None)
                    continue
                require(slack >= 0, f"majorant failed at {(h_value, u, v)}")
                row.append(slack.numerator)
                admissible_checks += 1
            rows.append(tuple(row))
        tables.append(tuple(rows))

    require(tables[0] == EXPECTED_SLACK_H0, "H=0 slack table changed")
    require(tables[1] == EXPECTED_SLACK_H1, "H=1 slack table changed")
    require(admissible_checks == 39, "admissible state count changed")
    return admissible_checks, tables[0], tables[1]


def format_slack_row(row: tuple[object, ...]) -> str:
    return "[" + ",".join("X" if value is None else str(value) for value in row) + "]"


def main() -> None:
    orbit_checks, basis_checks, boundary_checks, period_checks = (
        audit_seven_shift_independence()
    )
    reduced_count, negative_count, charge_rows = audit_three_charge_floor()
    subset_count, subset_state_checks = audit_five_subset_average()
    majorant_checks, table_h0, table_h1 = audit_pointwise_majorant()

    c3 = sum((row[2] for row in charge_rows[:3]), F(0))
    require(c3 == F(713, 16170), "top-three charge sum changed")
    hv_floor = F(6, 49) - c3
    require(hv_floor == F(181, 2310), "HV floor changed")

    mean_h = GUARD_MEASURE
    mean_u = 4 * DANGER_MEASURE
    mean_v = 3 * DANGER_MEASURE
    p_upper = (
        -F(65, 12) * mean_h
        + F(65, 42) * mean_u
        + F(13, 6) * mean_v
        - F(4, 3) * hv_floor
    )
    require(p_upper == F(3901, 24255), "integrated P upper bound changed")
    contradiction_gap = PAIR_FLOOR - p_upper
    require(contradiction_gap == F(107, 315315), "contradiction gap changed")
    require(contradiction_gap > 0, "pair-floor contradiction is not strict")

    print("MOD-7 FOURFOLD MOMENT OBSTRUCTION EXACT REFEREE")
    print(f"seven-shift open-cell orbit checks: {orbit_checks}")
    print(f"arbitrary-phi basis coefficient checks: {basis_checks}")
    print(f"null-boundary controls: {boundary_checks}")
    print(f"7-divisible period increment checks: {period_checks}")
    print("excluded 7|h control: orbit guard count 7 (not 2)")
    print(f"allowed reduced pairs with ab<=20: {reduced_count}")
    print(f"negative head rows: {negative_count}")
    print("top three negative charges:")
    for a, b, charge in charge_rows[:3]:
        print(f"  (a,b)=({a},{b}) -> {charge}")
    print(f"top-three charge sum C3: {c3}")
    print(f"mixed first-moment floor E[H V]: {hv_floor}")
    print(f"3-divisible+2-nondivisible five-subsets: {subset_count}")
    print(f"Boolean subset-average identity checks: {subset_state_checks}")
    print(f"admissible pointwise-majorant states: {majorant_checks}")
    print("scaled slack table (-420 L), H=0:")
    for u, row in enumerate(table_h0):
        print(f"  u{u} {format_slack_row(row)}")
    print("scaled slack table (-420 L), H=1:")
    for u, row in enumerate(table_h1):
        print(f"  u{u} {format_slack_row(row)}")
    print("forbidden uncovered state (H,u,v)=(0,0,0): scaled slack -910")
    print(f"THM-1234 lower bound E[P]: {PAIR_FLOOR}")
    print(f"integrated-majorant upper bound E[P]: {p_upper}")
    print(f"strict contradiction gap: {contradiction_gap}")
    print("PASS")


if __name__ == "__main__":
    main()
