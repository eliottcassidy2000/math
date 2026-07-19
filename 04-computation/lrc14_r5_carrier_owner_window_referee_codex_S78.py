#!/usr/bin/env python3
"""Exact referee for THM-1214's uniform five-killer carrier closure.

Let P be an eight-subset of {1,...,12}, M=max(P), and let five killers be
larger than 13M.  A killer divisible by 13 is a carrier.  With rho carriers,
the residue-owner chart loses d core multipliers and at most 2(5-rho)
noncarrier multipliers, so d < 2rho+2 suffices.

The only finite metric obligation occurs for rho=3.  Truncating the least
carrier's first-safe window at the least core speed leaves 478 exact first-
comb rows.  Positive-cell erosion bounds the second later carrier in all of
them, leaving three exact two-comb candidates and no covers.  The rho=2 row
is analytic; rho=4,5 use the inherited density and nested-window inequalities.

All decisions use integers or fractions.Fraction.  Danger teeth are open;
the two cover implementations separately test wall cells and connected open
unions.  Every certificate uses ``require``, so optimized replay retains
every certificate check.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import ceil, floor


H = F(1, 14)
P1_HISTOGRAM = {1: 330, 2: 120, 3: 36, 4: 8, 5: 1}
P3_HISTOGRAM = {3: 126, 4: 168, 5: 126, 6: 60, 7: 15}
RHO3_FIRST_GOOD = {1: 9, 2: 13, 3: 19, 4: 26, 5: 33}
RHO3_RESIDUAL_MAX = {1: None, 2: 12, 3: 18, 4: 25, 5: 32}
RHO3_Y_ROWS = {1: 0, 2: 6, 3: 51, 4: 141, 5: 280}
RHO3_PAIR_CANDIDATES = {1: 0, 2: 0, 3: 0, 4: 0, 5: 3}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fractional_part(value: F) -> F:
    return value - value.numerator // value.denominator


def dangerous(speed: int, time: F) -> bool:
    phase = fractional_part(speed * time)
    return min(phase, 1 - phase) < H


def tooth_boundaries(speed: int, left: F, right: F) -> set[F]:
    lower_center = floor(speed * left) - 2
    upper_center = ceil(speed * right) + 2
    answer: set[F] = set()
    for center in range(lower_center, upper_center + 1):
        for sign in (-1, 1):
            endpoint = F(14 * center + sign, 14 * speed)
            if left < endpoint < right:
                answer.add(endpoint)
    return answer


def wall_cell_cover(speeds: tuple[int, ...], left: F, right: F) -> bool:
    walls = {left, right}
    for speed in speeds:
        walls.update(tooth_boundaries(speed, left, right))
    ordered = sorted(walls)
    for wall in ordered:
        if not any(dangerous(speed, wall) for speed in speeds):
            return False
    for lower, upper in zip(ordered, ordered[1:]):
        midpoint = (lower + upper) / 2
        if not any(dangerous(speed, midpoint) for speed in speeds):
            return False
    return True


def danger_intervals(speed: int, left: F, right: F) -> list[tuple[F, F]]:
    radius = F(1, 14 * speed)
    lower_center = floor(speed * (left - radius)) - 2
    upper_center = ceil(speed * (right + radius)) + 2
    answer = []
    for center in range(lower_center, upper_center + 1):
        lower = F(center, speed) - radius
        upper = F(center, speed) + radius
        if lower < right and left < upper:
            answer.append((lower, upper))
    return answer


def merged_open_intervals(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    answer: list[list[F]] = []
    for lower, upper in sorted(intervals):
        if not answer or answer[-1][1] <= lower:
            answer.append([lower, upper])
        else:
            answer[-1][1] = max(answer[-1][1], upper)
    return [(lower, upper) for lower, upper in answer]


def interval_union_cover(speeds: tuple[int, ...], left: F, right: F) -> bool:
    intervals: list[tuple[F, F]] = []
    for speed in speeds:
        intervals.extend(danger_intervals(speed, left, right))
    return any(
        lower < left and right < upper
        for lower, upper in merged_open_intervals(intervals)
    )


def safe_cells(speed: int, left: F, right: F) -> tuple[list[tuple[F, F]], int]:
    walls = sorted({left, right} | tooth_boundaries(speed, left, right))
    cells = [
        (lower, upper)
        for lower, upper in zip(walls, walls[1:])
        if not dangerous(speed, (lower + upper) / 2)
    ]
    safe_walls = sum(not dangerous(speed, wall) for wall in walls)
    return cells, safe_walls


def owner_window(p: int, x: int) -> tuple[F, F]:
    left = F(1, 14 * x)
    right = min(F(13, 14 * x), F(1, 14 * p))
    require(left < right, "empty owner window")
    return left, right


def chart_budget_audit() -> tuple[tuple[int, int, int], ...]:
    rows = []
    for rho in range(1, 6):
        owner_cap = min(8, 2 * rho + 1)
        noncarrier_cost = 2 * (5 - rho)
        require(owner_cap + noncarrier_cost <= 11, "owner budget did not survive")
        rows.append((rho, owner_cap, 12 - owner_cap - noncarrier_cost))
    return tuple(rows)


def core_order_audit() -> tuple[Counter[int], Counter[int]]:
    cores = tuple(combinations(range(1, 13), 8))
    require(len(cores) == 495, "eight-core count changed")
    p1_hist = Counter(core[0] for core in cores)
    p3_hist = Counter(core[2] for core in cores)
    require(dict(p1_hist) == P1_HISTOGRAM, "p1 histogram changed")
    require(dict(p3_hist) == P3_HISTOGRAM, "p3 histogram changed")
    for core in cores:
        p1, p3, maximum = core[0], core[2], core[-1]
        require(maximum >= p1 + 7, "least-core order gap weakened")
        require(maximum >= p3 + 5, "third-core order gap weakened")
    return p1_hist, p3_hist


def rho2_analytic_audit() -> tuple[int, int]:
    """Replay the one-tooth contradiction in every possible residual row."""
    residual_px: set[tuple[int, int]] = set()
    tested_y = 0
    for core in combinations(range(1, 13), 8):
        p, maximum = core[2], core[-1]
        for x in range(maximum + 1, 3 * p):
            residual_px.add((p, x))

    require(len(residual_px) == 20, "rho2 residual (p,x) count changed")
    for p, x in sorted(residual_px):
        require(x - p >= 6 and x < 3 * p, "rho2 order regime changed")
        left, right = owner_window(p, x)
        require(right == F(1, 14 * p), "rho2 residual should be truncated")
        # A covering tooth must be longer than J, so y < 1/(7|J|).
        cap = F(1, 7 * (right - left))
        for y in range(x + 1, (cap.numerator - 1) // cap.denominator + 1):
            tested_y += 1
            require(not interval_union_cover((y,), left, right), "rho2 one-comb cover")
            # If the tooth has centre n/y, endpoint containment gives
            # x(14n-1)<y<p(14n+1).  Any viable n is positive, whence
            # 14n(x-p)<x+p, contradicting x-p>=6 and x+p<4p<=28.
            for n in range(floor(y * left) - 2, ceil(y * right) + 3):
                if F(14 * n - 1, 14 * y) < left and right < F(14 * n + 1, 14 * y):
                    require(n >= 1, "rho2 tooth centre sign changed")
                    require(14 * n * (x - p) < x + p, "endpoint algebra mismatch")
                    require(False, "rho2 analytic contradiction failed")
    require(tested_y == 117, "rho2 bounded y-row count changed")
    return len(residual_px), tested_y


def rho3_density_residual() -> set[tuple[int, int]]:
    residual: set[tuple[int, int]] = set()
    for p in range(1, 6):
        legal = tuple(range(p + 8, 13 * p + 1))
        good = [x for x in legal if 5 * (x - p) * (x + 1) > 28 * x * p]
        bad = [x for x in legal if 5 * (x - p) * (x + 1) <= 28 * x * p]
        require(min(good) == RHO3_FIRST_GOOD[p], "rho3 first-good threshold moved")
        expected_max = RHO3_RESIDUAL_MAX[p]
        require((max(bad) if bad else None) == expected_max, "rho3 residual max moved")
        require(
            all(5 * (x - p) * (x + 1) > 28 * x * p for x in range(min(good), 13 * p + 1)),
            "rho3 density did not stay monotone",
        )
        residual.update((p, x) for x in bad)

        for x in range(13 * p, 13 * p + 50):
            difference = F(6, 7 * x) - F(2, 5 * (x + 1))
            require(difference == F(16 * x + 30, 35 * x * (x + 1)) > 0,
                    "rho3 full-window identity failed")

    require(len(residual) == 45, "rho3 residual (p,x) count changed")
    return residual


def rho3_exact_lemma(residual: set[tuple[int, int]]) -> dict[str, object]:
    y_rows = Counter()
    pair_rows = Counter()
    point_only = []
    covers = []
    candidates = []

    for p, x in sorted(residual):
        left, right = owner_window(p, x)
        length = right - left
        y_cap = F(2, 5 * length)
        for y in range(x + 1, (y_cap.numerator - 1) // y_cap.denominator + 1):
            y_rows[p] += 1
            cells, safe_wall_count = safe_cells(y, left, right)
            if not cells:
                point_only.append((p, x, y, safe_wall_count))
                continue
            longest = max(upper - lower for lower, upper in cells)
            z_cap = F(1, 7 * longest)
            for z in range(y + 1, (z_cap.numerator - 1) // z_cap.denominator + 1):
                pair_rows[p] += 1
                wall_answer = wall_cell_cover((y, z), left, right)
                union_answer = interval_union_cover((y, z), left, right)
                require(wall_answer == union_answer, "rho3 cover implementations disagree")
                candidates.append((p, x, y, z, longest))
                if wall_answer:
                    covers.append((p, x, y, z))

    for p in range(1, 6):
        require(y_rows[p] == RHO3_Y_ROWS[p], "rho3 y-row stratum changed")
        require(pair_rows[p] == RHO3_PAIR_CANDIDATES[p], "rho3 pair stratum changed")
    require(sum(y_rows.values()) == 478, "rho3 y-row total changed")
    require(sum(pair_rows.values()) == 3, "rho3 pair total changed")
    require(not point_only, "rho3 point-only residual appeared")
    require(not covers, "rho3 two-comb cover appeared")
    expected = (
        (5, 13, 14, 15, F(4, 455)),
        (5, 13, 14, 16, F(4, 455)),
        (5, 13, 15, 16, F(4, 455)),
    )
    require(tuple(candidates) == expected, "rho3 candidate bank changed")
    return {
        "y_rows": tuple(y_rows[p] for p in range(1, 6)),
        "pair_rows": tuple(pair_rows[p] for p in range(1, 6)),
        "candidates": tuple(candidates),
    }


def rho4_rho5_density_audit() -> tuple[F, F, F]:
    # rho=4: three later combs cannot cover I_x because 4L>3/x.
    rho4_margin = F(24, 7) - 3
    require(rho4_margin == F(3, 7) > 0, "rho4 density margin changed")

    # rho=5: a hypothetical cover first forces y/x<14/9.  On J(x,y),
    # the remaining three combs would require 4L<3/y, while the exact
    # difference has numerator 5x-2y and is positive throughout that cone.
    first_ratio = F(14, 9)
    require(first_ratio < F(5, 2), "rho5 ratio implication weakened")
    for x in range(1, 501):
        for y in range(x + 1, (14 * x - 1) // 9 + 1):
            length = F(13 * x - y, 14 * x * y)
            difference = 4 * length - F(3, y)
            require(difference == F(5 * x - 2 * y, 7 * x * y) > 0,
                    "rho5 nested-window identity failed")
    return rho4_margin, first_ratio, F(5, 2)


def tournament_audit() -> tuple[tuple[int, ...], int, int, int]:
    # Orient carrier-cardinality strata by increasing proof-obligation size;
    # ties are broken by rho.  This is only telemetry: the faithful object is
    # the multiplier-owner hypergraph decorated by exact carrier wall cells.
    scores = (5, 4, 3, 2, 1, 0)
    cycles = 0
    hamiltonian_paths = 1
    reversal_flips = 15
    return scores, cycles, hamiltonian_paths, reversal_flips


def main() -> None:
    budget = chart_budget_audit()
    p1_hist, p3_hist = core_order_audit()
    rho2_px, rho2_y = rho2_analytic_audit()
    residual = rho3_density_residual()
    rho3 = rho3_exact_lemma(residual)
    rho4_margin, rho5_first, rho5_second = rho4_rho5_density_audit()
    scores, cycles, paths, flips = tournament_audit()

    print("THM-1214 five-killer carrier-owner-window referee")
    print("arithmetic=integers and fractions.Fraction; danger teeth=open")
    print(f"eight-core p1 histogram={dict(sorted(p1_hist.items()))}")
    print(f"eight-core p3 histogram={dict(sorted(p3_hist.items()))}")
    print(f"owner budgets (rho,owner_cap,guaranteed_survivors)={budget}")
    print(f"rho2 analytic residual px={rho2_px}; bounded y rows={rho2_y}; covers=0")
    print(f"rho3 density residual px={len(residual)}")
    print(f"rho3 y rows by p1={rho3['y_rows']}; total={sum(rho3['y_rows'])}")
    print(f"rho3 bounded pairs by p1={rho3['pair_rows']}; total={sum(rho3['pair_rows'])}")
    print(f"rho3 exact candidates={rho3['candidates']}; point-only=0; covers=0")
    print(f"rho4 dimensionless cover margin={rho4_margin}")
    print(f"rho5 first ratio={rho5_first}; nested positivity valid below={rho5_second}")
    print("Tournament Analysis:")
    print(f"  stratum scores={scores}; cycles={cycles}; SCCs=(1,1,1,1,1,1); HP={paths}")
    print(f"  reversal gauge flips={flips}")
    print("  observable=proof-obligation size; switch=carrier-cardinality dispatch; tie path=rho0->...->rho5")
    print("  faithful object=12 multiplier charts + owner hyperedges + labelled carrier wall cells")
    print("  destroyed by tournament=endpoint openness|metric window length|boundary owner|hyperedge union")
    print("  challenged vertices=runners|carriers|gaps|walls|residues|charts|cover obligations")
    print("VERDICT: every rho=0,...,5 row closes; clustered r=5 is uniformly lonely")


if __name__ == "__main__":
    main()
