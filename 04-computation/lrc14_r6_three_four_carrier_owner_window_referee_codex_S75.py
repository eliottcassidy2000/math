#!/usr/bin/env python3
"""Exact referee for THM-1169's three/four-carrier owner windows.

Let ``P`` be a seven-subset of ``{1,...,12}``, let ``M=max(P)``, and let
``p=P[1]`` be its second-smallest speed.  If the normalized carriers are
``x<...``, the first carrier is safe throughout

    [1/(14x), 13/(14x)].

For four carriers, THM-1155's interval-density inequality prevents the three
later carrier combs from covering this whole interval.  For three carriers we
truncate at the exact owner threshold ``1/(14p)``.  The same inequality leaves
only six tiny first-carrier ranges; this script exhausts their exact two-comb
open covers.

Every decision uses ``Fraction`` or integer arithmetic.  Strict danger arcs
remain open: the cover test checks every rational wall *and* every open cell.
Two independent cover implementations (wall cells and merged open intervals)
must agree.

Tournament / alternate-vertex audit.  The six ``p``-strata are telemetry
vertices, oriented by exact residual candidate count and tied by increasing
``p``.  Reversing the tie switch flips one edge.  The predicate-preserving
object is instead the labelled rational wall-cell complex on the owner window,
with the two later comb obligations attached.  Runner order, carrier order,
and the six-stratum tournament all destroy endpoint openness and owner time.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import ceil, floor


H = F(1, 14)
P_VALUES = tuple(range(2, 8))
EXPECTED_CORE_HISTOGRAM = {2: 252, 3: 252, 4: 168, 5: 84, 6: 30, 7: 6}
EXPECTED_FIRST_GOOD = {2: 13, 3: 19, 4: 26, 5: 33, 6: 39, 7: 46}
EXPECTED_RESIDUAL_MAX = {2: 12, 3: 18, 4: 25, 5: 32, 6: 38, 7: 45}
EXPECTED_PAIR_CANDIDATES = {2: 0, 3: 0, 4: 3, 5: 34, 6: 115, 7: 321}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fractional_part(value: F) -> F:
    return value - (value.numerator // value.denominator)


def dangerous(speed: int, time: F) -> bool:
    phase = fractional_part(speed * time)
    return min(phase, 1 - phase) < H


def tooth_boundaries(speed: int, left: F, right: F) -> set[F]:
    """All strict-danger endpoints lying in the open interval (left,right)."""
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
    """Exact open cover test by endpoints plus one point in every live cell."""
    walls = {left, right}
    for speed in speeds:
        walls.update(tooth_boundaries(speed, left, right))
    ordered = sorted(walls)

    # Open danger arcs do not own their endpoints.  A cover must therefore
    # protect every wall by a different arc (or by a nonincident tooth).
    for wall in ordered:
        if not any(dangerous(speed, wall) for speed in speeds):
            return False
    for lower, upper in zip(ordered, ordered[1:]):
        midpoint = (lower + upper) / 2
        if not any(dangerous(speed, midpoint) for speed in speeds):
            return False
    return True


def danger_intervals(speed: int, left: F, right: F) -> list[tuple[F, F]]:
    """Unclipped open teeth which meet [left,right]."""
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
    """Merge only strict overlaps; touching open intervals still miss a point."""
    answer: list[list[F]] = []
    for lower, upper in sorted(intervals):
        if not answer or answer[-1][1] <= lower:
            answer.append([lower, upper])
        else:
            answer[-1][1] = max(answer[-1][1], upper)
    return [(lower, upper) for lower, upper in answer]


def interval_union_cover(speeds: tuple[int, ...], left: F, right: F) -> bool:
    """Independent exact cover test using connected components of an open union."""
    intervals: list[tuple[F, F]] = []
    for speed in speeds:
        intervals.extend(danger_intervals(speed, left, right))
    return any(
        lower < left and right < upper
        for lower, upper in merged_open_intervals(intervals)
    )


def safe_cells(speed: int, left: F, right: F) -> tuple[list[tuple[F, F]], int]:
    """Positive cells and safe wall count of [left,right] after one open comb."""
    walls = sorted({left, right} | tooth_boundaries(speed, left, right))
    cells = [
        (lower, upper)
        for lower, upper in zip(walls, walls[1:])
        if not dangerous(speed, (lower + upper) / 2)
    ]
    safe_wall_count = sum(not dangerous(speed, wall) for wall in walls)
    return cells, safe_wall_count


def owner_window(p: int, x: int) -> tuple[F, F]:
    left = F(1, 14 * x)
    right = min(F(13, 14 * x), F(1, 14 * p))
    require(left < right, "owner window is empty")
    return left, right


def density_tail_audit() -> None:
    # Four carriers: if the three later combs covered the full first-safe
    # window, THM-1155 would force 6/(7x) <= 3/(4(x+1)).  The exact difference
    # is (3x+24)/(28x(x+1)), hence positive for every legal x.
    for x in range(8, 10001):
        difference = F(6, 7 * x) - F(3, 4 * (x + 1))
        require(
            difference == F(3 * x + 24, 28 * x * (x + 1)) > 0,
            "four-carrier density identity failed",
        )
        length = F(6, 7 * x)
        require(F(4, 1) / (3 * length) == F(14 * x, 9), "rho5 ratio cone moved")
        require(F(5, 1) / (2 * length) == F(35 * x, 12), "rho6 ratio cone moved")

    # Three carriers in the truncated owner window.  In the non-full-window
    # case, density closes exactly when
    #   5 (x-p) (x+1) > 28 x p.
    for p in P_VALUES:
        legal = range(p + 6, 13 * p + 1)
        first_good = min(
            x for x in legal if 5 * (x - p) * (x + 1) > 28 * x * p
        )
        bad = [x for x in legal if 5 * (x - p) * (x + 1) <= 28 * x * p]
        require(first_good == EXPECTED_FIRST_GOOD[p], "first-good threshold moved")
        require(max(bad) == EXPECTED_RESIDUAL_MAX[p], "residual maximum moved")
        require(
            all(5 * (x - p) * (x + 1) > 28 * x * p for x in range(first_good, 13 * p + 1)),
            "density inequality is not monotone after its recorded threshold",
        )

        # Once x>=13p, the owner cutoff lies beyond the first-safe interval.
        # The full interval beats the two-comb bound by a positive identity.
        for x in range(13 * p, 13 * p + 50):
            difference = F(6, 7 * x) - F(2, 5 * (x + 1))
            require(
                difference == F(16 * x + 30, 35 * x * (x + 1)) > 0,
                "three-carrier full-window identity failed",
            )


def finite_two_comb_audit() -> tuple[dict[int, int], dict[int, int], int, int]:
    candidate_counts: dict[int, int] = {}
    y_row_counts: dict[int, int] = {}
    point_only_rows = 0
    empty_after_y_rows = 0

    for p in P_VALUES:
        pair_candidates = 0
        y_rows = 0
        residual_max = EXPECTED_RESIDUAL_MAX[p]

        # p is the second member of a seven-subset, so its maximum is at least
        # p+5 and the least normalized carrier x is at least p+6.
        for x in range(p + 6, residual_max + 1):
            left, right = owner_window(p, x)
            length = right - left

            # A two-comb cover gives 5L <= 1/y+1/z < 2/y, hence
            # y < 2/(5L).  ceil(bound)-1 is the exact integer strict cap.
            y_cap = ceil(F(2, 5) / length) - 1
            for y in range(x + 1, y_cap + 1):
                y_rows += 1
                cells, safe_wall_count = safe_cells(y, left, right)
                if not cells:
                    # A point-only residual would require a separate congruence
                    # audit for unbounded z.  An empty residual would mean that
                    # D_y alone already covers the owner window.  Neither case
                    # occurs in the complete box, and both are counted rather
                    # than silently discarded.
                    if safe_wall_count:
                        point_only_rows += 1
                    else:
                        empty_after_y_rows += 1
                    continue

                longest = max(upper - lower for lower, upper in cells)
                # A closed interval of length ell cannot lie in one open z-tooth
                # of length 1/(7z) unless z < 1/(7ell).
                z_cap = ceil(F(1, 7) / longest) - 1
                for z in range(y + 1, z_cap + 1):
                    pair_candidates += 1
                    first = wall_cell_cover((y, z), left, right)
                    second = interval_union_cover((y, z), left, right)
                    require(first == second, "independent open-cover tests disagree")
                    require(not first, f"two later carriers cover the owner window: {(p,x,y,z)}")

        candidate_counts[p] = pair_candidates
        y_row_counts[p] = y_rows
        require(pair_candidates == EXPECTED_PAIR_CANDIDATES[p], "candidate count moved")

    require(point_only_rows == 0, "point-only residual needs an unbounded congruence audit")
    require(empty_after_y_rows == 0, "one later carrier covers the owner window")
    return candidate_counts, y_row_counts, point_only_rows, empty_after_y_rows


def tournament_audit(counts: dict[int, int]) -> tuple[tuple[int, ...], int, int]:
    vertices = list(P_VALUES)

    def orientation(tie_reverse: bool) -> set[tuple[int, int]]:
        edges = set()
        for left, right in combinations(vertices, 2):
            key_left = (counts[left], -left if tie_reverse else left)
            key_right = (counts[right], -right if tie_reverse else right)
            winner, loser = (left, right) if key_left > key_right else (right, left)
            edges.add((winner, loser))
        return edges

    up = orientation(False)
    down = orientation(True)
    scores = tuple(sorted(sum(source == v for source, _ in up) for v in vertices))
    flips = sum(((a, b) in up) != ((a, b) in down) for a, b in combinations(vertices, 2))
    cycles = sum(
        ((a, b) in up and (b, c) in up and (c, a) in up)
        or ((b, a) in up and (c, b) in up and (a, c) in up)
        for a, b, c in combinations(vertices, 3)
    )
    require(scores == tuple(range(6)), "residual-count tournament is not transitive")
    require(flips == 1 and cycles == 0, "tournament fingerprint changed")
    return scores, flips, cycles


def main() -> None:
    core_histogram = Counter(core[1] for core in combinations(range(1, 13), 7))
    require(len(list(combinations(range(1, 13), 7))) == 792, "core count changed")
    require(dict(core_histogram) == EXPECTED_CORE_HISTOGRAM, "p2 histogram changed")

    density_tail_audit()
    counts, y_rows, point_rows, empty_rows = finite_two_comb_audit()
    scores, flips, cycles = tournament_audit(counts)

    print("THM-1169 three/four-carrier owner-window closure")
    print("arithmetic=fractions.Fraction and integers only")
    print(f"seven_core_count=792")
    print(f"second_speed_histogram={dict(core_histogram)}")
    print("rho4=uniform: 6/(7x) > 3/(4(x+1)); three later combs cannot cover")
    print("rho3_owner_window=[1/(14x), min(13/(14x),1/(14p2))]")
    print(f"rho3_first_density_good={EXPECTED_FIRST_GOOD}")
    print(f"rho3_residual_x_max={EXPECTED_RESIDUAL_MAX}")
    print(f"rho3_finite_y_rows={y_rows}")
    print(f"rho3_exact_pair_candidates={counts}; total={sum(counts.values())}")
    print(f"rho3_point_only_rows={point_rows}")
    print(f"rho3_empty_after_first_comb_rows={empty_rows}")
    print("rho3_exact_two_comb_covers=0")
    print("boundary_readout=u=c/(14s), cM<13s; u<=1/(14p2) gives d<=5")
    print("rho5_density_residual=second carrier y < 14x/9")
    print("rho6_density_residual=second carrier y < 35x/12")
    print("Tournament Analysis:")
    print(f"  vertices=p2 residual strata 2..7; scores={scores}")
    print(f"  directed_cycles={cycles}; SCCs=(1,1,1,1,1,1); HP=1; tie_gauge_flips={flips}")
    print("  pairwise_observable=larger exact finite residual candidate count")
    print("  switch/gauge=reverse the p2 tie order; directed HP=7>6>5>4>3>2")
    print("  faithful_object=labelled rational wall cells plus two comb obligations")
    print("  destroyed_by_tournament=endpoint openness|owner time|comb incidence|metric length")
    print("  challenged_vertices=runners|carriers|ratios|walls|cells|owner obligations|cover obligations")
    print("VERDICT: rho=3,4 close uniformly; rho=5,6 remain only in the recorded close-pair cones")


if __name__ == "__main__":
    main()
