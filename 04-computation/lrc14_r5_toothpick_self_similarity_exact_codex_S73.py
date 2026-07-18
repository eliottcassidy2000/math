#!/usr/bin/env python3
"""Exact referee for THM-1127's fixed four-comb toothpick ray.

The theorem is deliberately local to

    P = (1,2,4,5,7,9,11,12),
    (k1,k2,k3,k4) = (K,K+4,K+5,K+9).

It does not close the uniform r=5 clustered stratum.  The script uses exact
``fractions.Fraction`` endpoint subtraction in the original time coordinate
and, independently, the fixed torus polygon

    H_t = {x in R/Z : ||x+a*t|| >= 1/14 for a in {0,4,5,9}}.

The identity ``x = K*t (mod 1)`` makes the remainder the slope-K section of
this one polygon.  This is the toothpick self-similarity carrier.

Tournament Analysis audit.  Candidate vertices include runners, combs, core
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, and proof obligations.  The useful vertices here are labelled
rational polygon walls.  Exact coordinate order orients their pair relation;
coincident walls are the ties and the coalesced sorted wall word is the tie
Hamiltonian path.  The resulting tournament is transitive.  Its naked order
forgets wall slopes, affine owners, strip adjacency, lengths, and the core
mask, so the faithful carrier is the labelled wall arrangement plus the
slope-K line, not a runner tournament.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import lcm


CORE = (1, 2, 4, 5, 7, 9, 11, 12)
OFFSETS = (0, 4, 5, 9)
H = F(1, 14)
SCAN_LO = 150
SCAN_HI = 1500
FORMULA_START = 771
TORUS_PERIOD = 194_040


@dataclass(frozen=True)
class Boundary:
    value: F
    owner: str


Interval = tuple[Boundary, Boundary]
Region = list[Interval]
PlainInterval = tuple[F, F]


def interval_length(interval: Interval) -> F:
    return interval[1].value - interval[0].value


def remove_danger(region: Region, speed: int) -> Region:
    """Subtract the closed 1/14-danger comb, retaining exact endpoint owners."""

    result: Region = []
    denominator = 14 * speed
    for left, right in region:
        cursor = left
        j_lo = int(left.value * speed) - 1
        j_hi = int(right.value * speed) + 1
        for tooth in range(j_lo, j_hi + 1):
            bad_left = F(14 * tooth - 1, denominator)
            bad_right = F(14 * tooth + 1, denominator)
            if bad_right <= left.value or right.value <= bad_left:
                continue
            clipped_left = max(left.value, bad_left)
            clipped_right = min(right.value, bad_right)
            if cursor.value < clipped_left:
                result.append(
                    (
                        cursor,
                        Boundary(clipped_left, f"D_{speed}:tooth_{tooth}:left"),
                    )
                )
            if cursor.value < clipped_right:
                cursor = Boundary(
                    clipped_right, f"D_{speed}:tooth_{tooth}:right"
                )
        if cursor.value < right.value:
            result.append((cursor, right))
    return result


def safe_region(speeds: tuple[int, ...]) -> Region:
    region = [(Boundary(F(0), "unit:0"), Boundary(F(1), "unit:1"))]
    for speed in speeds:
        region = remove_danger(region, speed)
    return region


def four_comb_remainder(core_safe: Region, k: int) -> Region:
    region = core_safe
    for offset in OFFSETS:
        region = remove_danger(region, k + offset)
    return region


def merge_intervals(intervals: list[PlainInterval]) -> list[PlainInterval]:
    merged: list[PlainInterval] = []
    for left, right in sorted(intervals):
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return merged


def torus_safe_section(t: F) -> list[PlainInterval]:
    """Return H_t exactly as intervals in the x-circle cut at x=0."""

    bad: list[PlainInterval] = []
    for offset in OFFSETS:
        centre = (-offset * t) % 1
        left = centre - H
        right = centre + H
        if left < 0:
            bad.extend(((F(0), right), (left + 1, F(1))))
        elif right > 1:
            bad.extend(((F(0), right - 1), (left, F(1))))
        else:
            bad.append((left, right))
    merged = merge_intervals(bad)
    safe: list[PlainInterval] = []
    cursor = F(0)
    for left, right in merged:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        safe.append((cursor, F(1)))
    # The offset-0 danger interval always covers the circle cut, so no two
    # surviving pieces at 0 and 1 need to be rejoined.
    return safe


def plain_measure(intervals: list[PlainInterval]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def core_plain(core_safe: Region) -> list[PlainInterval]:
    return [(left.value, right.value) for left, right in core_safe]


def inside_core(t: F, core: list[PlainInterval]) -> bool:
    return any(left < t < right for left, right in core)


def torus_breakpoints() -> set[F]:
    """Walls where two endpoints -a*t +/- 1/14 coincide modulo one."""

    points = {F(0), F(1)}
    for a in OFFSETS:
        for b in OFFSETS:
            if a >= b:
                continue
            difference = b - a
            for shift in (-F(1, 7), F(0), F(1, 7)):
                for integer in range(-1, difference + 2):
                    t = (F(integer) + shift) / difference
                    if 0 <= t <= 1:
                        points.add(t)
    return points


def longest_formula(k: int) -> F:
    m = (97 * k - 107) // 126
    return F(14 * m - 2 * k + 47, 14 * (k + 4) * (k + 5))


def containment_margin(k: int) -> int:
    """Nonnegative iff the rightmost full dominant strip lies in its core cell."""

    m = (97 * k - 107) // 126
    return 56 * m - 43 * k + 13


def main() -> None:
    core_safe = safe_region(CORE)
    core = core_plain(core_safe)
    core_measure = sum((interval_length(interval) for interval in core_safe), F(0))
    core_endpoint_period = 1
    for left, right in core:
        core_endpoint_period = lcm(
            core_endpoint_period, left.denominator, right.denominator
        )

    assert len(core_safe) == 20
    assert core_measure == F(549, 2695)
    assert core_endpoint_period == TORUS_PERIOD

    raw_walls = torus_breakpoints()
    raw_wall_period = 1
    for point in raw_walls:
        raw_wall_period = lcm(raw_wall_period, point.denominator)
    all_walls = set(raw_walls)
    for left, right in core:
        all_walls.update((left, right))
    wall_period = 1
    for point in all_walls:
        wall_period = lcm(wall_period, point.denominator)
    ordered_walls = sorted(all_walls)

    limit_measure = F(0)
    count_slope = F(0)
    core_cells: list[tuple[F, F, int, F]] = []
    maximum_gap = F(0)
    maximum_points: set[F] = set()
    dominant_cells = {
        (F(29, 126), F(13, 56)),
        (F(43, 56), F(97, 126)),
    }
    outside_gap = F(0)

    for left, right in zip(ordered_walls, ordered_walls[1:]):
        midpoint = (left + right) / 2
        if not inside_core(midpoint, core):
            continue
        middle_section = torus_safe_section(midpoint)
        component_count = len(middle_section)
        middle_measure = plain_measure(middle_section)
        left_measure = plain_measure(torus_safe_section(left))
        right_measure = plain_measure(torus_safe_section(right))
        assert left_measure + right_measure == 2 * middle_measure
        limit_measure += (right - left) * (left_measure + right_measure) / 2
        count_slope += (right - left) * component_count
        core_cells.append((left, right, component_count, middle_measure))

        for point in (left, right):
            gap = max(
                (b - a for a, b in torus_safe_section(point)), default=F(0)
            )
            if gap > maximum_gap:
                maximum_gap = gap
                maximum_points = {point}
            elif gap == maximum_gap:
                maximum_points.add(point)
            if (left, right) not in dominant_cells:
                outside_gap = max(outside_gap, gap)

    assert len(raw_walls) == 53
    assert raw_wall_period == 1260
    assert len(ordered_walls) == 91
    assert wall_period == TORUS_PERIOD
    assert len(core_cells) == 28
    assert limit_measure == F(59_151_097, 627_525_360)
    assert count_slope == F(15_163, 24_255)
    assert maximum_gap == F(79, 126)
    assert maximum_points == {F(29, 126), F(97, 126)}
    assert outside_gap == F(67, 154)

    mass_limit = count_slope / (6 * limit_measure)
    component_limit = F(1, 1) / (3 * maximum_gap)
    count_increment = TORUS_PERIOD * count_slope
    assert mass_limit == F(65_382_856, 59_151_097)
    assert component_limit == F(42, 79)
    assert count_increment == 121_304
    assert component_limit < mass_limit

    # The exact strip formula.  The margin has the useful residue recurrence
    # B(K+13) = B(K)+1, except at a remainder wrap where it jumps by 57.
    base_margins = tuple(containment_margin(k) for k in range(771, 784))
    assert base_margins == (12, 25, 38, 51, 8, 21, 34, 47, 4, 17, 30, 43, 0)
    assert containment_margin(770) == -1
    assert min(base_margins) == 0
    for k in range(771, 771 + 13 * 126):
        jump = containment_margin(k + 13) - containment_margin(k)
        assert jump in (1, 57)

    # Outside the two dominant cells every vertical gap is at most 67/154.
    # A labelled wall has slope of absolute value at most 9, hence a section
    # component there is at most (67/154)/(K-9).  The following increasing
    # quadratic proves that the dominant formula beats that bound from 771 on.
    def dominance_polynomial(k: int) -> int:
        return (
            79 * 154 * k * (k - 9)
            - 67 * 126 * (k + 5) * (k + 5)
        )

    assert dominance_polynomial(FORMULA_START) > 0
    assert (
        dominance_polynomial(FORMULA_START + 1)
        - dominance_polynomial(FORMULA_START)
        > 0
    )

    rows: list[tuple[int, int, F, F, F, F, F, F]] = []
    formula_mismatches: list[int] = []
    r_at_least_one: list[int] = []
    mass_branch: list[int] = []
    for k in range(SCAN_LO, SCAN_HI + 1):
        region = four_comb_remainder(core_safe, k)
        count = len(region)
        measure = sum((interval_length(interval) for interval in region), F(0))
        longest = max(interval_length(interval) for interval in region)
        threshold_mass = F(count, 1) / (6 * measure)
        threshold_component = F(1, 1) / (3 * longest)
        threshold = min(threshold_mass, threshold_component)
        ratio = threshold / (k + 9)
        rows.append(
            (
                k,
                count,
                measure,
                longest,
                threshold_mass,
                threshold_component,
                threshold,
                ratio,
            )
        )
        if longest != longest_formula(k):
            formula_mismatches.append(k)
        if ratio >= 1:
            r_at_least_one.append(k)
        if threshold_mass <= threshold_component:
            mass_branch.append(k)

    assert formula_mismatches
    assert max(formula_mismatches) == 770
    assert all(
        row[3] == longest_formula(row[0])
        for row in rows
        if row[0] >= FORMULA_START
    )
    assert r_at_least_one == [
        151,
        153,
        155,
        156,
        157,
        158,
        160,
        164,
        168,
        169,
        170,
        172,
        173,
        181,
        186,
        190,
        195,
        203,
        207,
        220,
        225,
        238,
        242,
        255,
        259,
        272,
        294,
        307,
        311,
        324,
        363,
    ]
    clustered_r_at_least_one = [k for k in r_at_least_one if k >= 157]
    assert len(clustered_r_at_least_one) == 27
    assert clustered_r_at_least_one[-1] == 363
    assert mass_branch == [207]

    maximum_ratio_row = max(rows, key=lambda row: row[7])
    clustered_maximum_ratio_row = max(
        (row for row in rows if row[0] >= 157), key=lambda row: row[7]
    )
    sharp_target_minimum_row = min(
        (row for row in rows if 157 <= row[0] <= 770),
        key=lambda row: 7 * (row[0] + 9) * row[3],
    )
    assert maximum_ratio_row[0] == 155
    assert maximum_ratio_row[7] == F(7, 6)
    assert clustered_maximum_ratio_row[0] == 173
    assert clustered_maximum_ratio_row[7] == F(2422, 2085)
    assert sharp_target_minimum_row[0] == 173
    assert (
        7 * (sharp_target_minimum_row[0] + 9) * sharp_target_minimum_row[3]
        == F(695, 346)
    )
    row_363 = next(row for row in rows if row[0] == 363)
    row_364 = next(row for row in rows if row[0] == 364)
    row_770 = next(row for row in rows if row[0] == 770)
    row_771 = next(row for row in rows if row[0] == 771)
    assert row_363[7] == F(242, 237)
    assert row_364[7] == F(633_696, 759_055)
    assert row_770[3] == F(17, 21_672)
    assert longest_formula(770) == F(6781, 8_397_900)
    assert row_771[3] == longest_formula(771)

    # The exact formula itself is stronger than the old 1/(3*k4) target.
    for k in range(FORMULA_START, FORMULA_START + 126):
        assert 3 * (k + 9) * longest_formula(k) > 1
        assert 7 * (k + 9) * longest_formula(k) > 1

    # High-slope spot checks of the exact affine count recurrence.  The proof
    # is the rational-cell translation in THM-1127; these checks exercise the
    # original endpoint carrier without constructing a 194040-row phase table.
    recurrence_checks: list[tuple[int, int, int]] = []
    for k in (150, 200, 771, 1500):
        low = len(four_comb_remainder(core_safe, k))
        high = len(four_comb_remainder(core_safe, k + TORUS_PERIOD))
        assert high - low == count_increment
        recurrence_checks.append((k, low, high))

    print("THM-1127 fixed-ray toothpick self-similarity: exact Fraction referee")
    print("scope=fixed core P=[1,2,4,5,7,9,11,12] and offsets [0,4,5,9] only")
    print("uniform_r5_claim=False")
    print(f"core_components={len(core_safe)}")
    print(f"core_measure={core_measure}")
    print(f"torus_raw_walls={len(raw_walls)}")
    print(f"torus_core_cells={len(core_cells)}")
    print(f"torus_period_Q={TORUS_PERIOD}")
    print(f"limit_mu={limit_measure}")
    print(f"limit_N_over_K={count_slope}")
    print(f"exact_N_increment_over_Q={count_increment}")
    print(f"limit_KL={maximum_gap}")
    print("limit_KL_maximizers=t=29/126,97/126")
    print(f"next_vertical_gap_outside_dominant_cells={outside_gap}")
    print(f"limit_Tmass_over_K={mass_limit}")
    print(f"limit_Tcomp_over_K={component_limit}")
    print(f"limit_R={component_limit}")
    print(
        "longest_formula_for_all_K_ge_771="
        "m=floor((97K-107)/126); "
        "L=(14m-2K+47)/(14(K+4)(K+5))"
    )
    print(f"containment_margin_K770={containment_margin(770)}")
    print(f"containment_margin_base_K771_to_783={base_margins}")
    print("containment_recurrence=B(K+13)-B(K) in {1,57}")
    print(f"scan_range={SCAN_LO}..{SCAN_HI}")
    print(f"scan_rows={len(rows)}")
    print("clustered_range=K>=157 because 13*max(P)=156")
    print("exact_finite_bank=K=157..770")
    print("proved_formula_tail=all K>=771")
    print(
        "sharp_target_min_157_to_770="
        "7*(K+9)*L=695/346 at K=173"
    )
    print("sharp_target_all_legal_K=7*(K+9)*L(K)>1 for every K>=157")
    print(f"scan_formula_last_mismatch={max(formula_mismatches)}")
    print(f"scan_R_ge_1_count={len(r_at_least_one)}")
    print(f"scan_clustered_R_ge_1_count={len(clustered_r_at_least_one)}")
    print(f"scan_R_ge_1_last_K={r_at_least_one[-1]}")
    print(f"scan_R_max=7/6 at K={maximum_ratio_row[0]}")
    print(
        "scan_clustered_K_ge_157_R_max="
        f"{clustered_maximum_ratio_row[7]} at K={clustered_maximum_ratio_row[0]}"
    )
    print(f"scan_mass_branch_K={mass_branch}")
    print(f"K363_R={row_363[7]}")
    print(f"K364_R={row_364[7]}")
    print("combined_exact_conclusion=R(K)<1 for every K>=364 on this fixed ray")
    print(f"K770_actual_L={row_770[3]}")
    print(f"K770_formula_candidate={longest_formula(770)}")
    print(f"K771_L={row_771[3]}")
    for k, low, high in recurrence_checks:
        print(
            f"N_recurrence_spot_K={k}:N(K)={low}:"
            f"N(K+Q)={high}:difference={high-low}"
        )
    print("functional_form:")
    print("  N(K+194040)=N(K)+121304 (affine quasipolynomial)")
    print(
        "  mu(K)=period-194040 quasirational; on each phase it is a rational "
        "sum of affine-wall intersection terms"
    )
    print("  L(K)=period-126 rational formula above for K>=771")
    print("  Tmass and R are quasirational; R tends to 42/79<1")
    print("tournament_analysis:")
    print(
        "  alternate_vertices=runners|combs|core_sections|section_boundaries|"
        "wall_crossings|residues|cover_arcs|Fourier_modes|proof_obligations"
    )
    print("  chosen_vertices=labelled_rational_polygon_walls")
    print("  pairwise_observable=exact rational wall-coordinate order")
    print("  switch_gauge=orient u->v iff t(u)<t(v)")
    print("  tie_hamiltonian_path=coalesce equal walls, then sorted wall word")
    print(f"  wall_vertices={len(ordered_walls)}")
    print(f"  score_histogram={{0..{len(ordered_walls)-1}: multiplicity 1}}")
    print("  directed_cycles=0")
    print(f"  strongly_connected_components={len(ordered_walls)} singleton")
    print("  hamiltonian_path_count=1")
    print("  edge_flip_locus=exact wall coincidence")
    print(
        "  preserved_by_labelled_polygon=exact remainder, components, mass, "
        "longest gap, owners, slope recurrence"
    )
    print(
        "  destroyed_by_order_only=wall slopes, affine owners, strip "
        "adjacency, metric lengths, core mask"
    )
    print(
        "  challenged_assumption=tournament vertices need not be runners; "
        "the faithful object is a Kakeya-like slope line through a fixed "
        "labelled torus polygon"
    )
    print("certificate=PASS")


if __name__ == "__main__":
    main()
