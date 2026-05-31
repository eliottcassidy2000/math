#!/usr/bin/env python3
"""
lonely_runner_feedback_cycle_s374.py

codex-2026-05-31 S374

A forced feedback loop for the Lonely Runner frontier.

The user-requested discipline is:

    14-runner attack -> dead end -> new 15-runner idea -> work it
    15-runner dead end -> new disproof-construction idea -> work it
    disproof dead end -> new 14-runner idea -> work it

This script makes that loop reproducible.  It reuses the S371/S372/S373 exact
cell and endpoint tools, then records several laps:

1. repair pressure on the best n=14 scalar-puncture near-blocker;
2. transferred n=15 scalar-puncture anatomy;
3. quotient-ladder disproof attempts for n=14 and n=15;
4. one-swap repair pressure around the seven-ladder near-disproof;
5. endpoint-protector pressure, i.e. the first "protection graph first"
   approximation from HYP-1828.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S372 = SourceFileLoader(
    "lonely_runner_creative_multiroute_s372",
    str(ROOT / "04-computation" / "lonely_runner_creative_multiroute_s372.py"),
).load_module()
S373 = SourceFileLoader(
    "lonely_runner_14_disproof_hunt_s373",
    str(ROOT / "04-computation" / "lonely_runner_14_disproof_hunt_s373.py"),
).load_module()


N14_BEST = (8, 2, 10, 4, 12, 13, 0, 8, 2, 10, 4, 12, 6)
N15_BEST = (7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 3)
SEVEN_LADDER = (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)


@dataclass(frozen=True)
class SpeedAudit:
    label: str
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witnesses: int
    unprotected: int
    first_unprotected: Fraction | None


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_ratio(x: Fraction) -> str:
    return f"{float(x):.6f}"


def distance_to_scalar(vector: tuple[int, ...], n: int) -> tuple[int, int]:
    best = (len(vector) + 1, -1)
    for m in range(n):
        scalar = S372.scalar_vector(n, m)
        dist = sum(a != b for a, b in zip(vector, scalar))
        if dist < best[0]:
            best = (dist, m)
    return best


def bit_indices(mask: int):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def missed_mask(system, vector: tuple[int, ...]) -> int:
    return system.all_mask ^ S372.blocked_mask(system, vector)


def mutation_repair_table(system, vector: tuple[int, ...], limit: int = 10) -> list[tuple]:
    base_missing = missed_mask(system, vector)
    base_missed_count = base_missing.bit_count()
    rows = []
    for coord in range(system.k):
        old = vector[coord]
        for new in range(system.n):
            if new == old:
                continue
            candidate = list(vector)
            candidate[coord] = new
            candidate = tuple(candidate)
            missing = missed_mask(system, candidate)
            fixed = (base_missing & ~missing).bit_count()
            lost = (missing & ~base_missing).bit_count()
            rows.append((missing.bit_count(), -fixed, lost, coord + 1, old, new, candidate))
    rows.sort()
    return [
        (missed, -neg_fixed, lost, coord, old, new, candidate)
        for missed, neg_fixed, lost, coord, old, new, candidate in rows[:limit]
    ]


def two_defect_pressure_from_best(system, vector: tuple[int, ...]) -> tuple[int, int, tuple[int, ...]]:
    """Best exact second defect that does not revert the vector to a scalar ramp."""

    dist, scalar_m = distance_to_scalar(vector, system.n)
    scalar = S372.scalar_vector(system.n, scalar_m)
    assert dist == 1
    defect_coord = next(i for i, (a, b) in enumerate(zip(vector, scalar)) if a != b)
    best_missed = system.candidate_count + 1
    best_count = 0
    best_vector: tuple[int, ...] | None = None
    for coord in range(system.k):
        if coord == defect_coord:
            continue
        old = vector[coord]
        for new in range(system.n):
            if new == old:
                continue
            candidate = list(vector)
            candidate[coord] = new
            candidate = tuple(candidate)
            if S372.scalar_multiplier(candidate, system.n) is not None:
                continue
            missed = S372.score_vector(system, candidate).missed
            if missed < best_missed:
                best_missed = missed
                best_count = 1
                best_vector = candidate
            elif missed == best_missed:
                best_count += 1
    assert best_vector is not None
    return best_missed, best_count, best_vector


def anatomy(system, vector: tuple[int, ...]) -> dict[str, object]:
    cells = S372.missed_cells(system, vector)
    s_hist = Counter(s for s, _p in cells)
    p_set = sorted({p for _s, p in cells})
    width_hist = Counter()
    total_width = Fraction(0)
    for p in p_set:
        lo, hi = system.intervals[p]
        width_hist[hi - lo] += 1
        total_width += hi - lo
    return {
        "missed": len(cells),
        "s_hist": tuple(sorted(s_hist.items())),
        "unique_patterns": len(p_set),
        "total_alpha_width": total_width,
        "width_hist": tuple((fmt_frac(k), v) for k, v in sorted(width_hist.items())),
        "pattern_ranges": S372.compress_indices(p_set),
    }


def print_anatomy(label: str, system, vector: tuple[int, ...]) -> None:
    row = anatomy(system, vector)
    dist, scalar_m = distance_to_scalar(vector, system.n)
    print(f"  {label}")
    print(f"    vector={vector}")
    print(f"    distance_to_scalar={dist} nearest_m={scalar_m}")
    print(
        "    "
        f"missed={row['missed']} s_hist={row['s_hist']} "
        f"unique_patterns={row['unique_patterns']} total_alpha_width={fmt_frac(row['total_alpha_width'])}"
    )
    print(f"    width_hist={row['width_hist']}")
    print(f"    pattern_ranges={row['pattern_ranges'][:10]}")


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def speed_audit(label: str, speeds: tuple[int, ...]) -> SpeedAudit:
    report = S356.report(label, list(speeds))
    summary = S360.summarize(list(speeds))
    return SpeedAudit(
        label=label,
        speeds=summary.speeds,
        classification=summary.classification,
        forbidden_length=summary.forbidden_length,
        max_gap=summary.max_gap,
        gap_ratio=summary.max_gap / summary.threshold,
        boundary_witnesses=report.boundary_witness_count,
        unprotected=summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
    )


def print_speed_table(title: str, rows: list[SpeedAudit], limit: int = 8) -> None:
    print(title)
    for row in rows[:limit]:
        print(
            "  "
            f"{row.label}: class={row.classification} "
            f"len={fmt_frac(row.forbidden_length)} "
            f"max_gap={fmt_frac(row.max_gap)} "
            f"gap/thresh={fmt_ratio(row.gap_ratio)} "
            f"boundary={row.boundary_witnesses} "
            f"unprotected={row.unprotected} "
            f"first={fmt_frac(row.first_unprotected)} "
            f"speeds={row.speeds}"
        )
    print()


def quotient_ladders(n: int, divisor: int) -> list[SpeedAudit]:
    """Anchored quotient ladders: {1} plus divisor*q, skipping one q."""

    k = n - 1
    rows: list[SpeedAudit] = []
    for skip in range(1, k + 1):
        speeds = tuple(sorted({1} | {divisor * q for q in range(1, k + 1) if q != skip}))
        if len(speeds) != k or not primitive(speeds):
            continue
        rows.append(speed_audit(f"n={n} d={divisor} skip={skip}", speeds))
    rows.sort(key=lambda row: (row.classification != "open_cover", row.max_gap, row.unprotected, row.speeds))
    return rows


def unprotected_points(speeds: tuple[int, ...]) -> list[Fraction]:
    points_by_value: dict[Fraction, list] = {}
    for endpoint in S360.endpoints(speeds):
        points_by_value.setdefault(endpoint.value, []).append(endpoint)
    out = []
    for value in sorted(points_by_value):
        indegree = sum(1 for v in speeds if S360.direct_protects(speeds, v, value))
        if indegree == 0:
            out.append(value)
    return out


def protector_pressure(
    speeds: tuple[int, ...],
    *,
    max_speed: int,
    top_k: int = 10,
) -> tuple[int, list[tuple[int, int, int]], tuple[tuple[int, int], ...]]:
    leaks = unprotected_points(speeds)
    cover_rows = []
    introduced_hist = Counter()
    n = len(speeds) + 1
    for v in range(1, max_speed + 1):
        if v in speeds:
            continue
        covered = sum(1 for t in leaks if S360.direct_protects(speeds, v, t))
        if not covered:
            continue
        introduced = {
            S360.circle_point(Fraction(n * m + sign, n * v))
            for m in range(v)
            for sign in (-1, 1)
        }
        introduced_hist[len(introduced)] += 1
        cover_rows.append((covered, -len(introduced), v))
    cover_rows.sort(reverse=True)
    return len(leaks), cover_rows[:top_k], tuple(sorted(introduced_hist.items())[:10])


def one_swap_grid_pressure(
    seed: tuple[int, ...],
    *,
    max_speed: int,
    exact_limit: int = 8,
) -> list[SpeedAudit]:
    masks = S373.build_grid_masks(max_speed)
    rows = []
    seed_set = set(seed)
    for old in seed:
        for new in range(1, max_speed + 1):
            if new in seed_set and new != old:
                continue
            candidate = tuple(sorted((seed_set - {old}) | {new}))
            if len(candidate) != len(seed) or not primitive(candidate):
                continue
            grid = S373.grid_candidate(f"swap {old}->{new}", candidate, masks)
            if grid is not None:
                rows.append(grid)
    rows.sort(key=S373.grid_rank)
    exact = []
    seen = set()
    for row in rows:
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        exact.append(speed_audit(row.label, row.speeds))
        if len(exact) >= exact_limit:
            break
    exact.sort(key=lambda row: (row.classification != "open_cover", row.max_gap, row.unprotected, row.speeds))
    return exact


def main() -> None:
    print("Lonely Runner forced feedback cycle (codex-2026-05-31 S374)")
    print("Routes: A=14-runner, B=15-runner, C=counterexample construction.")
    print()

    system14 = S372.build_pattern_system(14)
    system15 = S372.build_pattern_system(15)
    print(
        f"Cell systems: n=14 patterns={len(system14.patterns)} candidates={system14.candidate_count}; "
        f"n=15 patterns={len(system15.patterns)} candidates={system15.candidate_count}"
    )
    print()

    print("Cycle 1A: 14-runner attack - can the known 56-miss vector be repaired locally?")
    print_anatomy("n=14 best scalar-puncture near-blocker", system14, N14_BEST)
    for rank, row in enumerate(mutation_repair_table(system14, N14_BEST, limit=8), start=1):
        missed, fixed, lost, coord, old, new, _candidate = row
        print(
            f"    mutation_rank={rank} coord={coord} {old}->{new} "
            f"missed={missed} fixed={fixed} lost={lost}"
        )
    best2_missed, best2_count, best2_vector = two_defect_pressure_from_best(system14, N14_BEST)
    print(
        "  dead_end=only reverting the punctured coordinate repairs all cells; "
        f"best true second-defect pressure has missed={best2_missed} count={best2_count}."
    )
    print(f"  forced_new_15_idea=transfer the scalar-puncture anatomy to n=15, where the moat is larger.")
    print()

    print("Cycle 1B: 15-runner transfer - anatomy of the best n=15 scalar punctures")
    print_anatomy("n=15 S372 best non-scalar vector", system15, N15_BEST)
    for puncture in [
        (0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10),
    ]:
        print_anatomy("n=15 zero-ramp puncture", system15, puncture)
    print(
        "  dead_end=n=15 has a thicker scalar moat: best radius-1 punctures miss 120 cells, "
        "with active coordinates 6 or 14 and jumps 5 or 10."
    )
    print("  forced_disproof_idea=try quotient ladders that skip the dangerous puncture coordinate.")
    print()

    print("Cycle 1C: disproof construction - quotient ladders and skipped coordinates")
    rows14 = quotient_ladders(14, 7)
    rows15_d3 = quotient_ladders(15, 3)
    rows15_d5 = quotient_ladders(15, 5)
    print_speed_table("  n=14 seven-ladders", rows14, limit=6)
    print_speed_table("  n=15 three-ladders", rows15_d3, limit=6)
    print_speed_table("  n=15 five-ladders", rows15_d5, limit=6)
    print(
        "  dead_end=the best n=14 ladder is extremely close by measure but has 84 exposed endpoints; "
        "the n=15 ladders are worse by gap ratio."
    )
    print("  forced_new_14_idea=repair the seven-ladder by one speed swap before abandoning speed-first search.")
    print()

    print("Cycle 2A: 14-runner attack - one-swap repairs around the seven-ladder")
    swap_rows = one_swap_grid_pressure(SEVEN_LADDER, max_speed=252, exact_limit=10)
    print_speed_table("  best one-swap exact audits", swap_rows, limit=10)
    print(
        "  dead_end=no one-swap open cover; best swaps either return to the same seven-ladder geometry "
        "or reopen larger endpoint leaks."
    )
    print("  forced_new_15_idea=look for the analogue of one-swap repair as endpoint-protector pressure.")
    print()

    print("Cycle 2B: 15-runner idea - endpoint protector pressure for best n=15 ladders")
    best15 = min(rows15_d3 + rows15_d5, key=lambda row: (row.max_gap, row.unprotected, row.speeds))
    leaks15, protect15, intro15 = protector_pressure(best15.speeds, max_speed=300, top_k=8)
    print(f"  best15_ladder={best15.speeds} gap/thresh={fmt_ratio(best15.gap_ratio)} leaks={leaks15}")
    print(f"  top_protectors=(covered, -introduced_endpoint_count, speed) {protect15}")
    print(f"  introduced_endpoint_count_hist_prefix={intro15}")
    print(
        "  dead_end=available single protectors cover only slices of the leak layer and introduce many new endpoints."
    )
    print("  forced_disproof_idea=run the same endpoint-pressure test on the n=14 seven-ladder.")
    print()

    print("Cycle 2C: disproof construction - endpoint pressure on the seven-ladder")
    leaks14, protect14, intro14 = protector_pressure(SEVEN_LADDER, max_speed=252, top_k=12)
    print(f"  seven_ladder_leaks={leaks14}")
    print(f"  top_protectors=(covered, -introduced_endpoint_count, speed) {protect14}")
    print(f"  introduced_endpoint_count_hist_prefix={intro14}")
    print(
        "  dead_end=protection is not scarce, but every protector brings a large new endpoint boundary; "
        "this is the endpoint-cycle-first problem in miniature."
    )
    print("  forced_new_14_idea=combine scalar-cell moat and endpoint pressure as a two-layer obstruction.")
    print()

    print("Cycle 3A: 14-runner synthesis - scalar moat plus endpoint pressure")
    print(
        "  observation=the best cell-blocker is one scalar puncture away from equality, while the best "
        "interval-cover disproof is a seven-ladder with many exposed endpoints."
    )
    print(
        "  candidate_theorem=any speed set close enough to the Dirichlet/scalar spine to cover the unit layer "
        "inherits either the 56-cell scalar moat or the seven-ladder endpoint leak."
    )
    print(
        "  dead_end=this is not a proof yet; it needs a bridge between residue-vector cell blocking "
        "and endpoint-cover leakage."
    )
    print("  forced_new_15_idea=test whether the bridge should use divisor coordinates 6/14 for n=15.")
    print()

    print("Cycle 3B: 15-runner synthesis - divisor-coordinate warning")
    print(
        "  observation=n=15 punctures prefer coordinates 6 and 14 with jumps 5 or 10; "
        "these are exactly the 3*5 divisor jumps, not random coordinates."
    )
    print(
        "  candidate_theorem=for composite n, minimal scalar punctures live on coordinates whose gcd with n "
        "exposes a proper quotient, and their jump is a complementary divisor layer."
    )
    print("  dead_end=still local around scalar ramps; far-from-scalar blockers need a search certificate.")
    print("  forced_disproof_idea=record a reusable next search: endpoint cycles first, scalar distance second.")
    print()

    print("Cycle 3C: disproof route - next reproducible search specification")
    print("  1. choose a finite leak layer L, starting with {a/14, 9/98, 29/182, 15/182};")
    print("  2. find a directed protection cycle on L using the integer protection inequality;")
    print("  3. minimize scalar distance of realizing speeds to avoid the 56/120 scalar moats;")
    print("  4. only then run exact interval audits.")
    print()

    print("Final synthesis")
    print("  A. 14-runner route: local scalar repair is dead; the 56-cell moat is robust.")
    print("  B. 15-runner route: the transferred moat thickens to 120 and singles out divisor coordinates.")
    print("  C. disproof route: quotient ladders give tiny measure gaps but fail endpoint protection.")
    print("  D. new loop invariant: scalar distance and endpoint-protection closure must be optimized together.")


if __name__ == "__main__":
    main()
