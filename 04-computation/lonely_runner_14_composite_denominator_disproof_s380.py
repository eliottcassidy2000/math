#!/usr/bin/env python3
"""
lonely_runner_14_composite_denominator_disproof_s380.py

codex-2026-05-31 S380

Attempt a fourteen-runner disproof by leaning into the composite denominator
anomaly at n=14.  Previous sessions found two relevant facts:

* THM-360 forces any primitive open-cover counterexample to protect the unit
  endpoints a/14, hence to contain a 14-gate in the speed set.
* The scalar-puncture moat and quotient-ladder near misses expose denominator
  98/182 leak layers.  The best known speed-side near miss is the seven-ladder
  with a tiny positive gap but many exposed endpoints.

This script tries the disproof direction directly.  It builds a target layer
from unit endpoints, known 98/182 leak points, and the largest exact gaps of
the seven-ladder/gate seeds.  It then searches 13-speed sets that maximize
grid cover while protecting these composite-denominator targets, before
auditing the best candidates exactly by interval union and endpoint protection.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import ceil, floor, gcd
from pathlib import Path
import random


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()


K = 13
N = 14
MAX_SPEED = 420
# 35672 = 28 * lcm(14, 98, 182).  This samples the composite leak layers
# exactly while keeping bitset search cheap.
GRID = 35672
ONE = Fraction(1, 1)


@dataclass(frozen=True)
class Target:
    label: str
    point: Fraction
    weight: int


@dataclass(frozen=True)
class GridRow:
    label: str
    speeds: tuple[int, ...]
    uncovered: int
    longest_gap_grid: int
    target_weight: int
    target_uncovered: int
    composite_bonus: int


@dataclass(frozen=True)
class ExactRow:
    label: str
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witness_count: int
    unprotected_count: int
    first_unprotected: Fraction | None
    endpoint_q: int
    max_gap_ends: tuple[Fraction, Fraction] | None
    gap_endpoint_denominators: tuple[tuple[int, int], ...]
    residue_mod_14: tuple[tuple[int, int], ...]
    gcd_mod_14: tuple[tuple[int, int], ...]


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_float(x: Fraction) -> str:
    return f"{float(x):.6f}"


def circle(x: Fraction) -> Fraction:
    return x % ONE


def circular_distance(x: Fraction) -> Fraction:
    r = circle(x)
    return min(r, ONE - r)


def protects(speed: int, point: Fraction) -> bool:
    return circular_distance(speed * point) < Fraction(1, N)


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def canonical(raw: tuple[int, ...] | list[int]) -> tuple[int, ...] | None:
    speeds = tuple(sorted(set(int(v) for v in raw if v > 0)))
    if len(speeds) != K or not primitive(speeds):
        return None
    return speeds


def forbidden_mask(speed: int) -> int:
    mask = 0
    radius = Fraction(1, N * speed)
    for center_idx in range(speed):
        center = Fraction(center_idx, speed)
        for lo, hi in S356.split_circle_interval(center - radius, center + radius):
            # Grid points are j/GRID.  Endpoints are intentionally excluded,
            # except that wrapped intervals contain the circle point 0.
            if lo == 0:
                mask |= 1
            start = floor(lo * GRID) + 1
            end = ceil(hi * GRID) - 1
            if hi == 1:
                end = GRID - 1
            if 0 <= start <= end:
                length = end - start + 1
                mask |= ((1 << length) - 1) << start
    return mask


def build_masks(max_speed: int) -> dict[int, int]:
    return {speed: forbidden_mask(speed) for speed in range(1, max_speed + 1)}


def combined_mask(speeds: tuple[int, ...], masks: dict[int, int]) -> int:
    out = 0
    for speed in speeds:
        out |= masks[speed]
    return out


def longest_zero_run(mask: int) -> int:
    uncovered = [((mask >> j) & 1) == 0 for j in range(GRID)]
    if not any(uncovered):
        return 0
    doubled = uncovered + uncovered
    best = 0
    cur = 0
    for value in doubled:
        if value:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return min(best, GRID)


def exact_gaps(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals = S356.forbidden_intervals(speeds)
    components = S356.merge_intervals(intervals)
    return sorted(S356.circular_gaps(components), key=lambda gap: gap[1] - gap[0], reverse=True)


def orbit(point: Fraction, *, step_den: int = N) -> set[Fraction]:
    out = set()
    for j in range(step_den):
        shifted = circle(point + Fraction(j, step_den))
        out.add(shifted)
        out.add(circle(ONE - shifted))
    return out


def add_target(store: dict[Fraction, Target], point: Fraction, label: str, weight: int) -> None:
    point = circle(point)
    old = store.get(point)
    if old is None or weight > old.weight:
        store[point] = Target(label=label, point=point, weight=weight)


def seed_sets() -> list[tuple[str, tuple[int, ...]]]:
    seven = (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)
    gate = (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14)
    seeds = [
        ("initial segment", tuple(range(1, N))),
        ("single 14-gate", gate),
        ("seven-ladder", seven),
        ("even quotient ladder", tuple(sorted({1} | {2 * q for q in range(1, N) if q != 7}))),
    ]
    for skip in range(1, N):
        speeds = tuple(sorted({1} | {7 * q for q in range(1, N) if q != skip}))
        seeds.append((f"7-ladder skip={skip}", speeds))
    return [(label, speeds) for label, speeds in seeds if canonical(speeds) is not None]


def build_targets() -> list[Target]:
    targets: dict[Fraction, Target] = {}

    for a in range(1, N):
        point = Fraction(a, N)
        if gcd(a, N) == 1:
            add_target(targets, point, "unit a/14", 9)
        else:
            add_target(targets, point, "nonunit a/14", 4)

    known = [
        (Fraction(9, 98), "known 98 leak", 13),
        (Fraction(29, 182), "known 182 leak", 13),
        (Fraction(15, 182), "known 182 leak", 13),
        (Fraction(9, 56), "scalar-moat endpoint", 10),
        (Fraction(29, 182), "scalar-moat endpoint", 10),
    ]
    for point, label, weight in known:
        for image in orbit(point):
            add_target(targets, image, label, weight)

    for label, speeds in seed_sets():
        if "seven-ladder" not in label and label != "single 14-gate":
            continue
        for idx, gap in enumerate(exact_gaps(speeds)[:18], start=1):
            lo, hi = gap
            mid = circle((lo + hi) / 2)
            add_target(targets, lo, f"{label} gap{idx} left", 8)
            add_target(targets, hi, f"{label} gap{idx} right", 8)
            add_target(targets, mid, f"{label} gap{idx} mid", 12)

    return sorted(targets.values(), key=lambda item: (item.point.denominator, item.point.numerator))


def target_masks(targets: list[Target]) -> tuple[dict[int, int], dict[int, int]]:
    masks: dict[int, int] = {}
    weights: dict[int, int] = {}
    for speed in range(1, MAX_SPEED + 1):
        mask = 0
        weight = 0
        for idx, target in enumerate(targets):
            if protects(speed, target.point):
                mask |= 1 << idx
                weight += target.weight
        masks[speed] = mask
        weights[speed] = weight
    return masks, weights


def target_score(speeds: tuple[int, ...], masks: dict[int, int], total_mask: int) -> tuple[int, int]:
    mask = 0
    for speed in speeds:
        mask |= masks[speed]
    covered_weight = 0
    # Weight is recomputed by caller through the global targets list; this
    # helper returns only count-style coverage and uncovered target count.
    uncovered_count = (total_mask & ~mask).bit_count()
    return mask, uncovered_count


def composite_bonus(speeds: tuple[int, ...]) -> int:
    bonus = 0
    for speed in speeds:
        if speed % 14 == 0:
            bonus += 7
        elif speed % 7 == 0:
            bonus += 4
        elif speed % 2 == 0:
            bonus += 1
    return bonus


def row_for(
    label: str,
    raw: tuple[int, ...] | list[int],
    masks: dict[int, int],
    t_masks: dict[int, int],
    t_weights: dict[int, int],
    targets: list[Target],
    *,
    with_longest: bool = False,
) -> GridRow | None:
    speeds = canonical(raw)
    if speeds is None:
        return None
    grid_mask = combined_mask(speeds, masks)
    target_mask = 0
    for speed in speeds:
        target_mask |= t_masks[speed]
    target_weight = sum(target.weight for idx, target in enumerate(targets) if (target_mask >> idx) & 1)
    return GridRow(
        label=label,
        speeds=speeds,
        uncovered=GRID - grid_mask.bit_count(),
        longest_gap_grid=longest_zero_run(grid_mask) if with_longest else -1,
        target_weight=target_weight,
        target_uncovered=len(targets) - target_mask.bit_count(),
        composite_bonus=composite_bonus(speeds),
    )


def grid_rank(row: GridRow) -> tuple[int, int, int, int, int, tuple[int, ...]]:
    longest = row.longest_gap_grid if row.longest_gap_grid >= 0 else 0
    return (
        row.uncovered,
        longest,
        row.target_uncovered,
        -row.target_weight,
        -row.composite_bonus,
        row.speeds,
    )


def keep(rows: list[GridRow], row: GridRow | None, limit: int) -> None:
    if row is None:
        return
    rows.append(row)
    rows.sort(key=grid_rank)
    del rows[limit:]


def target_pool(t_weights: dict[int, int]) -> tuple[int, ...]:
    ranked = sorted(
        range(1, MAX_SPEED + 1),
        key=lambda speed: (
            -t_weights[speed],
            -(speed % 14 == 0),
            -(speed % 7 == 0),
            speed,
        ),
    )
    keepers = set(ranked[:90])
    keepers.update(range(1, 2 * N + 1))
    keepers.update(7 * q for q in range(1, MAX_SPEED // 7 + 1))
    keepers.update(14 * q for q in range(1, MAX_SPEED // 14 + 1))
    return tuple(sorted(v for v in keepers if 1 <= v <= MAX_SPEED))


def greedy_rows(
    masks: dict[int, int],
    t_masks: dict[int, int],
    t_weights: dict[int, int],
    targets: list[Target],
    pool: tuple[int, ...],
) -> list[GridRow]:
    rows: list[GridRow] = []
    forced_sets = [
        (14,),
        (1, 14),
        (7, 14),
        (1, 7, 14),
        (14, 28),
        (14, 98),
        (1, 14, 98),
    ]
    for forced in forced_sets:
        speeds = list(forced)
        grid_mask = combined_mask(tuple(speeds), masks)
        target_mask = 0
        for speed in speeds:
            target_mask |= t_masks[speed]
        while len(speeds) < K:
            best = None
            best_score = None
            for speed in pool:
                if speed in speeds:
                    continue
                new_grid = grid_mask | masks[speed]
                new_target = target_mask | t_masks[speed]
                score = (
                    new_grid.bit_count() - grid_mask.bit_count(),
                    sum(target.weight for idx, target in enumerate(targets) if ((new_target & ~target_mask) >> idx) & 1),
                    composite_bonus(tuple(speeds + [speed])),
                    -speed,
                )
                if best_score is None or score > best_score:
                    best_score = score
                    best = speed
            assert best is not None
            speeds.append(best)
            grid_mask |= masks[best]
            target_mask |= t_masks[best]
        keep(rows, row_for(f"target-grid greedy forced={forced}", speeds, masks, t_masks, t_weights, targets), 40)
    return rows


def splice_rows(
    masks: dict[int, int],
    t_masks: dict[int, int],
    t_weights: dict[int, int],
    targets: list[Target],
    pool: tuple[int, ...],
) -> list[GridRow]:
    rows: list[GridRow] = []
    bases = [
        ("seven-ladder", (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)),
        ("single 14-gate", (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14)),
    ]
    splice_pool = tuple(sorted(pool, key=lambda speed: (-t_weights[speed], speed))[:64])
    for base_label, base in bases:
        base_set = set(base)
        for remove_size in (1, 2):
            for removed in combinations(base, remove_size):
                stem = base_set - set(removed)
                for added in combinations(splice_pool, remove_size):
                    if stem & set(added):
                        continue
                    label = f"splice {base_label} remove={removed} add={added}"
                    keep(rows, row_for(label, tuple(sorted(stem | set(added))), masks, t_masks, t_weights, targets), 90)
    return rows


def mutation_rows(
    seeds: list[GridRow],
    masks: dict[int, int],
    t_masks: dict[int, int],
    t_weights: dict[int, int],
    targets: list[Target],
    pool: tuple[int, ...],
) -> list[GridRow]:
    rng = random.Random(14378)
    rows: list[GridRow] = []
    seed_speeds = [row.speeds for row in sorted(seeds, key=grid_rank)[:44]]
    seed_speeds.extend(speeds for _, speeds in seed_sets())
    seen: set[tuple[int, ...]] = set()
    for restart, seed in enumerate(seed_speeds[:52]):
        if seed in seen:
            continue
        seen.add(seed)
        current = row_for(f"mutation seed={restart}", seed, masks, t_masks, t_weights, targets)
        if current is None:
            continue
        best = current
        for step in range(520):
            old = rng.choice(current.speeds)
            new = rng.choice(pool)
            if new in current.speeds:
                continue
            trial = tuple(sorted((set(current.speeds) - {old}) | {new}))
            row = row_for(f"mutation restart={restart} step={step}", trial, masks, t_masks, t_weights, targets)
            if row is None:
                continue
            improves = grid_rank(row) < grid_rank(current)
            anneal = row.uncovered <= current.uncovered + 4 and rng.random() < 0.006
            if improves or anneal:
                current = row
                if grid_rank(row) < grid_rank(best):
                    best = row
        keep(rows, best, 80)
    return rows


def exact_row(row: GridRow) -> ExactRow:
    report = S356.report(row.label, list(row.speeds))
    summary = S360.summarize(list(row.speeds))
    gaps = exact_gaps(summary.speeds)
    max_gap = gaps[0] if gaps else None
    den_hist: Counter[int] = Counter()
    for lo, hi in gaps:
        den_hist[circle(lo).denominator] += 1
        den_hist[circle(hi).denominator] += 1
    return ExactRow(
        label=row.label,
        speeds=summary.speeds,
        classification=summary.classification,
        forbidden_length=summary.forbidden_length,
        max_gap=summary.max_gap,
        gap_ratio=summary.max_gap / summary.threshold,
        boundary_witness_count=report.boundary_witness_count,
        unprotected_count=summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
        endpoint_q=summary.boundary_modulus,
        max_gap_ends=(circle(max_gap[0]), circle(max_gap[1])) if max_gap else None,
        gap_endpoint_denominators=tuple(sorted(den_hist.items())),
        residue_mod_14=tuple(sorted(Counter(speed % N for speed in summary.speeds).items())),
        gcd_mod_14=tuple(sorted(Counter(gcd(speed, N) for speed in summary.speeds).items())),
    )


def exact_rank(row: ExactRow) -> tuple[int, int, Fraction, int, tuple[int, ...]]:
    class_rank = {"open_cover": 0, "boundary_only": 1, "positive_gap": 2}.get(row.classification, 3)
    return (
        class_rank,
        row.unprotected_count if row.classification != "positive_gap" else 10**9,
        row.max_gap,
        row.boundary_witness_count,
        row.speeds,
    )


def exact_audit(rows: list[GridRow], limit: int) -> list[ExactRow]:
    out: list[ExactRow] = []
    seen: set[tuple[int, ...]] = set()
    for row in sorted(rows, key=grid_rank):
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        out.append(exact_row(row))
        if len(out) >= limit:
            break
    return sorted(out, key=exact_rank)


def print_target_summary(targets: list[Target], t_weights: dict[int, int]) -> None:
    den_hist = Counter(target.point.denominator for target in targets)
    print("Composite-denominator target layer")
    print(f"  target_count={len(targets)}")
    print(f"  denominator_hist={dict(sorted(den_hist.items()))}")
    print("  top target-protecting speeds:")
    for speed in sorted(range(1, MAX_SPEED + 1), key=lambda v: (-t_weights[v], v))[:12]:
        print(
            f"    speed={speed:3d} weight={t_weights[speed]:4d} "
            f"mod14={speed % 14:2d} gcd14={gcd(speed, 14):2d}"
        )
    print()


def print_grid_table(title: str, rows: list[GridRow], limit: int = 10) -> None:
    print(title)
    for row in sorted(rows, key=grid_rank)[:limit]:
        print(
            "  "
            f"{row.label}: grid_uncovered={row.uncovered} longest={row.longest_gap_grid} "
            f"target_uncovered={row.target_uncovered} target_weight={row.target_weight} "
            f"composite_bonus={row.composite_bonus}"
        )
        print(f"    speeds={row.speeds}")
    print()


def print_exact_table(title: str, rows: list[ExactRow], limit: int = 10) -> None:
    print(title)
    if not rows:
        print("  (none)")
        print()
        return
    for row in rows[:limit]:
        ends = "-"
        if row.max_gap_ends is not None:
            ends = f"{fmt_frac(row.max_gap_ends[0])}->{fmt_frac(row.max_gap_ends[1])}"
        print(
            "  "
            f"{row.label}: class={row.classification} "
            f"forbidden={fmt_frac(row.forbidden_length)} "
            f"max_gap={fmt_frac(row.max_gap)} gap/th={fmt_float(row.gap_ratio)} "
            f"boundary_wit={row.boundary_witness_count} unprot={row.unprotected_count} "
            f"first_unprot={fmt_frac(row.first_unprotected)} Q={row.endpoint_q}"
        )
        print(f"    max_gap_ends={ends} gap_denoms={row.gap_endpoint_denominators[:8]}")
        print(f"    mod14={row.residue_mod_14} gcd14={row.gcd_mod_14}")
        print(f"    speeds={row.speeds}")
    print()


def main() -> None:
    print("Fourteen-runner composite-denominator disproof attempt (S380)")
    print(f"Convention: k={K}, threshold=1/{N}, speed search <= {MAX_SPEED}, grid={GRID}.")
    print("A disproof would be an exact open cover: full forbidden length and no unprotected endpoint.")
    print()

    print("Building grid masks...")
    masks = build_masks(MAX_SPEED)
    targets = build_targets()
    t_masks, t_weights = target_masks(targets)
    pool = target_pool(t_weights)
    print(f"  pool_size={len(pool)}")
    print_target_summary(targets, t_weights)

    seed_rows = [
        row_for(label, speeds, masks, t_masks, t_weights, targets)
        for label, speeds in seed_sets()
    ]
    seed_rows = [row for row in seed_rows if row is not None]
    print_grid_table("Seed composite/gate families on the target grid", seed_rows, limit=10)
    print_exact_table("Exact audit of seed families", exact_audit(seed_rows, limit=12), limit=8)

    greedy = greedy_rows(masks, t_masks, t_weights, targets, pool)
    print_grid_table("Target-grid greedy constructions", greedy, limit=8)
    print_exact_table("Exact audit of greedy constructions", exact_audit(greedy, limit=14), limit=8)

    spliced = splice_rows(masks, t_masks, t_weights, targets, pool)
    print_grid_table("Composite leak splices around the seven-ladder/gate seeds", spliced, limit=10)
    print_exact_table("Exact audit of splice constructions", exact_audit(spliced, limit=24), limit=10)

    mutation_input = seed_rows + greedy + spliced
    mutated = mutation_rows(mutation_input, masks, t_masks, t_weights, targets, pool)
    print_grid_table("Composite-denominator mutation search", mutated, limit=10)
    print_exact_table("Exact audit of mutation constructions", exact_audit(mutated, limit=30), limit=10)

    all_rows = seed_rows + greedy + spliced + mutated
    audited = exact_audit(all_rows, limit=90)
    print_exact_table("Best exact audited candidates overall", audited, limit=14)

    open_rows = [row for row in audited if row.classification == "open_cover"]
    boundary_rows = [row for row in audited if row.classification == "boundary_only"]
    positive_rows = [row for row in audited if row.classification == "positive_gap"]
    print("Disproof ledger")
    print(f"  open_cover_candidates={len(open_rows)}")
    print(f"  boundary_only_candidates={len(boundary_rows)}")
    print(f"  positive_gap_candidates={len(positive_rows)}")
    if positive_rows:
        best = min(positive_rows, key=lambda row: (row.max_gap, row.boundary_witness_count))
        print(
            "  best_positive_gap="
            f"{best.label}; gap/th={fmt_float(best.gap_ratio)}; "
            f"max_gap={fmt_frac(best.max_gap)}; ends={best.max_gap_ends}; speeds={best.speeds}"
        )
    if boundary_rows:
        best_boundary = min(boundary_rows, key=lambda row: (row.unprotected_count, row.boundary_witness_count))
        print(
            "  best_boundary_only="
            f"{best_boundary.label}; unprotected={best_boundary.unprotected_count}; "
            f"boundary_witnesses={best_boundary.boundary_witness_count}; speeds={best_boundary.speeds}"
        )
    print()

    print("Interpretation")
    print("  The composite-denominator route did not produce a fourteen-runner counterexample.")
    print("  It did stress the right anomaly: protecting the unit layer with 14-gates and")
    print("  targeting the 98/182 moat leaks keeps returning positive gaps instead of an")
    print("  open cover.  The best audited positive gaps still have endpoint denominators")
    print("  tied to the 2*7 and 2*7*13 layers, so the anomaly is real, but it behaves")
    print("  like an endpoint-transfer obstruction rather than a disproof architecture.")
    print("  Next disproof attempt should solve the endpoint protection graph first, with")
    print("  variables constrained to protect the 98/182 leak orbit, then realize those")
    print("  protection edges by speeds.")


if __name__ == "__main__":
    main()
