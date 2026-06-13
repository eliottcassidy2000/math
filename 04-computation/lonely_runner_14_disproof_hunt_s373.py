#!/usr/bin/env python3
"""
lonely_runner_14_disproof_hunt_s373.py

codex-2026-05-31 S373

Attempted disproof constructions for the fourteen-runner Lonely Runner case.

Repository convention: fourteen total runners means k=13 nonzero speeds and
loneliness threshold 1/14.  A disproof would be a primitive 13-speed set whose
forbidden arcs

    ||v t|| < 1/14

form a full open cover of R/Z.

This script is not a proof search in the positive direction.  It deliberately
tries to manufacture counterexamples using several construction grammars:

1. gate replacements of the Dirichlet-equality initial segment,
2. multi-gate overloads,
3. CRT residue lifts, including the S367 half-turn near-blocker residues,
4. greedy set-cover pressure on a fine grid,
5. random mutation with a mandatory 14-gate.

Every near-miss reported at the end is then audited by the exact Fraction
interval and endpoint-protection classifiers from S356/S360.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
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
GRID = 12012
MAX_SPEED = 252


@dataclass(frozen=True)
class GridCandidate:
    label: str
    speeds: tuple[int, ...]
    covered: int
    uncovered: int
    longest_uncovered_run: int


@dataclass(frozen=True)
class ExactCandidate:
    label: str
    speeds: tuple[int, ...]
    grid_uncovered: int
    grid_longest: int
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witness_count: int
    unprotected_count: int
    first_unprotected: Fraction | None
    endpoint_q: int


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_ratio(x: Fraction) -> str:
    return f"{float(x):.6f}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def canonical_speeds(raw: list[int] | tuple[int, ...]) -> tuple[int, ...] | None:
    if len(set(raw)) != K:
        return None
    speeds = tuple(sorted(raw))
    if not primitive(speeds):
        return None
    normalized = S356.normalize_speed_set(list(speeds))
    if len(normalized) != K:
        return None
    return speeds


def grid_forbidden_mask(speed: int) -> int:
    mask = 0
    for j in range(GRID):
        residue = (speed * j) % GRID
        distance = min(residue, GRID - residue)
        if N * distance < GRID:
            mask |= 1 << j
    return mask


def build_grid_masks(max_speed: int) -> dict[int, int]:
    return {v: grid_forbidden_mask(v) for v in range(1, max_speed + 1)}


def combined_mask(speeds: tuple[int, ...], masks: dict[int, int]) -> int:
    mask = 0
    for v in speeds:
        mask |= masks[v]
    return mask


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
            if cur > best:
                best = cur
        else:
            cur = 0
    return min(best, GRID)


def grid_candidate(
    label: str,
    raw: list[int] | tuple[int, ...],
    masks: dict[int, int],
    *,
    with_longest: bool = False,
) -> GridCandidate | None:
    speeds = canonical_speeds(raw)
    if speeds is None:
        return None
    mask = combined_mask(speeds, masks)
    covered = mask.bit_count()
    return GridCandidate(
        label=label,
        speeds=speeds,
        covered=covered,
        uncovered=GRID - covered,
        longest_uncovered_run=longest_zero_run(mask) if with_longest else -1,
    )


def exact_candidate(row: GridCandidate) -> ExactCandidate:
    report = S356.report(row.label, list(row.speeds))
    summary = S360.summarize(list(row.speeds))
    return ExactCandidate(
        label=row.label,
        speeds=summary.speeds,
        grid_uncovered=row.uncovered,
        grid_longest=row.longest_uncovered_run,
        classification=summary.classification,
        forbidden_length=summary.forbidden_length,
        max_gap=summary.max_gap,
        gap_ratio=summary.max_gap / summary.threshold,
        boundary_witness_count=report.boundary_witness_count,
        unprotected_count=summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
        endpoint_q=summary.boundary_modulus,
    )


def exact_rank(row: ExactCandidate) -> tuple:
    class_rank = {
        "open_cover": 0,
        "boundary_only": 1,
        "positive_gap": 2,
    }.get(row.classification, 3)
    return (
        class_rank,
        row.unprotected_count if row.classification != "positive_gap" else 10**9,
        row.max_gap,
        row.boundary_witness_count,
        row.grid_uncovered,
        row.speeds,
    )


def grid_rank(row: GridCandidate) -> tuple:
    longest = row.longest_uncovered_run if row.longest_uncovered_run >= 0 else 0
    return (row.uncovered, longest, row.speeds)


def keep_best(
    rows: list[GridCandidate],
    row: GridCandidate | None,
    *,
    limit: int,
) -> None:
    if row is None:
        return
    rows.append(row)
    rows.sort(key=grid_rank)
    del rows[limit:]


def exact_audit(
    rows: list[GridCandidate],
    *,
    limit: int,
    masks: dict[int, int],
) -> list[ExactCandidate]:
    seen: set[tuple[int, ...]] = set()
    out: list[ExactCandidate] = []
    for row in sorted(rows, key=grid_rank):
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        if row.longest_uncovered_run < 0:
            refreshed = grid_candidate(row.label, row.speeds, masks, with_longest=True)
            assert refreshed is not None
            row = refreshed
        out.append(exact_candidate(row))
        if len(out) >= limit:
            break
    return sorted(out, key=exact_rank)


def gate_replacements(masks: dict[int, int]) -> list[GridCandidate]:
    base = set(range(1, N))
    rows: list[GridCandidate] = []
    for removed in range(1, N):
        for q in range(1, 19):
            gate = N * q
            speeds = sorted((base - {removed}) | {gate})
            label = f"single-gate remove={removed} add={gate}"
            keep_best(rows, grid_candidate(label, speeds, masks), limit=40)
    return rows


def multi_gate_replacements(masks: dict[int, int]) -> list[GridCandidate]:
    base = set(range(1, N))
    rows: list[GridCandidate] = []
    for gate_count, q_max in ((2, 16), (3, 12)):
        gates = [N * q for q in range(1, q_max + 1)]
        for removed in combinations(range(1, N), gate_count):
            for added in combinations(gates, gate_count):
                speeds = sorted((base - set(removed)) | set(added))
                label = (
                    f"{gate_count}-gate remove={','.join(map(str, removed))} "
                    f"add={','.join(map(str, added))}"
                )
                keep_best(rows, grid_candidate(label, speeds, masks), limit=50)
    return rows


def lift_residue_pattern(
    residues: tuple[int, ...],
    rng: random.Random,
    *,
    height: int,
) -> tuple[int, ...] | None:
    speeds: list[int] = []
    used: set[int] = set()
    for residue in residues:
        for _ in range(40):
            q = rng.randrange(1, height + 1) if residue == 0 else rng.randrange(height + 1)
            speed = N * q if residue == 0 else residue + N * q
            if 1 <= speed <= MAX_SPEED and speed not in used:
                speeds.append(speed)
                used.add(speed)
                break
        else:
            return None
    return tuple(sorted(speeds))


def crt_lift_trials(masks: dict[int, int]) -> list[GridCandidate]:
    patterns: list[tuple[str, tuple[int, ...]]] = [
        ("residue-initial", tuple(range(1, N))),
        ("residue-gated-no-seven", (0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13)),
        ("s367-raw-halfturn", (8, 2, 10, 4, 12, 13, 0, 8, 2, 10, 4, 12, 6)),
        ("mostly-gates-one-seven", (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7)),
        ("unit-skeleton-plus-even", (0, 1, 3, 5, 9, 11, 13, 2, 4, 6, 8, 10, 12)),
        ("two-sheeted-42", (0, 0, 7, 7, 1, 3, 5, 9, 11, 13, 2, 6, 10)),
    ]
    rng = random.Random(37142)
    rows: list[GridCandidate] = []
    for name, residues in patterns:
        for trial in range(900):
            speeds = lift_residue_pattern(residues, rng, height=13)
            if speeds is None or not primitive(speeds):
                continue
            label = f"crt-lift {name} trial={trial}"
            keep_best(rows, grid_candidate(label, speeds, masks), limit=70)
    return rows


def greedy_cover(
    masks: dict[int, int],
    *,
    pool: tuple[int, ...],
    forced: tuple[int, ...],
    label: str,
) -> GridCandidate | None:
    speeds = list(forced)
    mask = combined_mask(tuple(speeds), masks) if speeds else 0
    while len(speeds) < K:
        best_v = None
        best_count = -1
        for v in pool:
            if v in speeds:
                continue
            count = (mask | masks[v]).bit_count()
            if count > best_count:
                best_count = count
                best_v = v
        assert best_v is not None
        speeds.append(best_v)
        mask |= masks[best_v]

    # One-pass replacement cleanup on the same grid objective.
    changed = True
    while changed:
        changed = False
        current = combined_mask(tuple(speeds), masks).bit_count()
        for idx, old in enumerate(list(speeds)):
            if old in forced:
                continue
            for v in pool:
                if v in speeds:
                    continue
                candidate = list(speeds)
                candidate[idx] = v
                count = combined_mask(tuple(candidate), masks).bit_count()
                if count > current:
                    speeds = candidate
                    current = count
                    changed = True
    return grid_candidate(label, speeds, masks)


def greedy_families(masks: dict[int, int]) -> list[GridCandidate]:
    pool = tuple(range(1, MAX_SPEED + 1))
    rows = [
        greedy_cover(masks, pool=pool, forced=(1,), label="greedy anchored 1"),
        greedy_cover(masks, pool=pool, forced=(1, N), label="greedy anchored 1 forced 14"),
        greedy_cover(masks, pool=pool, forced=(1, 2 * N), label="greedy anchored 1 forced 28"),
        greedy_cover(masks, pool=pool, forced=(1, N, 2 * N), label="greedy anchored 1 forced 14,28"),
    ]
    return [row for row in rows if row is not None]


def mutate_seed(
    seed: tuple[int, ...],
    masks: dict[int, int],
    rng: random.Random,
    *,
    pool: tuple[int, ...],
    steps: int,
    require_gate: bool,
) -> GridCandidate:
    current = tuple(sorted(seed))
    current_mask = combined_mask(current, masks)
    current_score = current_mask.bit_count()
    best = grid_candidate("mutation seed", current, masks)
    assert best is not None

    for step in range(steps):
        idx = rng.randrange(K)
        old = current[idx]
        v = rng.choice(pool)
        if v in current:
            continue
        candidate = list(current)
        candidate[idx] = v
        candidate = sorted(candidate)
        if require_gate and not any(x % N == 0 for x in candidate):
            continue
        if canonical_speeds(candidate) is None or not primitive(tuple(candidate)):
            continue
        candidate_tuple = tuple(candidate)
        candidate_mask = combined_mask(candidate_tuple, masks)
        candidate_score = candidate_mask.bit_count()
        temperature_accept = rng.random() < 0.002 and candidate_score + 6 >= current_score
        if candidate_score > current_score or temperature_accept:
            current = candidate_tuple
            current_mask = candidate_mask
            current_score = candidate_score
            row = grid_candidate(f"mutation step={step}", current, masks)
            assert row is not None
            if grid_rank(row) < grid_rank(best):
                best = row
    return best


def mutation_search(masks: dict[int, int], seeds: list[GridCandidate]) -> list[GridCandidate]:
    rng = random.Random(37199)
    pool = tuple(range(1, MAX_SPEED + 1))
    seed_sets = [row.speeds for row in seeds]
    seed_sets.extend(
        [
            tuple(range(1, N)),
            tuple(sorted([1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14])),
            tuple(sorted([1, 3, 5, 7, 9, 11, 13, 14, 28, 42, 56, 70, 84])),
        ]
    )
    rows: list[GridCandidate] = []
    for idx, seed in enumerate(seed_sets[:34]):
        if canonical_speeds(seed) is None:
            continue
        row = mutate_seed(
            seed,
            masks,
            rng,
            pool=pool,
            steps=420,
            require_gate=True,
        )
        keep_best(rows, GridCandidate(f"mutation restart={idx}", row.speeds, row.covered, row.uncovered, row.longest_uncovered_run), limit=50)
    return rows


def print_exact_table(title: str, rows: list[ExactCandidate], *, max_rows: int = 8) -> None:
    print(title)
    if not rows:
        print("  (no rows)")
        print()
        return
    for row in rows[:max_rows]:
        print(
            "  "
            f"{row.label}: class={row.classification} "
            f"grid_uncovered={row.grid_uncovered} "
            f"grid_longest={row.grid_longest} "
            f"forbidden_len={fmt_frac(row.forbidden_length)} "
            f"max_gap={fmt_frac(row.max_gap)} "
            f"gap/thresh={fmt_ratio(row.gap_ratio)} "
            f"boundary_witnesses={row.boundary_witness_count} "
            f"unprotected={row.unprotected_count} "
            f"first_unprotected={fmt_frac(row.first_unprotected)} "
            f"Q={row.endpoint_q}"
        )
        print(f"    speeds={row.speeds}")
    print()


def main() -> None:
    print("Lonely Runner fourteen-runner disproof hunt (codex-2026-05-31 S373)")
    print(f"Convention: k={K} moving speeds, threshold=1/{N}, grid={GRID}.")
    print("Exact audits use S356 interval union and S360 endpoint protection.")
    print()

    print(f"Building grid masks for speeds 1..{MAX_SPEED} ...")
    masks = build_grid_masks(MAX_SPEED)
    print("Grid masks ready.")
    print()

    baseline_row = grid_candidate("initial segment", tuple(range(1, N)), masks, with_longest=True)
    assert baseline_row is not None
    baseline = exact_candidate(baseline_row)
    print_exact_table("Baseline tight equality object", [baseline], max_rows=1)

    family_grid_rows: dict[str, list[GridCandidate]] = {
        "Single 14-gate replacements": gate_replacements(masks),
        "Multi-gate overloads": multi_gate_replacements(masks),
        "CRT residue lifts": crt_lift_trials(masks),
        "Greedy grid set covers": greedy_families(masks),
    }

    all_grid_rows: list[GridCandidate] = [baseline_row]
    for title, rows in family_grid_rows.items():
        all_grid_rows.extend(rows)
        audited = exact_audit(rows, limit=14, masks=masks)
        print_exact_table(title, audited)

    mutation_rows = mutation_search(masks, sorted(all_grid_rows, key=grid_rank)[:34])
    all_grid_rows.extend(mutation_rows)
    print_exact_table("Mandatory-gate mutation search", exact_audit(mutation_rows, limit=16, masks=masks))

    final_exact = exact_audit(all_grid_rows, limit=70, masks=masks)
    open_covers = [row for row in final_exact if row.classification == "open_cover"]
    boundary_rows = [row for row in final_exact if row.classification == "boundary_only"]
    positive_rows = [row for row in final_exact if row.classification == "positive_gap"]

    print_exact_table("Best exact audited candidates overall", final_exact, max_rows=12)

    print("Construction ledger")
    print(f"  open_cover_candidates={len(open_covers)}")
    print(f"  audited_boundary_only={len(boundary_rows)}")
    print(f"  audited_positive_gap={len(positive_rows)}")
    if positive_rows:
        best_positive = min(positive_rows, key=lambda row: (row.max_gap, row.grid_uncovered))
        print(
            "  best_positive_gap="
            f"{best_positive.label} gap/thresh={fmt_ratio(best_positive.gap_ratio)} "
            f"max_gap={fmt_frac(best_positive.max_gap)} speeds={best_positive.speeds}"
        )
    if boundary_rows:
        best_boundary = min(boundary_rows, key=lambda row: (row.unprotected_count, row.boundary_witness_count))
        print(
            "  best_boundary_only="
            f"{best_boundary.label} boundary_witnesses={best_boundary.boundary_witness_count} "
            f"unprotected={best_boundary.unprotected_count} speeds={best_boundary.speeds}"
        )
    print()

    print("Creative dead ends")
    print("  1. The initial segment already has full measure, but exposes the unit layer at a/14.")
    print("  2. A single 14-gate protects that unit layer and immediately opens positive gaps.")
    print("  3. Adding more gates worsens the leak; gate overload is not a viable cover architecture.")
    print("  4. S367 residue near-blocker lifts do not become interval counterexamples; repeated residues leak.")
    print("  5. Greedy and mutation pressure find dense covers, but exact audits still expose gaps or boundary witnesses.")
    print()

    print("New disproof idea worth keeping")
    print("  Instead of choosing speeds first, solve for a directed endpoint-protection cycle on the")
    print("  unit layer plus one higher quotient layer, then realize that cycle by speeds.  The")
    print("  speed-first searches repeatedly pay for unit protection by creating higher-denominator")
    print("  leaks, so the next counterexample attempt should prescribe the protection graph first.")


if __name__ == "__main__":
    main()
