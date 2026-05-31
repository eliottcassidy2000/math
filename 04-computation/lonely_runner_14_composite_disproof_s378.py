#!/usr/bin/env python3
"""
lonely_runner_14_composite_disproof_s378.py

codex-2026-05-31 S378

Targeted counterexample hunt for the 14-runner Lonely Runner frontier.

The S376 recursive atlas flagged n=14 as a composite anomaly:

* threshold denominator 14 has largest proper divisor 7;
* the scalar-puncture moat is unusually thin in density, with the best
  puncture at coordinate 6 and delta 7;
* the largest-proper-divisor quotient ladder has a tiny exact gap ratio
  0.005411, but many unprotected endpoints.

This script tries to turn that anomaly into a disproof.  It starts from the
seven-ladder near-miss and searches for endpoint repairs using speeds tuned
to the 2*7 structure: the missing 42 gate, other multiples of 7 and 14, and
top direct protectors of the exposed endpoint layer.

Success would be an exact S360 classification "open_cover".  The logged
near-misses are still useful if they explain why the composite anomaly does
not currently yield a counterexample.
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
S373 = SourceFileLoader(
    "lonely_runner_14_disproof_hunt_s373",
    str(ROOT / "04-computation" / "lonely_runner_14_disproof_hunt_s373.py"),
).load_module()


N = 14
K = 13
MAX_SPEED = 252
SEVEN_LADDER = (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)
MISSING_HALF_GATE = 42
BEAM_DEPTH = 1
BEAM_WIDTH = 40
EXACT_AUDITS_PER_STEP = 12
ANCHOR_TOP_PROTECTORS = 4
RUN_TOP_PROTECTOR_REPLACEMENTS = False
RUN_ANCHORED_COMBINATIONS = False


@dataclass(frozen=True)
class ProtectorRow:
    speed: int
    covers: int
    introduced: int
    residue: int
    quotient: int | None
    score: Fraction


@dataclass(frozen=True)
class AuditRow:
    label: str
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witnesses: int
    unprotected: int
    first_unprotected: Fraction | None
    peel_depth: int
    core_endpoints: int


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_ratio(x: Fraction) -> str:
    return f"{float(x):.6f}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def unprotected_points(speeds: tuple[int, ...]) -> list[Fraction]:
    points: dict[Fraction, list] = {}
    for endpoint in S360.endpoints(speeds):
        points.setdefault(endpoint.value, []).append(endpoint)

    out: list[Fraction] = []
    for value in sorted(points):
        indegree = sum(1 for v in speeds if S360.direct_protects(speeds, v, value))
        if indegree == 0:
            out.append(value)
    return out


def introduced_endpoint_values(speed: int) -> set[Fraction]:
    return {
        S360.circle_point(Fraction(N * m + sign, N * speed))
        for m in range(speed)
        for sign in (-1, 1)
    }


def protector_rows(seed: tuple[int, ...], max_speed: int) -> list[ProtectorRow]:
    leaks = unprotected_points(seed)
    rows: list[ProtectorRow] = []
    for speed in range(1, max_speed + 1):
        if speed in seed:
            continue
        covers = sum(1 for t in leaks if S360.direct_protects(seed, speed, t))
        if covers == 0:
            continue
        introduced = len(introduced_endpoint_values(speed))
        # Favor leak coverage, then lightly reward 7/14 compatibility and
        # penalize the number of new endpoint values.
        divisor_bonus = Fraction(8 if speed % 14 == 0 else 4 if speed % 7 == 0 else 0)
        score = Fraction(covers, 1) + divisor_bonus - Fraction(introduced, 28)
        rows.append(
            ProtectorRow(
                speed=speed,
                covers=covers,
                introduced=introduced,
                residue=speed % N,
                quotient=(speed // 7 if speed % 7 == 0 else None),
                score=score,
            )
        )
    return sorted(rows, key=lambda r: (-r.covers, -r.score, r.introduced, r.speed))


def audit(label: str, speeds: tuple[int, ...]) -> AuditRow:
    report = S356.report(label, list(speeds))
    summary = S360.summarize(list(speeds))
    return AuditRow(
        label=label,
        speeds=summary.speeds,
        classification=summary.classification,
        forbidden_length=summary.forbidden_length,
        max_gap=summary.max_gap,
        gap_ratio=summary.max_gap / summary.threshold,
        boundary_witnesses=report.boundary_witness_count,
        unprotected=summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
        peel_depth=-1,
        core_endpoints=-1,
    )


def audit_rank(row: AuditRow) -> tuple:
    class_rank = {"open_cover": 0, "boundary_only": 1, "positive_gap": 2}.get(
        row.classification, 3
    )
    return (
        class_rank,
        row.max_gap,
        row.unprotected,
        row.boundary_witnesses,
        row.speeds,
    )


def print_audits(title: str, rows: list[AuditRow], limit: int = 10) -> None:
    print(title)
    for row in sorted(rows, key=audit_rank)[:limit]:
        print(
            "  "
            f"{row.label}: class={row.classification} "
            f"len={fmt_frac(row.forbidden_length)} "
            f"gap={fmt_frac(row.max_gap)} "
            f"gap/th={fmt_ratio(row.gap_ratio)} "
            f"boundary={row.boundary_witnesses} "
            f"unprotected={row.unprotected} "
            f"first={fmt_frac(row.first_unprotected)} "
            f"peel={'-' if row.peel_depth < 0 else row.peel_depth} "
            f"coreE={'-' if row.core_endpoints < 0 else row.core_endpoints} "
            f"speeds={row.speeds}"
        )
    print()


def exact_half_gate_replacements() -> list[AuditRow]:
    rows: list[AuditRow] = []
    for removed in SEVEN_LADDER:
        speeds = tuple(sorted((set(SEVEN_LADDER) - {removed}) | {MISSING_HALF_GATE}))
        if len(speeds) != K or not primitive(speeds):
            continue
        rows.append(audit(f"insert42 remove={removed}", speeds))
    return rows


def exact_top_protector_replacements(
    protectors: list[ProtectorRow],
    top_k: int = 28,
) -> list[AuditRow]:
    rows: list[AuditRow] = []
    for protector in protectors[:top_k]:
        if protector.speed in SEVEN_LADDER:
            continue
        for removed in SEVEN_LADDER:
            speeds = tuple(sorted((set(SEVEN_LADDER) - {removed}) | {protector.speed}))
            if len(speeds) != K or not primitive(speeds):
                continue
            rows.append(audit(f"top-protector {protector.speed} remove={removed}", speeds))
    unique: dict[tuple[int, ...], AuditRow] = {}
    for row in rows:
        unique.setdefault(row.speeds, row)
    return sorted(unique.values(), key=audit_rank)


def candidate_pool(protectors: list[ProtectorRow]) -> tuple[int, ...]:
    pool = set(SEVEN_LADDER)
    pool.add(MISSING_HALF_GATE)
    pool.update(range(1, N))
    pool.update(7 * q for q in range(1, 37) if 7 * q <= MAX_SPEED)
    pool.update(N * q for q in range(1, 19) if N * q <= MAX_SPEED)
    pool.update(row.speed for row in protectors[:36])
    pool.update(row.speed for row in protectors if row.speed % 7 == 0 and row.covers >= 12)
    return tuple(sorted(v for v in pool if 1 <= v <= MAX_SPEED))


def grid_rank(row) -> tuple:
    longest = row.longest_uncovered_run if row.longest_uncovered_run >= 0 else 0
    return (row.uncovered, longest, row.speeds)


def beam_search(pool: tuple[int, ...], depth: int = 4, width: int = 180):
    masks = S373.build_grid_masks(max(pool))
    seed = S373.grid_candidate("seven-ladder", SEVEN_LADDER, masks, with_longest=True)
    assert seed is not None
    beam = [seed]
    exact_rows: list[AuditRow] = [audit("seven-ladder", SEVEN_LADDER)]
    seen_exact = {SEVEN_LADDER}

    for step in range(1, depth + 1):
        next_rows = []
        seen_grid: set[tuple[int, ...]] = set()
        for row in beam:
            current = set(row.speeds)
            for old in row.speeds:
                base = current - {old}
                for new in pool:
                    if new in base:
                        continue
                    speeds = tuple(sorted(base | {new}))
                    if speeds in seen_grid or len(speeds) != K:
                        continue
                    seen_grid.add(speeds)
                    label = f"beam{step} replace={old}->{new}"
                    grid = S373.grid_candidate(label, speeds, masks)
                    if grid is not None:
                        next_rows.append(grid)
        next_rows.sort(key=grid_rank)
        beam = []
        used: set[tuple[int, ...]] = set()
        for row in next_rows:
            if row.speeds in used:
                continue
            used.add(row.speeds)
            beam.append(row)
            if len(beam) >= width:
                break

        for row in beam[:EXACT_AUDITS_PER_STEP]:
            if row.speeds in seen_exact:
                continue
            seen_exact.add(row.speeds)
            exact_rows.append(audit(row.label, row.speeds))

        best = min(exact_rows, key=audit_rank)
        print(
            f"  beam_step={step} grid_best_uncovered={beam[0].uncovered} "
            f"exact_best={best.classification} gap/th={fmt_ratio(best.gap_ratio)} "
            f"unprotected={best.unprotected}"
        )

    return sorted(exact_rows, key=audit_rank), beam


def anchored_combinations(pool: tuple[int, ...], protectors: list[ProtectorRow]) -> list[AuditRow]:
    """Exact audit small hand-built combinations around 42 and the top protectors."""

    rows: list[AuditRow] = []
    top = [row.speed for row in protectors[:ANCHOR_TOP_PROTECTORS]]
    must_include_sets = [
        (42,),
        (42, 98),
        (42, 84),
        (42, 98, 196),
        (42, 70, 98),
    ]

    for required in must_include_sets:
        required_set = set(required)
        remove_count = len(required_set - set(SEVEN_LADDER))
        for removed in combinations(SEVEN_LADDER, remove_count):
            speeds = tuple(sorted((set(SEVEN_LADDER) - set(removed)) | required_set))
            if len(speeds) == K and primitive(speeds):
                rows.append(audit(f"forced={required} removed={removed}", speeds))

    for added in combinations(top, 2):
        add_set = set(added) | {42}
        remove_count = len(add_set - set(SEVEN_LADDER))
        if remove_count > 4:
            continue
        for removed in combinations(SEVEN_LADDER, remove_count):
            speeds = tuple(sorted((set(SEVEN_LADDER) - set(removed)) | add_set))
            if len(speeds) == K and primitive(speeds):
                rows.append(audit(f"top-protectors add={tuple(sorted(add_set))}", speeds))

    # Remove exact duplicates while preserving the best label seen first.
    unique: dict[tuple[int, ...], AuditRow] = {}
    for row in rows:
        unique.setdefault(row.speeds, row)
    return sorted(unique.values(), key=audit_rank)


def hand_anchored_audits() -> list[AuditRow]:
    """Small exact multi-protector probes distilled from the 42/7 layer."""

    rows: list[AuditRow] = []
    specs = [
        ((42, 98), (35, 49)),
        ((42, 98), (70, 84)),
        ((42, 196), (70, 84)),
        ((42, 119), (35, 70)),
        ((42, 119), (70, 84)),
        ((42, 119, 196), (21, 35, 84)),
        ((42, 119, 196), (21, 70, 84)),
        ((42, 98, 196), (21, 70, 84)),
    ]
    for add, remove in specs:
        speeds = tuple(sorted((set(SEVEN_LADDER) - set(remove)) | set(add)))
        if len(speeds) == K and primitive(speeds):
            rows.append(audit(f"hand add={add} remove={remove}", speeds))
    return sorted(rows, key=audit_rank)


def leak_anatomy(speeds: tuple[int, ...]) -> None:
    leaks = unprotected_points(speeds)
    denom_hist = Counter(t.denominator for t in leaks)
    residue_mod_14 = Counter(int((t * N) % N) if (t * N).denominator == 1 else -1 for t in leaks)
    residue_mod_98 = Counter(int((t * 98) % 98) if (t * 98).denominator == 1 else -1 for t in leaks)
    print("Seven-ladder exposed endpoint anatomy")
    print(f"  leaks={len(leaks)}")
    print(f"  denominator_hist={tuple(sorted(denom_hist.items()))}")
    print(f"  residue_mod_14_hist={tuple(sorted(residue_mod_14.items()))}")
    print(f"  residue_mod_98_hist_prefix={tuple(sorted(residue_mod_98.items())[:14])}")
    print(f"  first_leaks={[fmt_frac(t) for t in leaks[:12]]}")
    print()


def main() -> None:
    print("14-runner composite-anomaly disproof attempt (codex-2026-05-31 S378)")
    print("Target: k=13 speeds, threshold=1/14, exact open-cover classification.")
    print()

    base = audit("seven-ladder", SEVEN_LADDER)
    print_audits("Baseline half-divisor quotient ladder", [base], limit=1)
    leak_anatomy(SEVEN_LADDER)

    protectors = protector_rows(SEVEN_LADDER, MAX_SPEED)
    print("Top direct protectors of the seven-ladder leak layer")
    for row in protectors[:18]:
        print(
            "  "
            f"speed={row.speed:3d} covers={row.covers:2d} "
            f"introduced={row.introduced:3d} residue14={row.residue:2d} "
            f"q7={row.quotient if row.quotient is not None else '-':>2} "
            f"score={float(row.score):7.3f}"
        )
    print()

    half_gate_repairs = exact_half_gate_replacements()
    print_audits(
        "Exact replacement of one ladder speed by the missing half-gate 42",
        half_gate_repairs,
        limit=13,
    )

    pool = candidate_pool(protectors)
    print(
        f"Composite repair pool: size={len(pool)} max={max(pool)} "
        f"multiples_of_7={sum(1 for v in pool if v % 7 == 0)} "
        f"multiples_of_14={sum(1 for v in pool if v % 14 == 0)}"
    )
    print()

    if RUN_TOP_PROTECTOR_REPLACEMENTS:
        single_repairs = exact_top_protector_replacements(protectors, top_k=12)
        print_audits(
            "Exact one-speed replacement by top composite protectors",
            single_repairs,
            limit=14,
        )
    else:
        single_repairs = []
        print("Exact one-speed replacement by top composite protectors")
        print("  skipped in the stored bounded run; enable RUN_TOP_PROTECTOR_REPLACEMENTS for the larger pass.")
        print()

    anchored = hand_anchored_audits()
    print_audits("Hand-anchored exact audits forcing 42 plus high-cover protectors", anchored, limit=14)

    best_overall = sorted(half_gate_repairs + single_repairs + anchored + [base], key=audit_rank)[0]
    print("Final disproof status")
    if best_overall.classification == "open_cover":
        print("  FOUND_OPEN_COVER=YES")
    else:
        print("  FOUND_OPEN_COVER=NO")
    print(
        "  best="
        f"{best_overall.label} class={best_overall.classification} "
        f"gap/th={fmt_ratio(best_overall.gap_ratio)} "
        f"unprotected={best_overall.unprotected} "
        f"first={fmt_frac(best_overall.first_unprotected)} "
        f"speeds={best_overall.speeds}"
    )
    print()
    print("Interpretation")
    print(
        "  The 14=2*7 anomaly really does create a near-disproof channel: the "
        "largest-proper-divisor ladder is extremely close in measure."
    )
    print(
        "  But the same composite repair layer is endpoint-expensive.  The missing "
        "42 gate protects many old leaks, yet every replacement keeps a positive "
        "gap or creates a large new boundary layer."
    )
    print(
        "  In this search, the anomaly behaves more like evidence for an endpoint "
        "closure obstruction than like a counterexample recipe."
    )


if __name__ == "__main__":
    main()
