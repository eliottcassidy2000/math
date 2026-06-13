#!/usr/bin/env python3
"""
lonely_runner_composite_gate_feedback_s364.py

codex-2026-05-31 S364

Three-track exploratory loop requested by the user:

  A. attack the fourteen-runner frontier: reduced k=13, n=14;
  B. when that leaks, jump to a creative fifteen-runner idea: k=14, n=15;
  C. when that leaks, search for possible disproof constructions, then cycle.

The script is intentionally a "research notebook with teeth": cheap exact
interval scoring ranks many candidates, then the endpoint/core machinery is
run only on the most dangerous survivors.
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
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class ScoredSet:
    speeds: tuple[int, ...]
    n: int
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    classification: str


@dataclass(frozen=True)
class DeepSet:
    score: ScoredSet
    endpoints: int
    unprotected: int
    unit_skeleton: bool
    peel_depth: int
    core_endpoints: int
    first_layer_modulus: int


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def score_set(speeds: tuple[int, ...]) -> ScoredSet:
    row = S356.report("feedback-loop", list(speeds))
    if row.forbidden_length < ONE:
        classification = "positive_gap"
    elif row.boundary_witness_count:
        classification = "boundary_only"
    else:
        classification = "open_cover"
    return ScoredSet(
        speeds=tuple(row.speeds),
        n=len(row.speeds) + 1,
        forbidden_length=row.forbidden_length,
        max_gap=row.max_gap,
        gap_ratio=row.max_gap / row.threshold if row.threshold else Fraction(0),
        classification=classification,
    )


def deepen(score: ScoredSet) -> DeepSet:
    summary = S362.summarize(list(score.speeds))
    first_layer_modulus = 1
    if summary.peel_layers:
        first_layer_modulus = summary.peel_layers[0].removed_subgroup_modulus
    return DeepSet(
        score=score,
        endpoints=summary.endpoint_count,
        unprotected=summary.unprotected_count,
        unit_skeleton=summary.unit_skeleton,
        peel_depth=len(summary.peel_layers),
        core_endpoints=summary.core_endpoint_count,
        first_layer_modulus=first_layer_modulus,
    )


def gate_count(speeds: tuple[int, ...], n: int) -> int:
    return sum(1 for v in speeds if v % n == 0)


def gate_profile_string(speeds: tuple[int, ...], n: int) -> str:
    factors = [p for p in range(2, n + 1) if n % p == 0 and all(p % q for q in range(2, p))]
    pieces = [f"{n}:{gate_count(speeds, n)}"]
    for p in factors:
        pieces.append(f"{p}:{sum(1 for v in speeds if v % p == 0)}")
    units = sum(1 for v in speeds if gcd(v, n) == 1)
    pieces.append(f"units:{units}")
    return " ".join(pieces)


def print_deep(label: str, deep: DeepSet) -> None:
    score = deep.score
    print(f"[{label}]")
    print(f"  speeds={score.speeds}")
    print(f"  gates {gate_profile_string(score.speeds, score.n)}")
    print(
        "  "
        f"class={score.classification} "
        f"forbidden={fmt_frac(score.forbidden_length)} "
        f"max_gap={fmt_frac(score.max_gap)} "
        f"gap/thresh={float(score.gap_ratio):.6f}"
    )
    print(
        "  "
        f"endpoints={deep.endpoints} unprotected={deep.unprotected} "
        f"unit={deep.unit_skeleton} peel={deep.peel_depth} "
        f"core_E={deep.core_endpoints} first_mod={deep.first_layer_modulus}"
    )
    print()


def best_by_gap(candidates: list[tuple[int, ...]], keep: int = 6) -> list[ScoredSet]:
    best: list[ScoredSet] = []
    for combo in candidates:
        if not primitive(combo):
            continue
        score = score_set(combo)
        best.append(score)
        best.sort(key=lambda s: (s.classification != "open_cover", s.gap_ratio, s.speeds))
        best = best[:keep]
    return best


def gated_box(k: int, max_speed: int) -> list[ScoredSet]:
    n = k + 1
    candidates = [
        combo
        for combo in combinations(range(1, max_speed + 1), k)
        if any(v % n == 0 for v in combo)
    ]
    return best_by_gap(candidates)


def mutation_walk(
    seed: tuple[int, ...],
    max_speed: int,
    steps: int,
    rng: random.Random,
) -> list[ScoredSet]:
    k = len(seed)
    n = k + 1
    current = tuple(sorted(seed))
    current_score = score_set(current)
    best = [current_score]
    for temperature_step in range(steps):
        drop = rng.choice(current)
        add = rng.randrange(1, max_speed + 1)
        if add in current:
            continue
        proposal = tuple(sorted((set(current) - {drop}) | {add}))
        if len(proposal) != k or not primitive(proposal):
            continue
        if not any(v % n == 0 for v in proposal):
            continue
        proposal_score = score_set(proposal)
        better = proposal_score.gap_ratio <= current_score.gap_ratio
        lucky = rng.random() < max(0.01, 0.10 * (1 - temperature_step / steps))
        if better or lucky:
            current = proposal
            current_score = proposal_score
        best.append(proposal_score)
        best.sort(key=lambda s: (s.classification != "open_cover", s.gap_ratio, s.speeds))
        best = best[:8]
    return best


def disproof_pressure(k: int, max_speed: int, trials: int, rng: random.Random) -> list[ScoredSet]:
    """Try to maximize coverage by sampling dense gated sets.

    A true disproof candidate would be an open_cover.  The ranking keeps small
    max gaps, because open covers have gap_ratio=0 and no boundary witnesses.
    """

    n = k + 1
    multiples = [v for v in range(1, max_speed + 1) if v % n == 0]
    best: list[ScoredSet] = []
    seen: set[tuple[int, ...]] = set()
    for _ in range(trials):
        forced = rng.choice(multiples)
        pool = [v for v in range(1, max_speed + 1) if v != forced]
        combo = tuple(sorted([forced] + rng.sample(pool, k - 1)))
        if combo in seen or not primitive(combo):
            continue
        seen.add(combo)
        score = score_set(combo)
        best.append(score)
        best.sort(
            key=lambda s: (
                s.classification != "open_cover",
                s.classification != "boundary_only",
                s.gap_ratio,
                -gate_count(s.speeds, n),
                s.speeds,
            )
        )
        best = best[:8]
    return best


def cycle_fourteen(rng: random.Random) -> list[DeepSet]:
    print("Cycle A: fourteen-runner attack (k=13, n=14)")
    seeds = [
        tuple(range(1, 14)),
        (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14),
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14),
    ]
    scores = [score_set(seed) for seed in seeds]
    scores.extend(gated_box(13, 16))
    scores.extend(mutation_walk(seeds[1], 56, 90, rng))
    scores.sort(key=lambda s: (s.classification != "open_cover", s.gap_ratio, s.speeds))
    deeps = [deepen(score) for score in scores[:6]]
    for i, deep in enumerate(deeps, 1):
        print_deep(f"14.{i}", deep)
    print("Dead-end note: every 14-gated survivor still has positive gap or empty core.\n")
    return deeps


def cycle_fifteen(rng: random.Random) -> list[DeepSet]:
    print("Cycle B: fifteen-runner idea (k=14, n=15)")
    seeds = [
        tuple(range(1, 15)),
        tuple(range(1, 14)) + (15,),
        (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15),
        (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 30),
    ]
    scores = [score_set(tuple(sorted(seed))) for seed in seeds]
    scores.extend(gated_box(14, 17))
    scores.extend(mutation_walk(tuple(range(1, 14)) + (15,), 60, 90, rng))
    scores.sort(key=lambda s: (s.classification != "open_cover", s.gap_ratio, s.speeds))
    deeps = [deepen(score) for score in scores[:6]]
    for i, deep in enumerate(deeps, 1):
        print_deep(f"15.{i}", deep)
    print("Dead-end note: the 15=3*5 gate also leaks; no all-protected core found.\n")
    return deeps


def cycle_disproof(rng: random.Random) -> list[DeepSet]:
    print("Cycle C: disproof-construction pressure")
    scores: list[ScoredSet] = []
    scores.extend(disproof_pressure(13, 90, 120, rng))
    scores.extend(disproof_pressure(14, 90, 120, rng))
    scores.sort(key=lambda s: (s.classification != "open_cover", s.gap_ratio, s.speeds))
    deeps = [deepen(score) for score in scores[:8]]
    for i, deep in enumerate(deeps, 1):
        print_deep(f"X.{i}", deep)
    print("Dead-end note: random high-coverage gated searches did not produce an open cover.\n")
    return deeps


def main() -> None:
    rng = random.Random(364)
    print("Lonely Runner 14/15/disproof feedback loop (codex-2026-05-31 S364)")
    print("Exact interval scoring; endpoint/core deepening on dangerous survivors.\n")

    all_deeps: list[DeepSet] = []
    for lap in range(1, 3):
        print(f"=== Feedback lap {lap} ===")
        all_deeps.extend(cycle_fourteen(rng))
        all_deeps.extend(cycle_fifteen(rng))
        all_deeps.extend(cycle_disproof(rng))

    best = sorted(all_deeps, key=lambda d: (d.score.gap_ratio, -d.peel_depth, d.score.speeds))[:10]
    print("Cross-cycle hardest survivors")
    for i, deep in enumerate(best, 1):
        print_deep(f"hard.{i}", deep)

    print("Synthesis")
    print("  14-runner route: the n=14 gate is necessary but seems to create")
    print("  leaks; the hardest local example has a long peel rather than a core.")
    print("  15-runner route: n=15=3*5 suggests the same CRT-gate descent, now")
    print("  with two prime channels; exploratory gated sets also leak.")
    print("  Disproof route: sampling tries to maximize coverage under the same")
    print("  gates, but no open cover or nonempty protection core appears.")
    print("  New common target: prove composite-gate protection cores project to")
    print("  smaller Bohr-boundary cores along the n-divisible speed channel.")


if __name__ == "__main__":
    main()
