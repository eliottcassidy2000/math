#!/usr/bin/env python3
"""
lonely_runner_feedback_loop_s364.py

codex-2026-05-31 S364

Feedback-loop attack on the Lonely Runner frontier:

* Route A: fourteen total runners, reduced k=13, n=14.
* Route B: fifteen total runners, reduced k=14, n=15.
* Route C: attempted near-counterexample constructions.

The script is deliberately exploratory.  It packages the session's repeated
dead-end rule into a reproducible ledger:

1. Attack the n=14 micro-staircase ansatz.
2. When the first obstruction is structural, push the idea to n=15.
3. When that also stalls, try to manufacture a counterexample-shaped speed set.
4. Return with the new obstruction family excised.

The central finite object is the cell arrangement of

    floor(n * {i * alpha}),   i=1,...,n-1.

For a residue vector v mod n and a candidate time cell s/n + alpha, a coordinate
is safe when (s*v_i + floor(n*{i*alpha})) mod n is neither 0 nor n-1.  A vector
that blocks every cell is an obstruction to a naive cellwise tight-lift lemma.
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

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class PatternSystem:
    n: int
    k: int
    patterns: tuple[tuple[int, ...], ...]
    masks: tuple[tuple[int, ...], ...]
    all_mask: int
    candidate_count: int


@dataclass(frozen=True)
class SearchResult:
    vector: tuple[int, ...]
    score: int
    candidate_count: int
    missed: int
    scalar_multiplier: int | None


@dataclass(frozen=True)
class ScalarFamily:
    full_multipliers: tuple[int, ...]
    gcd_histogram: tuple[tuple[int, int], ...]
    unit_example: SearchResult
    nonunit_example: SearchResult


@dataclass(frozen=True)
class GapTry:
    label: str
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    witness: Fraction | None
    boundary_witness_count: int


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def cell_pattern(n: int, k: int, alpha: Fraction) -> tuple[int, ...]:
    return tuple(int((n * ((i * alpha) % ONE)) // ONE) for i in range(1, k + 1))


def cell_patterns(n: int, k: int) -> tuple[tuple[int, ...], ...]:
    breaks = {Fraction(0), ONE}
    for i in range(1, k + 1):
        for a in range(n * i + 1):
            breaks.add(Fraction(a, n * i))
    ordered = sorted(breaks)
    patterns: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for lo, hi in zip(ordered, ordered[1:]):
        if lo == hi:
            continue
        pattern = cell_pattern(n, k, (lo + hi) / 2)
        if pattern not in seen:
            seen.add(pattern)
            patterns.append(pattern)
    return tuple(patterns)


def build_pattern_system(n: int) -> PatternSystem:
    k = n - 1
    patterns = cell_patterns(n, k)
    candidate_count = n * len(patterns)
    all_mask = (1 << candidate_count) - 1

    masks: list[list[int]] = [[0 for _ in range(n)] for _ in range(k)]
    bit = 1
    for s in range(n):
        for pattern in patterns:
            for i, bin_value in enumerate(pattern):
                for residue in range(n):
                    if (s * residue + bin_value) % n in (0, n - 1):
                        masks[i][residue] |= bit
            bit <<= 1
    return PatternSystem(
        n=n,
        k=k,
        patterns=patterns,
        masks=tuple(tuple(row) for row in masks),
        all_mask=all_mask,
        candidate_count=candidate_count,
    )


def scalar_multiplier(vector: tuple[int, ...], n: int) -> int | None:
    for m in range(n):
        if all(value == (m * i) % n for i, value in enumerate(vector, start=1)):
            return m
    return None


def blocked_mask(system: PatternSystem, vector: tuple[int, ...]) -> int:
    mask = 0
    for i, residue in enumerate(vector):
        mask |= system.masks[i][residue]
    return mask


def score_vector(system: PatternSystem, vector: tuple[int, ...]) -> SearchResult:
    mask = blocked_mask(system, vector)
    score = mask.bit_count()
    return SearchResult(
        vector=vector,
        score=score,
        candidate_count=system.candidate_count,
        missed=system.candidate_count - score,
        scalar_multiplier=scalar_multiplier(vector, system.n),
    )


def best_scalar_ramp(system: PatternSystem) -> SearchResult:
    return scalar_family(system).nonunit_example


def scalar_family(system: PatternSystem) -> ScalarFamily:
    results: dict[int, SearchResult] = {}
    for m in range(system.n):
        vector = tuple((m * i) % system.n for i in range(1, system.k + 1))
        results[m] = score_vector(system, vector)

    full = tuple(m for m, result in results.items() if result.missed == 0)
    gcd_counts: dict[int, int] = {}
    for m in full:
        g = gcd(m, system.n)
        gcd_counts[g] = gcd_counts.get(g, 0) + 1

    unit_m = next(m for m in full if gcd(m, system.n) == 1)
    nonunit_m = next(m for m in full if 1 < gcd(m, system.n) < system.n)
    return ScalarFamily(
        full_multipliers=full,
        gcd_histogram=tuple(sorted(gcd_counts.items())),
        unit_example=results[unit_m],
        nonunit_example=results[nonunit_m],
    )


def print_scalar_family(n: int, family: ScalarFamily) -> None:
    print(
        "  "
        f"n={n} full_scalar_multipliers={family.full_multipliers} "
        f"gcd_histogram={family.gcd_histogram}"
    )
    print_result(f"n{n}_scalar_unit_example", family.unit_example)
    print_result(f"n{n}_scalar_nonunit_example", family.nonunit_example)


def improve_vector(
    system: PatternSystem,
    vector: tuple[int, ...],
    *,
    exclude_scalar: bool,
    rng: random.Random,
    sweeps: int,
) -> SearchResult:
    current = tuple(vector)
    current_result = score_vector(system, current)
    if exclude_scalar and current_result.scalar_multiplier is not None:
        current_result = SearchResult(current, -1, system.candidate_count, system.candidate_count + 1, current_result.scalar_multiplier)

    for _ in range(sweeps):
        changed = False
        order = list(range(system.k))
        rng.shuffle(order)
        for idx in order:
            best_local = current_result
            best_value = current[idx]
            for residue in range(system.n):
                if residue == current[idx]:
                    continue
                candidate = list(current)
                candidate[idx] = residue
                candidate_tuple = tuple(candidate)
                result = score_vector(system, candidate_tuple)
                if exclude_scalar and result.scalar_multiplier is not None:
                    continue
                if result.score > best_local.score:
                    best_local = result
                    best_value = residue
            if best_local.score > current_result.score:
                tmp = list(current)
                tmp[idx] = best_value
                current = tuple(tmp)
                current_result = best_local
                changed = True
        if not changed:
            break
    return current_result


def local_search(
    system: PatternSystem,
    *,
    exclude_scalar: bool,
    seed_vectors: list[tuple[int, ...]],
    restarts: int,
    sweeps: int,
    seed: int,
) -> SearchResult:
    rng = random.Random(seed)
    starts = list(seed_vectors)
    for _ in range(restarts):
        starts.append(tuple(rng.randrange(system.n) for _ in range(system.k)))

    best: SearchResult | None = None
    for vector in starts:
        result = improve_vector(
            system,
            vector,
            exclude_scalar=exclude_scalar,
            rng=rng,
            sweeps=sweeps,
        )
        if best is None or result.score > best.score:
            best = result
    assert best is not None
    return best


def gate_status(speeds: tuple[int, ...], n: int) -> str:
    forced = [v for v in speeds if v % n == 0]
    return f"mult_{n}={len(forced)}"


def classify_speed_set(label: str, speeds: list[int]) -> GapTry:
    row = S356.report(label, speeds)
    if row.forbidden_length < ONE:
        classification = "positive_gap"
    elif row.boundary_witness_count:
        classification = "boundary_only"
    else:
        classification = "open_cover_candidate"
    return GapTry(
        label=label,
        speeds=tuple(row.speeds),
        classification=classification,
        forbidden_length=row.forbidden_length,
        max_gap=row.max_gap,
        gap_ratio=row.max_gap / row.threshold if row.threshold else Fraction(0),
        witness=row.witness,
        boundary_witness_count=row.boundary_witness_count,
    )


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def random_gated_near_counterexamples(
    *,
    k: int,
    n: int,
    max_speed: int,
    trials: int,
    seed: int,
) -> list[GapTry]:
    rng = random.Random(seed)
    multiples = [v for v in range(1, max_speed + 1) if v % n == 0]
    population = list(range(1, max_speed + 1))
    seen: set[tuple[int, ...]] = set()
    rows: list[GapTry] = []
    for trial in range(trials):
        forced = rng.choice(multiples)
        rest = [v for v in population if v != forced]
        combo = tuple(sorted([forced] + rng.sample(rest, k - 1)))
        if combo in seen or not primitive(combo):
            continue
        seen.add(combo)
        rows.append(classify_speed_set(f"random n={n} trial={trial}", list(combo)))
    rows.sort(
        key=lambda row: (
            0 if row.classification == "open_cover_candidate" else 1,
            0 if row.classification == "boundary_only" else 1,
            row.gap_ratio,
            row.speeds,
        )
    )
    return rows[:8]


def deterministic_gap_tries(k: int, n: int) -> list[GapTry]:
    families: list[tuple[str, list[int]]] = [
        (f"initial segment k={k}", list(range(1, k + 1))),
        (f"single {n}-gate", list(range(1, k)) + [n]),
        (f"double {n}-gate", list(range(1, k - 1)) + [n, 2 * n]),
        (f"triple {n}-gate", list(range(1, k - 2)) + [n, 2 * n, 3 * n]),
        (f"remove middle add gate", [v for v in range(1, k + 2) if v != (k + 1) // 2][: k - 1] + [n]),
    ]
    if n == 14:
        families.extend(
            [
                ("factor ladder 2x7", [2, 4, 6, 7, 8, 10, 12, 14, 21, 28, 35, 42, 49]),
                ("near doubled initials", [2 * i for i in range(1, 7)] + [1, 3, 5, 7, 9, 11, 14]),
            ]
        )
    if n == 15:
        families.extend(
            [
                ("factor ladder 3x5", [3, 5, 6, 9, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50]),
                ("near tripled initials", [3 * i for i in range(1, 8)] + [1, 2, 4, 7, 8, 11, 15]),
            ]
        )
    return [classify_speed_set(label, speeds) for label, speeds in families]


def print_result(prefix: str, result: SearchResult) -> None:
    scalar = (
        "non-scalar"
        if result.scalar_multiplier is None
        else f"scalar ramp m={result.scalar_multiplier}, gcd(m,n)={gcd(result.scalar_multiplier, len(result.vector) + 1)}"
    )
    print(f"  {prefix}_vector={result.vector}")
    print(
        "  "
        f"{prefix}_covered={result.score}/{result.candidate_count} "
        f"missed={result.missed} {scalar}"
    )


def print_gap_rows(title: str, rows: list[GapTry], n: int) -> None:
    print(title)
    for row in rows:
        print(
            "  "
            f"{row.label}: class={row.classification} {gate_status(row.speeds, n)} "
            f"forbidden_len={fmt_frac(row.forbidden_length)} "
            f"max_gap={fmt_frac(row.max_gap)} "
            f"gap/thresh={float(row.gap_ratio):.6f} "
            f"boundary_witnesses={row.boundary_witness_count} "
            f"speeds={row.speeds}"
        )
    print()


def main() -> None:
    print("Lonely Runner 14/15/counterexample feedback loop (codex-2026-05-31 S364)")
    print("External inspiration checked during session:")
    print("  Rosenfeld arXiv:2509.14111 (8 runners), Trakulthongchai arXiv:2511.22427 (9/10),")
    print("  Sungkawichai-Trakulthongchai arXiv:2604.23906 (11/12/13),")
    print("  and Jensen arXiv:2605.27941 (mixed thresholds/Fourier formulas).")
    print()

    systems = {n: build_pattern_system(n) for n in (14, 15)}
    for n, system in systems.items():
        print(f"Cell system n={n}: patterns={len(system.patterns)} candidates={system.candidate_count}")
    print()

    print("Cycle 1A: attack fourteen runners with the naive full-cell tight lift")
    scalar_family14 = scalar_family(systems[14])
    scalar14 = scalar_family14.nonunit_example
    print_scalar_family(14, scalar_family14)
    print("  dead_end=the whole affine scalar-ramp line blocks every cell.")
    print("  lesson=split off Dirichlet-equality ramps before testing a micro-staircase lemma.")
    print()

    print("Cycle 1B: push that obstruction shape to fifteen runners")
    scalar_family15 = scalar_family(systems[15])
    scalar15 = scalar_family15.nonunit_example
    print_scalar_family(15, scalar_family15)
    print("  new_idea=unit ramps are reindexed initial segments; nonunit ramps are the quotient/descent cases.")
    print()

    print("Cycle 1C: try to turn gates into counterexample-shaped covers")
    print_gap_rows("  Deterministic 14-runner construction attempts", deterministic_gap_tries(13, 14), 14)
    print_gap_rows("  Deterministic 15-runner construction attempts", deterministic_gap_tries(14, 15), 15)

    print("Cycle 2A: return to fourteen runners, excising scalar ramps")
    seeds14 = [
        scalar14.vector,
        (8, 2, 10, 4, 12, 13, 0, 8, 2, 10, 4, 12, 6),
        (7, 4, 9, 6, 7, 8, 5, 0, 1, 12, 13, 12, 7),
    ]
    nonscalar14 = local_search(
        systems[14],
        exclude_scalar=True,
        seed_vectors=seeds14,
        restarts=70,
        sweeps=30,
        seed=36414,
    )
    print_result("n14_nonscalar_best", nonscalar14)
    print("  dead_end=not a full blocker after scalar excision; missed cells become the next certificate targets.")
    print()

    print("Cycle 2B: perform the same scalar-excision test for fifteen runners")
    seeds15 = [
        scalar15.vector,
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10),
        (3, 6, 9, 12, 0, 3, 6, 9, 12, 5, 3, 6, 9, 12),
    ]
    nonscalar15 = local_search(
        systems[15],
        exclude_scalar=True,
        seed_vectors=seeds15,
        restarts=70,
        sweeps=30,
        seed=36415,
    )
    print_result("n15_nonscalar_best", nonscalar15)
    print("  dead_end=again no full non-scalar blocker in this search; the residue gaps look structured.")
    print()

    print("Cycle 2C: randomized near-counterexample attempts after the residue lesson")
    print_gap_rows(
        "  Best random 14-gated attempts",
        random_gated_near_counterexamples(k=13, n=14, max_speed=90, trials=90, seed=364140),
        14,
    )
    print_gap_rows(
        "  Best random 15-gated attempts",
        random_gated_near_counterexamples(k=14, n=15, max_speed=105, trials=90, seed=364150),
        15,
    )

    print("Synthesis")
    print("  1. The naive n=14 micro-staircase lemma is false unless the affine scalar-ramp family is split off.")
    print("  2. The same scalar line appears immediately at n=15; nonunit members are quotient/descent cases.")
    print("  3. After scalar excision, deterministic local search found only near-blockers, not full blockers.")
    print("  4. Gated speed-set construction attempts did not produce an open-cover candidate; all sampled hard cases")
    print("     retained either a boundary witness or a positive open gap.")
    print("  5. Next proof move: classify scalar ramps as quotient/descent cases, then certify the missed cells of")
    print("     the best non-scalar blockers by a mixed-threshold or endpoint-pressure argument.")


if __name__ == "__main__":
    main()
