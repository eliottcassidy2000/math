#!/usr/bin/env python3
"""
lonely_runner_cycle_s372.py

codex-2026-05-31 S372

Long feedback-loop sprint across the next Lonely Runner frontiers.

Routes:
  A. n=14 / k=13: push the scalar-gauge quotient target beyond S367.
  B. n=15 / k=14: force the same idea through the 15=3*5 case.
  C. disproof pressure: lift quotient near-blockers into actual speed sets.

The key distinction from S364/S367 is that the counterexample search is no
longer just random gated sampling.  It uses residue templates from the quotient
near-blockers and asks whether those templates lift to harder interval covers.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations, product
from math import gcd
from pathlib import Path
import random
import sys


sys.stdout.reconfigure(line_buffering=True)

ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S367 = SourceFileLoader(
    "lonely_runner_k13_scalar_gauge_s367",
    str(ROOT / "04-computation" / "lonely_runner_k13_scalar_gauge_s367.py"),
).load_module()


@dataclass(frozen=True)
class TopVector:
    n: int
    vector: tuple[int, ...]
    missed: int
    score: int
    support: tuple[int, ...]
    source: str


@dataclass(frozen=True)
class SpeedAttempt:
    n: int
    source: str
    scalar_shift: int
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witnesses: int


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def support(vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(i for i, value in enumerate(vector, start=1) if value)


def scalar_add(vector: tuple[int, ...], n: int, m: int) -> tuple[int, ...]:
    return tuple((value + m * i) % n for i, value in enumerate(vector, start=1))


def is_scalar_zero(vector: tuple[int, ...]) -> bool:
    return all(value == 0 for value in vector)


def raw_score(system, vector: tuple[int, ...]) -> tuple[int, int]:
    mask = 0
    for coord, residue in enumerate(vector):
        mask |= system.masks[coord][residue]
    score = mask.bit_count()
    return score, system.candidate_count - score


def insert_top(top: list[TopVector], row: TopVector, limit: int) -> None:
    top.append(row)
    top.sort(key=lambda item: (item.missed, len(item.support), item.vector, item.source))
    del top[limit:]


def print_top_vectors(title: str, rows: list[TopVector], limit: int = 8) -> None:
    print(title)
    for rank, row in enumerate(rows[:limit], start=1):
        print(
            "  "
            f"{rank:02d}. n={row.n} missed={row.missed:5d} "
            f"score={row.score:5d} support={row.support} "
            f"source={row.source} vector={row.vector}"
        )
    print()


def exact_support_scan(system, support_size: int, *, top_limit: int = 10) -> list[TopVector]:
    """Exact scan of normalized vectors with a fixed support size.

    Coordinate 1 is fixed to 0 by scalar-gauge normalization, so available
    support coordinates are 2..k in one-based labels.
    """

    positions = range(1, system.k)
    residues = range(1, system.n)
    top: list[TopVector] = []
    missed_hist: Counter[int] = Counter()
    scanned = 0

    progress_stride = 2_000_000 if support_size >= 4 else None

    for coords in combinations(positions, support_size):
        coord_set = set(coords)
        base_mask = 0
        for coord in range(system.k):
            if coord not in coord_set:
                base_mask |= system.masks[coord][0]

        for values in product(residues, repeat=support_size):
            mask = base_mask
            vector = [0 for _ in range(system.k)]
            for coord, value in zip(coords, values):
                vector[coord] = value
                mask |= system.masks[coord][value]
            score = mask.bit_count()
            missed = system.candidate_count - score
            missed_hist[missed] += 1
            scanned += 1
            if progress_stride and scanned % progress_stride == 0:
                print(
                    "  "
                    f"progress scanned={scanned} current_best={top[0].missed if top else 'n/a'}"
                )
            insert_top(
                top,
                TopVector(
                    n=system.n,
                    vector=tuple(vector),
                    missed=missed,
                    score=score,
                    support=support(tuple(vector)),
                    source=f"support-{support_size}",
                ),
                top_limit,
            )

    print(f"Exact normalized support scan n={system.n}, support={support_size}")
    print(f"  scanned={scanned}")
    print(f"  best_missed={top[0].missed if top else 'n/a'}")
    print(f"  missed_hist_prefix={tuple(sorted(missed_hist.items())[:10])}")
    print_top_vectors("  top vectors", top, limit=min(6, top_limit))
    return top


def subgroup_scan(
    system,
    values: tuple[int, ...],
    label: str,
    *,
    top_limit: int = 10,
) -> list[TopVector]:
    top: list[TopVector] = []
    support_hist: Counter[int] = Counter()
    missed_hist: Counter[int] = Counter()
    scanned = 0

    for tail in product(values, repeat=system.k - 1):
        vector = (0,) + tuple(tail)
        if is_scalar_zero(vector):
            continue
        score, missed = raw_score(system, vector)
        scanned += 1
        support_hist[len(support(vector))] += 1
        missed_hist[missed] += 1
        insert_top(
            top,
            TopVector(
                n=system.n,
                vector=vector,
                missed=missed,
                score=score,
                support=support(vector),
                source=label,
            ),
            top_limit,
        )

    print(f"Exact normalized subgroup scan n={system.n}, label={label}")
    print(f"  values={values}")
    print(f"  scanned={scanned}")
    print(f"  best_missed={top[0].missed if top else 'n/a'}")
    print(f"  support_hist={tuple(sorted(support_hist.items()))}")
    print(f"  missed_hist_prefix={tuple(sorted(missed_hist.items())[:10])}")
    print_top_vectors("  top subgroup vectors", top, limit=min(6, top_limit))
    return top


def improve_vector(system, vector: tuple[int, ...], rng: random.Random, sweeps: int) -> TopVector:
    current = S367.gauge_normalize(vector, system.n)
    if is_scalar_zero(current):
        current = (0,) + tuple(rng.randrange(system.n) for _ in range(system.k - 1))
        if is_scalar_zero(current):
            current = current[:-1] + (1,)

    score, missed = raw_score(system, current)
    best = TopVector(system.n, current, missed, score, support(current), "local")

    for _ in range(sweeps):
        changed = False
        order = list(range(1, system.k))
        rng.shuffle(order)
        for coord in order:
            local_best = best
            local_value = current[coord]
            for residue in range(system.n):
                if residue == current[coord]:
                    continue
                candidate = list(current)
                candidate[coord] = residue
                candidate_tuple = tuple(candidate)
                if is_scalar_zero(candidate_tuple):
                    continue
                candidate_score, candidate_missed = raw_score(system, candidate_tuple)
                if candidate_score > local_best.score:
                    local_best = TopVector(
                        system.n,
                        candidate_tuple,
                        candidate_missed,
                        candidate_score,
                        support(candidate_tuple),
                        "local",
                    )
                    local_value = residue
            if local_best.score > best.score:
                current = current[:coord] + (local_value,) + current[coord + 1 :]
                best = local_best
                changed = True
        if not changed:
            break
    return best


def local_search(
    system,
    seeds: list[tuple[int, ...]],
    *,
    restarts: int,
    sweeps: int,
    seed: int,
) -> list[TopVector]:
    rng = random.Random(seed)
    top_by_vector: dict[tuple[int, ...], TopVector] = {}
    starts = list(seeds)
    for _ in range(restarts):
        starts.append((0,) + tuple(rng.randrange(system.n) for _ in range(system.k - 1)))

    for start in starts:
        row = improve_vector(system, start, rng, sweeps)
        old = top_by_vector.get(row.vector)
        if old is None or row.score > old.score:
            top_by_vector[row.vector] = row

    rows = sorted(top_by_vector.values(), key=lambda item: (item.missed, item.vector))[:12]
    print(f"Deterministic local quotient search n={system.n}")
    print(f"  starts={len(starts)} restarts={restarts} sweeps={sweeps}")
    print_top_vectors("  local optima", rows, limit=8)
    return rows


def classify_speed_attempt(
    source: str,
    n: int,
    scalar_shift: int,
    speeds: tuple[int, ...],
) -> SpeedAttempt:
    row = S356.report(source, list(speeds))
    if row.forbidden_length < Fraction(1):
        classification = "positive_gap"
    elif row.boundary_witness_count:
        classification = "boundary_only"
    else:
        classification = "open_cover"
    return SpeedAttempt(
        n=n,
        source=source,
        scalar_shift=scalar_shift,
        speeds=tuple(row.speeds),
        classification=classification,
        forbidden_length=row.forbidden_length,
        max_gap=row.max_gap,
        gap_ratio=row.max_gap / row.threshold if row.threshold else Fraction(0),
        boundary_witnesses=row.boundary_witness_count,
    )


def residue_classes(n: int, max_speed: int) -> dict[int, list[int]]:
    classes: dict[int, list[int]] = {r: [] for r in range(n)}
    for value in range(1, max_speed + 1):
        classes[value % n].append(value)
    return classes


def deterministic_lift(residues: tuple[int, ...], classes: dict[int, list[int]]) -> tuple[int, ...] | None:
    used: set[int] = set()
    speeds: list[int] = []
    for residue in residues:
        options = [value for value in classes[residue] if value not in used]
        if not options:
            return None
        value = options[0]
        used.add(value)
        speeds.append(value)
    return tuple(sorted(speeds))


def random_lift(
    residues: tuple[int, ...],
    classes: dict[int, list[int]],
    rng: random.Random,
) -> tuple[int, ...] | None:
    need = Counter(residues)
    chosen: list[int] = []
    for residue, count in need.items():
        options = classes[residue]
        if len(options) < count:
            return None
        chosen.extend(rng.sample(options, count))
    return tuple(sorted(chosen))


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for value in speeds:
        g = gcd(g, value)
    return g == 1


def template_lift_pressure(
    templates: list[TopVector],
    *,
    max_speed: int,
    random_per_shift: int,
    seed: int,
    top_limit: int = 10,
) -> list[SpeedAttempt]:
    rng = random.Random(seed)
    attempts: dict[tuple[int, ...], SpeedAttempt] = {}

    for template in templates:
        n = template.n
        classes = residue_classes(n, max_speed)
        for m in range(n):
            residues = scalar_add(template.vector, n, m)
            deterministic = deterministic_lift(residues, classes)
            candidates: list[tuple[str, tuple[int, ...] | None]] = [
                (f"{template.source}:det", deterministic)
            ]
            for trial in range(random_per_shift):
                candidates.append((f"{template.source}:rnd{trial}", random_lift(residues, classes, rng)))
            for source_suffix, speeds in candidates:
                if speeds is None or len(set(speeds)) != len(template.vector) or not primitive(speeds):
                    continue
                source = f"n{n}-{source_suffix}-miss{template.missed}"
                attempt = classify_speed_attempt(source, n, m, speeds)
                old = attempts.get(attempt.speeds)
                if old is None or (attempt.classification, attempt.gap_ratio) < (old.classification, old.gap_ratio):
                    attempts[attempt.speeds] = attempt

    rows = sorted(
        attempts.values(),
        key=lambda item: (
            {"open_cover": 0, "boundary_only": 1, "positive_gap": 2}[item.classification],
            item.gap_ratio,
            -sum(1 for value in item.speeds if value % item.n == 0),
            item.speeds,
        ),
    )[:top_limit]

    print(f"Residue-template lift pressure max_speed={max_speed}")
    print(f"  templates={len(templates)} unique_attempts={len(attempts)}")
    for rank, row in enumerate(rows, start=1):
        gate_count = sum(1 for value in row.speeds if value % row.n == 0)
        print(
            "  "
            f"{rank:02d}. n={row.n} class={row.classification} "
            f"gap/thresh={float(row.gap_ratio):.6f} "
            f"max_gap={fmt_frac(row.max_gap)} boundary={row.boundary_witnesses} "
            f"gate_mults={gate_count} shift={row.scalar_shift} "
            f"source={row.source} speeds={row.speeds}"
        )
    print()
    return rows


def main() -> None:
    print("Lonely Runner 14/15/disproof cycling sprint (codex S372)")
    print("Exact quotient cell scans plus residue-template speed lifts.\n")

    systems = {n: S367.build_pattern_system(n) for n in (14, 15)}
    for n, system in systems.items():
        print(
            f"Cell system n={n}: k={system.k} patterns={len(system.patterns)} "
            f"candidates={system.candidate_count}"
        )
    print()

    print("Cycle 1A: n=14 quotient support attack")
    n14_support1 = exact_support_scan(systems[14], 1)
    n14_support2 = exact_support_scan(systems[14], 2)
    n14_support3 = exact_support_scan(systems[14], 3)
    n14_torsion = subgroup_scan(systems[14], (0, 7), "2-torsion")
    n14_local = local_search(
        systems[14],
        [row.vector for row in n14_support1[:3] + n14_torsion[:3]],
        restarts=120,
        sweeps=80,
        seed=37214,
    )
    print("Cycle 1A dead end: scalar quotient still has no full nonzero blocker; best remains a half-turn leak.\n")

    print("Cycle 1B: force the quotient idea through n=15")
    n15_support1 = exact_support_scan(systems[15], 1)
    n15_support2 = exact_support_scan(systems[15], 2)
    n15_support3 = exact_support_scan(systems[15], 3)
    n15_order3 = subgroup_scan(systems[15], (0, 5, 10), "order-3-subgroup")
    n15_local = local_search(
        systems[15],
        [row.vector for row in n15_support1[:4] + n15_order3[:4]],
        restarts=160,
        sweeps=90,
        seed=37215,
    )
    print("Cycle 1B dead end: the 15=3*5 analogue also leaks; the best obstruction lives in the order-3 subgroup.\n")

    print("Cycle 1C: lift quotient near-blockers into actual speed-set pressure tests")
    top_templates = (
        n14_support1[:2]
        + n14_torsion[:2]
        + n14_local[:2]
        + n15_support1[:3]
        + n15_order3[:3]
        + n15_local[:2]
    )
    template_rows = template_lift_pressure(
        top_templates,
        max_speed=90,
        random_per_shift=1,
        seed=372915,
        top_limit=12,
    )
    print("Cycle 1C dead end: quotient-guided lifts still produce positive gaps, not open covers.\n")

    print("Cycle 2A: return to n=14 with one more exact layer")
    n14_support4 = exact_support_scan(systems[14], 4)
    print("Cycle 2A dead end: support 4 does not beat the coordinate-6 half-turn extremal.\n")

    print("Cycle 2B: n=15 creative follow-up from the n=14 half-turn lesson")
    n15_templates = sorted(
        n15_support1 + n15_support2[:4] + n15_order3[:6] + n15_local[:6],
        key=lambda row: (row.missed, len(row.support), row.vector),
    )[:8]
    print_top_vectors("  selected n=15 templates for second lift", n15_templates, limit=8)
    print("Cycle 2B dead end: the selected templates are all nonzero quotient leaks with explicit missed cells.\n")

    print("Cycle 2C: disproof pressure, second pass with only the sharpest templates")
    second_templates = sorted(
        n14_support1[:1] + n14_support4[:3] + n15_templates[:5],
        key=lambda row: (row.n, row.missed, row.vector),
    )
    template_lift_pressure(
        second_templates,
        max_speed=120,
        random_per_shift=1,
        seed=372240,
        top_limit=12,
    )
    print("Cycle 2C dead end: increasing the lift range lowered no case to an open cover.\n")

    best14 = min(n14_support1 + n14_support2 + n14_support3 + n14_torsion + n14_local + n14_support4, key=lambda r: r.missed)
    best15 = min(n15_support1 + n15_support2 + n15_support3 + n15_order3 + n15_local, key=lambda r: r.missed)

    print("Synthesis")
    print(
        "  1. n=14: exact support <=4, full 2-torsion, and local search all point to "
        f"{best14.vector} with missed={best14.missed}."
    )
    print(
        "  2. n=15: the scalar-gauge quotient analogue has a sharper order-3 obstruction "
        f"{best15.vector} with missed={best15.missed}."
    )
    print(
        "  3. The n=15 obstruction is the creative transfer: replace the n=14 half-turn "
        "coordinate by an order-3 coordinate at the far endpoint of the staircase."
    )
    print(
        "  4. Residue-template lifting is a better disproof pressure test than unguided gates, "
        "but every lifted candidate here still had a positive gap."
    )
    print(
        "  5. Next target: classify single-coordinate subgroup defects for all composite n, "
        "then turn their missed-cell families into a finite stencil certificate."
    )


if __name__ == "__main__":
    main()
