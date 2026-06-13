#!/usr/bin/env python3
"""
lonely_runner_creative_multiroute_s372.py

codex-2026-05-31 S372

Creative multi-route pressure on the fourteen-runner frontier.  The session
tries several small, exact, and deliberately different attacks:

1. scalar-ramp moat atlas for many n;
2. exact one- and two-puncture searches around the scalar ramps for n=14,15;
3. anatomy of the witness cells opened by the best n=14 near-blocker;
4. greedy non-scalar blocker search from unrelated starts;
5. gated speed-set pressure with endpoint-peel summaries.

The main new clue is that the S364 "best non-scalar blockers" are not generic:
they are scalar ramps with one coordinate punctured by a divisor-layer jump.
For n=14 every one-puncture scalar deviation opens at least 56 witness cells,
and every two-puncture deviation opens at least 112 in the exact full-cell
system.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations, product
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
class PatternSystem:
    n: int
    k: int
    patterns: tuple[tuple[int, ...], ...]
    intervals: tuple[tuple[Fraction, Fraction], ...]
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
class PunctureSummary:
    n: int
    pattern_count: int
    radius: int
    checked: int
    best_missed: int
    best_count: int
    positions: tuple[int, ...]
    deltas: tuple[int, ...]
    examples: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class SpeedSummary:
    label: str
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witnesses: int
    unprotected: int
    peel_depth: int
    core_endpoints: int


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def cell_pattern(n: int, k: int, alpha: Fraction) -> tuple[int, ...]:
    return tuple(int((n * ((i * alpha) % ONE)) // ONE) for i in range(1, k + 1))


def build_pattern_system(n: int) -> PatternSystem:
    k = n - 1
    breaks = {Fraction(0), ONE}
    for i in range(1, k + 1):
        for a in range(n * i + 1):
            breaks.add(Fraction(a, n * i))

    patterns: list[tuple[int, ...]] = []
    intervals: list[tuple[Fraction, Fraction]] = []
    seen: set[tuple[int, ...]] = set()
    ordered = sorted(breaks)
    for lo, hi in zip(ordered, ordered[1:]):
        if lo == hi:
            continue
        pattern = cell_pattern(n, k, (lo + hi) / 2)
        if pattern not in seen:
            seen.add(pattern)
            patterns.append(pattern)
            intervals.append((lo, hi))

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
        patterns=tuple(patterns),
        intervals=tuple(intervals),
        masks=tuple(tuple(row) for row in masks),
        all_mask=all_mask,
        candidate_count=candidate_count,
    )


def scalar_vector(n: int, m: int) -> tuple[int, ...]:
    return tuple((m * i) % n for i in range(1, n))


def scalar_multiplier(vector: tuple[int, ...], n: int) -> int | None:
    for m in range(n):
        if vector == scalar_vector(n, m):
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


def scalar_audit(system: PatternSystem) -> tuple[bool, tuple[tuple[int, int], ...]]:
    full = []
    for m in range(system.n):
        result = score_vector(system, scalar_vector(system.n, m))
        if result.missed == 0:
            full.append(m)
    hist = Counter(gcd(m, system.n) for m in full)
    return len(full) == system.n, tuple(sorted(hist.items()))


def puncture_vectors(n: int, radius: int):
    """Return (vector, positions, deltas) for scalar-ramp punctures."""

    for m in range(n):
        base = list(scalar_vector(n, m))
        for positions in combinations(range(n - 1), radius):
            choices = []
            for pos in positions:
                choices.append([r for r in range(n) if r != base[pos]])
            for residues in product(*choices):
                vector = list(base)
                deltas = []
                for pos, residue in zip(positions, residues):
                    deltas.append((residue - base[pos]) % n)
                    vector[pos] = residue
                yield (tuple(vector), tuple(pos + 1 for pos in positions), tuple(deltas))


def puncture_summary(system: PatternSystem, radius: int) -> PunctureSummary:
    best_missed = system.candidate_count + 1
    best: list[tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]] = []
    checked = 0
    for vector, positions, deltas in puncture_vectors(system.n, radius):
        checked += 1
        missed = score_vector(system, vector).missed
        if missed < best_missed:
            best_missed = missed
            best = [(vector, positions, deltas)]
        elif missed == best_missed:
            best.append((vector, positions, deltas))

    return PunctureSummary(
        n=system.n,
        pattern_count=len(system.patterns),
        radius=radius,
        checked=checked,
        best_missed=best_missed,
        best_count=len(best),
        positions=tuple(sorted({pos for _v, positions, _d in best for pos in positions})),
        deltas=tuple(sorted({delta for _v, _p, deltas in best for delta in deltas})),
        examples=tuple(vector for vector, _positions, _deltas in best[:5]),
    )


def scalar_puncture_atlas(n_min: int, n_max: int) -> list[PunctureSummary]:
    out = []
    for n in range(n_min, n_max + 1):
        out.append(puncture_summary(build_pattern_system(n), radius=1))
    return out


def missed_cells(system: PatternSystem, vector: tuple[int, ...]) -> list[tuple[int, int]]:
    mask = blocked_mask(system, vector)
    out = []
    per_s = len(system.patterns)
    for candidate in range(system.candidate_count):
        if not ((mask >> candidate) & 1):
            out.append((candidate // per_s, candidate % per_s))
    return out


def compress_indices(values: list[int]) -> list[tuple[int, int]]:
    if not values:
        return []
    values = sorted(values)
    out = []
    lo = hi = values[0]
    for value in values[1:]:
        if value == hi + 1:
            hi = value
        else:
            out.append((lo, hi))
            lo = hi = value
    out.append((lo, hi))
    return out


def print_missed_cell_anatomy(
    label: str, system: PatternSystem, vector: tuple[int, ...]
) -> None:
    cells = missed_cells(system, vector)
    by_s: dict[int, list[int]] = defaultdict(list)
    for s, pattern_index in cells:
        by_s[s].append(pattern_index)

    unique_patterns = sorted({pattern_index for _s, pattern_index in cells})
    total_alpha_width = sum(
        system.intervals[idx][1] - system.intervals[idx][0] for idx in unique_patterns
    )

    print(label)
    print(f"  vector={vector}")
    print(f"  missed_cells={len(cells)}")
    print(f"  s_hist={tuple(sorted(Counter(s for s, _idx in cells).items()))}")
    print(f"  unique_patterns={len(unique_patterns)} total_alpha_width={fmt_frac(total_alpha_width)}")
    for s in sorted(by_s):
        print(f"  s={s:2d} pattern_ranges={compress_indices(by_s[s])}")
    print("  alpha_cells=")
    for idx in unique_patterns[:16]:
        lo, hi = system.intervals[idx]
        print(f"    p={idx:4d} [{fmt_frac(lo)}, {fmt_frac(hi)}) width={fmt_frac(hi - lo)}")
    print()


def improve_vector(
    system: PatternSystem,
    vector: tuple[int, ...],
    rng: random.Random,
    sweeps: int,
) -> SearchResult:
    current = vector
    current_result = score_vector(system, current)
    if current_result.scalar_multiplier is not None:
        current_result = SearchResult(current, -1, system.candidate_count, system.candidate_count + 1, current_result.scalar_multiplier)

    for _ in range(sweeps):
        changed = False
        order = list(range(system.k))
        rng.shuffle(order)
        for idx in order:
            best_local = current_result
            best_residue = current[idx]
            for residue in range(system.n):
                if residue == current[idx]:
                    continue
                candidate = list(current)
                candidate[idx] = residue
                candidate_tuple = tuple(candidate)
                result = score_vector(system, candidate_tuple)
                if result.scalar_multiplier is not None:
                    continue
                if result.score > best_local.score:
                    best_local = result
                    best_residue = residue
            if best_local.score > current_result.score:
                current = tuple(
                    best_residue if i == idx else value for i, value in enumerate(current)
                )
                current_result = best_local
                changed = True
        if not changed:
            break
    return current_result


def greedy_search(system: PatternSystem, restarts: int, seed: int) -> SearchResult:
    rng = random.Random(seed)
    starts = [
        scalar_vector(system.n, 0),
        scalar_vector(system.n, 1),
        scalar_vector(system.n, system.n // 2),
    ]
    for _ in range(restarts):
        starts.append(tuple(rng.randrange(system.n) for _ in range(system.k)))

    best: SearchResult | None = None
    for vector in starts:
        result = improve_vector(system, vector, rng, sweeps=18)
        if best is None or result.score > best.score:
            best = result
    assert best is not None
    return best


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def classify_speed_set(label: str, speeds: tuple[int, ...]) -> SpeedSummary:
    row = S356.report(label, list(speeds))
    if row.forbidden_length < ONE:
        classification = "positive_gap"
    elif row.boundary_witness_count:
        classification = "boundary_only"
    else:
        classification = "open_cover_candidate"
    descent = S362.summarize(list(row.speeds))
    threshold = Fraction(1, len(row.speeds) + 1)
    return SpeedSummary(
        label=label,
        speeds=tuple(row.speeds),
        classification=classification,
        forbidden_length=row.forbidden_length,
        max_gap=row.max_gap,
        gap_ratio=row.max_gap / threshold,
        boundary_witnesses=row.boundary_witness_count,
        unprotected=descent.unprotected_count,
        peel_depth=len(descent.peel_layers),
        core_endpoints=descent.core_endpoint_count,
    )


def random_gate_pressure(k: int, n: int, max_speed: int, trials: int, seed: int) -> list[SpeedSummary]:
    rng = random.Random(seed)
    population = list(range(1, max_speed + 1))
    gates = [v for v in population if v % n == 0]
    raw = []
    seen: set[tuple[int, ...]] = set()

    for trial in range(trials):
        forced = rng.choice(gates)
        rest = [v for v in population if v != forced]
        speeds = tuple(sorted([forced] + rng.sample(rest, k - 1)))
        if speeds in seen or not primitive(speeds):
            continue
        seen.add(speeds)
        row = S356.report(f"random gate {n} #{trial}", list(speeds))
        if row.forbidden_length < ONE:
            rank = 2
        elif row.boundary_witness_count:
            rank = 1
        else:
            rank = 0
        threshold = Fraction(1, k + 1)
        raw.append((rank, row.max_gap / threshold, tuple(row.speeds), f"random gate {n} #{trial}"))
        raw.sort(key=lambda item: (item[0], item[1], item[2]))
        raw = raw[:5]

    return [classify_speed_set(label, speeds) for _rank, _ratio, speeds, label in raw]


def print_speed_summary(title: str, rows: list[SpeedSummary]) -> None:
    print(title)
    for row in rows:
        gate = len(row.speeds) + 1
        gate_count = sum(1 for speed in row.speeds if speed % gate == 0)
        print(
            "  "
            f"{row.label}: class={row.classification} mult_{gate}={gate_count} "
            f"gap/thresh={float(row.gap_ratio):.6f} max_gap={fmt_frac(row.max_gap)} "
            f"boundary={row.boundary_witnesses} unprotected={row.unprotected} "
            f"peel={row.peel_depth} core_E={row.core_endpoints} speeds={row.speeds}"
        )
    print()


def main() -> None:
    print("Lonely Runner creative multi-route sprint (codex-2026-05-31 S372)")
    print("Reduced frontier: k=13 moving speeds, n=14, threshold=1/14.\n")

    systems = {n: build_pattern_system(n) for n in (14, 15)}
    for n, system in systems.items():
        all_scalar, hist = scalar_audit(system)
        print(
            f"Cell system n={n}: patterns={len(system.patterns)} "
            f"candidates={system.candidate_count} all_scalar_full={all_scalar} "
            f"gcd_hist={hist}"
        )
    print()

    print("Route 1: scalar-puncture moat atlas")
    print("  n  patterns  best_h1_missed  best_count  positions  deltas")
    for summary in scalar_puncture_atlas(4, 22):
        print(
            "  "
            f"{summary.n:2d} {summary.pattern_count:8d} "
            f"{summary.best_missed:15d} {summary.best_count:11d} "
            f"{summary.positions!s:16s} {summary.deltas}"
        )
    print()

    print("Route 2: exact radius-1 and radius-2 punctures around scalar ramps")
    for n in (14, 15):
        for radius in (1, 2):
            summary = puncture_summary(systems[n], radius)
            print(
                "  "
                f"n={n} radius={radius} checked={summary.checked} "
                f"best_missed={summary.best_missed} best_count={summary.best_count} "
                f"positions={summary.positions} deltas={summary.deltas}"
            )
            for vector in summary.examples[:3]:
                print(f"    example={vector}")
    print()

    print("Route 3: greedy non-scalar blocker searches")
    for n in (14, 15):
        result = greedy_search(systems[n], restarts=90, seed=36700 + n)
        distance_to_scalar = min(
            sum(a != b for a, b in zip(result.vector, scalar_vector(n, m))) for m in range(n)
        )
        print(
            "  "
            f"n={n} best={result.score}/{result.candidate_count} "
            f"missed={result.missed} distance_to_scalar={distance_to_scalar} "
            f"vector={result.vector}"
        )
    print()

    print("Route 4: anatomy of the n=14 one-puncture witness moat")
    best14 = (0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0)
    transported14 = (8, 2, 10, 4, 12, 13, 0, 8, 2, 10, 4, 12, 6)
    print_missed_cell_anatomy("  zero-ramp puncture at i=6 by +7", systems[14], best14)
    print_missed_cell_anatomy("  S364 transported best near-blocker", systems[14], transported14)

    print("Route 5: speed-set gate pressure")
    hand14 = [
        ("initial segment", tuple(range(1, 14))),
        ("single 14 gate", tuple(list(range(1, 13)) + [14])),
        ("drop six add 14", tuple([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14])),
        ("factor 2x7 ladder", (2, 4, 6, 7, 8, 10, 12, 14, 21, 28, 35, 42, 49)),
    ]
    hand15 = [
        ("initial segment", tuple(range(1, 15))),
        ("single 15 gate", tuple(list(range(1, 14)) + [15])),
        ("drop six add 15", tuple([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15])),
        ("factor 3x5 ladder", (3, 5, 6, 9, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50)),
    ]
    print_speed_summary("  Hand-built n=14 families", [classify_speed_set(label, speeds) for label, speeds in hand14])
    print_speed_summary("  Hand-built n=15 families", [classify_speed_set(label, speeds) for label, speeds in hand15])
    print_speed_summary(
        "  Random n=14 gate pressure",
        random_gate_pressure(k=13, n=14, max_speed=84, trials=45, seed=36714),
    )
    print_speed_summary(
        "  Random n=15 gate pressure",
        random_gate_pressure(k=14, n=15, max_speed=90, trials=45, seed=36715),
    )

    print("Synthesis")
    print("  1. The old near-blockers are scalar ramps with a single divisor jump, not generic obstructions.")
    print("  2. In the exact n=14 full-cell system, radius-1 scalar punctures leave at least 56 witnesses;")
    print("     radius-2 punctures leave at least 112.  The scalar line has a visible moat.")
    print("  3. The n=14 moat is organized as seven odd s-layers times eight alpha-cells, with")
    print("     total alpha width 7/858 for the zero-ramp representative.")
    print("  4. Random and hand-built gated speed sets again leaked by positive gaps or boundary witnesses.")
    print("  5. Next proof target: prove the scalar-puncture moat symbolically, then combine it with")
    print("     endpoint descent so any putative blocker must be far from every scalar ramp and still core-free.")


if __name__ == "__main__":
    main()
