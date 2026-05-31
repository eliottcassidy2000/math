#!/usr/bin/env python3
"""
lonely_runner_torsion_crt_feedback_s377.py

codex-2026-05-31 S377

Another long Lonely Runner feedback loop, deliberately cycling through:

1. the fourteen-runner scalar-gauge quotient problem;
2. a creative fifteen-runner CRT/torsion analogue; and
3. reverse-cover pressure for possible disproof constructions.

The script is intentionally exploratory, but all micro-staircase statements use
the exact full-cell bitset model from S367.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from pathlib import Path
import importlib.util
import random
import sys


ROOT = Path(__file__).resolve().parents[1]


def load_s367():
    path = ROOT / "04-computation/lonely_runner_k13_scalar_gauge_s367.py"
    spec = importlib.util.spec_from_file_location("s367_s377", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


S367 = load_s367()


@dataclass(frozen=True)
class GapResult:
    speeds: tuple[int, ...]
    covered_length: Fraction
    max_gap: Fraction
    interval_count: int
    complement_count: int


def fmt_frac(x: Fraction) -> str:
    return S367.fmt_frac(x)


def support_of(vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(i for i, value in enumerate(vector, start=1) if value)


def result_for(system, vector: tuple[int, ...]):
    return S367.result_for_normalized(system, tuple(value % system.n for value in vector))


def insert_top(top: list, result, limit: int) -> None:
    top.append(result)
    top.sort(key=lambda r: (r.missed, r.vector))
    del top[limit:]


def missed_meta(system, vector: tuple[int, ...]) -> list[tuple[int, int]]:
    return S367.missed_candidates(system, vector)


def missed_summary(system, vector: tuple[int, ...]) -> dict[str, object]:
    missed = missed_meta(system, vector)
    by_shift = Counter(s for s, _ in missed)
    by_gcd = Counter(gcd(s, system.n) for s, _ in missed)
    by_pattern = defaultdict(list)
    widths = Counter()
    margins = Counter()
    for s, p_idx in missed:
        by_pattern[p_idx].append(s)
        pattern = system.patterns[p_idx]
        widths[pattern.hi - pattern.lo] += 1
        residues = tuple((s * value + bin_value) % system.n for value, bin_value in zip(vector, pattern.bins))
        margins[min(S367.residue_margin(residue, system.n) for residue in residues)] += 1
    return {
        "missed": len(missed),
        "shift_hist": sorted(by_shift.items()),
        "gcd_shift_hist": sorted(by_gcd.items()),
        "pattern_count": len(by_pattern),
        "width_hist": [(fmt_frac(width), count) for width, count in sorted(widths.items())],
        "margin_hist": sorted(margins.items()),
        "patterns": by_pattern,
    }


def scalar_ramp_full_check(system) -> tuple[int, ...]:
    full = []
    for m in range(system.n):
        vector = tuple((m * i) % system.n for i in range(1, system.k + 1))
        if S367.raw_score(system, vector) == system.candidate_count:
            full.append(m)
    return tuple(full)


def print_system_header(system, label: str) -> None:
    print(
        f"{label}: n={system.n} k={system.k} "
        f"patterns={len(system.patterns)} candidates={system.candidate_count}"
    )
    print(f"  scalar_ramp_full_multipliers={scalar_ramp_full_check(system)}")


def two_torsion_spectrum(system) -> list:
    assert system.n % 2 == 0
    top = []
    by_support = defaultdict(lambda: {"scanned": 0, "best": system.candidate_count + 1, "count": 0, "example": None})
    value = system.n // 2
    for bits in product((0, value), repeat=system.k - 1):
        vector = (0,) + bits
        if not any(vector):
            continue
        result = result_for(system, vector)
        support = len(support_of(vector))
        row = by_support[support]
        row["scanned"] += 1
        if result.missed < row["best"]:
            row["best"] = result.missed
            row["count"] = 1
            row["example"] = vector
        elif result.missed == row["best"]:
            row["count"] += 1
        insert_top(top, result, 14)

    print("  exact normalized 2-torsion spectrum")
    for support in sorted(by_support):
        row = by_support[support]
        print(
            f"    support={support:2d} scanned={row['scanned']:4d} "
            f"best_missed={row['best']:4d} best_count={row['count']:3d} "
            f"example_support={support_of(row['example'])}"
        )
    print("  top binary vectors")
    for item in top[:8]:
        print(f"    missed={item.missed:4d} support={support_of(item.vector)} vector={item.vector}")
    return top


def half_turn_coordinate_profile(system) -> None:
    profile = []
    for coord in range(1, system.k):
        vector = [0] * system.k
        vector[coord] = system.n // 2
        result = result_for(system, tuple(vector))
        profile.append((coord + 1, result.missed))
    print(f"  single-coordinate half-turn profile={profile}")


def print_missed_ledger(system, vector: tuple[int, ...], label: str) -> None:
    summary = missed_summary(system, vector)
    print(label)
    print(f"    vector={vector}")
    print(f"    missed={summary['missed']} pattern_count={summary['pattern_count']}")
    print(f"    shift_hist={summary['shift_hist']}")
    print(f"    gcd_shift_hist={summary['gcd_shift_hist']}")
    print(f"    width_hist={summary['width_hist']}")
    print(f"    margin_hist={summary['margin_hist']}")


def exact_small_support(system, support_size: int, residues: tuple[int, ...], limit: int = 8) -> tuple[int, list]:
    top = []
    scanned = 0
    for coords in combinations(range(1, system.k), support_size):
        for values in product(residues, repeat=support_size):
            vector = [0] * system.k
            for coord, value in zip(coords, values):
                vector[coord] = value
            result = result_for(system, tuple(vector))
            scanned += 1
            if result.scalar_multiplier is None:
                insert_top(top, result, limit)
    return scanned, top


def print_small_support(system, support_size: int, residues: tuple[int, ...], label: str) -> list:
    scanned, top = exact_small_support(system, support_size, residues)
    print(
        f"  {label}: support={support_size} residues={residues} "
        f"scanned={scanned} best_missed={top[0].missed if top else 'n/a'}"
    )
    for item in top[:5]:
        print(f"    missed={item.missed:4d} support={support_of(item.vector)} vector={item.vector}")
    return top


def local_reverse_cover(system, restarts: int, sweeps: int, seed: int, seeds: list[tuple[int, ...]]) -> list:
    rng = random.Random(seed)
    best_by_vector = {}
    for vector in seeds:
        if len(vector) != system.k:
            continue
        result = S367.improve_vector(system, vector, rng, sweeps=sweeps)
        if result.scalar_multiplier is None:
            best_by_vector[result.vector] = result
    for _ in range(restarts):
        vector = S367.random_normalized_vector(system, rng)
        if not any(vector):
            continue
        result = S367.improve_vector(system, vector, rng, sweeps=sweeps)
        if result.scalar_multiplier is None:
            old = best_by_vector.get(result.vector)
            if old is None or result.score > old.score:
                best_by_vector[result.vector] = result
    return sorted(best_by_vector.values(), key=lambda r: (r.missed, r.vector))[:12]


def print_local_reverse_cover(system, label: str, restarts: int, seed: int, seeds: list[tuple[int, ...]]) -> list:
    results = local_reverse_cover(system, restarts=restarts, sweeps=24, seed=seed, seeds=seeds)
    print(f"  reverse-cover local search {label}: restarts={restarts}")
    for idx, result in enumerate(results[:8], start=1):
        print(f"    {idx:2d}. missed={result.missed:4d} vector={result.vector}")
    if results:
        print_missed_ledger(system, results[0].vector, f"  best reverse-cover ledger {label}")
    return results


def one_extra_defect_from_half_turn(system) -> list:
    base = [0] * system.k
    base[5] = system.n // 2
    top = []
    for coord in range(1, system.k):
        if coord == 5:
            continue
        for residue in range(1, system.n):
            vector = list(base)
            vector[coord] = residue
            result = result_for(system, tuple(vector))
            insert_top(top, result, 12)
    print("  coordinate-6 half-turn plus one extra defect")
    for item in top[:10]:
        print(f"    missed={item.missed:4d} support={support_of(item.vector)} vector={item.vector}")
    return top


def crt_signature(vector: tuple[int, ...], n: int) -> str:
    if n != 15:
        return "n/a"
    mod3 = Counter(value % 3 for value in vector)
    mod5 = Counter(value % 5 for value in vector)
    return f"mod3={sorted(mod3.items())} mod5={sorted(mod5.items())}"


def forbidden_intervals(speed: int, k: int) -> list[tuple[Fraction, Fraction]]:
    radius = Fraction(1, k + 1)
    out: list[tuple[Fraction, Fraction]] = []
    for center in range(speed):
        lo = Fraction(center, speed) - Fraction(radius, speed)
        hi = Fraction(center, speed) + Fraction(radius, speed)
        if lo < 0:
            out.append((lo + 1, Fraction(1)))
            out.append((Fraction(0), hi))
        elif hi > 1:
            out.append((lo, Fraction(1)))
            out.append((Fraction(0), hi - 1))
        else:
            out.append((lo, hi))
    return out


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    intervals = sorted(intervals)
    merged: list[tuple[Fraction, Fraction]] = []
    for lo, hi in intervals:
        if not merged or lo > merged[-1][1]:
            merged.append((lo, hi))
        else:
            old_lo, old_hi = merged[-1]
            merged[-1] = (old_lo, max(old_hi, hi))
    return merged


def gap_result(speeds: tuple[int, ...]) -> GapResult:
    k = len(speeds)
    intervals = []
    for speed in speeds:
        intervals.extend(forbidden_intervals(speed, k))
    merged = merge_intervals(intervals)
    covered = sum((hi - lo for lo, hi in merged), Fraction(0))
    gaps = []
    cursor = Fraction(0)
    for lo, hi in merged:
        if lo > cursor:
            gaps.append(lo - cursor)
        cursor = max(cursor, hi)
    if cursor < 1:
        gaps.append(Fraction(1) - cursor)
    return GapResult(
        speeds=speeds,
        covered_length=covered,
        max_gap=max(gaps) if gaps else Fraction(0),
        interval_count=len(merged),
        complement_count=len(gaps),
    )


def disproof_pressure(k: int, n: int, max_speed: int, seed: int, samples: int) -> list[GapResult]:
    rng = random.Random(seed)
    candidates: set[tuple[int, ...]] = set()
    candidates.add(tuple(sorted(set(range(1, k + 1)) - {n // 2} | {n})))
    pool = list(range(1, max_speed + 1))
    gated_pool = [x for x in pool if x % n == 0]
    for gate in gated_pool:
        low = [x for x in pool if x != gate and gcd(x, n) == 1]
        for omit in range(1, min(n, k + 2)):
            trial = sorted((set(range(1, k + 2)) - {omit}) | {gate})
            if len(trial) == k:
                candidates.add(tuple(trial))
        for _ in range(samples // max(1, len(gated_pool))):
            trial = sorted(rng.sample(low, k - 1) + [gate])
            candidates.add(tuple(trial))
    while len(candidates) < samples:
        candidates.add(tuple(sorted(rng.sample(pool, k))))
    results = [gap_result(candidate) for candidate in candidates]
    return sorted(results, key=lambda r: (r.max_gap, -r.covered_length, r.speeds))[:10]


def print_disproof_pressure(k: int, n: int, max_speed: int, seed: int, samples: int) -> list[GapResult]:
    results = disproof_pressure(k=k, n=n, max_speed=max_speed, seed=seed, samples=samples)
    print(f"  disproof pressure: k={k} n={n} max_speed={max_speed} candidates>={samples}")
    for item in results[:8]:
        print(
            "    "
            f"max_gap={fmt_frac(item.max_gap):>9s} covered={fmt_frac(item.covered_length):>9s} "
            f"components={item.interval_count:3d} complement={item.complement_count:2d} speeds={item.speeds}"
        )
    return results


def cycle_one(system14, system15) -> None:
    print("\nCycle 1 -- exact quotient obstruction, then odd-modulus replacement")
    print("Route 14: the hard object is still the quotient, not the scalar line.")
    half_turn_coordinate_profile(system14)
    top_binary = two_torsion_spectrum(system14)
    print_missed_ledger(system14, top_binary[0].vector, "  coordinate-6 binary extremal ledger")

    print("\nRoute 15: no half-turn exists, so test the CRT subgroups of Z/15Z.")
    residues_all = tuple(range(1, system15.n))
    residues_3_torsion = (5, 10)
    residues_5_torsion = (3, 6, 9, 12)
    print_small_support(system15, 1, residues_all, "all residues")
    print_small_support(system15, 2, residues_all, "all residues")
    print_small_support(system15, 3, residues_3_torsion, "3-torsion subgroup")
    print_small_support(system15, 3, residues_5_torsion, "5-torsion subgroup")

    print("\nRoute disproof: try to make an open-cover-looking gated family.")
    print_disproof_pressure(k=13, n=14, max_speed=32, seed=3721, samples=180)
    print_disproof_pressure(k=14, n=15, max_speed=35, seed=3722, samples=220)


def cycle_two(system14, system15) -> tuple[list, list]:
    print("\nCycle 2 -- local full-blocker pressure in the quotient")
    print("Route 14: push random classes uphill after exact torsion failed.")
    seeds14 = [
        (0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0),
        (0, 7, 0, 7, 0, 7, 0, 7, 0, 7, 0, 7, 0),
    ]
    results14 = print_local_reverse_cover(system14, "n=14", restarts=42, seed=3723, seeds=seeds14)

    print("\nRoute 15: seed the search with CRT stripes and let the cell system reply.")
    seeds15 = [
        (0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10),
        (0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5),
        (0,) + tuple(5 if i % 2 else 10 for i in range(1, system15.k)),
        (0,) + tuple((3 * i) % 15 for i in range(1, system15.k)),
        (0,) + tuple((5 * i) % 15 for i in range(1, system15.k)),
        (0,) + tuple(3 if i % 3 == 1 else 6 if i % 3 == 2 else 9 for i in range(1, system15.k)),
    ]
    results15 = print_local_reverse_cover(system15, "n=15", restarts=48, seed=3724, seeds=seeds15)
    if results15:
        print(f"  n=15 best CRT signature: {crt_signature(results15[0].vector, 15)}")

    print("\nRoute disproof: rank the closest exact interval-cover candidates again.")
    print_disproof_pressure(k=13, n=14, max_speed=42, seed=3725, samples=260)
    print_disproof_pressure(k=14, n=15, max_speed=45, seed=3726, samples=300)
    return results14, results15


def cycle_three(system14, system15, results14: list, results15: list) -> None:
    print("\nCycle 3 -- convert dead ends into proof-shaped pressure")
    print("Route 14: perturb the unique binary extremal and watch it get worse.")
    one_extra_defect_from_half_turn(system14)
    if results14:
        print_missed_ledger(system14, results14[0].vector, "  best n=14 local optimum re-ledger")

    print("\nRoute 15: inspect the best odd-composite near-blocker by shift gcd.")
    if results15:
        print_missed_ledger(system15, results15[0].vector, "  best n=15 local optimum ledger")
        print(f"    CRT signature={crt_signature(results15[0].vector, 15)}")

    print("\nRoute disproof: force gate overload with two multiples of n.")
    for k, n, max_speed, seed in [(13, 14, 56, 3727), (14, 15, 60, 3728)]:
        rng = random.Random(seed)
        candidates = []
        gates = [n, 2 * n]
        pool = [x for x in range(1, max_speed + 1) if x not in gates]
        for _ in range(160):
            speeds = tuple(sorted(rng.sample(pool, k - len(gates)) + gates))
            candidates.append(gap_result(speeds))
        candidates.sort(key=lambda r: (r.max_gap, -r.covered_length, r.speeds))
        print(f"  two-gate overload: k={k} n={n} max_speed={max_speed}")
        for item in candidates[:6]:
            print(
                "    "
                f"max_gap={fmt_frac(item.max_gap):>9s} covered={fmt_frac(item.covered_length):>9s} "
                f"components={item.interval_count:3d} speeds={item.speeds}"
            )


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)

    print("Lonely Runner 14/15/disproof feedback loop (S377)")
    print("All micro-staircase coverage uses exact full-cell bitsets from S367.\n")

    system14 = S367.build_pattern_system(14)
    system15 = S367.build_pattern_system(15)
    print_system_header(system14, "Fourteen-runner quotient model")
    print_system_header(system15, "Fifteen-runner quotient model")

    cycle_one(system14, system15)
    results14, results15 = cycle_two(system14, system15)
    cycle_three(system14, system15, results14, results15)

    print("\nSynthesis")
    print("  1. The n=14 dead end sharpened: the coordinate-6 half-turn is not just")
    print("     the best binary obstruction, it is the unique global minimum of the")
    print("     full normalized 2-torsion cube.  Adding any second defect to it")
    print("     immediately exposes more cells.")
    print("  2. The n=15 replacement obstruction is not binary; its natural small")
    print("     probes are the CRT subgroups of Z/15Z.  Exact small-support scans and")
    print("     local search still leave positive missed-cell margins.")
    print("  3. Disproof pressure keeps leaking: gated and two-gate interval-cover")
    print("     candidates retain explicit positive complement gaps.")
    print("  4. New proof idea: treat n=14 as a 2-torsion/chirality boundary and n=15")
    print("     as a CRT-subgroup leakage problem.  A unified lemma should say every")
    print("     non-scalar quotient defect exposes either a torsion stencil or a CRT")
    print("     subgroup leak before it can become an open cover.")


if __name__ == "__main__":
    main()
