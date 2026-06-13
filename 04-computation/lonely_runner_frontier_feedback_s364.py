#!/usr/bin/env python3
"""
lonely_runner_frontier_feedback_s364.py

opus-2026-05-31 S364

Feedback-loop research sprint for the next Lonely Runner frontiers.

The loop intentionally cycles through three lanes:

1. k=13 / n=14, the first public-open finite case.
2. k=14 / n=15, a deliberately creative next-case analogue.
3. counterexample pressure tests, trying to make open covers rather than proofs.

Each lane is designed to hit a concrete wall and then hand one new idea to the
next lane.  The computations stay exact for interval unions and endpoint-core
peeling, using the S356/S362 helpers.
"""

from __future__ import annotations

from collections import Counter
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
class CompactSummary:
    speeds: tuple[int, ...]
    n: int
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    endpoints: int
    unprotected: int
    peel_depth: int
    core_endpoints: int


@dataclass(frozen=True)
class CoarseVectorResult:
    n: int
    vector: tuple[int, ...]
    coarse_blocked: int
    coarse_total: int
    micro_denominator: int | None
    micro_s: int | None
    micro_r: int | None
    micro_bins: tuple[int, ...] | None
    micro_residues: tuple[int, ...] | None
    no_wrap_interval: tuple[Fraction, Fraction] | None


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def classify_report(row) -> str:
    if row.forbidden_length < ONE:
        return "positive_gap"
    if row.boundary_witness_count:
        return "boundary_only"
    return "open_cover"


def compact(raw_speeds: tuple[int, ...] | list[int]) -> CompactSummary:
    row = S356.report("frontier-feedback", list(raw_speeds))
    descent = S362.summarize(list(row.speeds))
    gap_ratio = row.max_gap / row.threshold if row.threshold else Fraction(0)
    return CompactSummary(
        speeds=tuple(row.speeds),
        n=len(row.speeds) + 1,
        classification=classify_report(row),
        forbidden_length=row.forbidden_length,
        max_gap=row.max_gap,
        gap_ratio=gap_ratio,
        endpoints=descent.endpoint_count,
        unprotected=descent.unprotected_count,
        peel_depth=len(descent.peel_layers),
        core_endpoints=descent.core_endpoint_count,
    )


def print_compact(label: str, summary: CompactSummary) -> None:
    print(f"[{label}]")
    print(f"  speeds={summary.speeds}")
    print(
        "  "
        f"n={summary.n} class={summary.classification} "
        f"forbidden_length={fmt_frac(summary.forbidden_length)} "
        f"max_gap={fmt_frac(summary.max_gap)} "
        f"gap/thresh={float(summary.gap_ratio):.6f}"
    )
    print(
        "  "
        f"endpoints={summary.endpoints} "
        f"unprotected={summary.unprotected} "
        f"peel_depth={summary.peel_depth} "
        f"core_E={summary.core_endpoints}"
    )


def divisors(n: int) -> list[int]:
    return [d for d in range(2, n + 1) if n % d == 0]


def gate_profile(speeds: tuple[int, ...], n: int) -> dict[int, int]:
    return {d: sum(1 for v in speeds if v % d == 0) for d in divisors(n)}


def gate_scan(k: int, max_speed: int, gate: int) -> list[CompactSummary]:
    total = 0
    gated = 0
    positive = 0
    full = 0
    boundary = 0
    open_cover = 0
    best_raw: list[tuple[Fraction, tuple[int, ...]]] = []
    full_examples: list[CompactSummary] = []
    gate_hist: Counter[tuple[tuple[int, int], ...]] = Counter()

    for combo in combinations(range(1, max_speed + 1), k):
        if not primitive(combo):
            continue
        total += 1
        if not any(v % gate == 0 for v in combo):
            continue
        gated += 1
        gate_hist[tuple(sorted(gate_profile(combo, gate).items()))] += 1
        row = S356.report("gate-scan", list(combo))
        cls = classify_report(row)
        if cls == "positive_gap":
            positive += 1
            ratio = row.max_gap / row.threshold
            best_raw.append((ratio, tuple(row.speeds)))
            best_raw.sort(key=lambda item: (item[0], item[1]))
            best_raw = best_raw[:6]
        else:
            full += 1
            summary = compact(tuple(row.speeds))
            full_examples.append(summary)
            boundary += int(cls == "boundary_only")
            open_cover += int(cls == "open_cover")

    print(f"Exact primitive gate scan k={k}, n={k + 1}, max_speed={max_speed}")
    print(f"  total_primitive={total}")
    print(f"  pass_{gate}_gate={gated}")
    print(f"  positive_gap={positive}")
    print(f"  full_measure={full}")
    print(f"  boundary_only={boundary}")
    print(f"  open_cover={open_cover}")
    print(f"  top_gate_profiles={gate_hist.most_common(5)}")
    best = [compact(combo) for _ratio, combo in best_raw]
    print("  tightest_positive=")
    for summary in best:
        print(
            "    "
            f"speeds={summary.speeds} "
            f"gap/thresh={float(summary.gap_ratio):.6f} "
            f"unprotected={summary.unprotected} "
            f"peel={summary.peel_depth} "
            f"core_E={summary.core_endpoints}"
        )
    if full_examples:
        print("  full_measure_examples=")
        for summary in full_examples[:5]:
            print(
                "    "
                f"speeds={summary.speeds} class={summary.classification} "
                f"unprotected={summary.unprotected} core_E={summary.core_endpoints}"
            )
    print()
    return best + full_examples[:3]


def gate_overload_family(k: int, n: int, max_gate_multiple: int = 6) -> list[CompactSummary]:
    """Try the obvious disproof construction: stuff the set with gate multiples."""

    out: list[CompactSummary] = []
    base_pool = [v for v in range(1, n) if v != n]
    for mult_count in range(1, min(max_gate_multiple, k) + 1):
        multiples = tuple(n * j for j in range(1, mult_count + 1))
        fillers = tuple(v for v in base_pool if v not in multiples)[: k - mult_count]
        speeds = tuple(sorted(multiples + fillers))
        if len(speeds) == k and primitive(speeds):
            out.append(compact(speeds))
    return out


def random_gate_search(
    k: int, n: int, max_speed: int, trials: int, seed: int
) -> list[CompactSummary]:
    rng = random.Random(seed)
    population = list(range(1, max_speed + 1))
    gate_values = [v for v in population if v % n == 0]
    best_raw: list[tuple[int, Fraction, tuple[int, ...], str]] = []
    seen: set[tuple[int, ...]] = set()

    for _ in range(trials):
        forced = rng.choice(gate_values)
        rest = [v for v in population if v != forced]
        combo = tuple(sorted([forced] + rng.sample(rest, k - 1)))
        if combo in seen or not primitive(combo):
            continue
        seen.add(combo)
        row = S356.report("random-gate-search", list(combo))
        cls = classify_report(row)
        ratio = row.max_gap / row.threshold if row.threshold else Fraction(0)
        # Full covers would sort first.  Otherwise we minimize gap ratio.
        class_rank = {"open_cover": 0, "boundary_only": 1, "positive_gap": 2}[cls]
        best_raw.append((class_rank, ratio, tuple(row.speeds), cls))
        best_raw.sort(key=lambda item: (item[0], item[1], item[2]))
        best_raw = best_raw[:8]

    best = [compact(combo) for _rank, _ratio, combo, _cls in best_raw]
    print(f"Random counterexample pressure k={k}, n={n}, max_speed={max_speed}, trials={trials}")
    for summary in best:
        print(
            "  "
            f"speeds={summary.speeds} "
            f"class={summary.classification} "
            f"gap/thresh={float(summary.gap_ratio):.6f} "
            f"gate_mults={gate_profile(summary.speeds, n)} "
            f"unprotected={summary.unprotected} "
            f"peel={summary.peel_depth} "
            f"core_E={summary.core_endpoints}"
        )
    print()
    return best


def peel_layers(label: str, speeds: tuple[int, ...]) -> None:
    summary = S362.summarize(list(speeds))
    print(f"{label}: full endpoint peel")
    print(
        "  "
        f"speeds={summary.speeds} class={summary.classification} "
        f"Q={summary.q} unprotected={summary.unprotected_count} "
        f"quotient_layer={summary.quotient_layer} "
        f"core_E={summary.core_endpoint_count}"
    )
    for layer in summary.peel_layers:
        print(
            "  "
            f"layer={layer.step:02d} "
            f"removed_E={layer.removed_endpoints:4d} "
            f"removed_I={layer.removed_intervals:4d} "
            f"remaining_E={layer.remaining_endpoints:4d} "
            f"remaining_I={layer.remaining_intervals:4d} "
            f"subgroup={layer.removed_subgroup_modulus}"
        )
    print()


def rbf(n: int, k: int, r: int, denominator: int) -> tuple[int, ...]:
    return tuple((n * ((i * r) % denominator)) // denominator for i in range(1, k + 1))


def residues_after_shift(
    vector: tuple[int, ...], n: int, s: int, bins: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple((s * value + bin_value) % n for value, bin_value in zip(vector, bins))


def good_residue(residue: int, n: int) -> bool:
    return residue not in (0, n - 1)


def blocked_count(vector: tuple[int, ...], n: int, denominator: int) -> int:
    blocked = 0
    for s in range(n):
        for r in range(denominator):
            bins = rbf(n, len(vector), r, denominator)
            residues = residues_after_shift(vector, n, s, bins)
            if any(not good_residue(residue, n) for residue in residues):
                blocked += 1
    return blocked


def no_wrap_interval(n: int, pattern: tuple[int, ...]) -> tuple[Fraction, Fraction] | None:
    lo = Fraction(0)
    hi = Fraction(1)
    for i, value in enumerate(pattern, start=1):
        lo = max(lo, Fraction(value, n * i))
        hi = min(hi, Fraction(value + 1, n * i))
    if lo < hi and hi <= Fraction(1, len(pattern)):
        return (lo, hi)
    return None


def resolve_on_denominators(
    vector: tuple[int, ...], n: int, denominators: list[int]
) -> tuple[int, int, tuple[int, ...], tuple[int, ...]] | None:
    k = len(vector)
    for denominator in denominators:
        for s in range(n):
            for r in range(denominator):
                bins = rbf(n, k, r, denominator)
                residues = residues_after_shift(vector, n, s, bins)
                if all(good_residue(residue, n) for residue in residues):
                    return denominator, s, r, bins, residues
    return None


def find_coarse_obstruction(n: int, seed: int, trials: int = 4000) -> CoarseVectorResult:
    rng = random.Random(seed)
    k = n - 1
    best_vector: tuple[int, ...] | None = None
    best_blocked = -1
    total = n * n

    structured_seeds: list[tuple[int, ...]] = []
    if n == 14:
        structured_seeds.append((7, 4, 9, 6, 7, 8, 5, 0, 1, 12, 13, 12, 7))
    if n == 15:
        structured_seeds.extend(
            [
                tuple((i * 7 + (i // 3)) % n for i in range(1, n)),
                tuple((i * 4 + (i // 5)) % n for i in range(1, n)),
                tuple((i * i + 2 * i) % n for i in range(1, n)),
            ]
        )

    for vector in structured_seeds:
        count = blocked_count(vector, n, n)
        if count > best_blocked:
            best_blocked = count
            best_vector = vector

    for _ in range(trials):
        vector = tuple(rng.randrange(n) for _ in range(k))
        count = blocked_count(vector, n, n)
        if count > best_blocked:
            best_blocked = count
            best_vector = vector
        if count == total:
            break

    assert best_vector is not None
    resolution = resolve_on_denominators(
        best_vector,
        n,
        [d for d in [17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 127, 149, 181, 211, 251, 307, 401] if d != n],
    )
    if resolution is None:
        return CoarseVectorResult(
            n=n,
            vector=best_vector,
            coarse_blocked=best_blocked,
            coarse_total=total,
            micro_denominator=None,
            micro_s=None,
            micro_r=None,
            micro_bins=None,
            micro_residues=None,
            no_wrap_interval=None,
        )
    denominator, s, r, bins, residues = resolution
    return CoarseVectorResult(
        n=n,
        vector=best_vector,
        coarse_blocked=best_blocked,
        coarse_total=total,
        micro_denominator=denominator,
        micro_s=s,
        micro_r=r,
        micro_bins=bins,
        micro_residues=residues,
        no_wrap_interval=no_wrap_interval(n, bins),
    )


def print_coarse_obstruction(label: str, result: CoarseVectorResult) -> None:
    print(label)
    print(f"  n={result.n} vector={result.vector}")
    print(f"  coarse_blocked={result.coarse_blocked}/{result.coarse_total}")
    if result.micro_denominator is None:
        print("  micro_resolution=not found in tested denominators")
    else:
        print(
            "  "
            f"micro_resolution=denom {result.micro_denominator}, "
            f"s={result.micro_s}, r={result.micro_r}"
        )
        print(f"  bins={result.micro_bins}")
        print(f"  residues={result.micro_residues}")
        if result.no_wrap_interval is not None:
            lo, hi = result.no_wrap_interval
            print(
                "  "
                f"no_wrap_cell=[{fmt_frac(lo)}, {fmt_frac(hi)}) "
                f"width={fmt_frac(hi - lo)}"
            )
    print()


def main() -> None:
    print("Lonely Runner frontier feedback loop (opus-2026-05-31 S364)")
    print("All interval/core computations use exact rational arithmetic.\n")

    print("Cycle 1: attack k=13 / n=14 by extending the CRT gate scan")
    # max_speed=17 is a strict extension beyond S363's max_speed=16, while
    # still small enough for exact Fraction unions in an interactive session.
    best14 = gate_scan(k=13, max_speed=17, gate=14)
    if best14:
        peel_layers("Dead-end audit for tightest n=14 gated positive example", best14[0].speeds)

    print("Cycle 2: forced creative jump to k=14 / n=15")
    examples15 = [
        ("initial segment k=14", tuple(range(1, 15))),
        ("first 15-gated window", tuple(range(1, 14)) + (15,)),
        ("double 15-gated window", tuple(range(1, 13)) + (15, 30)),
        ("3x5 pressure", (1, 2, 3, 4, 5, 6, 9, 10, 12, 15, 20, 25, 30, 45)),
    ]
    for label, speeds in examples15:
        print_compact(label, compact(tuple(sorted(speeds))))
    print()
    best15 = gate_scan(k=14, max_speed=18, gate=15)

    print("Cycle 3: disproof-construction pressure by overloading gate multiples")
    for k, n in [(13, 14), (14, 15)]:
        print(f"Gate-overload family k={k}, n={n}")
        for summary in gate_overload_family(k, n):
            print(
                "  "
                f"mults={gate_profile(summary.speeds, n)[n]} "
                f"class={summary.classification} "
                f"gap/thresh={float(summary.gap_ratio):.6f} "
                f"unprotected={summary.unprotected} "
                f"peel={summary.peel_depth} "
                f"speeds={summary.speeds}"
            )
        print()

    print("Cycle 4: back to n=14, micro-staircase obstruction instead of endpoint gates")
    print_coarse_obstruction(
        "n=14 coarse r/14 obstruction and prime-grid repair",
        find_coarse_obstruction(n=14, seed=36414, trials=200),
    )

    print("Cycle 5: forced n=15 analogue, looking for a 3x5 micro-staircase")
    print_coarse_obstruction(
        "n=15 coarse r/15 obstruction and prime-grid repair",
        find_coarse_obstruction(n=15, seed=36415, trials=6000),
    )

    print("Cycle 6: counterexample search after both staircase lanes hit finite certificates")
    random_gate_search(k=13, n=14, max_speed=80, trials=30, seed=3641314)
    random_gate_search(k=14, n=15, max_speed=90, trials=30, seed=3641415)

    print("Cycle 7: synthesis")
    print("  n=14 dead end: larger exact gated boxes still leaked; the long peel")
    print("  has no terminal endpoint core, so the route needs a structural descent.")
    print("  n=15 new idea: the same composite phenomenon appears with 15=3*5;")
    print("  it has a stronger CRT gate but likely needs a two-prime micro-staircase.")
    print("  disproof pressure: gate-overload and random searches made positive gaps,")
    print("  not open covers.  The best candidates are boundary-adjacent leaks,")
    print("  suggesting future attacks should rank peel depth and cell width together.")


if __name__ == "__main__":
    main()
