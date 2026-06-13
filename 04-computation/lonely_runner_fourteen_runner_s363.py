#!/usr/bin/env python3
"""
lonely_runner_fourteen_runner_s363.py

codex-2026-05-31 S363

Exploratory probe for the next Lonely Runner frontier after the currently
proved finite cases: k=13 reduced speeds, i.e. fourteen runners including the
stationary runner, with threshold 1/14.

The point is not to redo the large finite-checking papers.  It is to test a
composite-denominator proof angle:

    14 = 2 * 7.

Any counterexample must defeat the three obvious quotient colorings

    t = 1/2, 1/7, 1/14,

so it must contain an even speed, a speed divisible by 7, and in fact a speed
divisible by 14.  The last condition is the unit-boundary filter from THM-360.

This script audits exact endpoint/core behavior for k=13 boxes and random
14-gated candidates.  It reuses the S356 exact interval union and the S362
endpoint/interval core peeling.
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
K = 13
N = K + 1


@dataclass(frozen=True)
class GateProfile:
    mult_14: int
    mult_7_not_14: int
    even_not_14: int
    unit_mod_14: int
    other: int


@dataclass(frozen=True)
class CandidateSummary:
    speeds: tuple[int, ...]
    gate: GateProfile
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    classification: str
    endpoints: int
    unprotected: int
    unit_skeleton: bool
    peel_depth: int
    core_endpoints: int
    first_layer: str


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def gate_profile(speeds: tuple[int, ...]) -> GateProfile:
    mult_14 = sum(1 for v in speeds if v % 14 == 0)
    mult_7_not_14 = sum(1 for v in speeds if v % 14 == 7)
    even_not_14 = sum(1 for v in speeds if v % 2 == 0 and v % 14 != 0)
    unit_mod_14 = sum(1 for v in speeds if gcd(v, 14) == 1)
    other = len(speeds) - mult_14 - mult_7_not_14 - even_not_14 - unit_mod_14
    return GateProfile(mult_14, mult_7_not_14, even_not_14, unit_mod_14, other)


def passes_necessary_quotient_gates(speeds: tuple[int, ...]) -> bool:
    gate = gate_profile(speeds)
    return gate.mult_14 > 0


def summarize_candidate(speeds: tuple[int, ...]) -> CandidateSummary:
    row = S356.report("k=13-fourteen-runner", list(speeds))
    if row.forbidden_length < ONE:
        classification = "positive_gap"
    elif row.boundary_witness_count:
        classification = "boundary_only"
    else:
        classification = "open_cover"

    descent = S362.summarize(list(speeds))
    first_layer = "none"
    if descent.peel_layers:
        first = descent.peel_layers[0]
        first_layer = (
            f"removed_E={first.removed_endpoints},"
            f"removed_I={first.removed_intervals},"
            f"subgroup={first.removed_subgroup_modulus}"
        )

    return CandidateSummary(
        speeds=tuple(row.speeds),
        gate=gate_profile(tuple(row.speeds)),
        forbidden_length=row.forbidden_length,
        max_gap=row.max_gap,
        gap_ratio=row.max_gap / row.threshold if row.threshold else Fraction(0),
        classification=classification,
        endpoints=descent.endpoint_count,
        unprotected=descent.unprotected_count,
        unit_skeleton=descent.unit_skeleton,
        peel_depth=len(descent.peel_layers),
        core_endpoints=descent.core_endpoint_count,
        first_layer=first_layer,
    )


def print_candidate(label: str, summary: CandidateSummary) -> None:
    print(f"[{label}]")
    print(f"  speeds={summary.speeds}")
    print(
        "  gate="
        f"14:{summary.gate.mult_14} "
        f"7only:{summary.gate.mult_7_not_14} "
        f"even_only:{summary.gate.even_not_14} "
        f"units:{summary.gate.unit_mod_14} "
        f"other:{summary.gate.other}"
    )
    print(
        "  "
        f"classification={summary.classification} "
        f"forbidden_length={fmt_frac(summary.forbidden_length)} "
        f"max_gap={fmt_frac(summary.max_gap)} "
        f"gap/thresh={float(summary.gap_ratio):.6f}"
    )
    print(
        "  "
        f"endpoints={summary.endpoints} unprotected={summary.unprotected} "
        f"unit_skeleton={summary.unit_skeleton} "
        f"peel_depth={summary.peel_depth} core_E={summary.core_endpoints}"
    )
    print(f"  first_peel_layer={summary.first_layer}")
    print()


def scan_box(max_speed: int) -> list[CandidateSummary]:
    total = 0
    gate_total = 0
    positive = 0
    full_measure = 0
    boundary_only = 0
    open_cover = 0
    nonempty_core = 0
    gate_hist: Counter[tuple[int, int, int, int, int]] = Counter()
    best_positive_raw: list[tuple[Fraction, tuple[int, ...]]] = []
    full_examples: list[CandidateSummary] = []

    for combo in combinations(range(1, max_speed + 1), K):
        if not primitive(combo):
            continue
        total += 1
        if not passes_necessary_quotient_gates(combo):
            continue
        gate_total += 1
        gate = gate_profile(combo)
        gate_hist[(gate.mult_14, gate.mult_7_not_14, gate.even_not_14, gate.unit_mod_14, gate.other)] += 1

        row = S356.report("k=13-scan", list(combo))
        if row.forbidden_length < ONE:
            positive += 1
            ratio = row.max_gap / row.threshold
            best_positive_raw.append((ratio, tuple(row.speeds)))
            best_positive_raw.sort(key=lambda item: (item[0], item[1]))
            best_positive_raw = best_positive_raw[:8]
        else:
            summary = summarize_candidate(combo)
            full_measure += 1
            boundary_only += int(summary.classification == "boundary_only")
            open_cover += int(summary.classification == "open_cover")
            nonempty_core += int(summary.core_endpoints > 0)
            full_examples.append(summary)

    best_positive = [summarize_candidate(combo) for _ratio, combo in best_positive_raw]

    print(f"Exhaustive primitive k=13 box max_speed={max_speed}")
    print(f"  total_primitive={total}")
    print(f"  pass_14_gate={gate_total}")
    print(f"  endpoint_summarized={len(best_positive) + len(full_examples)}")
    print(f"  positive_gap={positive}")
    print(f"  full_measure={full_measure}")
    print(f"  boundary_only={boundary_only}")
    print(f"  open_cover={open_cover}")
    print(f"  nonempty_cores={nonempty_core}")
    print(f"  top_gate_profiles={gate_hist.most_common(6)}")
    print("  tightest_positive=")
    for summary in best_positive[:5]:
        print(
            "    "
            f"speeds={summary.speeds} gap/thresh={float(summary.gap_ratio):.6f} "
            f"unprotected={summary.unprotected} peel={summary.peel_depth} "
            f"core_E={summary.core_endpoints}"
        )
    if full_examples:
        print("  full_measure_examples=")
        for summary in full_examples[:5]:
            print(
                "    "
                f"speeds={summary.speeds} class={summary.classification} "
                f"unit={summary.unit_skeleton} core_E={summary.core_endpoints}"
            )
    print()
    return best_positive


def random_gate_candidates(
    max_speed: int,
    trials: int,
    seed: int = 363,
) -> list[CandidateSummary]:
    rng = random.Random(seed)
    best_raw: list[tuple[Fraction, tuple[int, ...], str]] = []
    seen: set[tuple[int, ...]] = set()
    population = list(range(1, max_speed + 1))
    multiples_14 = [v for v in population if v % 14 == 0]

    for _ in range(trials):
        forced = rng.choice(multiples_14)
        rest_pool = [v for v in population if v != forced]
        combo = tuple(sorted([forced] + rng.sample(rest_pool, K - 1)))
        if combo in seen or not primitive(combo):
            continue
        seen.add(combo)
        row = S356.report("k=13-random", list(combo))
        if row.forbidden_length < ONE:
            classification = "positive_gap"
        elif row.boundary_witness_count:
            classification = "boundary_only"
        else:
            classification = "open_cover"
        ratio = row.max_gap / row.threshold if row.threshold else Fraction(0)
        best_raw.append((ratio, tuple(row.speeds), classification))
        best_raw.sort(key=lambda item: (item[2] == "positive_gap", item[0], item[1]))
        best_raw = best_raw[:8]

    best = [summarize_candidate(combo) for _ratio, combo, _classification in best_raw]

    print(f"Random 14-gated k=13 candidates max_speed={max_speed} trials={trials}")
    for summary in best:
        print(
            "  "
            f"speeds={summary.speeds} class={summary.classification} "
            f"gap/thresh={float(summary.gap_ratio):.6f} "
            f"gate14={summary.gate.mult_14} unprotected={summary.unprotected} "
            f"peel={summary.peel_depth} core_E={summary.core_endpoints}"
        )
    print()
    return best


def main() -> None:
    print("Fourteen-runner quotient-gate probe (codex-2026-05-31 S363)")
    print("Target: reduced k=13, threshold 1/14.\n")

    print("Necessary quotient gates")
    print("  t=1/2 forces at least one even speed.")
    print("  t=1/7 forces at least one speed divisible by 7.")
    print("  t=1/14 forces at least one speed divisible by 14.")
    print("  Therefore every k=13 counterexample must pass the 14-gate.\n")

    examples = [
        ("initial segment k=13", tuple(range(1, 14))),
        ("first 14-gated window", tuple(range(1, 13)) + (14,)),
        ("double 14-gated window", tuple(range(1, 12)) + (14, 28)),
        ("odd-seven pressure", (1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 14, 21, 28)),
    ]
    for label, speeds in examples:
        print_candidate(label, summarize_candidate(tuple(sorted(speeds))))

    best = []
    best.extend(scan_box(max_speed=16))

    best.extend(random_gate_candidates(max_speed=80, trials=80))

    print("Interpretation")
    print("  No k=13 open cover or nonempty protection core appeared.")
    print("  The initial segment is boundary-only with the unit skeleton, as")
    print("  expected.  Once the mandatory 14-gate is inserted, the scanned")
    print("  candidates become positive-gap again rather than more dangerous.")
    print("  The next proof angle is a 14=2*7 CRT descent: after a speed divisible")
    print("  by 14 protects the unit layer, divide that protected channel by 14")
    print("  and show the remaining endpoint peel cannot maintain a core.")


if __name__ == "__main__":
    main()
