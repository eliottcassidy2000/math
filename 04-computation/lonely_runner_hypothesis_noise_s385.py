#!/usr/bin/env python3
"""
lonely_runner_hypothesis_noise_s385.py

codex-2026-05-31 S385

Creative hypothesis-generation pass for the Lonely Runner thread.

This script deliberately treats repo tangents as noise sources:

* endpoint-transfer/private pivots -> endpoint incidence rank;
* Fejer/Riesz kernel pressure -> critical-radius surplus and kernel centers;
* quotient transport/good-cut peeling -> endpoint peel debt;
* natural operation shadows -> additive/multiplicative speed closure;
* composite denominator anomaly -> debt export across leak layers.

The goal is not to prove LRC or run a huge search.  The goal is to generate
new models that are numerically anchored enough to be worth keeping.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, lcm
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
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()
S379 = SourceFileLoader(
    "lonely_runner_shape_questions_s379",
    str(ROOT / "04-computation" / "lonely_runner_shape_questions_s379.py"),
).load_module()


ONE = Fraction(1, 1)
PRODUCT_SUM_TARGETS = {
    4,
    6,
    8,
    9,
    10,
    12,
    14,
    15,
    16,
    18,
    20,
    21,
    22,
    24,
    25,
    26,
    27,
    28,
}


@dataclass(frozen=True)
class IncidenceFingerprint:
    label: str
    n: int
    speeds: tuple[int, ...]
    classification: str
    max_gap_ratio: Fraction
    critical_ratio: Fraction
    critical_time: Fraction
    critical_mode: str
    endpoint_count: int
    interval_count: int
    unprotected: int
    peel_depth: int
    core_endpoints: int
    protection_rank: int
    boundary_rank: int
    protection_nullity: int
    rank_pressure: Fraction
    min_protector_degree: int
    protector_hist: tuple[tuple[int, int], ...]
    private_endpoint_count: int
    private_interval_count: int
    quotient_layer: str
    additive_gates: int
    multiplicative_gates: int
    divisor_edges: int
    low_layer_debt: int


@dataclass(frozen=True)
class SwapAudit:
    base_label: str
    action: str
    speeds: tuple[int, ...]
    classification: str
    max_gap_ratio: Fraction
    forbidden_length: Fraction
    boundary_witnesses: int
    unprotected: int
    first_unprotected: Fraction | None
    q: int


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_float(x: Fraction) -> str:
    return f"{float(x):.6f}"


def circle(x: Fraction) -> Fraction:
    return x % ONE


def residue(q: int, t: Fraction) -> int:
    return int((circle(t) * q) % q)


def gf2_rank(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        x = row
        while x:
            pivot = x.bit_length() - 1
            if pivot in basis:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                break
    return len(basis)


def bit_rows_for_endpoint_relation(
    endpoints: list[Fraction],
    intervals: list,
    relation: dict,
) -> list[int]:
    index = {interval: i for i, interval in enumerate(intervals)}
    rows: list[int] = []
    for endpoint in endpoints:
        row = 0
        for interval in relation[endpoint]:
            row |= 1 << index[interval]
        rows.append(row)
    return rows


def operation_profile(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    speed_set = set(speeds)
    additive = 0
    multiplicative = 0
    for x, y in combinations(speeds, 2):
        if x + y in speed_set:
            additive += 1
        if x * y in speed_set:
            multiplicative += 1
    for x in speeds:
        if 2 * x in speed_set:
            additive += 1
        if x * x in speed_set:
            multiplicative += 1
    divisor_edges = sum(1 for x in speeds for z in speeds if x < z and z % x == 0)
    return additive, multiplicative, divisor_edges


def classify(row) -> str:
    if row.max_gap > 0:
        return "positive_gap"
    if row.boundary_witness_count:
        return "boundary_only"
    return "open_cover_candidate"


def low_layer_debt(n: int, q: int, endpoints: set[Fraction]) -> int:
    """Count unprotected endpoints in low denominator layers tied to n."""

    debt = 0
    for endpoint in endpoints:
        d = endpoint.denominator
        if d == n or (d % n == 0 and d <= n * n) or (n % 2 == 0 and d <= n * 13):
            debt += 1
    return debt


def endpoint_modulus(endpoints: list[Fraction]) -> int:
    q = 1
    for endpoint in endpoints:
        q = lcm(q, endpoint.denominator)
    return q


def loneliness_at(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(min(circle(v * t), ONE - circle(v * t)) for v in speeds)


def critical_probe(
    speeds: tuple[int, ...], report: S356.GapReport, exact_sum_limit: int = 360
) -> tuple[Fraction, Fraction, str]:
    """Use exact critical radii when cheap; otherwise probe threshold-gap centers."""

    if sum(speeds) <= exact_sum_limit:
        critical_radius, critical_time = S379.critical_loneliness(speeds)
        return critical_radius, critical_time, "exact"

    components = S356.merge_intervals(S356.forbidden_intervals(speeds))
    gaps = S356.circular_gaps(components)
    candidates = {S356.midpoint_on_circle(gap) for gap in gaps}
    if report.witness is not None:
        candidates.add(report.witness)
    if report.boundary_witness is not None:
        candidates.add(report.boundary_witness)
    if not candidates:
        return report.threshold, Fraction(0), "covered_probe"

    best_time = min(candidates)
    best_radius = Fraction(-1)
    for candidate in candidates:
        radius = loneliness_at(speeds, candidate)
        if radius > best_radius or (radius == best_radius and candidate < best_time):
            best_radius = radius
            best_time = candidate
    return best_radius, best_time, "gap_probe"


def incidence_fingerprint(label: str, raw_speeds: tuple[int, ...]) -> IncidenceFingerprint:
    report = S356.report(label, list(raw_speeds))
    speeds = report.speeds
    n = len(speeds) + 1
    endpoints_set, intervals_set, owners, protectors, boundary = S362.build_endpoint_system(speeds)
    endpoints = sorted(endpoints_set)
    intervals = sorted(intervals_set)
    critical_radius, critical_time, critical_mode = critical_probe(speeds, report)
    q = endpoint_modulus(endpoints)
    unprotected = {endpoint for endpoint in endpoints if not protectors[endpoint]}
    peel_layers, core_endpoints, _core_intervals = S362.peel_protection_core(
        q, endpoints_set, intervals_set, protectors, boundary
    )

    protection_rows = bit_rows_for_endpoint_relation(endpoints, intervals, protectors)
    owner_rows = bit_rows_for_endpoint_relation(endpoints, intervals, owners)
    protection_rank = gf2_rank(protection_rows)
    owner_rank = gf2_rank(owner_rows)
    protected_count = len(endpoints) - len(unprotected)
    protection_nullity = protected_count - protection_rank
    rank_pressure = Fraction(protection_rank, max(1, protected_count))

    protector_degrees = Counter(len(protectors[endpoint]) for endpoint in endpoints)
    private_endpoint_count = protector_degrees.get(1, 0)
    protected_by_interval: Counter = Counter()
    for endpoint in endpoints:
        for interval in protectors[endpoint]:
            protected_by_interval[interval] += 1
    private_interval_count = sum(1 for count in protected_by_interval.values() if count == 1)
    additive, multiplicative, divisor_edges = operation_profile(speeds)

    return IncidenceFingerprint(
        label=label,
        n=n,
        speeds=speeds,
        classification=classify(report),
        max_gap_ratio=report.max_gap / report.threshold,
        critical_ratio=critical_radius / report.threshold,
        critical_time=critical_time,
        critical_mode=critical_mode,
        endpoint_count=len(endpoints),
        interval_count=len(intervals),
        unprotected=len(unprotected),
        peel_depth=len(peel_layers),
        core_endpoints=len(core_endpoints),
        protection_rank=protection_rank,
        boundary_rank=owner_rank,
        protection_nullity=protection_nullity,
        rank_pressure=rank_pressure,
        min_protector_degree=min(protector_degrees) if protector_degrees else 0,
        protector_hist=tuple(sorted(protector_degrees.items())),
        private_endpoint_count=private_endpoint_count,
        private_interval_count=private_interval_count,
        quotient_layer=S362.quotient_layer(n, q, unprotected),
        additive_gates=additive,
        multiplicative_gates=multiplicative,
        divisor_edges=divisor_edges,
        low_layer_debt=low_layer_debt(n, q, unprotected),
    )


def sample_sets() -> list[tuple[str, tuple[int, ...]]]:
    return [
        ("initial n=8", tuple(range(1, 8))),
        ("sporadic tight n=8A", (1, 4, 5, 6, 7, 11, 13)),
        ("initial n=14", tuple(range(1, 14))),
        (
            "n14 seven-ladder",
            (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91),
        ),
        (
            "n14 14-ladder debt",
            (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182),
        ),
        ("n14 single-gate", tuple(list(range(1, 13)) + [14])),
        ("initial n=15", tuple(range(1, 15))),
        (
            "n15 3x5 ladder",
            (3, 5, 6, 9, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50),
        ),
        (
            "n15 mixed gates",
            (1, 3, 5, 6, 9, 10, 14, 15, 20, 25, 30, 45, 60, 75),
        ),
    ]


def semantic_pool(n: int) -> tuple[int, ...]:
    pool = set(range(1, 2 * n + 1))
    pool.update(PRODUCT_SUM_TARGETS)
    pool.update(n * q for q in range(1, n))
    for d in range(2, n):
        if n % d == 0:
            pool.update(d * q for q in range(1, n))
    return tuple(sorted(v for v in pool if v > 0))


def semantic_priority(
    old: int, new: int, speeds: tuple[int, ...], n: int
) -> tuple[int, int, int, int]:
    if new <= 2 * n:
        bucket = 0
    elif new in PRODUCT_SUM_TARGETS:
        bucket = 1
    elif new % n == 0:
        bucket = 2
    elif any(n % d == 0 and new % d == 0 for d in range(2, n)):
        bucket = 3
    else:
        bucket = 4

    speed_set = set(speeds)
    relation_bonus = int(
        not any(new == x + y or new == x * y for x in speed_set for y in speed_set)
    )
    return (bucket, relation_bonus, abs(new - old), old)


def swap_audit(
    base_label: str,
    speeds: tuple[int, ...],
    limit: int = 12,
    candidate_cap: int = 84,
) -> list[SwapAudit]:
    n = len(speeds) + 1
    rows: list[SwapAudit] = []
    seen: set[tuple[int, ...]] = set()
    pool = semantic_pool(n)
    candidates: list[tuple[tuple[int, int, int, int], int, int, tuple[int, ...]]] = []
    for old in speeds:
        for new in pool:
            if new in speeds or new == old:
                continue
            candidate = tuple(sorted((set(speeds) - {old}) | {new}))
            if len(candidate) != len(speeds):
                continue
            g = 0
            for v in candidate:
                g = gcd(g, v)
            if g != 1 or candidate in seen:
                continue
            seen.add(candidate)
            candidates.append(
                (semantic_priority(old, new, speeds, n), old, new, candidate)
            )

    for _priority, old, new, candidate in sorted(candidates)[:candidate_cap]:
        report = S356.report(f"{base_label}: {old}->{new}", list(candidate))
        rows.append(
            SwapAudit(
                base_label=base_label,
                action=f"{old}->{new}",
                speeds=report.speeds,
                classification=classify(report),
                max_gap_ratio=report.max_gap / report.threshold,
                forbidden_length=report.forbidden_length,
                boundary_witnesses=report.boundary_witness_count,
                unprotected=-1,
                first_unprotected=None,
                q=report.boundary_modulus,
            )
        )

    rows.sort(
        key=lambda row: (
            {"open_cover_candidate": 0, "boundary_only": 1, "positive_gap": 2}[row.classification],
            row.max_gap_ratio,
            row.boundary_witnesses,
            row.action,
        )
    )
    top = rows[:limit]
    hydrated: list[SwapAudit] = []
    for row in top:
        summary = S360.summarize(list(row.speeds))
        hydrated.append(
            SwapAudit(
                base_label=row.base_label,
                action=row.action,
                speeds=row.speeds,
                classification=row.classification,
                max_gap_ratio=row.max_gap_ratio,
                forbidden_length=row.forbidden_length,
                boundary_witnesses=row.boundary_witnesses,
                unprotected=summary.unprotected_count,
                first_unprotected=summary.first_unprotected,
                q=summary.boundary_modulus,
            )
        )
    return hydrated


def print_noise_sources() -> None:
    print("Lonely Runner creative hypothesis noise pass (S385)")
    print()
    print("Noise sources used as model generators")
    print("  projection defects        -> compare support shadows with rank/residue")
    print("  endpoint transfer         -> private pivots and all-protected cores")
    print("  Fejer/Riesz kernels       -> critical-radius surplus and kernel centers")
    print("  quotient transport        -> peel layers as exported debt")
    print("  natural operation modes   -> additive closure vs divisor closure")
    print("  composite denominator 14  -> protect low layers, watch higher leaks")
    print()


def print_fingerprint_table(fps: list[IncidenceFingerprint]) -> None:
    print("Endpoint-incidence fingerprints")
    header = (
        "label                  n class          gap/th  crit/th mode        endpts unprot "
        "peel coreE rank/null rankP privateE debt divE add mul quotient"
    )
    print(header)
    print("-" * len(header))
    for fp in fps:
        print(
            f"{fp.label:<22} {fp.n:2d} {fp.classification:<14} "
            f"{fmt_float(fp.max_gap_ratio):>7s} {fmt_float(fp.critical_ratio):>8s} "
            f"{fp.critical_mode:<10s} "
            f"{fp.endpoint_count:6d} {fp.unprotected:7d} {fp.peel_depth:4d} "
            f"{fp.core_endpoints:5d} "
            f"{fp.protection_rank:4d}/{fp.protection_nullity:<4d} "
            f"{fmt_float(fp.rank_pressure):>5s} {fp.private_endpoint_count:8d} "
            f"{fp.low_layer_debt:4d} {fp.divisor_edges:4d} "
            f"{fp.additive_gates:3d} {fp.multiplicative_gates:3d} {fp.quotient_layer}"
        )
    print()


def print_semantic_swaps() -> None:
    print("Semantic one-swap provocations")
    print("  Pool = low speeds, product-sum targets, divisor ladders, and n-gates.")
    print()
    bases = [
        (
            "n14 seven-ladder",
            (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91),
        ),
        (
            "n14 14-ladder debt",
            (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182),
        ),
        ("n14 single-gate", tuple(list(range(1, 13)) + [14])),
        (
            "n15 mixed gates",
            (1, 3, 5, 6, 9, 10, 14, 15, 20, 25, 30, 45, 60, 75),
        ),
    ]
    for label, speeds in bases:
        print(f"[{label}]")
        for row in swap_audit(label, speeds, limit=6):
            print(
                "  "
                f"{row.action:<9s} class={row.classification:<12s} "
                f"gap/th={fmt_float(row.max_gap_ratio):>8s} "
                f"forbidden={fmt_frac(row.forbidden_length):>12s} "
                f"bdyW={row.boundary_witnesses:4d} "
                f"unprot={row.unprotected:4d} first={fmt_frac(row.first_unprotected)}"
            )
            print(f"    speeds={row.speeds}")
        print()


def print_models(fps: list[IncidenceFingerprint]) -> None:
    print("New models and reframes generated")
    print()
    tiny_gap = min(
        [fp for fp in fps if fp.classification == "positive_gap"],
        key=lambda fp: fp.max_gap_ratio,
    )
    high_rank = max(fps, key=lambda fp: fp.rank_pressure)
    high_debt = max(fps, key=lambda fp: fp.low_layer_debt)
    print("1. Endpoint-incidence matroid model")
    print(
        "   Treat the protection relation as a binary matrix.  A counterexample "
        "would need no leaf and an unusually self-supporting protection matroid. "
        f"In the sample, the strongest rank pressure is {high_rank.label} "
        f"with rank pressure {fmt_float(high_rank.rank_pressure)}, but it still "
        f"has {high_rank.unprotected} exposed endpoints."
    )
    print()
    print("2. Morse landscape model")
    print(
        "   Visible gap width is a bad proxy for danger.  The tiny-gap example "
        f"{tiny_gap.label} has gap/th={fmt_float(tiny_gap.max_gap_ratio)} but "
        f"critical/th={fmt_float(tiny_gap.critical_ratio)} at "
        f"t={fmt_frac(tiny_gap.critical_time)}.  A tiny gap can be a steep "
        "safe valley, not a near-counterexample plateau."
    )
    print()
    print("3. Debt export model")
    print(
        "   Protecting low quotient layers exports unprotected endpoints upward. "
        f"The largest low-layer debt in the sample is {high_debt.label}: "
        f"{high_debt.low_layer_debt} low-layer exposed endpoints, despite "
        f"{high_debt.divisor_edges} divisor edges."
    )
    print()
    print("4. Two-shadow operation model")
    print(
        "   Addition supplies the Dirichlet/scalar equality spine; multiplication "
        "supplies divisor gates and endpoint protection channels.  Product-sum "
        "targets should be viewed as interference coordinates between these "
        "two shadows, not as numerology."
    )
    print()
    print("5. Abstract-arc mirage model")
    print(
        "   Search separately in abstract circular-arc endpoint systems and in "
        "integer-speed realizations.  Some all-protected arc systems may be "
        "topologically possible but arithmetically unrealizable."
    )
    print()


def print_new_hypotheses() -> None:
    print("Hypotheses to preserve")
    print("  HYP-1842: endpoint protection is controlled by a private-pivot matroid.")
    print("  HYP-1843: near-disproof danger is critical-radius surplus, not gap width.")
    print("  HYP-1844: protecting quotient layers obeys a debt-export law.")
    print()


def main() -> None:
    print_noise_sources()
    fps: list[IncidenceFingerprint] = []
    for label, speeds in sample_sets():
        print(f"  fingerprinting {label}...", flush=True)
        fps.append(incidence_fingerprint(label, speeds))
    print()
    print_fingerprint_table(fps)
    print_semantic_swaps()
    print_models(fps)
    print_new_hypotheses()


if __name__ == "__main__":
    main()
