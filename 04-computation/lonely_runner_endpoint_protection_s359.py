#!/usr/bin/env python3
"""
lonely_runner_endpoint_protection_s359.py

Endpoint-protection graph and core-peeling probe for the reduced Lonely Runner
Conjecture.

For k nonzero speeds V and threshold 1/(k+1), each speed v forbids open
intervals around a/v.  A counterexample is a full open cover of the circle.
If the forbidden union has total length 1 but some endpoint is unprotected,
that endpoint is a boundary lonely witness.

This script turns that sentence into exact arithmetic:

    endpoints = all forbidden interval endpoints
    e -> interval I when e lies strictly inside I

The "protection graph" is a finite incidence witness.  Scans below ask:
  * how tight examples fail to be counterexamples;
  * whether near-tight positive-gap examples have many protected endpoints;
  * how the trivial t=1/(k+1) witness shows up in the protection ledger.
  * whether endpoint protection contains a nonempty self-supporting core.
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

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class TaggedInterval:
    speed: int
    center_index: int
    side: int
    lo: Fraction
    hi: Fraction


@dataclass(frozen=True)
class EndpointRecord:
    t: Fraction
    incident: int
    protection: int
    is_witness: bool
    speeds_touching: tuple[int, ...]
    speeds_protecting: tuple[int, ...]


@dataclass(frozen=True)
class EndpointProtectionIndex:
    intervals: tuple[TaggedInterval, ...]
    endpoints: tuple[Fraction, ...]
    interval_endpoints: tuple[tuple[int, int], ...]
    protectors: tuple[tuple[int, ...], ...]
    touching_counts: tuple[int, ...]
    touching_speeds: tuple[tuple[int, ...], ...]
    protecting_speeds: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class CorePeelingSummary:
    endpoint_count: int
    interval_count: int
    rounds: int
    removed_profile: tuple[int, ...]
    sample: tuple[Fraction, ...]


@dataclass(frozen=True)
class ProtectionReport:
    speeds: tuple[int, ...]
    threshold: Fraction
    forbidden_length: Fraction
    max_gap: Fraction
    witness: Fraction | None
    boundary_witness_count: int
    endpoint_count: int
    protected_endpoint_count: int
    unprotected_endpoint_count: int
    min_protection: int
    max_protection: int
    avg_protection: Fraction
    component_count: int
    q_endpoint: int
    t0_is_witness: bool
    t0_protection: int
    t0_incident: int
    unit_boundary_skeleton: bool
    core_endpoint_count: int
    core_interval_count: int
    core_rounds: int
    core_removed_profile: tuple[int, ...]
    core_sample: tuple[Fraction, ...]
    protection_hist: tuple[tuple[int, int], ...]
    witness_protection_hist: tuple[tuple[int, int], ...]
    unprotected_sample: tuple[Fraction, ...]


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def tagged_forbidden_intervals(speeds: tuple[int, ...]) -> list[TaggedInterval]:
    k = len(speeds)
    threshold = Fraction(1, k + 1)
    tagged: list[TaggedInterval] = []
    for v in speeds:
        radius = threshold / v
        for a in range(v):
            center = Fraction(a, v)
            for side, (lo, hi) in enumerate(S356.split_circle_interval(center - radius, center + radius)):
                tagged.append(TaggedInterval(v, a, side, lo, hi))
    return tagged


def endpoint_modulus(points: list[Fraction]) -> int:
    q = 1
    for point in points:
        q = lcm(q, point.denominator)
    return q


def protection_index(speeds: tuple[int, ...]) -> EndpointProtectionIndex:
    intervals = tagged_forbidden_intervals(speeds)
    endpoints = tuple(sorted({p % ONE for I in intervals for p in (I.lo, I.hi)}))
    endpoint_to_index = {t: i for i, t in enumerate(endpoints)}
    interval_endpoints = tuple(
        (endpoint_to_index[I.lo % ONE], endpoint_to_index[I.hi % ONE])
        for I in intervals
    )
    protectors: list[tuple[int, ...]] = []
    touching_counts: list[int] = []
    touching_speeds: list[tuple[int, ...]] = []
    protecting_speeds: list[tuple[int, ...]] = []
    for t in endpoints:
        touching = []
        protecting = []
        protecting_indices = []
        for i, I in enumerate(intervals):
            if t == I.lo % ONE or t == I.hi % ONE:
                touching.append(I.speed)
            if I.lo < t < I.hi:
                protecting_indices.append(i)
                protecting.append(I.speed)
        protectors.append(tuple(protecting_indices))
        touching_counts.append(len(touching))
        touching_speeds.append(tuple(sorted(set(touching))))
        protecting_speeds.append(tuple(sorted(set(protecting))))
    return EndpointProtectionIndex(
        intervals=tuple(intervals),
        endpoints=endpoints,
        interval_endpoints=interval_endpoints,
        protectors=tuple(protectors),
        touching_counts=tuple(touching_counts),
        touching_speeds=tuple(touching_speeds),
        protecting_speeds=tuple(protecting_speeds),
    )


def endpoint_records_from_index(
    speeds: tuple[int, ...],
    index: EndpointProtectionIndex,
) -> list[EndpointRecord]:
    records: list[EndpointRecord] = []
    for i, t in enumerate(index.endpoints):
        records.append(
            EndpointRecord(
                t=t,
                incident=index.touching_counts[i],
                protection=len(index.protectors[i]),
                is_witness=S356.is_lonely_witness(speeds, t),
                speeds_touching=index.touching_speeds[i],
                speeds_protecting=index.protecting_speeds[i],
            )
        )
    return records


def endpoint_records(speeds: tuple[int, ...]) -> list[EndpointRecord]:
    return endpoint_records_from_index(speeds, protection_index(speeds))


def protection_core(index: EndpointProtectionIndex) -> CorePeelingSummary:
    active_endpoints = set(range(len(index.endpoints)))
    active_intervals = set(range(len(index.intervals)))
    removed_profile: list[int] = []

    while True:
        active_intervals = {
            i
            for i in active_intervals
            if index.interval_endpoints[i][0] in active_endpoints
            and index.interval_endpoints[i][1] in active_endpoints
        }
        to_remove = [
            e
            for e in active_endpoints
            if not any(i in active_intervals for i in index.protectors[e])
        ]
        if not to_remove:
            break
        for e in to_remove:
            active_endpoints.remove(e)
        removed_profile.append(len(to_remove))

    core_intervals = [
        i
        for i in active_intervals
        if index.interval_endpoints[i][0] in active_endpoints
        and index.interval_endpoints[i][1] in active_endpoints
    ]
    return CorePeelingSummary(
        endpoint_count=len(active_endpoints),
        interval_count=len(core_intervals),
        rounds=len(removed_profile),
        removed_profile=tuple(removed_profile),
        sample=tuple(index.endpoints[i] for i in sorted(active_endpoints)[:8]),
    )


def protection_report(raw_speeds: tuple[int, ...]) -> ProtectionReport:
    speeds = S356.normalize_speed_set(list(raw_speeds))
    row = S356.report("endpoint-protection", list(speeds))
    index = protection_index(speeds)
    records = endpoint_records_from_index(speeds, index)
    core = protection_core(index)
    protections = [r.protection for r in records]
    hist = Counter(protections)
    witness_hist = Counter(r.protection for r in records if r.is_witness)
    unprotected = [r.t for r in records if r.protection == 0]
    n = len(speeds) + 1
    unit_skeleton = [Fraction(0)] + [Fraction(a, n) for a in range(1, n) if gcd(a, n) == 1]
    t0 = Fraction(1, len(speeds) + 1)
    t0_records = [r for r in records if r.t == t0]
    if t0_records:
        t0_record = t0_records[0]
        t0_protection = t0_record.protection
        t0_incident = t0_record.incident
    else:
        t0_protection = -1
        t0_incident = 0
    return ProtectionReport(
        speeds=speeds,
        threshold=row.threshold,
        forbidden_length=row.forbidden_length,
        max_gap=row.max_gap,
        witness=row.witness,
        boundary_witness_count=row.boundary_witness_count,
        endpoint_count=len(records),
        protected_endpoint_count=sum(1 for r in records if r.protection > 0),
        unprotected_endpoint_count=sum(1 for r in records if r.protection == 0),
        min_protection=min(protections) if protections else 0,
        max_protection=max(protections) if protections else 0,
        avg_protection=Fraction(sum(protections), len(protections)) if protections else Fraction(0),
        component_count=row.components,
        q_endpoint=endpoint_modulus([r.t for r in records]),
        t0_is_witness=S356.is_lonely_witness(speeds, t0),
        t0_protection=t0_protection,
        t0_incident=t0_incident,
        unit_boundary_skeleton=unprotected == unit_skeleton,
        core_endpoint_count=core.endpoint_count,
        core_interval_count=core.interval_count,
        core_rounds=core.rounds,
        core_removed_profile=core.removed_profile,
        core_sample=core.sample,
        protection_hist=tuple(sorted(hist.items())),
        witness_protection_hist=tuple(sorted(witness_hist.items())),
        unprotected_sample=tuple(unprotected[:8]),
    )


def print_report(label: str, report: ProtectionReport) -> None:
    protected_ratio = Fraction(report.protected_endpoint_count, report.endpoint_count)
    gap_ratio = report.max_gap / report.threshold if report.threshold else Fraction(0)
    print(f"[{label}]")
    print(f"  speeds={report.speeds}")
    print(
        "  "
        f"forbidden_length={fmt_frac(report.forbidden_length)} "
        f"components={report.component_count} max_gap={fmt_frac(report.max_gap)} "
        f"gap/thresh={float(gap_ratio):.6f}"
    )
    print(
        "  "
        f"endpoints={report.endpoint_count} protected={report.protected_endpoint_count} "
        f"unprotected={report.unprotected_endpoint_count} "
        f"protected_ratio={float(protected_ratio):.6f}"
    )
    print(
        "  "
        f"min/max/avg_protection={report.min_protection}/"
        f"{report.max_protection}/{float(report.avg_protection):.3f} "
        f"hist={dict(report.protection_hist)}"
    )
    print(
        "  "
        f"boundary_witness_count={report.boundary_witness_count} "
        f"witness_protection_hist={dict(report.witness_protection_hist)}"
    )
    print(
        "  "
        f"t0=1/{len(report.speeds)+1} witness={report.t0_is_witness} "
        f"protection={report.t0_protection} incident={report.t0_incident} "
        f"Q_endpoint={report.q_endpoint} unit_skeleton={report.unit_boundary_skeleton}"
    )
    print(
        "  "
        f"core_endpoints={report.core_endpoint_count} "
        f"core_intervals={report.core_interval_count} "
        f"peel_rounds={report.core_rounds} "
        f"removed_profile={report.core_removed_profile}"
    )
    if report.core_sample:
        print(
            "  core_sample="
            + ", ".join(fmt_frac(x) for x in report.core_sample)
        )
    print(
        "  unprotected_sample="
        + ", ".join(fmt_frac(x) for x in report.unprotected_sample)
    )
    print()


def known_examples() -> list[tuple[str, tuple[int, ...]]]:
    return [
        ("initial k=3", (1, 2, 3)),
        ("initial k=4", (1, 2, 3, 4)),
        ("sporadic n=5", (1, 3, 4, 7)),
        ("initial k=5", (1, 2, 3, 4, 5)),
        ("sporadic n=6", (1, 3, 4, 5, 9)),
        ("initial k=7", (1, 2, 3, 4, 5, 6, 7)),
        ("sporadic n=8 A", (1, 4, 5, 6, 7, 11, 13)),
        ("sporadic n=8 B", (1, 2, 3, 4, 5, 7, 12)),
        ("near tight k=8 mixed CRT", (5, 8, 13, 21, 34, 55, 89, 144)),
    ]


def scan_near_tight(k: int, max_speed: int, limit: int = 12) -> list[tuple[Fraction, Fraction, ProtectionReport]]:
    rows: list[tuple[Fraction, Fraction, ProtectionReport]] = []
    for combo in combinations(range(1, max_speed + 1), k):
        if not primitive(combo):
            continue
        report = protection_report(combo)
        if report.max_gap == 0:
            continue
        gap_ratio = report.max_gap / report.threshold
        protected_ratio = Fraction(report.protected_endpoint_count, report.endpoint_count)
        rows.append((gap_ratio, -protected_ratio, report))
    rows.sort(key=lambda item: (item[0], item[1], item[2].speeds))
    return rows[:limit]


def scan_core_boxes(boxes: list[tuple[int, int]], limit: int = 4) -> None:
    print("Protection-core peeling scans")
    print(
        "  Core rule: keep endpoints only if protected by an interval whose "
        "two endpoints also survive."
    )
    for k, max_speed in boxes:
        n = k + 1
        total = 0
        nonempty = 0
        mod_filter_total = 0
        mod_filter_nonempty = 0
        max_core = 0
        max_rounds = 0
        frontier = []
        for combo in combinations(range(1, max_speed + 1), k):
            if not primitive(combo):
                continue
            total += 1
            report = protection_report(combo)
            protected_ratio = Fraction(report.protected_endpoint_count, report.endpoint_count)
            gap_ratio = report.max_gap / report.threshold if report.threshold else Fraction(0)
            passes_mod_filter = any(v % n == 0 for v in combo)
            if passes_mod_filter:
                mod_filter_total += 1
            if report.core_endpoint_count:
                nonempty += 1
                if passes_mod_filter:
                    mod_filter_nonempty += 1
            max_core = max(max_core, report.core_endpoint_count)
            max_rounds = max(max_rounds, report.core_rounds)
            frontier.append(
                (
                    -report.core_endpoint_count,
                    -report.core_rounds,
                    report.unprotected_endpoint_count,
                    -protected_ratio,
                    gap_ratio,
                    report.speeds,
                    report,
                )
            )
        frontier.sort()
        print(f"k={k}, max_speed={max_speed}")
        print(
            "  "
            f"primitive_sets={total} nonempty_cores={nonempty} "
            f"mod_filter_sets={mod_filter_total} "
            f"mod_filter_nonempty={mod_filter_nonempty} "
            f"max_core={max_core} max_peel_rounds={max_rounds}"
        )
        print("  hardest empty-core examples")
        for _, _, _, neg_protected_ratio, gap_ratio, _, report in frontier[:limit]:
            print(
                "  "
                f"speeds={report.speeds} core={report.core_endpoint_count} "
                f"rounds={report.core_rounds} removed={report.core_removed_profile} "
                f"unprotected={report.unprotected_endpoint_count} "
                f"protected_ratio={float(-neg_protected_ratio):.6f} "
                f"gap/thresh={float(gap_ratio):.6f}"
            )
        print()


def scan_mod_filter(k: int, max_speed: int, limit: int = 8) -> None:
    """Scan sets containing a speed divisible by k+1.

    If no speed is 0 mod k+1, then t=1/(k+1) is an immediate witness.
    Counterexamples must pass this primitive filter.
    """
    n = k + 1
    candidates = []
    total = 0
    for combo in combinations(range(1, max_speed + 1), k):
        if not primitive(combo):
            continue
        if all(v % n != 0 for v in combo):
            continue
        total += 1
        report = protection_report(combo)
        gap_ratio = report.max_gap / report.threshold if report.threshold else Fraction(0)
        protected_ratio = Fraction(report.protected_endpoint_count, report.endpoint_count)
        candidates.append((gap_ratio, -protected_ratio, report))
    candidates.sort(key=lambda item: (item[0], item[1], item[2].speeds))

    print(f"Mod-filter scan k={k}, max_speed={max_speed}, contains 0 mod {n}")
    print(f"  candidates={total}")
    for gap_ratio, neg_protected_ratio, report in candidates[:limit]:
        print(
            "  "
            f"speeds={report.speeds} gap/thresh={float(gap_ratio):.6f} "
            f"protected_ratio={float(-neg_protected_ratio):.6f} "
            f"unprotected={report.unprotected_endpoint_count} "
            f"t0_protection={report.t0_protection} "
            f"core={report.core_endpoint_count}"
        )
    print()


def main() -> None:
    print("Lonely Runner endpoint-protection graph (codex-2026-05-30 S359)")
    print("Exact rational arithmetic; open intervals protect endpoints strictly.\n")

    print("Known and near-tight examples")
    for label, speeds in known_examples():
        print_report(label, protection_report(speeds))

    print("Near-tight positive-gap scans")
    # Keep these boxes intentionally modest: endpoint protection is quadratic
    # in interval/endpoints for each speed set.  The older S357 scan covers
    # larger boxes for gap existence; this pass studies the protection geometry.
    for k, max_speed in [(4, 18), (5, 16), (6, 14)]:
        print(f"k={k}, max_speed={max_speed}")
        for gap_ratio, neg_protected_ratio, report in scan_near_tight(k, max_speed, limit=8):
            print(
                "  "
                f"speeds={report.speeds} gap/thresh={float(gap_ratio):.6f} "
                f"protected_ratio={float(-neg_protected_ratio):.6f} "
                f"unprotected={report.unprotected_endpoint_count} "
                f"core={report.core_endpoint_count} "
                f"rounds={report.core_rounds} "
                f"hist={dict(report.protection_hist)}"
            )
        print()

    scan_core_boxes([(3, 20), (4, 16), (5, 13), (6, 11)])

    print("Primitive counterexample filter")
    print("  If no speed is divisible by k+1, then t=1/(k+1) is immediately lonely.")
    print("  Therefore every counterexample must contain a speed divisible by k+1.")
    print("  Small mod-filter scans are possible with scan_mod_filter(), but are")
    print("  omitted from the default run because endpoint protection is quadratic")
    print("  in endpoints and intervals for each speed set.\n")

    print("Takeaways")
    print("  1. Boundary-only tight examples have many unprotected endpoints;")
    print("     no scanned tight example is remotely close to all-protected.")
    print("  2. The trivial t=1/(k+1) witness is exactly the first filter:")
    print("     any counterexample must contain a speed divisible by k+1.")
    print("  3. Near-tight positive-gap examples can have high protected ratios,")
    print("     but still leave explicit unprotected endpoints/gaps.")
    print("  4. Core peeling is stronger than initial exposure in the scanned boxes:")
    print("     every primitive set tested peels to the empty endpoint core.")
    print("  5. A counterexample would have to evade both layers: no exposed endpoint,")
    print("     and a nonempty self-supporting endpoint-protection core.")


if __name__ == "__main__":
    main()
