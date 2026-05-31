#!/usr/bin/env python3
"""
lonely_runner_shape_questions_s379.py

codex-2026-05-31 S379

Conceptual probe for the Lonely Runner thread.

The point is not to search another large box.  The point is to compare
representations:

  speed list                arithmetic data
  pulled-back intervals      a circular cover
  endpoints                  a finite incidence/protection graph
  interval overlaps          a circular-arc nerve shadow
  max-min landscape          the actual loneliness height

The script asks which representation feels closest to "the shape that
determines" the instance.
"""

from __future__ import annotations

from collections import Counter, deque
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
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class ShapeReport:
    label: str
    speeds: tuple[int, ...]
    n: int
    threshold: Fraction
    critical_radius: Fraction
    critical_ratio: Fraction
    critical_time: Fraction
    classification: str
    forbidden_length: Fraction
    components: int
    gap_count: int
    max_gap_ratio: Fraction
    boundary_witnesses: int
    boundary_modulus: int
    unit_skeleton: bool
    quotient_layer: str
    endpoint_count: int
    unprotected: int
    protection_pressure: Fraction
    peel_depth: int
    core_endpoints: int
    owner_hist: tuple[tuple[int, int], ...]
    protector_hist: tuple[tuple[int, int], ...]
    interval_count: int
    overlap_edges: int
    overlap_components: int
    overlap_density: Fraction
    additive_gates: int
    multiplicative_gates: int
    divisor_edges: int


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_float(x: Fraction) -> str:
    return f"{float(x):.6f}"


def circle_point(x: Fraction) -> Fraction:
    return x % ONE


def dist_to_integer(x: Fraction) -> Fraction:
    r = circle_point(x)
    return min(r, ONE - r)


def loneliness_at(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(dist_to_integer(v * t) for v in speeds)


def critical_loneliness(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    """Return exact max_t min_v ||v t|| and one maximizing time.

    The lower envelope of the piecewise-linear distance functions changes at
    half-integer breakpoints and at equalities between signed linear branches.
    Enumerating those candidate times is small for the examples in this probe
    and gives an exact answer over Fractions.
    """

    candidates: set[Fraction] = {Fraction(0), ONE}
    for v in speeds:
        for a in range(v + 1):
            candidates.add(Fraction(a, v))
        for a in range(v):
            candidates.add(Fraction(2 * a + 1, 2 * v))

    branches: list[tuple[int, int, int]] = []
    for v in speeds:
        for a in range(v + 1):
            branches.append((1, v, a))
            branches.append((-1, v, a))

    for s1, v1, a1 in branches:
        for s2, v2, a2 in branches:
            denom = s1 * v1 - s2 * v2
            if denom == 0:
                continue
            numer = s1 * a1 - s2 * a2
            t = Fraction(numer, denom)
            if 0 <= t <= 1:
                candidates.add(t)

    best_t = Fraction(0)
    best_r = Fraction(-1)
    for t in candidates:
        r = loneliness_at(speeds, t)
        if r > best_r or (r == best_r and t < best_t):
            best_r = r
            best_t = circle_point(t)
    return best_r, best_t


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

    divisor_edges = sum(
        1 for x in speeds for z in speeds if x < z and z % x == 0
    )
    return additive, multiplicative, divisor_edges


def interval_overlap_stats(
    speeds: tuple[int, ...],
) -> tuple[int, int, Fraction]:
    intervals = list(S362.interval_keys(speeds))
    m = len(intervals)
    if m <= 1:
        return 0, m, Fraction(0)

    parent = list(range(m))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    edges = 0
    for i, left in enumerate(intervals):
        c1 = S362.interval_center(left)
        r1 = S362.interval_radius(speeds, left)
        for j in range(i + 1, m):
            right = intervals[j]
            c2 = S362.interval_center(right)
            r2 = S362.interval_radius(speeds, right)
            if S362.circular_distance(c1, c2) < r1 + r2:
                edges += 1
                union(i, j)

    components = len({find(i) for i in range(m)})
    possible = m * (m - 1) // 2
    return edges, components, Fraction(edges, possible)


def classify(row) -> str:
    if row.max_gap > 0:
        return "positive_gap"
    if row.boundary_witness_count:
        return "boundary_only"
    return "open_cover_candidate"


def shape_report(label: str, raw_speeds: tuple[int, ...]) -> ShapeReport:
    base = S356.report(label, list(raw_speeds))
    descent = S362.summarize(list(base.speeds))
    endpoints, intervals, owners, protectors, _boundary = S362.build_endpoint_system(
        base.speeds
    )
    critical_radius, critical_time = critical_loneliness(base.speeds)
    additive, multiplicative, divisor_edges = operation_profile(base.speeds)
    overlap_edges, overlap_components, overlap_density = interval_overlap_stats(
        base.speeds
    )
    endpoint_count = len(endpoints)
    unprotected = descent.unprotected_count
    threshold = base.threshold
    return ShapeReport(
        label=label,
        speeds=base.speeds,
        n=len(base.speeds) + 1,
        threshold=threshold,
        critical_radius=critical_radius,
        critical_ratio=critical_radius / threshold,
        critical_time=critical_time,
        classification=classify(base),
        forbidden_length=base.forbidden_length,
        components=base.components,
        gap_count=base.gap_count,
        max_gap_ratio=base.max_gap / threshold if threshold else Fraction(0),
        boundary_witnesses=base.boundary_witness_count,
        boundary_modulus=base.boundary_modulus,
        unit_skeleton=descent.unit_skeleton,
        quotient_layer=descent.quotient_layer,
        endpoint_count=endpoint_count,
        unprotected=unprotected,
        protection_pressure=Fraction(endpoint_count - unprotected, endpoint_count),
        peel_depth=len(descent.peel_layers),
        core_endpoints=descent.core_endpoint_count,
        owner_hist=tuple(sorted(Counter(len(owners[e]) for e in endpoints).items())),
        protector_hist=tuple(
            sorted(Counter(len(protectors[e]) for e in endpoints).items())
        ),
        interval_count=len(intervals),
        overlap_edges=overlap_edges,
        overlap_components=overlap_components,
        overlap_density=overlap_density,
        additive_gates=additive,
        multiplicative_gates=multiplicative,
        divisor_edges=divisor_edges,
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
            "n14 single-gate",
            tuple(list(range(1, 13)) + [14]),
        ),
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


def print_representation_map() -> None:
    print("WHAT DO THE OBJECTS REPRESENT?")
    print("=" * 88)
    rows = [
        (
            "speed v",
            "a character t -> v*t on R/Z; one frequency of the forbidden wave",
        ),
        (
            "speed set V",
            "a family of characters; equivalently a pullback machine from torus to time circle",
        ),
        (
            "threshold 1/(k+1)",
            "Dirichlet equality radius; the point where volume count alone is tight",
        ),
        (
            "forbidden interval",
            "one runner's unsafe fiber pulled back to the common time circle",
        ),
        (
            "endpoint",
            "a boundary event where one runner is exactly at distance 1/(k+1)",
        ),
        (
            "protected endpoint",
            "a boundary event hidden inside another runner's unsafe interior",
        ),
        (
            "boundary-only example",
            "the cover has no open gap, but its open intervals fail at finite endpoints",
        ),
        (
            "counterexample",
            "a full open cover whose every endpoint is protected",
        ),
        (
            "scalar ramp",
            "Dirichlet equality spine / gauge orbit in the micro-staircase quotient",
        ),
        (
            "torsion defect",
            "a quotient leak through a proper subgroup of Z/nZ",
        ),
    ]
    for name, meaning in rows:
        print(f"  {name:<22} -> {meaning}")
    print()


def print_shape_table(reports: list[ShapeReport]) -> None:
    print("REPRESENTATION COMPARISON TABLE")
    print("=" * 88)
    header = (
        "label                  n class          crit/th maxgap/th bdyW "
        "unprot peel overlap dens divE addG mulG quotient"
    )
    print(header)
    print("-" * len(header))
    for row in reports:
        print(
            f"{row.label:<22} {row.n:>2} {row.classification:<14} "
            f"{fmt_float(row.critical_ratio):>7} "
            f"{fmt_float(row.max_gap_ratio):>9} {row.boundary_witnesses:>4} "
            f"{row.unprotected:>7} {row.peel_depth:>4} "
            f"{row.overlap_components:>7} {fmt_float(row.overlap_density):>7} "
            f"{row.divisor_edges:>4} {row.additive_gates:>4} "
            f"{row.multiplicative_gates:>4} {row.quotient_layer}"
        )
    print()


def print_endpoint_fingerprints(reports: list[ShapeReport]) -> None:
    print("ENDPOINT SHAPE FINGERPRINTS")
    print("=" * 88)
    print(
        "The endpoint graph remembers something the speed list hides: which "
        "boundary failures are isolated, multiply-owned, or deeply protected."
    )
    print()
    for row in reports:
        print(f"[{row.label}]")
        print(
            "  "
            f"Q={row.boundary_modulus} endpoints={row.endpoint_count} "
            f"intervals={row.interval_count} components={row.components} "
            f"critical={fmt_frac(row.critical_radius)} at t={fmt_frac(row.critical_time)}"
        )
        print(
            "  "
            f"unit_skeleton={row.unit_skeleton} "
            f"protection_pressure={fmt_float(row.protection_pressure)} "
            f"owner_hist={row.owner_hist} protector_hist={row.protector_hist}"
        )
    print()


def print_small_answers(reports: list[ShapeReport]) -> None:
    print("SMALL ANSWERS FROM THIS PASS")
    print("=" * 88)

    tight = [r for r in reports if r.classification == "boundary_only"]
    near = [r for r in reports if r.classification == "positive_gap"]

    print("1. The speed list is not the fundamental shape.")
    print(
        "   The same arithmetic-looking operations can produce very different "
        "endpoint fingerprints.  The cover/endpoint incidence graph is closer "
        "to the object that decides."
    )
    print()

    print("2. Tight examples in this sample are endpoint, not interval, failures.")
    for row in tight:
        print(
            f"   {row.label}: critical/th={fmt_float(row.critical_ratio)}, "
            f"boundary_witnesses={row.boundary_witnesses}, "
            f"unit_skeleton={row.unit_skeleton}, quotient={row.quotient_layer}"
        )
    print()

    print("3. Near-disproofs are not nearly tight in the same representation.")
    for row in near:
        print(
            f"   {row.label}: gap/th={fmt_float(row.max_gap_ratio)}, "
            f"critical/th={fmt_float(row.critical_ratio)}, "
            f"unprotected={row.unprotected}, quotient={row.quotient_layer}"
        )
    print()

    high_pressure = max(reports, key=lambda r: r.protection_pressure)
    low_gap = min(near, key=lambda r: r.max_gap_ratio)
    print("4. Protection pressure and visible gap are different axes.")
    print(
        f"   Highest protection pressure here: {high_pressure.label} "
        f"({fmt_float(high_pressure.protection_pressure)})."
    )
    print(
        f"   Smallest positive visible gap here: {low_gap.label} "
        f"(gap/th={fmt_float(low_gap.max_gap_ratio)}), but it still has "
        f"{low_gap.unprotected} unprotected endpoints."
    )
    print()

    print("5. The likely fundamental object is a finite circular-arc incidence system.")
    print(
        "   LRC asks whether an open circular-arc cover can be full while its "
        "boundary incidence graph has no leaf after the right protection peel. "
        "The speed set is a generator of that system, not the system itself."
    )
    print()


def print_new_questions() -> None:
    print("NEW QUESTIONS GENERATED")
    print("=" * 88)
    questions = [
        "Can every primitive LRC endpoint-protection graph be peeled by a potential that sees only cyclic order and owner/protector degrees?",
        "Is the scalar-ramp family best understood as the complete/additive shadow, with all real difficulty living in the residue after quotienting it?",
        "Can critical radius, endpoint peel depth, and repair deficit be unified as one Morse function on a circular-arc metagraph?",
        "Which endpoint fingerprints are realizable by integer speed sets, and which are abstract circular-arc mirages?",
        "Does every tiny-gap near-disproof fail because it protects a low quotient layer while increasing unprotected endpoints in a higher quotient layer?",
        "For n=14, can the eight alpha stencils be characterized without mentioning alpha, just as the minimal leaf types of the endpoint incidence graph?",
        "Is product-sum target status a coordinate-level shadow of the same additive/multiplicative split visible in endpoint protection?",
        "What is the right notion of homology for an open circular-arc cover with protected endpoints: nerve homology, leaf-peeling core, or repair graph cycles?",
    ]
    for i, question in enumerate(questions, 1):
        print(f"{i}. {question}")
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 88)
    print(
        "The most fundamental representation seems to be neither speeds nor "
        "times alone, but the boundary of the pulled-back forbidden cover.  "
        "Speeds generate maps; intervals are unsafe fibers; endpoints are the "
        "places where the proof can fail; protection is the local mechanism "
        "that could turn a tight example into a counterexample."
    )
    print(
        "So the central shape is a finite circular-arc incidence system with a "
        "peeling dynamics.  Scalar ramps are the equality spine of this system. "
        "Torsion and product-sum gates mark where quotient layers leak.  Repair "
        "deficit measures the cost of moving inside the residue metagraph."
    )
    print(
        "A good next proof question is not 'where is the lonely time?' but "
        "'why must the endpoint-protection graph have a leaf after the scalar "
        "equality spine is quotiented?'"
    )


def main() -> None:
    print("Lonely Runner shape questions (codex-2026-05-31 S379)")
    print()
    reports = [shape_report(label, speeds) for label, speeds in sample_sets()]
    print_representation_map()
    print_shape_table(reports)
    print_endpoint_fingerprints(reports)
    print_small_answers(reports)
    print_new_questions()
    print_synthesis()


if __name__ == "__main__":
    main()
