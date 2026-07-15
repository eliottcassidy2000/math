#!/usr/bin/env python3
"""Primitive scale-one Hamming-six AP-pin reduction.

This exact scout has three deliberately separated jobs.

1. Independently audit the algebraic part and terminal certificates of the
   nonprimitive contraction in the THM-815 C.1 addendum.
2. Classify the 923 primitive-core missing-label rows by antipodal AP-grid
   pins, retained-core divisor obligations, and the height-zero Kakeya ledger.
3. Use a height-independent one-sided germ handoff lemma to reduce the 20
   rows made from three full antipodal pairs.  Six rows close locally; eight
   force one of ten fixed-coordinate slices, all of which close by an exact
   row-wise THM-815 recursion.  Six rows remain open.

The exact cover carrier is the residual strict-safe interval union together
with the labelled remaining comb bank.  The antipodal count is a useful
stratification, not a quotient of that carrier.
"""

from __future__ import annotations

import hashlib
import itertools
from collections import Counter
from fractions import Fraction as F
from math import gcd
from pathlib import Path
from struct import pack


P13 = 13
BASE = tuple(range(1, P13))
ODD_ROW = (1, 3, 5, 7, 9, 11)

CONTRACTION_SOURCE = Path(
    "04-computation/lrc13_hamming_six_nonprimitive_contraction_scout_codex_S11.cpp"
)
CONTRACTION_OUTPUT = Path(
    "05-knowledge/results/lrc13_hamming_six_nonprimitive_contraction_scout_codex_S11.out"
)
CONTRACTION_SOURCE_SHA256 = (
    "ee57510a4796e23da1408b383af1146478f067a5e0e98d2ad52220e55e2e8bf1"
)
CONTRACTION_OUTPUT_SHA256 = (
    "aa87c107ff7e5fdc6cb6a3803ecc0751ac41c267dd59540e4bd6dd764f7769ed"
)
SLICE_TRACE_SHA256 = "f8c8465455fdbf1f21aec0438a8894c11f56accf53579dc84658a703aa5c227e"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_bands(speed: int) -> tuple[tuple[F, F], ...]:
    return tuple(
        (F(13 * k + 1, 13 * speed), F(13 * (k + 1) - 1, 13 * speed))
        for k in range(speed)
    )


def danger_bands(speed: int) -> tuple[tuple[F, F], ...]:
    answer = []
    for k in range(speed + 1):
        lo = max(F(0), F(13 * k - 1, 13 * speed))
        hi = min(F(1), F(13 * k + 1, 13 * speed))
        if lo < hi:
            answer.append((lo, hi))
    return tuple(answer)


SAFE_CACHE: dict[int, tuple[tuple[F, F], ...]] = {}
DANGER_CACHE: dict[int, tuple[tuple[F, F], ...]] = {}


def safe_for(speed: int) -> tuple[tuple[F, F], ...]:
    if speed not in SAFE_CACHE:
        SAFE_CACHE[speed] = safe_bands(speed)
    return SAFE_CACHE[speed]


def danger_for(speed: int) -> tuple[tuple[F, F], ...]:
    if speed not in DANGER_CACHE:
        DANGER_CACHE[speed] = danger_bands(speed)
    return DANGER_CACHE[speed]


def meet(
    left: tuple[tuple[F, F], ...], right: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    answer = []
    i = 0
    j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            answer.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(answer)


def merged(intervals: list[tuple[F, F]]) -> tuple[tuple[F, F], ...]:
    answer: list[tuple[F, F]] = []
    for lo, hi in sorted(intervals):
        if answer and lo <= answer[-1][1]:
            answer[-1] = (answer[-1][0], max(answer[-1][1], hi))
        else:
            answer.append((lo, hi))
    return tuple(answer)


def measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((hi - lo for lo, hi in intervals), F(0))


def safe_set(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    answer: tuple[tuple[F, F], ...] = ((F(0), F(1)),)
    for speed in speeds:
        answer = meet(answer, safe_for(speed))
    return answer


def safe_by_closed_danger(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    teeth: list[tuple[F, F]] = []
    for speed in speeds:
        teeth.extend(danger_for(speed))
    covered = merged(teeth)
    answer = []
    cursor = F(0)
    for lo, hi in covered:
        if cursor < lo:
            answer.append((cursor, lo))
        cursor = max(cursor, hi)
    if cursor < 1:
        answer.append((cursor, F(1)))
    return tuple(answer)


def longest(components: tuple[tuple[F, F], ...]) -> F:
    if not components:
        raise RuntimeError("longest component requested from an empty set")
    return max(hi - lo for lo, hi in components)


def discrepancy_cap(remaining: int, component_length: F) -> int:
    return (
        22 * remaining * component_length.denominator
        // (13 * (13 - 2 * remaining) * component_length.numerator)
    )


def retained(missing: tuple[int, ...]) -> tuple[int, ...]:
    missing_set = set(missing)
    return tuple(label for label in BASE if label not in missing_set)


def audit_nonprimitive_contraction() -> None:
    assert sha256(CONTRACTION_SOURCE) == CONTRACTION_SOURCE_SHA256
    assert sha256(CONTRACTION_OUTPUT) == CONTRACTION_OUTPUT_SHA256

    candidates = []
    for missing in itertools.combinations(BASE, 6):
        core = retained(missing)
        core_gcd = 0
        for label in core:
            core_gcd = gcd(core_gcd, label)
        if core_gcd > 1:
            candidates.append((missing, core, core_gcd))
    assert candidates == [(ODD_ROW, (2, 4, 6, 8, 10, 12), 2)]

    # The retained label 2 makes the total gcd divide two.  Exact parity
    # enumeration confirms that all six odd-label heights must be odd.
    nonprimitive_patterns = []
    for parity in itertools.product((0, 1), repeat=6):
        heights = tuple(2 if bit == 0 else 1 for bit in parity)
        packet = list((2, 4, 6, 8, 10, 12))
        packet.extend(r + 13 * h for r, h in zip(ODD_ROW, heights))
        packet_gcd = 0
        for speed in packet:
            packet_gcd = gcd(packet_gcd, speed)
        if packet_gcd > 1:
            nonprimitive_patterns.append(parity)
    assert nonprimitive_patterns == [(1, 1, 1, 1, 1, 1)]

    # Independently reconstruct the two terminal safe sets and their dead
    # next-label arithmetic from the frozen C++ certificate.
    leaves = (
        (((7, 33), (10, 36), (9, 48), (8, 60), (11, 63)), 12, 64, 58, F(17, 3120)),
        (((7, 33), (10, 36), (9, 48), (8, 60), (12, 64)), 11, 76, 54, F(47, 8580)),
    )
    for path, unused, least_next, component_count, expected_longest in leaves:
        packet = (1, 2, 3, 4, 5, 6) + tuple(speed for _, speed in path)
        open_safe = safe_set(packet)
        danger_complement = safe_by_closed_danger(packet)
        assert open_safe == danger_complement
        assert len(open_safe) == component_count
        assert longest(open_safe) == expected_longest
        assert discrepancy_cap(1, expected_longest) == 28
        candidate = unused + 13
        while candidate <= path[-1][1]:
            candidate += 13
        assert candidate == least_next and candidate > 28


def antipodal_count(missing: tuple[int, ...]) -> int:
    missing_set = set(missing)
    return sum(r in missing_set and 13 - r in missing_set for r in range(1, 7))


def ap_grid_interior(core: tuple[int, ...]) -> tuple[int, ...]:
    answer = []
    for a in range(1, 13):
        minimum = min(min((a * p) % 13, 13 - (a * p) % 13) for p in core)
        if minimum > 1:
            answer.append(a)
    return tuple(answer)


def expected_ap_grid(missing: tuple[int, ...]) -> tuple[int, ...]:
    missing_set = set(missing)
    answer = set()
    for r in range(1, 7):
        if r in missing_set and 13 - r in missing_set:
            a = pow(r, -1, 13)
            answer.update((a, 13 - a))
    return tuple(sorted(answer))


def clean_moduli(core: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(q for q in range(2, 13) if all(p % q for p in core))


def kakeya_row(
    missing: tuple[int, ...], core_components: tuple[tuple[F, F], ...]
) -> tuple[int, int, bool]:
    zero_debt = 0
    unique_full = 0
    every_zero_is_degree_one = True
    for component in core_components:
        component_length = component[1] - component[0]
        pieces: list[tuple[F, F]] = []
        sum_measure = F(0)
        touching = []
        full = []
        for label in missing:
            intersection = meet((component,), danger_for(label))
            intersection_measure = measure(intersection)
            if intersection_measure:
                touching.append(label)
                pieces.extend(intersection)
            if intersection_measure == component_length:
                full.append(label)
            sum_measure += intersection_measure
        union_measure = measure(merged(pieces))
        assert union_measure == component_length  # [12] is tight at 1/13.
        if sum_measure == union_measure:
            zero_debt += 1
            if len(full) == 1:
                unique_full += 1
            if len(touching) != 1 or full != touching:
                every_zero_is_degree_one = False
    return zero_debt, unique_full, every_zero_is_degree_one


def primitive_row_atlas() -> dict[str, object]:
    f_histogram: Counter[int] = Counter()
    ap_histogram: Counter[int] = Counter()
    q_size_histogram: Counter[int] = Counter()
    q_signatures: Counter[tuple[int, ...]] = Counter()
    zero_histogram: Counter[int] = Counter()
    unique_histogram: Counter[int] = Counter()
    no_unique_rows = []
    no_zero_rows = []
    f0_zero_counts = []
    root_lengths = set()
    root_caps = []
    rows = {}

    for missing in itertools.combinations(BASE, 6):
        if missing == ODD_ROW:
            continue
        core = retained(missing)
        core_gcd = 0
        for label in core:
            core_gcd = gcd(core_gcd, label)
        assert core_gcd == 1

        f_value = antipodal_count(missing)
        ap_points = ap_grid_interior(core)
        assert ap_points == expected_ap_grid(missing)
        assert len(ap_points) == 2 * f_value

        components = safe_set(core)
        assert components == safe_by_closed_danger(core)
        root_length = longest(components)
        root_cap = discrepancy_cap(6, root_length)
        zero_debt, unique_full, degree_one = kakeya_row(missing, components)
        q_values = clean_moduli(core)

        f_histogram[f_value] += 1
        ap_histogram[len(ap_points)] += 1
        q_size_histogram[len(q_values)] += 1
        q_signatures[q_values] += 1
        zero_histogram[zero_debt] += 1
        unique_histogram[unique_full] += 1
        root_lengths.add(root_length)
        root_caps.append(root_cap)
        if not zero_debt:
            no_zero_rows.append(missing)
        if not unique_full:
            no_unique_rows.append(missing)
        if f_value == 0:
            assert degree_one
            f0_zero_counts.append(zero_debt)

        rows[missing] = {
            "core": core,
            "f": f_value,
            "q": q_values,
            "root_length": root_length,
            "root_cap": root_cap,
            "zero": zero_debt,
            "unique": unique_full,
        }

    assert len(rows) == 923
    assert f_histogram == {0: 63, 1: 480, 2: 360, 3: 20}
    assert ap_histogram == {0: 63, 2: 480, 4: 360, 6: 20}
    assert q_size_histogram == {0: 1, 1: 26, 2: 137, 3: 274, 4: 285, 5: 160, 6: 40}
    assert len(q_signatures) == 182
    assert sum(count for pins, count in unique_histogram.items() if pins > 0) == 916
    assert sum(count for pins, count in zero_histogram.items() if pins > 0) == 919
    assert min(f0_zero_counts) == 4 and max(f0_zero_counts) == 22

    expected_no_unique = [
        (2, 3, 6, 7, 11, 12),
        (2, 6, 7, 8, 11, 12),
        (3, 4, 5, 9, 10, 11),
        (4, 5, 7, 8, 9, 10),
        (4, 5, 8, 9, 11, 12),
        (5, 7, 8, 9, 10, 11),
        (5, 7, 8, 9, 11, 12),
    ]
    expected_no_zero = [
        (4, 5, 7, 8, 9, 10),
        (4, 5, 8, 9, 11, 12),
        (5, 7, 8, 9, 10, 11),
        (5, 7, 8, 9, 11, 12),
    ]
    assert no_unique_rows == expected_no_unique
    assert no_zero_rows == expected_no_zero

    # Clean-modulus data is a necessary owner ledger, not a safe quotient.
    liar_a = (1, 3, 4, 7, 8, 9)
    liar_b = (1, 2, 3, 7, 8, 9)
    assert rows[liar_a]["q"] == rows[liar_b]["q"] == (7, 8, 9)
    assert rows[liar_a]["root_length"] == F(3, 65)
    assert rows[liar_b]["root_length"] == F(3, 52)
    assert rows[liar_a]["root_cap"] == 220
    assert rows[liar_b]["root_cap"] == 176

    # At t=1/q, q<=12, a replacement is dangerous iff it is divisible by q.
    for q in range(2, 13):
        for speed in range(1, 1001):
            residue = speed % q
            distance = F(min(residue, q - residue), q)
            assert (distance <= F(1, 13)) == (residue == 0)

    return {
        "rows": rows,
        "f_histogram": f_histogram,
        "ap_histogram": ap_histogram,
        "q_size_histogram": q_size_histogram,
        "q_signature_count": len(q_signatures),
        "zero_histogram": zero_histogram,
        "unique_histogram": unique_histogram,
        "no_unique": tuple(no_unique_rows),
        "no_zero": tuple(no_zero_rows),
        "root_length_count": len(root_lengths),
        "root_cap_range": (min(root_caps), max(root_caps)),
    }


def signed_residue(value: int) -> int:
    residue = value % 13
    return residue if residue <= 6 else residue - 13


def left_start_coefficient(z: int) -> int:
    """First leftward danger start from phase z/13, excluding z=0."""
    assert -6 <= z <= 6 and z != 0
    return z - 1 if z > 0 else 12 + z


def canonical_cycles(
    graph: dict[int, dict[int, F]]
) -> tuple[tuple[tuple[int, ...], F], ...]:
    seen = set()
    answer = []

    def visit(start: int, path: tuple[int, ...]) -> None:
        for nxt in graph[path[-1]]:
            if nxt == start and len(path) >= 2:
                canonical = min(path[i:] + path[:i] for i in range(len(path)))
                if canonical not in seen:
                    seen.add(canonical)
                    product = F(1)
                    for left, right in zip(canonical, canonical[1:] + canonical[:1]):
                        product *= graph[left][right]
                    answer.append((canonical, product))
            elif nxt not in path:
                visit(start, path + (nxt,))

    for vertex in sorted(graph):
        visit(vertex, (vertex,))
    return tuple(sorted(answer))


LOCAL_CLOSED_ROWS = (
    (1, 2, 3, 10, 11, 12),
    (1, 3, 4, 9, 10, 12),
    (1, 3, 5, 8, 10, 12),
    (1, 4, 6, 7, 9, 12),
    (2, 3, 4, 9, 10, 11),
    (2, 3, 5, 8, 10, 11),
)

SLICE_RESTRICTIONS = {
    (1, 2, 5, 8, 11, 12): ((5, 18),),
    (1, 4, 5, 8, 9, 12): ((5, 18),),
    (1, 5, 6, 7, 8, 12): ((5, 18),),
    (2, 3, 6, 7, 10, 11): ((6, 19),),
    (2, 4, 6, 7, 9, 11): ((6, 19),),
    (2, 5, 6, 7, 8, 11): ((5, 18), (6, 19)),
    (3, 4, 6, 7, 9, 10): ((6, 19),),
    (4, 5, 6, 7, 8, 9): ((5, 18), (6, 19)),
}

OPEN_F3_ROWS = (
    (1, 2, 4, 9, 11, 12),
    (1, 2, 6, 7, 11, 12),
    (1, 3, 6, 7, 10, 12),
    (2, 4, 5, 8, 9, 11),
    (3, 4, 5, 8, 9, 10),
    (3, 5, 6, 7, 8, 10),
)


def f3_handoff_classification() -> dict[str, object]:
    min_product_histogram: Counter[F] = Counter()
    closed = []
    restricted = {}
    unresolved = []
    direct_rows = 0

    for pair_choices in itertools.combinations(range(1, 7), 3):
        missing = tuple(sorted(pair_choices + tuple(13 - r for r in pair_choices)))
        core = retained(missing)
        graph: dict[int, dict[int, F]] = {}
        exit_caps: dict[int, F] = {}
        exit_options: list[tuple[int, int]] = []

        for owner in missing:
            a = pow(owner, -1, 13)
            graph[owner] = {}
            for provider in missing:
                if provider == owner:
                    continue
                z = signed_residue(a * provider)
                graph[owner][provider] = F(left_start_coefficient(z), 2)
            core_extent = min(
                F(left_start_coefficient(signed_residue(a * label)), label)
                for label in core
            )
            exit_caps[owner] = F(2, 1) / core_extent
            speed = owner + 13
            while speed <= exit_caps[owner]:
                exit_options.append((owner, speed))
                speed += 13

        cycles = canonical_cycles(graph)
        assert len(cycles) == 409
        minimum_product = min(product for _, product in cycles)
        min_product_histogram[minimum_product] += 1
        nonexpanding = [(cycle, product) for cycle, product in cycles if product <= 1]

        noncontracting = minimum_product > 1
        if minimum_product == 1:
            # Equality around a product-one cycle would be necessary.  Here
            # every such edge has weight one, forcing distinct-residue speeds
            # to be equal, which is impossible.
            assert len(nonexpanding) == 2
            assert all(
                all(graph[left][right] == 1 for left, right in zip(cycle, cycle[1:] + cycle[:1]))
                for cycle, _ in nonexpanding
            )
            noncontracting = True
        if minimum_product < 1:
            assert minimum_product == F(1, 16) and len(nonexpanding) == 7

        if noncontracting:
            if exit_options:
                restricted[missing] = tuple(sorted(exit_options))
            else:
                closed.append(missing)
        else:
            unresolved.append(missing)

    assert min_product_histogram == {F(3, 2): 12, F(1): 2, F(1, 16): 6}
    assert tuple(closed) == LOCAL_CLOSED_ROWS
    assert restricted == SLICE_RESTRICTIONS
    assert tuple(unresolved) == OPEN_F3_ROWS

    # Bounded direct interval evidence for the local theorem: all height-1/2
    # rows that violate its necessary exit disjunction remain strict-safe.
    for missing in LOCAL_CLOSED_ROWS:
        core = retained(missing)
        for heights in itertools.product((1, 2), repeat=6):
            packet = core + tuple(r + 13 * h for r, h in zip(missing, heights))
            assert safe_set(packet)
            direct_rows += 1
    for missing, options in SLICE_RESTRICTIONS.items():
        core = retained(missing)
        for heights in itertools.product((1, 2), repeat=6):
            lifts = tuple(r + 13 * h for r, h in zip(missing, heights))
            if any((label, speed) in options for label, speed in zip(missing, lifts)):
                continue
            assert safe_set(core + lifts)
            direct_rows += 1
    assert direct_rows == 608

    return {
        "min_product_histogram": min_product_histogram,
        "closed": tuple(closed),
        "restricted": restricted,
        "unresolved": tuple(unresolved),
        "direct_rows": direct_rows,
    }


SLICE_EXPECTED = {
    ((1, 2, 5, 8, 11, 12), (5, 18)): ((1, 28, 301, 119, 11, 0), (28, 301, 119, 11, 0), F(19, 585)),
    ((1, 4, 5, 8, 9, 12), (5, 18)): ((1, 22, 196, 119, 0, 0), (22, 196, 119, 0, 0), F(101, 2574)),
    ((1, 5, 6, 7, 8, 12), (5, 18)): ((1, 21, 194, 95, 1, 0), (21, 194, 95, 1, 0), F(5, 117)),
    ((2, 3, 6, 7, 10, 11), (6, 19)): ((1, 20, 170, 130, 5, 0), (20, 170, 130, 5, 0), F(21, 494)),
    ((2, 4, 6, 7, 9, 11), (6, 19)): ((1, 20, 166, 66, 0, 0), (20, 166, 66, 0, 0), F(11, 247)),
    ((3, 4, 6, 7, 9, 10), (6, 19)): ((1, 30, 334, 166, 4, 0), (30, 334, 166, 4, 0), F(87, 2717)),
    ((4, 5, 6, 7, 8, 9), (5, 18)): ((1, 20, 171, 196, 14, 0), (20, 171, 196, 14, 0), F(7, 156)),
    ((4, 5, 6, 7, 8, 9), (6, 19)): ((1, 20, 172, 80, 3, 0), (20, 172, 80, 3, 0), F(11, 247)),
    ((2, 5, 6, 7, 8, 11), (5, 18)): ((1, 21, 192, 115, 3, 0), (21, 192, 115, 3, 0), F(19, 468)),
    ((2, 5, 6, 7, 8, 11), (6, 19)): ((1, 26, 280, 150, 8, 0), (26, 280, 150, 8, 0), F(4, 117)),
}


def fixed_slice_replay() -> dict[str, object]:
    aggregate_nodes = [0] * 6
    aggregate_edges = [0] * 5
    aggregate_covering = 0
    independent_component_checks = 0
    digest = hashlib.sha256()
    rows = []

    def digest_int(value: int) -> None:
        digest.update(pack(">q", value))

    for (missing, fixed), (expected_nodes, expected_edges, expected_root) in SLICE_EXPECTED.items():
        core = retained(missing)
        labels = tuple(label for label in missing if label != fixed[0])
        fixed_core = core + (fixed[1],)
        root = safe_set(fixed_core)
        assert root == safe_by_closed_danger(fixed_core)
        assert longest(root) == expected_root

        nodes = [0] * 6
        edges = [0] * 5
        covering = 0
        path: list[tuple[int, int]] = []

        def recurse(
            components: tuple[tuple[F, F], ...],
            used: tuple[int, ...],
            previous: int,
        ) -> None:
            nonlocal covering, independent_component_checks
            depth = len(used)
            nodes[depth] += 1
            current_speeds = fixed_core + tuple(speed for _, speed in path)
            assert components == safe_by_closed_danger(current_speeds)
            independent_component_checks += 1

            digest_int(depth)
            digest_int(len(path))
            for label, speed in path:
                digest_int(label)
                digest_int(speed)
            digest_int(len(components))
            for lo, hi in components:
                for value in (lo.numerator, lo.denominator, hi.numerator, hi.denominator):
                    digest_int(value)

            if not components:
                covering += 1
                digest_int(-2)
                return
            if depth == 5:
                digest_int(-1)
                return

            remaining = 5 - depth
            cap = discrepancy_cap(remaining, longest(components))
            digest_int(cap)
            for label in labels:
                if label in used:
                    continue
                speed = label + 13
                if speed <= previous:
                    speed += 13 * ((previous - speed) // 13 + 1)
                while speed <= cap:
                    edges[depth] += 1
                    path.append((label, speed))
                    recurse(meet(components, safe_for(speed)), used + (label,), speed)
                    path.pop()
                    speed += 13

        recurse(root, tuple(), 0)
        assert tuple(nodes) == expected_nodes
        assert tuple(edges) == expected_edges
        assert covering == 0 and nodes[5] == 0
        for depth in range(6):
            aggregate_nodes[depth] += nodes[depth]
        for depth in range(5):
            aggregate_edges[depth] += edges[depth]
        aggregate_covering += covering
        rows.append((missing, fixed, tuple(nodes), tuple(edges), expected_root))

    assert aggregate_nodes == [10, 228, 2176, 1236, 49, 0]
    assert aggregate_edges == [228, 2176, 1236, 49, 0]
    assert sum(aggregate_nodes) == 3699
    assert independent_component_checks == 3699
    assert aggregate_covering == 0
    assert digest.hexdigest() == SLICE_TRACE_SHA256

    return {
        "rows": tuple(rows),
        "nodes": tuple(aggregate_nodes),
        "edges": tuple(aggregate_edges),
        "covering": aggregate_covering,
        "component_checks": independent_component_checks,
        "trace": digest.hexdigest(),
    }


def row_text(row: tuple[int, ...]) -> str:
    return "{" + ",".join(map(str, row)) + "}"


def main() -> None:
    audit_nonprimitive_contraction()
    atlas = primitive_row_atlas()
    handoff = f3_handoff_classification()
    slices = fixed_slice_replay()

    print("=" * 96)
    print("LRC13 PRIMITIVE-CORE SCALE-ONE HAMMING-SIX AP-PIN REDUCTION")
    print("=" * 96)
    print("NONPRIMITIVE CONTRACTION AUDIT")
    print("retained_core_gcd>1 rows=1/924 missing={1,3,5,7,9,11} core={2,4,6,8,10,12}")
    print("height_parity criterion=all six odd; gcd exactly 2")
    print("divide_by_two={1,...,6} union {6+i+13k_i}")
    print("independent deepest-leaf endpoint reconstructions=2/2 cap=28 outgoing=0")
    print("compiled canonical C++ replay checked separately byte-for-byte against stored output")
    print()

    print("PRIMITIVE-CORE 923-ROW ATLAS")
    print("scope warning: primitive packet != primitive retained core")
    print("the exceptional odd-label row with mixed height parity lies outside this atlas")
    print(f"antipodal_full_pair_f={dict(sorted(atlas['f_histogram'].items()))}")
    print(f"interior_AP_gridpoints={dict(sorted(atlas['ap_histogram'].items()))}")
    print("f interpretation: 2f points a/13 lie in E_P; their one-sided germ owners are +/-a^(-1)")
    print(f"root_longest_values={atlas['root_length_count']} root_first_cap_range={atlas['root_cap_range']}")
    print(f"clean_modulus_Q_size={dict(sorted(atlas['q_size_histogram'].items()))}")
    print(f"clean_modulus_signatures={atlas['q_signature_count']}")
    print("clean-modulus owner law: q in Q(P) => some lift is divisible by q")
    print("Q-liar: Q={7,8,9} but L=3/65,cap=220 versus L=3/52,cap=176")
    print("Kakeya base ledger: rows with zero-debt component=919; with unique full owner=916")
    print("f=0 primitive transversals=63; every zero-debt component is degree-one, count range=4..22")
    print("no_unique_full_owner_rows=" + ";".join(row_text(row) for row in atlas["no_unique"]))
    print("no_zero_debt_rows=" + ";".join(row_text(row) for row in atlas["no_zero"]))
    print()

    print("F=3 HEIGHT-INDEPENDENT GERM HANDOFF")
    print("left coefficient c(z)=z-1 for z>0, 12+z for z<0")
    print("owner r either exits: u_r<=2/min_(p in P)c(a*p)/p, or hands to s:")
    print("u_s >= c(a*s)u_r/2, where a=r^(-1) mod 13")
    print(f"minimum_cycle_product_histogram={dict(sorted(handoff['min_product_histogram'].items()))}")
    print("bounded direct height-{1,2} local checks=" + str(handoff["direct_rows"]))
    print("boundary_closed_6=" + ";".join(row_text(row) for row in handoff["closed"]))
    print("forced_slice_rows_8:")
    for row, options in handoff["restricted"].items():
        print("  " + row_text(row) + " => " + " or ".join(f"u_{r}={u}" for r, u in options))
    print("cycle_open_6=" + ";".join(row_text(row) for row in handoff["unresolved"]))
    print()

    print("TEN FIXED-COORDINATE SLICE REPLAYS (8 ROWS)")
    for missing, fixed, nodes, edges, root_length in slices["rows"]:
        print(
            f"  R={row_text(missing)} fixed=u_{fixed[0]}={fixed[1]} "
            f"rootL={root_length} nodes={','.join(map(str,nodes))} "
            f"edges={','.join(map(str,edges))} covering=0"
        )
    print(f"aggregate_depth_nodes={','.join(map(str,slices['nodes']))} total=3699")
    print(f"aggregate_edges={','.join(map(str,slices['edges']))}")
    print(f"closed-danger component crosschecks={slices['component_checks']}")
    print(f"slice_trace_sha256={slices['trace']}")
    print("slice verdict: all 10 slices empty; therefore all 8 forced-slice rows are loose")
    print("FINAL REDUCTION: 14 of 20 primitive-core f=3 rows closed; row frontier 923 -> 909")
    print("SCOPE: the other 903 f<=2 rows and the six listed product-1/16 f=3 rows remain open")
    print("the exceptional odd-label row with primitive mixed-parity packets is separate and open")
    print()

    print("REFLECTION / ASSUMPTION CHALLENGE / TOURNAMENT ANALYSIS")
    print("f is invariant under r->13-r, every unit dilation mod 13, and R-complement")
    print("f preserves antipodal AP-cusp count and forced germ ownership only")
    print("f destroys pair identities, signs, core endpoints/L, divisor obligations, and Kakeya overlap")
    print("these actions need not preserve the primitive-core domain or exact LRC geometry")
    print("they are not common integer dilations of the normalized scale-one packet")
    print("pairwise observable=provider start c(a*s)/(13u_s) versus owner reach 2/(13u_r)")
    print("switch=edge r->s when the exact handoff inequality holds; gauge=oriented AP left germ")
    print("tie Hamiltonian path=increasing (speed,label), used only to enumerate remaining lifts")
    print("the handoff carrier is a weighted digraph, not a tournament; cycle products are proof-bearing")
    print("faithful global carrier=residual endpoint union x labelled remaining comb bank")
    print()

    print("PROVENANCE")
    print(f"contraction_source_sha256={sha256(CONTRACTION_SOURCE)}")
    print(f"contraction_output_sha256={sha256(CONTRACTION_OUTPUT)}")
    print(f"script_sha256={sha256(Path(__file__))}")


if __name__ == "__main__":
    main()
