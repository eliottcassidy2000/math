#!/usr/bin/env python3
"""Exact continuation-congruence counterexamples for the THM-822 codecs.

The replay separates three operations which must not be conflated:

* a labelled replacement-vertex insertion in the H0/H1 handoff relation;
* monotone intersection of an exact strict residual with one new safe comb;
* deletion/replacement/transport of a comb already used to form a residual.

All handoff tests use THM-806/822's half-open convention

    -u_owner < z*u_owner - 2*u_provider - 13*m*u_owner <= u_owner.

All residual endpoints and maximin candidates are exact ``Fraction`` values.
No sampled time grid or floating-point comparison occurs.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
from typing import Iterable, Sequence


P = 13
THRESHOLD = F(1, P)
EXPECTED_PAYLOAD_SHA256 = "b5c04163a21209f82231228b6fcb8c4eb5beb96fe5d254157253d21f21d14f27"
EXPECTED_CERTIFICATE_SHA256 = "b9580c2c4793555ae76ae2c5bea5a2cc0b22a1f6bc2900feff6acf1735ce63d2"

Labels = tuple[int, ...]
Speeds = tuple[int, ...]
Interval = tuple[F, F]
IntervalWord = tuple[Interval, ...]
Edge = tuple[int, int]
DecoratedEdge = tuple[int, int, int]


def circle_norm(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def handoff_relation(labels: Labels, speeds: Speeds) -> tuple[DecoratedEdge, ...]:
    """Return every live directed edge with its unique integer centre."""

    assert len(labels) == len(speeds)
    assert all(speed % P == label for label, speed in zip(labels, speeds))
    answer: list[DecoratedEdge] = []
    for i, provider in enumerate(labels):
        for j, owner in enumerate(labels):
            if i == j:
                continue
            z = provider * pow(owner, -1, P) % P
            found: list[int] = []
            for m in range(-100, 101):
                scaled = (
                    z * speeds[j]
                    - 2 * speeds[i]
                    - P * m * speeds[j]
                )
                if -speeds[j] < scaled <= speeds[j]:
                    found.append(z - P * m)
            assert len(found) <= 1
            if found:
                answer.append((provider, owner, found[0]))
    return tuple(answer)


def edge_support(relation: Sequence[DecoratedEdge]) -> frozenset[Edge]:
    return frozenset((provider, owner) for provider, owner, _ in relation)


def complete_tournament(
    labels: Labels, relation: Sequence[DecoratedEdge]
) -> frozenset[Edge]:
    """Complete silent pairs in increasing-label order."""

    live = edge_support(relation)
    answer: set[Edge] = set()
    for left, right in combinations(labels, 2):
        forward = (left, right) in live
        reverse = (right, left) in live
        assert not (forward and reverse)
        if forward:
            answer.add((left, right))
        elif reverse:
            answer.add((right, left))
        else:
            answer.add((left, right))
    return frozenset(answer)


def tournament_fingerprint(
    labels: Labels, edges: frozenset[Edge]
) -> tuple[object, ...]:
    scores = Counter(
        sum((vertex, other) in edges for other in labels if other != vertex)
        for vertex in labels
    )
    triangles = 0
    for a, b, c in combinations(labels, 3):
        triangles += (
            ((a, b) in edges and (b, c) in edges and (c, a) in edges)
            or ((b, a) in edges and (c, b) in edges and (a, c) in edges)
        )

    reach = {
        (a, b): a == b or (a, b) in edges
        for a in labels
        for b in labels
    }
    for pivot in labels:
        for a in labels:
            for b in labels:
                reach[a, b] = reach[a, b] or (
                    reach[a, pivot] and reach[pivot, b]
                )
    unseen = set(labels)
    scc: list[int] = []
    while unseen:
        seed = min(unseen)
        block = {
            vertex
            for vertex in unseen
            if reach[seed, vertex] and reach[vertex, seed]
        }
        scc.append(len(block))
        unseen -= block

    paths = [
        path
        for path in permutations(labels)
        if all((path[i], path[i + 1]) in edges for i in range(len(path) - 1))
    ]
    return (
        tuple(sorted(scores.items())),
        triangles,
        tuple(sorted(scc, reverse=True)),
        len(paths),
        paths[0],
    )


def safe_bands(speed: int) -> IntervalWord:
    return tuple(
        (
            F(P * cell + 1, P * speed),
            F(P * (cell + 1) - 1, P * speed),
        )
        for cell in range(speed)
    )


def intersect_unions(left: IntervalWord, right: IntervalWord) -> IntervalWord:
    answer: list[Interval] = []
    i = j = 0
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


def strict_residual(speeds: Iterable[int]) -> IntervalWord:
    current: IntervalWord = ((F(0), F(1)),)
    for speed in sorted(speeds):
        current = intersect_unions(current, safe_bands(speed))
        if not current:
            break
    return current


def add_comb(residual: IntervalWord, speed: int) -> IntervalWord:
    """The owner-free Markov action for a pure monotone insertion."""

    return intersect_unions(residual, safe_bands(speed))


def exact_maximin(speeds: Speeds) -> tuple[F, tuple[F, ...]]:
    """Evaluate max_t min_w ||wt|| at every exact affine candidate."""

    denominators = {2 * speed for speed in speeds}
    for left, right in combinations(speeds, 2):
        denominators.add(left + right)
        if left != right:
            denominators.add(abs(left - right))
    candidates = {
        F(numerator, denominator)
        for denominator in denominators
        if denominator
        for numerator in range(denominator + 1)
    }
    values = {
        time: min(circle_norm(speed * time) for speed in speeds)
        for time in candidates
    }
    maximum = max(values.values())
    witnesses = tuple(sorted(time for time, value in values.items() if value == maximum))
    return maximum, witnesses


def format_fraction(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def format_relation(relation: Sequence[DecoratedEdge]) -> str:
    return "{" + ",".join(
        f"{provider}->{owner}:k={centre}"
        for provider, owner, centre in relation
    ) + "}"


def format_intervals(intervals: IntervalWord) -> str:
    return "{" + ",".join(
        f"({format_fraction(lo)},{format_fraction(hi)})"
        for lo, hi in intervals
    ) + "}"


def format_fingerprint(fingerprint: tuple[object, ...]) -> str:
    score, triangles, scc, path_count, first_path = fingerprint
    return (
        f"score={score} triangles={triangles} SCC={scc} "
        f"Hamiltonian_paths={path_count} first_path={first_path}"
    )


def main() -> None:
    base_labels = (1, 2, 3, 4, 5)
    child_labels = (1, 2, 3, 4, 5, 6)
    low_speeds = (14, 15, 16, 17, 18)
    high_speeds = (27, 28, 29, 30, 31)
    incoming = 19

    low_base = handoff_relation(base_labels, low_speeds)
    high_base = handoff_relation(base_labels, high_speeds)
    expected_base = ((2, 1, 2), (3, 1, 3), (4, 2, 2))
    assert low_base == high_base == expected_base

    low_child = handoff_relation(child_labels, (*low_speeds, incoming))
    high_child = handoff_relation(child_labels, (*high_speeds, incoming))
    low_incident = tuple(edge for edge in low_child if 6 in edge[:2])
    high_incident = tuple(edge for edge in high_child if 6 in edge[:2])
    assert low_incident == ((6, 2, 3), (6, 3, 2))
    assert high_incident == ((5, 6, 3), (6, 3, 2))
    assert low_child != high_child

    base_tournament = complete_tournament(base_labels, low_base)
    low_tournament = complete_tournament(child_labels, low_child)
    high_tournament = complete_tournament(child_labels, high_child)
    base_fingerprint = tournament_fingerprint(base_labels, base_tournament)
    low_fingerprint = tournament_fingerprint(child_labels, low_tournament)
    high_fingerprint = tournament_fingerprint(child_labels, high_tournament)
    completed_edge_flips = sum(
        ((left, right) in low_tournament) != ((left, right) in high_tournament)
        for left, right in combinations(child_labels, 2)
    )
    live_support_delta = edge_support(low_child) ^ edge_support(high_child)
    assert completed_edge_flips == 1
    assert live_support_delta == frozenset({(6, 2), (5, 6)})

    source_s = tuple(range(1, 14))
    source_t = (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 20)
    deleted_s = tuple(speed for speed in source_s if speed != 13)
    deleted_t = tuple(speed for speed in source_t if speed != 13)

    residual_s = strict_residual(source_s)
    residual_t = strict_residual(source_t)
    residual_deleted_s = strict_residual(deleted_s)
    residual_deleted_t = strict_residual(deleted_t)
    assert residual_s == residual_t == ()
    assert residual_deleted_s == ()
    assert residual_deleted_t == ((F(79, 260), F(4, 13)),
                                  (F(9, 13), F(181, 260)))

    maximin_s = exact_maximin(source_s)
    maximin_t = exact_maximin(source_t)
    maximin_deleted_s = exact_maximin(deleted_s)
    maximin_deleted_t = exact_maximin(deleted_t)
    assert maximin_s[0] == F(1, 14)
    assert maximin_t[0] == F(2, 27)
    assert maximin_deleted_s[0] == F(1, 13)
    assert maximin_deleted_t[0] == F(2, 23)
    assert maximin_s[0] != maximin_t[0]

    # Equality of exact residuals is a congruence for every pure insertion:
    # intersection with a common set preserves equality.  Speed 17 is a
    # concrete replay row; the identity itself is operation-independent.
    addition_probe = 17
    added_s = add_comb(residual_s, addition_probe)
    added_t = add_comb(residual_t, addition_probe)
    assert added_s == added_t == ()
    assert add_comb(residual_s, addition_probe) == strict_residual(
        (*source_s, addition_probe)
    )
    assert add_comb(residual_t, addition_probe) == strict_residual(
        (*source_t, addition_probe)
    )

    truth_table = (
        ("H0", "insert labelled (6,19)", True, False, "FAIL"),
        ("H1", "insert labelled (6,19)", True, False, "FAIL"),
        ("exact residual E", "add common safe comb", True, True, "PASS_BY_INTERSECTION"),
        ("exact residual E", "delete common comb 13", True, False, "FAIL"),
        ("exact residual E", "evaluate exact M", True, False, "FAIL"),
        ("E + labelled tooth bank", "delete/recompute", False, True, "SUFFICIENT_CARRIER"),
    )

    payload = (
        ("base", low_base),
        ("low_child", low_child),
        ("high_child", high_child),
        ("S", source_s, residual_s, maximin_s),
        ("T", source_t, residual_t, maximin_t),
        ("S_delete13", deleted_s, residual_deleted_s, maximin_deleted_s),
        ("T_delete13", deleted_t, residual_deleted_t, maximin_deleted_t),
    )
    payload_sha = sha256(repr(payload).encode()).hexdigest()
    certificate = (
        payload,
        truth_table,
        base_fingerprint,
        low_fingerprint,
        high_fingerprint,
        tuple(sorted(live_support_delta)),
        completed_edge_flips,
    )
    certificate_sha = sha256(repr(certificate).encode()).hexdigest()
    if EXPECTED_PAYLOAD_SHA256 != "PENDING":
        assert payload_sha == EXPECTED_PAYLOAD_SHA256
    if EXPECTED_CERTIFICATE_SHA256 != "PENDING":
        assert certificate_sha == EXPECTED_CERTIFICATE_SHA256

    source_sha = sha256(Path(__file__).read_bytes()).hexdigest()

    print("LRC13 HAMMING-FIVE CONTINUATION CONGRUENCE - EXACT REPLAY")
    print("scope=THM-822 bounded liar plus exact common-deletion witness")
    print("arithmetic=integer_half_open_handoff+Fraction_residuals+Fraction_maximin")
    print("claim=pure_addition_only; deletion_replacement_transport_require_tooth_provenance")
    print()
    print("H0_H1_INSERTION_COUNTEREXAMPLE")
    print(f"base_labels={base_labels}")
    print(f"low_speeds={low_speeds} high_speeds={high_speeds}")
    print(f"shared_base_relation={format_relation(low_base)}")
    print("half_open_rule=-u_owner<z*u_owner-2*u_provider-13*m*u_owner<=u_owner")
    print(f"incoming=(label=6,speed={incoming})")
    print(f"low_incident={format_relation(low_incident)}")
    print(f"high_incident={format_relation(high_incident)}")
    print("H0_base_equal=true H1_base_equal=true child_equal=false")
    print()
    print("STRICT_RESIDUAL_COMMON_DELETION_COUNTEREXAMPLE")
    print(f"S={source_s} E_S={format_intervals(residual_s)} M_S={maximin_s[0]}")
    print(f"M_S_witnesses={maximin_s[1]}")
    print(f"T={source_t} E_T={format_intervals(residual_t)} M_T={maximin_t[0]}")
    print(f"M_T_witnesses={maximin_t[1]}")
    print(f"delete=13 E_Sminus={format_intervals(residual_deleted_s)} M_Sminus={maximin_deleted_s[0]}")
    print(f"M_Sminus_witnesses={maximin_deleted_s[1]}")
    print(f"delete=13 E_Tminus={format_intervals(residual_deleted_t)} M_Tminus={maximin_deleted_t[0]}")
    print(f"M_Tminus_witnesses={maximin_deleted_t[1]}")
    print("base_endpoint_words_equal=true deleted_endpoint_words_equal=false")
    print()
    print("KERNEL_CONGRUENCE_TRUTH_TABLE")
    print("columns=observation|operation|base_equal|outputs_equal|verdict")
    for row in truth_table:
        print("|".join(map(str, row)))
    print("identity=E(S_union_{u})=E(S)_intersection_G_u")
    print("endpoint_word_preserves=one_step_exact_residual_geometry_and_threshold_nonemptiness")
    print("endpoint_word_destroys=exact_M+deleted_comb_identity+owner_ancestry+transport_action")
    print()
    print("TOURNAMENT_ANALYSIS")
    print("vertices=replacement_labels")
    print("pair_observable=directed_half_open_handoff")
    print("switch=low_lift_branch_vs_high_lift_branch_after_common_insertion")
    print(f"tie_Hamiltonian_path={child_labels} rule=increasing_label_on_silent_pairs")
    print(f"base {format_fingerprint(base_fingerprint)}")
    print(f"low_child {format_fingerprint(low_fingerprint)}")
    print(f"high_child {format_fingerprint(high_fingerprint)}")
    print(f"live_support_symmetric_difference={tuple(sorted(live_support_delta))}")
    print(f"completed_tournament_edge_flips={completed_edge_flips}")
    print("tournament_preserves=completed_pairwise_orientation_for_this_gauge")
    print("tournament_destroys=live_vs_tie_provenance+integer_centres+residual_geometry+LRC_threshold_truth")
    print()
    print(f"payload_sha256={payload_sha}")
    print(f"certificate_sha256={certificate_sha}")
    print(f"source_sha256={source_sha}")
    print("PASS: both continuation counterexamples and the pure-addition identity are exact")


if __name__ == "__main__":
    main()
