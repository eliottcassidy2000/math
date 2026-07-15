#!/usr/bin/env python3
"""Exact replay for THM-806, scale-one Hamming-three collar closure.

For distinct labels ``r,s,t`` and positive lift heights, put

    A = ([12] - {r,s,t}) union {r+13*i,s+13*j,t+13*k}.

The replay has three logically separate layers.

1.  At the inverse thirteenth belonging to a missing owner, the nine-speed
    core is strictly safe on a universal left collar of length ``1/156``.
    If the owner's replacement exceeds 24, its aligned danger tooth ends
    inside that collar.  Tightness would force a cross-colour left handoff.
    The exact half-open phase rule is ``(-1,1]`` after multiplying phases by
    13.  A finite directed graph with one incoming handoff at every owner has
    a 2- or 3-cycle; elementary residue/ratio bands rule both out.  Hence one
    replacement is at most 24.

2.  Choose the least such replacement ``x`` and order the other two ``v<w``.
    The settled 10-speed and 11-speed Lonely Runner bounds, combined with the
    periodic-danger measure lemma, give ``v<=381`` and ``w<=12*v``.

3.  Every row of that finite box is checked with exact ``Fraction`` safe
    components.  If ``(lo,hi)`` is a component for the eleven-speed core,
    with centre ``c`` and half-width ``h``, the last speed ``w`` covers it iff

        ||w*c|| + w*h <= 1/13.

    The hot loop evaluates this identity after clearing denominators.  There
    is no floating point and no sampled-circle step.

Tournament Analysis / assumption challenge
------------------------------------------
The theorem-bearing local vertices are owner-collar exit obligations and
replacement colours; the exact pair observable is half-open handoff
eligibility.  The finite closure then changes carrier to strict-safe
components versus last-speed teeth.  Runner vertices, residue vertices, and a
bare tournament all forget either the speed ratio, the boundary orientation,
or component width.  For telemetry only, the residue tournament declares the
subunit ``provider/owner=2`` handoff as a live edge, uses increasing residue on
silent pairs, and reverses the live arrows as its switch.  Its fingerprints
are printed below; the proof uses the decorated incidence relations.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from struct import pack
from typing import Iterable, Sequence


DELTA = F(1, 13)
BASE = tuple(range(1, 13))


def fmt(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def intersect_interval_unions(
    left: Sequence[tuple[F, F]], right: Sequence[tuple[F, F]]
) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return out


def strict_safe_components(speeds: Iterable[int]) -> tuple[tuple[F, F], ...]:
    """Exact components of ``{t: min_v ||v*t||>1/13}`` in ``(0,1)``."""

    current: list[tuple[F, F]] = [(F(0), F(1))]
    for speed in sorted(speeds):
        bands = [
            (F(13 * k + 1, 13 * speed), F(13 * (k + 1) - 1, 13 * speed))
            for k in range(speed)
        ]
        current = intersect_interval_unions(current, bands)
        if not current:
            break
    return tuple(current)


# A packed component is ``(centre_numerator, halfwidth_numerator, denominator)``:
# centre=cn/den and halfwidth=hn/den.  It avoids Fraction allocation in the
# 5.7-million-row final containment loop.
PackedComponent = tuple[int, int, int]


def pack_components(components: Sequence[tuple[F, F]]) -> tuple[PackedComponent, ...]:
    packed: list[PackedComponent] = []
    for lo, hi in components:
        ln, ld = lo.numerator, lo.denominator
        rn, rd = hi.numerator, hi.denominator
        denominator = 2 * ld * rd
        packed.append((ln * rd + rn * ld, rn * ld - ln * rd, denominator))
    return tuple(packed)


def first_containment_failure(
    components: Sequence[PackedComponent], speed: int
) -> tuple[int, int, int, int, int, int] | None:
    """First component not contained in ``D_speed(1/13)``.

    Returns ``(index, positive_surplus, cn, hn, den, lhs)`` where

        lhs = 13*(distance_numerator(speed*centre) + speed*halfwidth_numerator)

    and containment is exactly ``lhs<=den``.
    """

    for index, (centre_num, halfwidth_num, denominator) in enumerate(components):
        residue = (speed * centre_num) % denominator
        distance_num = min(residue, denominator - residue)
        lhs = 13 * (distance_num + speed * halfwidth_num)
        if lhs > denominator:
            return (
                index,
                lhs - denominator,
                centre_num,
                halfwidth_num,
                denominator,
                lhs,
            )
    return None


def scaled_left_phase(provider: int, owner: int) -> F:
    """The phase at the owner's own-tooth left exit, multiplied by 13."""

    owner_label = owner % 13
    provider_label = provider % 13
    residue = provider_label * pow(owner_label, -1, 13) % 13
    phase = F(residue, 1) - 2 * F(provider, owner)
    # Canonical representative modulo 13 in (-13/2,13/2].
    quotient = (2 * phase.numerator + 13 * phase.denominator) // (
        26 * phase.denominator
    )
    phase -= 13 * quotient
    while phase <= F(-13, 2):
        phase += 13
    while phase > F(13, 2):
        phase -= 13
    return phase


def left_handoff(provider: int, owner: int) -> bool:
    phase = scaled_left_phase(provider, owner)
    return -1 < phase <= 1


def collar_handoff_audit() -> None:
    collar = F(1, 156)
    assert F(2, 13 * 24) == collar
    assert F(2, 13 * 25) < collar

    # The symbolic cycle exclusions used in the proof.
    assert pow(2, -1, 13) == 7
    assert pow(4, -1, 13) == 10
    middle_bands = {2, 3, 4}
    products = {a * b % 13 for a in middle_bands for b in middle_bands}
    assert products == {3, 4, 6, 8, 9, 12}
    assert 7 not in products

    # Exact stress test of the half-open predicate.  This is not the infinite
    # proof, but checks every labelled triple through lift height twelve.
    rows = 0
    all_large_rows = 0
    all_indegree_rows = 0
    for labels in combinations(BASE, 3):
        for i in range(1, 13):
            for j in range(1, 13):
                for k in range(1, 13):
                    speeds = (
                        labels[0] + 13 * i,
                        labels[1] + 13 * j,
                        labels[2] + 13 * k,
                    )
                    rows += 1
                    if min(speeds) <= 24:
                        continue
                    all_large_rows += 1
                    indegree = [0, 0, 0]
                    for provider in range(3):
                        for owner in range(3):
                            if provider != owner and left_handoff(
                                speeds[provider], speeds[owner]
                            ):
                                indegree[owner] += 1
                    if all(indegree):
                        all_indegree_rows += 1
    assert rows == 380_160
    assert all_indegree_rows == 0

    print("ORIENTED_COLLAR_HANDOFF")
    print("core_left_collar=1/156 own_tooth_length=2/(13u) all_u>24_forces_cross_handoff")
    print("left_handoff_phase=(-1,1] after multiplying by 13")
    print("subunit_handoff_residue=2 inverse_residue=7 three_cycle_middle_bands={2,3,4}")
    print(
        f"height12_stress_rows={rows} all_replacements_gt24={all_large_rows} "
        f"all_three_indegrees_positive={all_indegree_rows}"
    )
    print("conclusion=hypothetical tight proper triple lift has a replacement in [14,24]")
    print()


def analytic_box_audit() -> None:
    ten_gap = F(1, 11) - DELTA
    ten_radius = ten_gap / 24
    ten_length = 2 * ten_radius
    required_reciprocal_sum = F(9, 2) * ten_length
    assert ten_gap == F(2, 143)
    assert ten_radius == F(1, 1716)
    assert ten_length == F(1, 858)
    assert required_reciprocal_sum == F(3, 572)
    assert F(2, 382) < required_reciprocal_sum
    assert int(F(2, 1) / required_reciprocal_sum) == 381

    eleven_gap = F(1, 12) - DELTA
    assert eleven_gap == F(1, 156)
    # An eleven-core whose maximum is v has a delta-safe interval of length
    # 1/(78v).  One w-tooth has length 2/(13w), hence w<=12v.

    print("ANALYTIC_FINITE_BOX")
    print(
        f"ten_core_gap={fmt(ten_gap)} radius={fmt(ten_radius)} "
        f"safe_interval_length={fmt(ten_length)}"
    )
    print("two-comb cover requires 1/y+1/z>=3/572; with v=min(y,z), v<=381")
    print("eleven_core_gap=1/156 interval_length=1/(78v); one-tooth cover forces w<=12v")
    print()


def row_count_formula() -> int:
    """Independent closed-form count of the canonical finite box."""

    total = 0
    for labels in combinations(BASE, 3):
        for anchor_label in labels:
            x = anchor_label + 13
            if x > 24:
                continue
            others = tuple(label for label in labels if label != anchor_label)
            for v_label, w_label in (others, tuple(reversed(others))):
                for v in range(v_label + 13, 382, 13):
                    if v <= x:
                        continue
                    first_height = max(1, (v - w_label) // 13 + 1)
                    last_height = (12 * v - w_label) // 13
                    total += max(0, last_height - first_height + 1)
    return total


def finite_component_closure_audit() -> None:
    expected_rows = row_count_formula()
    assert expected_rows == 5_713_539

    rows = 0
    tight_rows = 0
    component_cores = 0
    anchor_histogram: Counter[int] = Counter()
    component_count_histogram: Counter[int] = Counter()
    failure_index_histogram: Counter[int] = Counter()
    digest = sha256()
    closest: tuple[int, int, tuple[int, ...]] | None = None

    for labels in combinations(BASE, 3):
        core = tuple(value for value in BASE if value not in labels)
        for anchor_label in labels:
            x = anchor_label + 13
            if x > 24:
                continue
            others = tuple(label for label in labels if label != anchor_label)
            for v_label, w_label in (others, tuple(reversed(others))):
                cache: dict[int, tuple[PackedComponent, ...]] = {}
                for v in range(v_label + 13, 382, 13):
                    if v <= x:
                        continue
                    components = cache.get(v)
                    if components is None:
                        exact_components = strict_safe_components((*core, x, v))
                        assert exact_components
                        components = pack_components(exact_components)
                        cache[v] = components
                        component_cores += 1
                        component_count_histogram[len(components)] += 1

                    first_height = max(1, (v - w_label) // 13 + 1)
                    last_height = (12 * v - w_label) // 13
                    for height in range(first_height, last_height + 1):
                        w = w_label + 13 * height
                        assert x < v < w <= 12 * v
                        rows += 1
                        anchor_histogram[x] += 1
                        failure = first_containment_failure(components, w)
                        if failure is None:
                            tight_rows += 1
                            digest.update(pack(">6q", *labels, x, v, w))
                            continue

                        index, surplus, cn, hn, denominator, lhs = failure
                        failure_index_histogram[index] += 1
                        digest.update(
                            pack(
                                ">12q",
                                *labels,
                                x,
                                v,
                                w,
                                index,
                                surplus,
                                cn,
                                hn,
                                denominator,
                                lhs,
                            )
                        )
                        record = (surplus, denominator, (*labels, x, v, w, index))
                        if closest is None or surplus * closest[1] < closest[0] * denominator:
                            closest = record

    assert rows == expected_rows
    assert tight_rows == 0
    assert closest is not None

    print("FINITE_COMPONENT_CLOSURE")
    print(
        "canonicalization=unique least replacement x<=24; remaining labels assigned "
        "to ordered x<v<w"
    )
    print(f"exact_rows={rows} distinct_eleven_cores={component_cores} tight_rows={tight_rows}")
    print(f"anchor_row_histogram={dict(sorted(anchor_histogram.items()))}")
    print(f"component_count_histogram={dict(sorted(component_count_histogram.items()))}")
    print(f"failure_component_index_histogram={dict(sorted(failure_index_histogram.items()))}")
    print(
        f"closest_first_certificate_surplus={closest[0]}/{closest[1]} "
        f"row_and_component={closest[2]}"
    )
    print(f"certificate_digest={digest.hexdigest()}")
    print("result=every canonical row leaves a nonempty strict-safe subinterval")
    print()


def distance_num(speed: int, numerator: int, denominator: int) -> int:
    residue = speed * numerator % denominator
    return min(residue, denominator - residue)


def exact_maximin(speeds: Sequence[int]) -> tuple[F, F]:
    denominators = {2 * speed for speed in speeds}
    for left, right in combinations(speeds, 2):
        denominators.add(left + right)
        denominators.add(abs(left - right))
    denominators.discard(0)

    best = F(0)
    witness = F(0)
    for denominator in denominators:
        for numerator in range(1, denominator):
            value = F(
                min(distance_num(speed, numerator, denominator) for speed in speeds),
                denominator,
            )
            if value > best:
                best = value
                witness = F(numerator, denominator)
    return best, witness


def height_one_crosscheck() -> None:
    rows = 0
    minimum: F | None = None
    minimizers: list[tuple[int, int, int]] = []
    for labels in combinations(BASE, 3):
        speeds = tuple(
            sorted((set(BASE) - set(labels)) | {label + 13 for label in labels})
        )
        value, _ = exact_maximin(speeds)
        assert value > DELTA
        rows += 1
        if minimum is None or value < minimum:
            minimum = value
            minimizers = [labels]
        elif value == minimum:
            minimizers.append(labels)
    assert rows == 220
    assert minimum == F(2, 21)
    assert minimizers == [(2, 5, 6), (4, 5, 6)]

    print("HEIGHT_ONE_INDEPENDENT_CROSSCHECK")
    print(
        f"rows={rows} exact_piecewise_maximin_minimum={fmt(minimum)} "
        f"missing_label_minimizers={minimizers}"
    )
    print()


def build_residue_tournament(reverse_live: bool) -> tuple[list[list[bool]], int]:
    adjacency = [[False] * 12 for _ in range(12)]
    live_edges = 0
    for low in BASE:
        for high in range(low + 1, 13):
            high_hits_low = high * pow(low, -1, 13) % 13 == 2
            low_hits_high = low * pow(high, -1, 13) % 13 == 2
            assert not (high_hits_low and low_hits_high)
            if high_hits_low or low_hits_high:
                live_edges += 1
                provider, owner = (high, low) if high_hits_low else (low, high)
                if reverse_live:
                    provider, owner = owner, provider
                adjacency[provider - 1][owner - 1] = True
            else:
                adjacency[low - 1][high - 1] = True
    return adjacency, live_edges


def tournament_fingerprint(
    adjacency: Sequence[Sequence[bool]],
) -> tuple[dict[int, int], int, list[int], int, list[int]]:
    n = len(adjacency)
    score_histogram = dict(sorted(Counter(sum(row) for row in adjacency).items()))
    triangles = 0
    for i, j, k in combinations(range(n), 3):
        triangles += bool(
            (adjacency[i][j] and adjacency[j][k] and adjacency[k][i])
            or (adjacency[i][k] and adjacency[k][j] and adjacency[j][i])
        )

    reach = [[i == j or adjacency[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    unseen = set(range(n))
    scc_sizes: list[int] = []
    while unseen:
        root = min(unseen)
        component = {v for v in unseen if reach[root][v] and reach[v][root]}
        scc_sizes.append(len(component))
        unseen -= component
    scc_sizes.sort(reverse=True)

    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if count == 0:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adjacency[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count

    path: list[int] = []
    for vertex in range(n):
        position = 0
        while position < len(path) and adjacency[path[position]][vertex]:
            position += 1
        path.insert(position, vertex)
    assert all(adjacency[path[i]][path[i + 1]] for i in range(n - 1))
    return (
        score_histogram,
        triangles,
        scc_sizes,
        sum(dp[-1]),
        [vertex + 1 for vertex in path],
    )


def tournament_audit() -> None:
    forward, live = build_residue_tournament(False)
    switched, switched_live = build_residue_tournament(True)
    assert live == switched_live == 12
    edge_flips = sum(
        forward[i][j] != switched[i][j]
        for i in range(12)
        for j in range(i + 1, 12)
    )
    assert edge_flips == 12
    first = tournament_fingerprint(forward)
    second = tournament_fingerprint(switched)

    print("TOURNAMENT_TELEMETRY")
    print("vertices=residue_owners observable=subunit_ratio_2_handoff tie_gauge=increasing_residue")
    print(f"switch=reverse_live_handoffs live_edges={live} edge_flips={edge_flips}")
    print(
        f"forward_score_histogram={first[0]} directed_triangles={first[1]} "
        f"scc_sizes={first[2]} hamiltonian_paths={first[3]} tie_path={first[4]}"
    )
    print(
        f"switched_score_histogram={second[0]} directed_triangles={second[1]} "
        f"scc_sizes={second[2]} hamiltonian_paths={second[3]} tie_path={second[4]}"
    )
    print("interpretation=tournament_is_telemetry; proof_retains_ratios_half_open_flags_and_widths")
    print()


def main() -> None:
    print("LRC13 SCALE-ONE HAMMING-THREE COLLAR CLOSURE — EXACT REPLAY")
    print("arithmetic=integer+fractions.Fraction floating_point=none")
    print("scope=all distinct labels and all three positive lift heights")
    print()

    collar_handoff_audit()
    analytic_box_audit()
    finite_component_closure_audit()
    height_one_crosscheck()
    tournament_audit()

    print("PASS: every proper labelled scale-one triple lift has M>1/13.")


if __name__ == "__main__":
    main()
