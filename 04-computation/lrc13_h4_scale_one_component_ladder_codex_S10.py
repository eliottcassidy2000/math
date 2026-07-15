#!/usr/bin/env python3
"""Exact replay for THM-815: scale-one Hamming-four closure.

For four distinct labels ``R`` in ``{1,...,12}`` and positive lift heights,
the packet is

    ({1,...,12} - R) union {r + 13*h_r : r in R}.

The proof is a recursive safe-component ladder.  Exact component censuses and
the sharp one-period danger discrepancy successively bound the ordered lifts
``x < v < w < z`` by ``x<=105``, ``v<=118``, ``w<=88``, then ``w<=83``.
An eleven-core census gives ``z<=50``.  The remaining 35,640 rows are closed
by exact component containment.  Fractions and integer cross-products are
used throughout; there is no floating point or sampled-circle decision.

The same discrepancy has positive cover deficit whenever at most six danger
combs remain.  The replay records the first exact recursion step at Hamming
radii five and six, while making no claim that either finite tree is empty.

Tournament Analysis is faithful only at the earlier collar-exit layer.  The
replay also checks the unique four-cycle band word ``(2,2,2,5)`` and records
the tournament fingerprint of a loose row realizing it.  The theorem-bearing
finite carrier is instead bipartite: safe components versus remaining danger
combs, decorated by component length.  Antisymmetrizing that carrier would
destroy the covering predicate.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import prod
from typing import Iterable, Sequence


DELTA = F(1, 13)
BASE = tuple(range(1, 13))


def fmt(q: F) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def safe_bands(speed: int) -> tuple[tuple[F, F], ...]:
    return tuple(
        (F(13 * k + 1, 13 * speed), F(13 * (k + 1) - 1, 13 * speed))
        for k in range(speed)
    )


BANDS = {speed: safe_bands(speed) for speed in range(1, 119)}


def intersect_interval_unions(
    left: Sequence[tuple[F, F]], right: Sequence[tuple[F, F]]
) -> tuple[tuple[F, F], ...]:
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
    return tuple(out)


def strict_safe_components(speeds: Iterable[int]) -> tuple[tuple[F, F], ...]:
    current: tuple[tuple[F, F], ...] = ((F(0), F(1)),)
    for speed in sorted(speeds):
        current = intersect_interval_unions(current, BANDS[speed])
        if not current:
            break
    return current


def longest_component(components: Sequence[tuple[F, F]]) -> tuple[F, F]:
    assert components
    return max(components, key=lambda interval: interval[1] - interval[0])


def fraction_record(q: F) -> tuple[int, int]:
    return q.numerator, q.denominator


@dataclass
class Census:
    name: str
    rows: int = 0
    minimum: F | None = None
    minimizers: list[tuple[object, ...]] = field(default_factory=list)
    digest: object = field(default_factory=sha256)

    def add(self, key: tuple[object, ...], components: Sequence[tuple[F, F]]) -> None:
        lo, hi = longest_component(components)
        length = hi - lo
        record = key + (
            len(components),
            lo.numerator,
            lo.denominator,
            hi.numerator,
            hi.denominator,
        )
        self.digest.update((repr(record) + "\n").encode())
        self.rows += 1
        witness = key + ((lo, hi), len(components))
        if self.minimum is None or length < self.minimum:
            self.minimum = length
            self.minimizers = [witness]
        elif length == self.minimum:
            self.minimizers.append(witness)

    def line(self) -> str:
        assert self.minimum is not None
        witnesses: object = (
            tuple(self.minimizers) if len(self.minimizers) <= 5 else self.minimizers[0]
        )
        return (
            f"{self.name}.rows={self.rows} minimum_longest_component={fmt(self.minimum)} "
            f"minimizer_count={len(self.minimizers)} minimizers={witnesses} "
            f"sha256={self.digest.hexdigest()}"
        )


def sharp_discrepancy_audit() -> None:
    # For a circle arc A of length a and a residual interval of length s,
    # max_s(min(a,s)-a*s)=a*(1-a).  Here a=2/13.
    a = 2 * DELTA
    error = a * (1 - a)
    assert a == F(2, 13) and error == F(22, 169)

    # If m danger combs cover an interval of length L, then
    # sum(1/u_i) >= 13*(13-2m)*L/22.
    assert F(13 * (13 - 8), 22) * F(1, 78) == F(65, 1716)
    assert int(F(4, 1) / F(65, 1716)) == 105
    assert int(F(22 * 3, 13 * 7) / F(7, 1144)) == 118
    assert int(F(22 * 2, 13 * 9) / F(23, 5434)) == 88
    assert int(F(22 * 2, 13 * 9) / F(1, 221)) == 83
    assert int(F(2, 13) / F(1, 325)) == 50


def higher_radius_first_step_audit() -> tuple[tuple[Census, int, int], ...]:
    """Record the first recursive cap at Hamming radii five and six.

    The rest of those finite trees is not enumerated here.  This census only
    verifies the exact initial component minima, the resulting discrepancy
    caps, and the number of labelled one-lift prefixes below each cap.
    """

    rows: list[tuple[Census, int, int]] = []
    expected = {
        5: (
            F(1, 52),
            146,
            40_590,
            "5943566f8db6dd00177782f3403691ef616d2be921a179ff00f2cfca9aea2814",
        ),
        6: (
            F(31, 1430),
            468,
            194_040,
            "1e338bd839aa68384bdcd66e4061bd31e918399263e5fb9ffe59e8e867d5930e",
        ),
    }
    base_set = set(BASE)
    for missing_count in (5, 6):
        census = Census(f"hamming_{missing_count}_initial_core")
        for labels in combinations(BASE, missing_count):
            retained = tuple(sorted(base_set - set(labels)))
            census.add((labels,), strict_safe_components(retained))

        minimum, cap, labelled_prefixes, digest = expected[missing_count]
        assert census.rows == len(tuple(combinations(BASE, missing_count)))
        assert census.minimum == minimum
        assert census.digest.hexdigest() == digest
        computed_cap = int(
            F(22 * missing_count, 13 * (13 - 2 * missing_count)) / minimum
        )
        assert computed_cap == cap
        computed_prefixes = sum(
            len(tuple(range(label + 13, cap + 1, 13)))
            for labels in combinations(BASE, missing_count)
            for label in labels
        )
        assert computed_prefixes == labelled_prefixes
        rows.append((census, cap, labelled_prefixes))
    return tuple(rows)


def ladder_censuses() -> tuple[Census, Census, Census, Census, Census, Census]:
    c8 = Census("eight_core")
    c9 = Census("nine_core_x_le_105")
    c10 = Census("ten_core_x_le_105_v_le_118")
    c10_fixed = Census("ten_core_fixed_x_lt_v_lt_88")
    c11 = Census("eleven_core_x_lt_v_lt_w_le_83")
    c11_fixed = Census("eleven_core_fixed_x_lt_v_lt_w_lt_51")

    base_set = set(BASE)
    for labels in combinations(BASE, 4):
        retained = tuple(sorted(base_set - set(labels)))
        components8 = strict_safe_components(retained)
        c8.add((labels,), components8)

        for x_label in labels:
            for x in range(x_label + 13, 106, 13):
                components9 = intersect_interval_unions(components8, BANDS[x])
                c9.add((labels, x_label, x), components9)

                for v_label in labels:
                    if v_label == x_label:
                        continue
                    for v in range(v_label + 13, 119, 13):
                        if v <= x:
                            continue
                        components10 = intersect_interval_unions(components9, BANDS[v])
                        key10 = (labels, x_label, x, v_label, v)
                        c10.add(key10, components10)
                        if v < 88:
                            c10_fixed.add(key10, components10)

                        if v >= 83:
                            continue
                        for w_label in labels:
                            if w_label in (x_label, v_label):
                                continue
                            for w in range(w_label + 13, 84, 13):
                                if w <= v:
                                    continue
                                components11 = intersect_interval_unions(
                                    components10, BANDS[w]
                                )
                                key11 = key10 + (w_label, w)
                                c11.add(key11, components11)
                                if w < 51:
                                    c11_fixed.add(key11, components11)

    assert c8.rows == 495 and c8.minimum == F(1, 78)
    assert c9.rows == 14_025 and c9.minimum == F(7, 1144)
    assert c10.rows == 191_070 and c10.minimum == F(23, 5434)
    assert c10_fixed.rows == 98_145 and c10_fixed.minimum == F(1, 221)
    assert c11.rows == 313_965 and c11.minimum == F(1, 325)
    assert c11_fixed.rows == 49_005 and c11_fixed.minimum == F(1, 325)
    expected_digests = (
        "c2c884d9c1eff1f011da1ad06c927af75e16f3f64d629fa1929f37a83b380de3",
        "5a99e2b06830e2159bc17be322c85ac3091ba95066745dde22b0ad88d6935ae7",
        "226c3e58e9d4502a7906e3f0e7ce5b9248ba5b34919ad9bb6ce39348d0b58eb8",
        "c6bc04110ed271bded8aacf5f315b1af727102dea3b31f63cbc4ade7ffb34846",
        "fe3e611b3059c2443b019e673d349a9d218a493db5ea0a78bec03159b26cf01c",
        "8049d45c937f05c41102428820e70d7ff63ca3fc2d2501785f507eeeabb33286",
    )
    assert tuple(census.digest.hexdigest() for census in (c8, c9, c10, c10_fixed, c11, c11_fixed)) == expected_digests
    return c8, c9, c10, c10_fixed, c11, c11_fixed


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
    """Return the first component not contained in D_speed(1/13)."""

    for index, (centre_num, halfwidth_num, denominator) in enumerate(components):
        residue = (speed * centre_num) % denominator
        distance_num = min(residue, denominator - residue)
        lhs = 13 * (distance_num + speed * halfwidth_num)
        if lhs > denominator:
            return index, lhs - denominator, centre_num, halfwidth_num, denominator, lhs
    return None


def final_containment_census() -> tuple[int, int, str, object, Counter[int], Counter[int]]:
    rows = tight_rows = 0
    digest = sha256()
    closest: tuple[int, int, object] | None = None
    component_histogram: Counter[int] = Counter()
    failure_index_histogram: Counter[int] = Counter()
    base_set = set(BASE)

    for labels in combinations(BASE, 4):
        retained = tuple(sorted(base_set - set(labels)))
        lift_banks = [tuple(range(label + 13, 51, 13)) for label in labels]
        cache: dict[tuple[int, int, int], tuple[PackedComponent, ...]] = {}
        for lifts in product(*lift_banks):
            x, v, w, z = sorted(lifts)
            assert x < v < w < z <= 50
            key = (x, v, w)
            packed = cache.get(key)
            if packed is None:
                components = strict_safe_components((*retained, x, v, w))
                packed = pack_components(components)
                cache[key] = packed
                component_histogram[len(packed)] += 1

            rows += 1
            failure = first_containment_failure(packed, z)
            if failure is None:
                tight_rows += 1
                digest.update((repr((labels, lifts, (x, v, w, z), "tight")) + "\n").encode())
                continue

            index, surplus, cn, hn, denominator, lhs = failure
            failure_index_histogram[index] += 1
            record = (
                labels,
                lifts,
                (x, v, w, z),
                index,
                surplus,
                cn,
                hn,
                denominator,
                lhs,
            )
            digest.update((repr(record) + "\n").encode())
            normalized = (surplus, denominator, record)
            if closest is None or surplus * closest[1] < closest[0] * denominator:
                closest = normalized

    assert rows == 35_640 and tight_rows == 0 and closest is not None
    assert digest.hexdigest() == "118e9413d8e9b4daf3a240b96a6f70d4760ae0771485cec44a1cdb3af8f704cf"
    return (
        rows,
        tight_rows,
        digest.hexdigest(),
        closest,
        component_histogram,
        failure_index_histogram,
    )


def norm(q: F) -> F:
    residue = q % 1
    return min(residue, 1 - residue)


def direct_endpoint_crosscheck() -> tuple[int, int, str, object]:
    """Independent full-packet endpoint-cell witness over the final box.

    This does not form an eleven-core and does not invoke the last-comb
    containment criterion.  All threshold endpoints for all twelve speeds are
    sorted at once; a midpoint is tested in each resulting open cell.
    """

    rows = failures = 0
    digest = sha256()
    smallest_margin: tuple[F, object] | None = None
    base_set = set(BASE)
    for labels in combinations(BASE, 4):
        retained = tuple(sorted(base_set - set(labels)))
        lift_banks = [tuple(range(label + 13, 51, 13)) for label in labels]
        for lifts in product(*lift_banks):
            speeds = tuple(sorted((*retained, *lifts)))
            endpoints = {F(0), F(1)}
            for speed in speeds:
                for lo, hi in BANDS[speed]:
                    endpoints.add(lo)
                    endpoints.add(hi)
            ordered = sorted(endpoints)

            witness: F | None = None
            clearance: F | None = None
            cell: tuple[F, F] | None = None
            for lo, hi in zip(ordered, ordered[1:]):
                midpoint = (lo + hi) / 2
                value = min(norm(speed * midpoint) for speed in speeds)
                if value > DELTA:
                    witness, clearance, cell = midpoint, value, (lo, hi)
                    break
            rows += 1
            if witness is None or clearance is None or cell is None:
                failures += 1
                digest.update((repr((labels, lifts, speeds, "failure")) + "\n").encode())
                continue
            record = (
                labels,
                lifts,
                speeds,
                witness.numerator,
                witness.denominator,
                clearance.numerator,
                clearance.denominator,
                cell[0].numerator,
                cell[0].denominator,
                cell[1].numerator,
                cell[1].denominator,
            )
            digest.update((repr(record) + "\n").encode())
            margin = clearance - DELTA
            if smallest_margin is None or margin < smallest_margin[0]:
                smallest_margin = (margin, record)

    assert rows == 35_640 and failures == 0 and smallest_margin is not None
    assert digest.hexdigest() == "82823bf934a438e4dcbcff2724f322da095ae23b54a333a55e91e8eb54face8c"
    return rows, failures, digest.hexdigest(), smallest_margin


def scaled_left_phase(provider: int, owner: int) -> F:
    provider_label = provider % 13
    owner_label = owner % 13
    residue = provider_label * pow(owner_label, -1, 13) % 13
    phase = F(residue) - 2 * F(provider, owner)
    while phase <= F(-13, 2):
        phase += 13
    while phase > F(13, 2):
        phase -= 13
    return phase


def left_handoff(provider: int, owner: int) -> bool:
    return -1 < scaled_left_phase(provider, owner) <= 1


def tournament_fingerprint(edge: Sequence[Sequence[bool]]) -> tuple[object, ...]:
    n = len(edge)
    scores = tuple(sum(row) for row in edge)
    triangles = sum(
        bool(
            edge[i][j] and edge[j][k] and edge[k][i]
            or edge[i][k] and edge[k][j] and edge[j][i]
        )
        for i, j, k in combinations(range(n), 3)
    )
    reach = [[i == j or edge[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or reach[i][k] and reach[k][j]
    unseen = set(range(n))
    scc: list[int] = []
    while unseen:
        root = min(unseen)
        component = {j for j in unseen if reach[root][j] and reach[j][root]}
        scc.append(len(component))
        unseen -= component
    hamiltonian_paths = sum(
        all(edge[path[i]][path[i + 1]] for i in range(n - 1))
        for path in permutations(range(n))
    )
    return scores, triangles, tuple(sorted(scc, reverse=True)), hamiltonian_paths


def exact_norm(speed: int, numerator: int, denominator: int) -> int:
    residue = speed * numerator % denominator
    return min(residue, denominator - residue)


def exact_maximin(speeds: Sequence[int]) -> tuple[F, F]:
    denominators = {2 * speed for speed in speeds}
    for a, b in combinations(speeds, 2):
        denominators.add(a + b)
        denominators.add(abs(a - b))
    denominators.discard(0)
    best = F(0)
    witness = F(0)
    for denominator in denominators:
        for numerator in range(1, denominator):
            value = F(
                min(exact_norm(speed, numerator, denominator) for speed in speeds),
                denominator,
            )
            if value > best:
                best, witness = value, F(numerator, denominator)
    return best, witness


def collar_cycle_and_tournament_audit() -> tuple[object, ...]:
    words = tuple(
        word
        for word in product(range(2, 18), repeat=4)
        if 2 in word
        and prod(a - 1 for a in word) <= 16 < prod(a + 1 for a in word)
        and prod(word) % 13 == 1
    )
    assert words == ((2, 2, 2, 5), (2, 2, 5, 2), (2, 5, 2, 2), (5, 2, 2, 2))

    labels = (1, 2, 4, 8)
    speeds = {1: 79, 2: 54, 4: 30, 8: 34}
    live = {
        (provider, owner)
        for provider in labels
        for owner in labels
        if provider != owner and left_handoff(speeds[provider], speeds[owner])
    }
    assert live == {(1, 8), (8, 4), (4, 2), (2, 1)}

    def gauge(reverse: bool) -> list[list[bool]]:
        edge = [[False] * 4 for _ in range(4)]
        for i, a in enumerate(labels):
            for j in range(i + 1, 4):
                b = labels[j]
                if (a, b) in live:
                    source, target = (b, a) if reverse else (a, b)
                elif (b, a) in live:
                    source, target = (a, b) if reverse else (b, a)
                else:
                    source, target = a, b
                edge[labels.index(source)][labels.index(target)] = True
        return edge

    forward = gauge(False)
    switched = gauge(True)
    edge_flips = sum(
        forward[i][j] != switched[i][j] for i in range(4) for j in range(i + 1, 4)
    )
    assert edge_flips == 4

    retained = tuple(sorted(set(BASE) - set(labels)))
    packet = tuple(sorted((*retained, *speeds.values())))
    maximum, witness = exact_maximin(packet)
    assert maximum == F(3, 19)
    return words, live, tournament_fingerprint(forward), tournament_fingerprint(switched), edge_flips, maximum, witness


def main() -> None:
    print("THM-815 SCALE-ONE HAMMING-FOUR COMPONENT LADDER — EXACT REPLAY")
    print("arithmetic=fractions.Fraction+integer_cross_products floating_point=none")
    print("scope=all four distinct labels and all four positive lift heights")
    print()

    sharp_discrepancy_audit()
    print("SHARP_INTERVAL_DISCREPANCY")
    print("bound=|I intersect D_u(1/13)| <= (2/13)|I| + 22/(169u)")
    print("derivation=max_s(min(2/13,s)-(2/13)s)=22/169")
    print()

    print("HAMMING_FIVE_SIX_FIRST_RECURSION_STEP")
    for census, cap, labelled_prefixes in higher_radius_first_step_audit():
        print(
            f"{census.line()} next_speed_cap={cap} "
            f"labelled_first_prefixes={labelled_prefixes}"
        )
    print("density_barrier=remaining_combs_7 gives 13-2m=-1; no discrepancy cap")
    print()

    censuses = ladder_censuses()
    print("SAFE_COMPONENT_LADDER")
    for census in censuses:
        print(census.line())
    print("bounds=x<=105, v<=118, w<=88, fixed_point_w<=83, z<=50")
    print()

    rows, tight, digest, closest, component_hist, failure_hist = final_containment_census()
    print("FINAL_COMPONENT_CONTAINMENT")
    print("domain=all proper lifts <=50; canonical numerical order x<v<w<z")
    print(f"rows={rows} tight_rows={tight}")
    print(f"distinct_core_component_histogram={dict(sorted(component_hist.items()))}")
    print(f"failure_component_index_histogram={dict(sorted(failure_hist.items()))}")
    print(f"closest_first_failure={closest}")
    print(f"sha256={digest}")
    print("result=every row leaves a nonempty strict-safe interval")
    print()

    endpoint = direct_endpoint_crosscheck()
    print("INDEPENDENT_FULL_PACKET_ENDPOINT_CROSSCHECK")
    print("method=sort_all_twelve_speed_threshold_endpoints_and_test_exact_cell_midpoints")
    print(f"rows={endpoint[0]} failures={endpoint[1]} sha256={endpoint[2]}")
    print(f"smallest_selected_witness_margin={endpoint[3]}")
    print()

    cycle = collar_cycle_and_tournament_audit()
    print("COLLAR_CYCLE_AND_TOURNAMENT_AUDIT")
    print(f"admissible_band_words={cycle[0]}")
    print(f"method_limit_live_handoffs={tuple(sorted(cycle[1]))}")
    print(f"forward_fingerprint={cycle[2]}")
    print(f"switched_fingerprint={cycle[3]} edge_flips={cycle[4]}")
    print(f"method_limit_M={cycle[5]} witness={cycle[6]}")
    print("carrier=owner_collar_exit_obligations; component ladder is bipartite, not tournament-exact")
    print()

    print("PASS: every proper scale-one Hamming-four lift has M>1/13.")


if __name__ == "__main__":
    main()
