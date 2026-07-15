#!/usr/bin/env python3
"""Exact compact verifier for THM-860.

The four independent finite checks are:

1. optimize the full primewise upper-set relative cuts on six valuations;
2. reconstruct the longest strict-safe root component for all C(12,6) cores;
3. enumerate every primitive-compatible non-all-one c=2 order/label
   common-sheet presentation and audit its signed-doubling graph; and
4. replay the explicit c=3 mixed-order common-sheet false positive and its
   exact lonely-runner maximum.

Only integer arithmetic and fractions.Fraction are used.  The theorem's
all-order relative-cut proof and finite-decidability conclusion are analytic;
this source freezes the small exact arithmetic on which its optimized
constants and boundary examples depend.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd, lcm
from pathlib import Path


LINES: list[str] = []


def emit(line: str = "") -> None:
    LINES.append(line)
    print(line)


def ceil_fraction(x: F) -> int:
    return (x.numerator + x.denominator - 1) // x.denominator


def rho(m: int) -> F:
    assert m >= 1
    return F((2 * m + 12) // 13, m)


def relative_cut_sum(orders: tuple[int, ...], mask: int) -> F:
    complement = [orders[i] for i in range(6) if not ((mask >> i) & 1)]
    complement_lcm = lcm(*complement) if complement else 1
    return sum(
        (
            rho(orders[i] // gcd(orders[i], complement_lcm))
            for i in range(6)
            if (mask >> i) & 1
        ),
        F(0),
    )


def all_relative_cuts_hold(orders: tuple[int, ...]) -> bool:
    return all(relative_cut_sum(orders, mask) >= 1 for mask in range(1, 1 << 6))


def upper_cuts_hold(p: int, valuations: tuple[int, ...]) -> bool:
    assert len(valuations) == 6
    assert tuple(sorted(valuations)) == valuations
    assert valuations[0] == 0
    assert valuations[-2] == valuations[-1]
    for k in range(5):
        if valuations[k] == valuations[k + 1]:
            continue
        contribution = sum(
            (rho(p ** (valuations[i] - valuations[k])) for i in range(k + 1, 6)),
            F(0),
        )
        if contribution < 1:
            return False
    return True


def valuation_optimization() -> None:
    emit("THM-860 PRIMITIVE H6 FINITE-RAMIFICATION EXACT VERIFIER")
    emit("SECTION valuation upper-set cuts")
    expected = {
        2: (8, 5, ((0, 2, 3, 4, 5, 5), (0, 3, 3, 4, 5, 5))),
        3: (4, 3, ((0, 1, 2, 3, 3, 3), (0, 2, 2, 3, 3, 3))),
        5: (1, 1, ((0, 1, 1, 1, 1, 1),)),
        7: (2, 1, ((0, 0, 1, 1, 1, 1), (0, 1, 1, 1, 1, 1))),
        11: (1, 0, ((0, 0, 0, 0, 0, 0),)),
    }
    for p, (coarse_cap, wanted_range, wanted_words) in expected.items():
        survivors: list[tuple[int, ...]] = []
        for tail in combinations_with_replacement(range(coarse_cap + 1), 5):
            word = (0,) + tail
            if word[-2] != word[-1]:
                continue
            if upper_cuts_hold(p, word):
                survivors.append(word)
        best = max(word[-1] for word in survivors)
        extremals = tuple(word for word in survivors if word[-1] == best)
        assert best == wanted_range
        assert extremals == wanted_words
        for word in extremals:
            orders = tuple(p**a for a in word)
            assert all_relative_cuts_hold(orders)
        emit(
            f"p={p} coarse_range={coarse_cap} optimized_range={best} "
            f"extremals={extremals}"
        )

    # The attenuation pivot is 77.  The all-six relative cut is sharper:
    # rho(72)=1/6, while rho(m)<1/6 for every allowed m>=73.  The finite
    # boundary checks below complement the analytic m>=79 estimate in THM-860.
    attenuation_pivot = 77
    all_six_pivot = 72
    assert rho(all_six_pivot) == F(1, 6)
    assert all(rho(m) < F(1, 6) for m in range(73, 78))
    assert 78 % 13 == 0

    quotient_multiplier = 2**5 * 3**3 * 5 * 7
    scale_cap = all_six_pivot * quotient_multiplier
    assert quotient_multiplier == 30_240
    assert attenuation_pivot * quotient_multiplier == 2_328_480
    assert scale_cap == 2_177_280
    emit(
        "quotient_divisor=2^5*3^3*5*7="
        f"{quotient_multiplier} attenuation_pivot={attenuation_pivot} "
        f"all_six_pivot={all_six_pivot} "
        f"scale_cap={all_six_pivot}*{quotient_multiplier}={scale_cap}"
    )
    emit(
        "rho_table="
        f"2:[{rho(2)},{rho(4)},{rho(8)},{rho(16)},{rho(32)},{rho(64)}] "
        f"3:[{rho(3)},{rho(9)},{rho(27)},{rho(81)}] "
        f"5:[{rho(5)},{rho(25)}] 7:[{rho(7)},{rho(49)}] 11:[{rho(11)}]"
    )


def danger_intervals(speed: int) -> tuple[tuple[F, F], ...]:
    radius = F(1, 13 * speed)
    intervals: list[tuple[F, F]] = []
    for k in range(speed):
        lo = F(k, speed) - radius
        hi = F(k, speed) + radius
        if lo < 0:
            intervals.append((F(0), hi))
            intervals.append((lo + 1, F(1)))
        elif hi > 1:
            intervals.append((lo, F(1)))
            intervals.append((F(0), hi - 1))
        else:
            intervals.append((lo, hi))
    return tuple(intervals)


DANGER = {speed: danger_intervals(speed) for speed in range(1, 13)}


def longest_strict_safe_component(speeds: tuple[int, ...]) -> F:
    intervals = sorted(interval for speed in speeds for interval in DANGER[speed])
    merged: list[list[F]] = []
    for lo, hi in intervals:
        if not merged or lo > merged[-1][1]:
            merged.append([lo, hi])
        elif hi > merged[-1][1]:
            merged[-1][1] = hi
    last = F(0)
    longest = F(0)
    for lo, hi in merged:
        longest = max(longest, lo - last)
        last = max(last, hi)
    return max(longest, F(1) - last)


def fnv64(data: bytes) -> int:
    value = 0xCBF29CE484222325
    for byte in data:
        value ^= byte
        value = (value * 0x100000001B3) & ((1 << 64) - 1)
    return value


def root_component_census() -> None:
    emit("SECTION 924 root components")
    rows: list[tuple[F, tuple[int, ...]]] = []
    payload: list[str] = []
    for core in combinations(range(1, 13), 6):
        length = longest_strict_safe_component(core)
        rows.append((length, core))
        payload.append(f"{','.join(map(str, core))}:{length.numerator}/{length.denominator}")
    rows.sort()
    minimum = rows[0][0]
    minimizers = tuple(core for length, core in rows if length == minimum)
    assert len(rows) == 924
    assert minimum == F(31, 1430)
    assert minimizers == ((1, 4, 7, 9, 10, 11),)
    normalized_cap = F(132, 13) / minimum
    per_ray = ceil_fraction(normalized_cap / 13)
    c2_per_ray = (normalized_cap - F(1, 2)) // 13 + 1
    assert normalized_cap == F(14520, 31)
    assert per_ray == 37
    assert c2_per_ray == 36
    root_hash = fnv64(("\n".join(payload) + "\n").encode())
    emit(
        f"roots={len(rows)} min_length={minimum} multiplicity={len(minimizers)} "
        f"core={minimizers[0]}"
    )
    emit(
        f"normalized_real_cap={normalized_cap} cap_over_13={normalized_cap / 13} "
        f"first_heights_per_ray<={per_ray} c2_first_heights_per_ray<={c2_per_ray}"
    )
    emit(f"root_length_payload_fnv64={root_hash:016x}")


def centered_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, order: int, unit: int) -> int:
    return next(
        u
        for u in range(13 * order)
        if u % 13 == (order * label) % 13 and u % order == unit % order
    )


def sheet_mask(label: int, order: int, unit: int, owner: int, scale: int) -> int:
    u = crt_base(label, order, unit)
    inverse_owner = pow(owner, -1, 13)
    mask = 0
    for sheet in range(scale):
        z = centered_residue(u * (inverse_owner + 13 * sheet), 13 * order)
        if -order < z <= order:
            mask |= 1 << sheet
    return mask


def common_sheet_c2(labels: tuple[int, ...], orders: tuple[int, ...]) -> bool:
    full = 0b11
    return all(
        reduce(
            int.__or__,
            (
                sheet_mask(label, order, 0 if order == 1 else 1, owner, 2)
                for label, order in zip(labels, orders)
            ),
            0,
        )
        == full
        for owner in labels
    )


def signed_doubling_edges(labels: tuple[int, ...]) -> set[tuple[int, int]]:
    return {
        (provider, owner)
        for provider in labels
        for owner in labels
        if provider != owner and owner * pow(provider, -1, 13) % 13 in (2, 11)
    }


def directed_cycles(
    labels: tuple[int, ...], edges: set[tuple[int, int]]
) -> set[tuple[int, ...]]:
    cycles: set[tuple[int, ...]] = set()

    def walk(start: int, path: tuple[int, ...]) -> None:
        for nxt in labels:
            if (path[-1], nxt) not in edges:
                continue
            if nxt == start and len(path) >= 2:
                rotations = tuple(path[i:] + path[:i] for i in range(len(path)))
                cycles.add(min(rotations))
            elif nxt not in path:
                walk(start, path + (nxt,))

    for start in labels:
        walk(start, (start,))
    return cycles


def c2_census() -> None:
    emit("SECTION c=2 direct masks and signed cycles")
    rows: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for labels in combinations(range(1, 13), 6):
        for orders in product((1, 2), repeat=6):
            if sum(order == 2 for order in orders) < 2:
                continue
            if common_sheet_c2(labels, orders):
                rows.append((labels, orders))
    assert len(rows) == 64
    assert {orders for _, orders in rows} == {(2, 2, 2, 2, 2, 2)}
    label_rows = tuple(labels for labels, _ in rows)

    unseen = set(label_rows)
    orbits: list[set[tuple[int, ...]]] = []
    while unseen:
        representative = min(unseen)
        orbit = {
            tuple(sorted(multiplier * label % 13 for label in representative))
            for multiplier in range(1, 13)
        }
        orbits.append(orbit)
        unseen -= orbit
    orbit_sizes = tuple(len(orbit) for orbit in orbits)
    orbit_representatives = tuple(min(orbit) for orbit in orbits)
    assert orbit_sizes == (12, 12, 12, 4, 12, 12)
    assert orbit_representatives == (
        (1, 2, 3, 4, 5, 6),
        (1, 2, 3, 4, 5, 7),
        (1, 2, 3, 4, 6, 8),
        (1, 2, 3, 5, 6, 9),
        (1, 2, 3, 5, 7, 9),
        (1, 2, 6, 8, 9, 10),
    )

    fingerprints: Counter[tuple[int, int, int, int, int]] = Counter()
    negative_histogram: Counter[int] = Counter()
    row_payload: list[str] = []
    for labels in label_rows:
        edges = signed_doubling_edges(labels)
        cycles = directed_cycles(labels, edges)
        hamilton_paths = sum(
            all((path[i], path[i + 1]) in edges for i in range(5))
            for path in permutations(labels)
        )
        assert len(edges) == 6
        assert len(cycles) == 1
        cycle = next(iter(cycles))
        assert len(cycle) == 6
        indegrees = Counter(owner for _, owner in edges)
        outdegrees = Counter(provider for provider, _ in edges)
        assert set(indegrees.values()) == {1}
        assert set(outdegrees.values()) == {1}
        ratios = tuple(
            cycle[(i + 1) % 6] * pow(cycle[i], -1, 13) % 13 for i in range(6)
        )
        negatives = ratios.count(11)
        assert negatives % 2 == 1
        negative_histogram[negatives] += 1
        fingerprints[(len(edges), len(cycles), hamilton_paths, 6, negatives % 2)] += 1
        row_payload.append(
            f"{','.join(map(str, labels))}:"
            f"{','.join(f'{a}>{b}' for a, b in sorted(edges))}"
        )

    assert fingerprints == Counter({(6, 1, 6, 6, 1): 64})
    assert 12 * 2**5 // 6 == 64
    rows_hash = fnv64(("\n".join(row_payload) + "\n").encode())
    emit(
        f"presentations={len(rows)} order_patterns={{(2,2,2,2,2,2):64}} "
        f"symbolic_count=12*32/6={12 * 2**5 // 6}"
    )
    emit(
        f"multiplicative_orbit_sizes={orbit_sizes} representatives={orbit_representatives}"
    )
    emit(
        "graph_fingerprint=(edges=6,SCC=[6],directed_cycles=1,"
        f"Hamiltonian_paths=6) count=64 negative_edge_histogram={dict(sorted(negative_histogram.items()))}"
    )
    emit(f"c2_row_edge_payload_fnv64={rows_hash:016x}")


def circle_norm(numerator: int, denominator: int) -> F:
    residue = numerator % denominator
    return F(min(residue, denominator - residue), denominator)


def exact_maximin(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    denominators = {2 * speed for speed in speeds}
    denominators |= {
        denominator
        for u, v in combinations(speeds, 2)
        for denominator in (u + v, abs(u - v))
        if denominator
    }
    best = F(0)
    witnesses: set[F] = set()
    for denominator in sorted(denominators):
        for numerator in range(denominator):
            value = min(circle_norm(speed * numerator, denominator) for speed in speeds)
            witness = F(numerator, denominator)
            if value > best:
                best = value
                witnesses = {witness}
            elif value == best:
                witnesses.add(witness)
    return best, tuple(sorted(witnesses))


def mask_to_word(mask: int, scale: int) -> str:
    sheets = "".join(str(sheet) for sheet in range(scale) if (mask >> sheet) & 1)
    return sheets if sheets else "-"


def c3_counterexample() -> None:
    emit("SECTION c=3 mixed common-sheet false positive")
    scale = 3
    replacements = {1: 16, 2: 45, 3: 48, 5: 28, 8: 37, 12: 10}
    providers = (1, 5, 8, 12, 2, 3)
    owners = tuple(sorted(replacements))
    rows: list[str] = []
    for owner in owners:
        words: list[str] = []
        union = 0
        for provider in providers:
            speed = replacements[provider]
            order = scale // gcd(scale, speed)
            u = order * speed // scale
            unit = u % order
            mask = sheet_mask(provider, order, unit, owner, scale)
            union |= mask
            words.append(mask_to_word(mask, scale))
        assert union == 0b111
        rows.append(f"owner={owner} masks={','.join(words)}")
    packet = tuple(
        sorted(
            [scale * r for r in range(1, 13) if r not in replacements]
            + list(replacements.values())
        )
    )
    assert packet == (10, 12, 16, 18, 21, 27, 28, 30, 33, 37, 45, 48)
    maximum, witnesses = exact_maximin(packet)
    assert maximum == F(5, 29)
    assert witnesses == (F(1, 58), F(57, 58))
    emit(
        f"scale={scale} labels={owners} packet={packet} orders=(3,1,1,3,3,3)"
    )
    for row in rows:
        emit(row)
    emit(f"exact_M={maximum} witnesses={witnesses} verdict=loose")


def main() -> None:
    valuation_optimization()
    root_component_census()
    c2_census()
    c3_counterexample()
    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()
    emit(f"source_sha256={source_digest}")
    payload = ("\n".join(LINES) + "\n").encode()
    emit(f"payload_fnv64={fnv64(payload):016x}")
    emit("VERDICT finite-ramification constants and boundary censuses reproduced exactly")


if __name__ == "__main__":
    main()
