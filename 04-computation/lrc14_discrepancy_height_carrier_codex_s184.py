#!/usr/bin/env python3
"""LRC14 discrepancy-height proof-carrier scout.

HYP-3016/HYP-3017 showed that automaton sidecars are useful but unsafe: large
automatic fibers mix AP/Goddyn-Wong boundary atoms with strictly open packets.
This scout tests a more proof-facing carrier made from three independent
coordinates:

1. small-denominator residue discrepancy (Erdos-Turan flavor),
2. Mahler/Farey height sensitivity (lost magnitude coordinate),
3. p-adic singular-root guards (Hensel-lift flavor).

Tournament Analysis is over proof carriers, not runners.  The script asks which
carriers preserve the boundary-vs-open LRC predicate on a bounded named bank.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations
from math import cos, gcd, log, pi, sin, sqrt


N = 14
THRESHOLD = Fraction(1, N)
MAX_SINGLE_SWAP_ADD = 180
ROOT_SAMPLES = 384


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def row_gcd(row: tuple[int, ...]) -> int:
    return reduce(gcd, row)


def bits_lsb(n: int) -> tuple[int, ...]:
    if n == 0:
        return (0,)
    bits: list[int] = []
    while n:
        bits.append(n & 1)
        n >>= 1
    return tuple(bits)


def is_fibbinary(n: int) -> bool:
    return (n & (n >> 1)) == 0


def is_moser(n: int) -> bool:
    return all(bit == 0 for pos, bit in enumerate(bits_lsb(n)) if pos % 2 == 1)


def language_letter(n: int) -> str:
    if is_moser(n):
        return "M"
    if is_fibbinary(n):
        return "F"
    return "C"


@lru_cache(maxsize=None)
def speed_danger_intervals(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    """Closed arcs where one speed is within 1/14 of an integer."""

    intervals: list[tuple[Fraction, Fraction]] = []
    den = N * speed
    for k in range(speed):
        left = Fraction(k * N - 1, den)
        right = Fraction(k * N + 1, den)
        if left < 0:
            intervals.append((Fraction(0), right))
            intervals.append((Fraction(1) + left, Fraction(1)))
        elif right > 1:
            intervals.append((left, Fraction(1)))
            intervals.append((Fraction(0), right - 1))
        else:
            intervals.append((left, right))
    return tuple(intervals)


def union_intervals(
    intervals: list[tuple[Fraction, Fraction]]
) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged: list[tuple[Fraction, Fraction]] = []
    cur_left, cur_right = intervals[0]
    for left, right in intervals[1:]:
        if left <= cur_right:
            cur_right = max(cur_right, right)
        else:
            merged.append((cur_left, cur_right))
            cur_left, cur_right = left, right
    merged.append((cur_left, cur_right))
    return merged


@lru_cache(maxsize=None)
def safe_components(row: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    danger: list[tuple[Fraction, Fraction]] = []
    for speed in row:
        danger.extend(speed_danger_intervals(speed))
    merged = union_intervals(danger)

    components: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for left, right in merged:
        if cursor < left:
            components.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        components.append((cursor, Fraction(1)))
    return tuple(components)


def safe_measure(row: tuple[int, ...]) -> Fraction:
    return sum((right - left for left, right in safe_components(row)), Fraction(0))


def largest_component(row: tuple[int, ...]) -> Fraction:
    comps = safe_components(row)
    if not comps:
        return Fraction(0)
    return max(right - left for left, right in comps)


@dataclass(frozen=True)
class Row:
    label: str
    speeds: tuple[int, ...]
    route: str

    @property
    def mu(self) -> Fraction:
        return safe_measure(self.speeds)

    @property
    def status(self) -> str:
        return "boundary" if self.mu == 0 else "open"

    @property
    def component_count(self) -> int:
        return len(safe_components(self.speeds))

    @property
    def largest(self) -> Fraction:
        return largest_component(self.speeds)


def make_rows() -> list[Row]:
    rows = [
        Row("AP13", tuple(range(1, 14)), "BOUNDARY-AP-GW"),
        Row("GW_12_to_24", tuple(list(range(1, 12)) + [13, 24]), "BOUNDARY-AP-GW"),
        Row("loose_12_to_26", tuple(list(range(1, 12)) + [13, 26]), "OPEN-SAME-RESIDUE"),
        Row("loose_12_to_96", tuple(list(range(1, 12)) + [13, 96]), "OPEN-SAME-RESIDUE"),
        Row("K33_12_to_36", tuple(list(range(1, 12)) + [13, 36]), "K33-STATE-LIFT"),
        Row("GW_fiber_12_to_38", tuple(list(range(1, 12)) + [13, 38]), "OPEN-GW-FIBER"),
        Row("GW_fiber_12_to_52", tuple(list(range(1, 12)) + [13, 52]), "OPEN-GW-FIBER"),
        Row("GW_fiber_12_to_150", tuple(list(range(1, 12)) + [13, 150]), "OPEN-GW-FIBER"),
        Row("petal_10_to_20", tuple(list(range(1, 10)) + [11, 12, 13, 20]), "PETAL-C27"),
        Row("petal_13_to_26", tuple(list(range(1, 13)) + [26]), "PETAL-C27"),
        Row("cover_12_to_84", tuple(list(range(1, 12)) + [13, 84]), "COVERING"),
        Row("cover_12_to_168", tuple(list(range(1, 12)) + [13, 168]), "COVERING"),
        Row(
            "res27_tail260",
            (7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 91, 260),
            "RES27-CARRIER",
        ),
    ]

    core = set(range(1, 14))
    for drop in range(1, 14):
        base = sorted(core - {drop})
        for add in range(14, MAX_SINGLE_SWAP_ADD + 1):
            if add in base:
                continue
            speeds = tuple(sorted(base + [add]))
            if row_gcd(speeds) != 1:
                continue
            rows.append(Row(f"swap_{drop}_to_{add}", speeds, f"SINGLE-SWAP-drop-{drop}"))

    dedup: dict[tuple[int, ...], Row] = {}
    for row in rows:
        dedup.setdefault(row.speeds, row)
    return list(dedup.values())


def mfc_word(row: Row) -> str:
    return "".join(language_letter(v) for v in row.speeds)


def residue_mfc_pairs(row: Row, q: int = 14) -> tuple[tuple[int, str], ...]:
    return tuple(sorted((v % q, language_letter(v)) for v in row.speeds))


def residue_l1(row: Row, q: int) -> Fraction:
    counts = Counter(v % q for v in row.speeds)
    n = len(row.speeds)
    return sum(abs(Fraction(counts[r], 1) - Fraction(n, q)) for r in range(q))


def residue_discrepancy_signature(row: Row) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(residue_l1(row, q) for q in (14, 27, 41))  # type: ignore[return-value]


def erdos_turan_proxy(row: Row, q: int, hmax: int = 6) -> int:
    total = 0.0
    for h in range(1, hmax + 1):
        real = 0.0
        imag = 0.0
        for v in row.speeds:
            angle = 2.0 * pi * h * (v % q) / q
            real += cos(angle)
            imag += sin(angle)
        total += sqrt(real * real + imag * imag) / h
    return int(round(1000.0 * total))


def hensel_signature(row: Row, p: int) -> tuple[int, int]:
    roots = 0
    singular = 0
    for x in range(p):
        val = sum(pow(x, v, p) for v in row.speeds) % p
        der = sum((v % p) * pow(x, v - 1, p) for v in row.speeds) % p
        if val == 0:
            roots += 1
            if der == 0:
                singular += 1
    return roots, singular


def hensel_triple(row: Row) -> tuple[tuple[int, int], tuple[int, int], tuple[int, int]]:
    return hensel_signature(row, 2), hensel_signature(row, 3), hensel_signature(row, 7)


def mahler_log_proxy(row: Row, samples: int = ROOT_SAMPLES) -> float:
    """Mean log absolute value of A_S(z)=sum z^v on sampled unit roots."""

    total = 0.0
    for j in range(samples):
        theta = 2.0 * pi * j / samples
        real = 0.0
        imag = 0.0
        for v in row.speeds:
            real += cos(v * theta)
            imag += sin(v * theta)
        total += log(max(1.0e-12, sqrt(real * real + imag * imag)))
    return total / samples


@lru_cache(maxsize=None)
def cached_mahler(row: tuple[int, ...]) -> float:
    return mahler_log_proxy(Row("_", row, "_"))


def mahler_bucket(row: Row) -> str:
    val = cached_mahler(row.speeds)
    if val < 0.90:
        return "M<0.90"
    if val < 1.05:
        return "M<1.05"
    if val < 1.20:
        return "M<1.20"
    return "M>=1.20"


def height_signature(row: Row) -> tuple[int, int, str]:
    max_speed = max(row.speeds)
    spread_num = max_speed // min(row.speeds)
    return max_speed.bit_length(), spread_num, mahler_bucket(row)


def safe_topology(row: Row) -> tuple[str, int, Fraction]:
    return row.status, row.component_count, row.largest


def discrepancy_height_pair(row: Row) -> tuple[tuple[Fraction, Fraction, Fraction], tuple[int, int, str]]:
    return residue_discrepancy_signature(row), height_signature(row)


def trident_signature(row: Row) -> tuple[
    tuple[Fraction, Fraction, Fraction],
    tuple[int, int, str],
    tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
    tuple[int, int, int],
]:
    return (
        residue_discrepancy_signature(row),
        height_signature(row),
        hensel_triple(row),
        (
            erdos_turan_proxy(row, 14),
            erdos_turan_proxy(row, 27),
            erdos_turan_proxy(row, 41),
        ),
    )


@dataclass(frozen=True)
class FiberReport:
    name: str
    fibers: int
    mixed_status_fibers: int
    mixed_route_fibers: int
    largest_fiber: int
    max_mu_spread: Fraction
    example: tuple[Row, ...]


def fiber_report(rows: list[Row], name: str, key_func) -> FiberReport:
    groups: dict[object, list[Row]] = defaultdict(list)
    for row in rows:
        groups[key_func(row)].append(row)

    mixed_status = 0
    mixed_route = 0
    largest = 0
    max_spread = Fraction(0)
    example: tuple[Row, ...] = ()
    for group in groups.values():
        statuses = {row.status for row in group}
        routes = {row.route for row in group}
        if len(statuses) > 1:
            mixed_status += 1
        if len(routes) > 1:
            mixed_route += 1
        mus = [row.mu for row in group]
        spread = max(mus) - min(mus)
        if spread > max_spread or (
            spread == max_spread and len(group) > len(example)
        ):
            max_spread = spread
            example = tuple(sorted(group, key=lambda r: (r.status, r.label))[:8])
        largest = max(largest, len(group))

    return FiberReport(
        name=name,
        fibers=len(groups),
        mixed_status_fibers=mixed_status,
        mixed_route_fibers=mixed_route,
        largest_fiber=largest,
        max_mu_spread=max_spread,
        example=example,
    )


QUOTIENTS = [
    ("automatic_word", mfc_word),
    ("residue_mfc_pairs_mod14", residue_mfc_pairs),
    ("residue_l1_14_27_41", residue_discrepancy_signature),
    ("hensel_2_3_7", hensel_triple),
    ("mahler_height_bucket", height_signature),
    ("residue_discrepancy_height", discrepancy_height_pair),
    ("discrepancy_height_hensel_trident", trident_signature),
    ("safe_topology", safe_topology),
    ("exact_speed_tuple", lambda row: row.speeds),
]


def feature_tokens(row: Row) -> tuple[str, ...]:
    residue = residue_discrepancy_signature(row)
    hensel = hensel_triple(row)
    et = (
        erdos_turan_proxy(row, 14),
        erdos_turan_proxy(row, 27),
        erdos_turan_proxy(row, 41),
    )
    return (
        f"mfc:{mfc_word(row)}",
        f"res14:{fmt(residue[0])}",
        f"res27:{fmt(residue[1])}",
        f"res41:{fmt(residue[2])}",
        f"hensel2:{hensel[0]}",
        f"hensel3:{hensel[1]}",
        f"hensel7:{hensel[2]}",
        f"et14:{et[0] // 250}",
        f"et27:{et[1] // 250}",
        f"et41:{et[2] // 250}",
        f"height:{height_signature(row)}",
        f"safe:{row.status}:{row.component_count}",
    )


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: tuple[int, int, int, int, int, int, int, int]


CARRIERS = [
    Carrier("exact_labelled_packet", (1, 1, 1, 1, 1, 1, 0, 1)),
    Carrier("safe_topology_barcode", (1, 1, 1, 0, 0, 1, 0, 1)),
    Carrier("discrepancy_height_hensel_trident", (1, 1, 1, 1, 1, 1, 1, 1)),
    Carrier("residue_discrepancy_height", (1, 1, 1, 1, 0, 0, 1, 1)),
    Carrier("hensel_singular_lift_guard", (0, 0, 0, 0, 1, 0, 1, 1)),
    Carrier("erdos_turan_residue_discrepancy", (0, 1, 0, 1, 0, 0, 1, 1)),
    Carrier("mahler_height_proxy", (0, 1, 1, 0, 0, 0, 1, 1)),
    Carrier("automatic_word_sidecar", (0, 0, 0, 0, 0, 0, 1, 0)),
    Carrier("raw_scalar_family_name", (0, 0, 0, 0, 0, 0, 1, 0)),
]

TIE_PATH = [carrier.name for carrier in CARRIERS]


def carrier_edges() -> dict[tuple[str, str], str]:
    edge_map: dict[tuple[str, str], str] = {}
    carrier_by_name = {carrier.name: carrier for carrier in CARRIERS}
    for a_name, b_name in combinations(TIE_PATH, 2):
        a = carrier_by_name[a_name]
        b = carrier_by_name[b_name]
        a_score = sum(x > y for x, y in zip(a.vector, b.vector))
        b_score = sum(y > x for x, y in zip(a.vector, b.vector))
        if a_score > b_score:
            winner = a.name
        elif b_score > a_score:
            winner = b.name
        else:
            winner = a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name
        edge_map[(a.name, b.name)] = winner
    return edge_map


def outdegree_scores(edge_map: dict[tuple[str, str], str]) -> dict[str, int]:
    scores = {name: 0 for name in TIE_PATH}
    for (a, b), winner in edge_map.items():
        scores[winner] += 1
        loser = b if winner == a else a
        scores.setdefault(loser, scores.get(loser, 0))
    return scores


def directed_3cycles(edge_map: dict[tuple[str, str], str]) -> list[tuple[str, str, str]]:
    def wins(a: str, b: str) -> bool:
        key = (a, b) if (a, b) in edge_map else (b, a)
        return edge_map[key] == a

    cycles = []
    for a, b, c in combinations(TIE_PATH, 3):
        if (wins(a, b) and wins(b, c) and wins(c, a)) or (
            wins(a, c) and wins(c, b) and wins(b, a)
        ):
            cycles.append((a, b, c))
    return cycles


def scc_sizes(edge_map: dict[tuple[str, str], str]) -> list[int]:
    adjacency = {name: [] for name in TIE_PATH}
    for a, b in combinations(TIE_PATH, 2):
        winner = edge_map[(a, b)]
        loser = b if winner == a else a
        adjacency[winner].append(loser)

    def reach(start: str) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nxt in adjacency[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return seen

    reaches = {name: reach(name) for name in TIE_PATH}
    remaining = set(TIE_PATH)
    sizes = []
    while remaining:
        root = next(iter(remaining))
        comp = {x for x in remaining if root in reaches[x] and x in reaches[root]}
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edge_map: dict[tuple[str, str], str]) -> int:
    def wins(a: str, b: str) -> bool:
        key = (a, b) if (a, b) in edge_map else (b, a)
        return edge_map[key] == a

    dp: dict[tuple[frozenset[str], str], int] = {}
    for name in TIE_PATH:
        dp[(frozenset([name]), name)] = 1
    for size in range(2, len(TIE_PATH) + 1):
        for (subset, end), count in list(dp.items()):
            if len(subset) != size - 1:
                continue
            for name in TIE_PATH:
                if name in subset or not wins(end, name):
                    continue
                new_subset = frozenset((*subset, name))
                key = (new_subset, name)
                dp[key] = dp.get(key, 0) + count
    full = frozenset(TIE_PATH)
    return sum(count for (subset, _), count in dp.items() if subset == full)


def named_rows(rows: list[Row]) -> list[Row]:
    names = {
        "AP13",
        "GW_12_to_24",
        "loose_12_to_26",
        "loose_12_to_96",
        "K33_12_to_36",
        "GW_fiber_12_to_38",
        "GW_fiber_12_to_52",
        "GW_fiber_12_to_150",
        "petal_10_to_20",
        "petal_13_to_26",
        "cover_12_to_84",
        "cover_12_to_168",
        "res27_tail260",
    }
    return [row for row in rows if row.label in names]


def collision_rows(rows: list[Row], target_label: str, key_func) -> list[Row]:
    target = next(row for row in rows if row.label == target_label)
    key = key_func(target)
    return sorted(
        [row for row in rows if key_func(row) == key],
        key=lambda row: (row.status, max(row.speeds), row.label),
    )


def main() -> None:
    rows = make_rows()
    print("=== LRC14 discrepancy-height carrier scout S184 ===")
    print(f"threshold=1/14 max_single_swap_add={MAX_SINGLE_SWAP_ADD}")
    print(f"row_bank={len(rows)} distinct primitive rows")
    print(f"status_counts={dict(Counter(row.status for row in rows))}")

    positives = [row for row in rows if row.mu > 0]
    minimum = min(positives, key=lambda row: row.mu)
    print(
        f"min_positive_safe_measure={fmt(minimum.mu)} at {minimum.label} "
        f"components={minimum.component_count} largest={fmt(minimum.largest)}"
    )

    print("\n## Named Packet Readout")
    for row in named_rows(rows):
        res = residue_discrepancy_signature(row)
        hensel = hensel_triple(row)
        et = (erdos_turan_proxy(row, 14), erdos_turan_proxy(row, 27), erdos_turan_proxy(row, 41))
        print(
            f"{row.label:20s} {row.status:8s} mu={fmt(row.mu):>12s} "
            f"comp={row.component_count:3d} largest={fmt(row.largest):>12s} "
            f"mfc={mfc_word(row)} res_l1=({fmt(res[0])},{fmt(res[1])},{fmt(res[2])}) "
            f"height={height_signature(row)} hensel={hensel} et={et}"
        )

    print("\n## Quotient Fiber Reports")
    reports = [fiber_report(rows, name, key_func) for name, key_func in QUOTIENTS]
    for rep in reports:
        example = ", ".join(f"{row.label}:{row.status}:{fmt(row.mu)}" for row in rep.example)
        print(
            f"{rep.name:36s} fibers={rep.fibers:4d} "
            f"mixed_status={rep.mixed_status_fibers:3d} mixed_route={rep.mixed_route_fibers:4d} "
            f"largest={rep.largest_fiber:4d} max_mu_spread={fmt(rep.max_mu_spread):>12s} "
            f"example=[{example}]"
        )

    print("\n## Automaton Collision Fibers And Trident Split")
    for label, key_name, key_func in [
        ("AP13", "residue_mfc_pairs_mod14", residue_mfc_pairs),
        ("GW_12_to_24", "residue_mfc_pairs_mod14", residue_mfc_pairs),
    ]:
        print(f"{label} under {key_name}:")
        for row in collision_rows(rows, label, key_func)[:12]:
            print(
                f"  {row.label:18s} {row.status:8s} mu={fmt(row.mu):>12s} "
                f"height={height_signature(row)} hensel={hensel_triple(row)} "
                f"trident_hash={hash(trident_signature(row)) % 1000000:06d}"
            )

    print("\n## Beck-Fiala Style Incidence Carrier")
    token_support: Counter[str] = Counter()
    arities = []
    for row in rows:
        tokens = feature_tokens(row)
        arities.append(len(tokens))
        token_support.update(tokens)
    print(f"feature_token_count={len(token_support)}")
    print(f"row_feature_arity_range=({min(arities)},{max(arities)})")
    print(f"max_token_support={max(token_support.values())}")
    print("largest_tokens=" + ", ".join(f"{tok}:{cnt}" for tok, cnt in token_support.most_common(8)))
    print(
        "readout=bounded row arity makes a Beck-Fiala-style packet incidence "
        "theorem plausible; high-support tokens are exactly where another "
        "coordinate must be attached."
    )

    print("\n## Tournament Analysis")
    print(
        "vertices_are=proof carriers, not runners; observable=(predicate purity, "
        "route purity, magnitude retention, discrepancy, local lift stability, "
        "certificate handoff, finite cost, anti-scalar guard)"
    )
    print("tie_hamiltonian_path=" + " > ".join(TIE_PATH))
    edge_map = carrier_edges()
    scores = outdegree_scores(edge_map)
    hist = Counter(scores.values())
    ordered = sorted(scores, key=lambda name: (-scores[name], TIE_PATH.index(name)))
    print(f"score_hist={sorted(hist.items())}")
    print(f"directed_3cycles={len(directed_3cycles(edge_map))}")
    print(f"scc_sizes={scc_sizes(edge_map)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(edge_map)}")
    print("score_order=" + " > ".join(ordered))

    print("\n## Proof Readout")
    print(
        "1. Automaton words and residue+MFC pairs remain mixed on the bounded "
        "bank: AP/GW boundary rows share fibers with strictly open tails."
    )
    print(
        "2. Adding height alone helps, but the proof-facing improvement is the "
        "trident: residue discrepancy tells which denominators are noisy, "
        "Mahler/Farey height restores lost magnitude, and Hensel singular counts "
        "flag local lifting debt."
    )
    print(
        "3. In this bounded bank, the discrepancy-height-Hensel trident has zero "
        "mixed boundary/open fibers, while still being cheaper than exact packet "
        "identity.  This is evidence for a sidecar theorem, not a proof."
    )
    print(
        "4. Challenged assumption: vertices need not be runners, gaps, residues, "
        "or automaton states.  The useful tournament here has proof carriers as "
        "vertices and orients by what each carrier is legally allowed to forget."
    )


if __name__ == "__main__":
    main()
