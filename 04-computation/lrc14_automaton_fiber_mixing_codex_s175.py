#!/usr/bin/env python3
"""LRC14 automaton quotient fiber-mixing stress test.

This extends the Moser/fibbinary/lacunary carrier work by asking a harder
question: which finite-state quotients mix exact LRC14 boundary atoms with
strictly open rows?

The row bank is deliberately modest and exact.  It contains named LRC14 rows
and the AP single-swap atlas {1,...,13}\{drop} U {add} through a bounded tail.
For each row we compute the threshold-1/14 safe components by rational interval
union, then group rows by sequence/automaton quotients.

Tournament Analysis is over quotient/proof carriers, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd


N = 14
THRESHOLD = Fraction(1, N)
MAX_SINGLE_SWAP_ADD = 180


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def bits_lsb(n: int) -> tuple[int, ...]:
    if n == 0:
        return (0,)
    out: list[int] = []
    while n:
        out.append(n & 1)
        n >>= 1
    return tuple(out)


def bit_length(n: int) -> int:
    return max(1, n.bit_length())


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


def terminal_state(n: int) -> str:
    """A deliberately coarse terminal-state quotient.

    It keeps only accepted/dead status for the two small automata.  This is the
    kind of quotient that looks useful but may forget the magnitude coordinate.
    """

    f = "Fok" if is_fibbinary(n) else "Fdead"
    m = "Mok" if is_moser(n) else "Mdead"
    return f"{f}/{m}"


def perfect_power_tag(n: int) -> str:
    if n <= 1:
        return "unit"
    max_exp = n.bit_length()
    for exp in range(max_exp, 1, -1):
        lo, hi = 1, int(round(n ** (1.0 / exp))) + 3
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            if mid**exp <= n:
                lo = mid
            else:
                hi = mid
        if lo > 1 and lo**exp == n:
            return f"{exp}pow"
    return "npow"


@lru_cache(maxsize=None)
def speed_danger_intervals(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    """Closed danger arcs for one speed, split inside [0,1]."""

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


def union_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    out: list[tuple[Fraction, Fraction]] = []
    cur_left, cur_right = intervals[0]
    for left, right in intervals[1:]:
        if left <= cur_right:
            cur_right = max(cur_right, right)
        else:
            out.append((cur_left, cur_right))
            cur_left, cur_right = left, right
    out.append((cur_left, cur_right))
    return out


@lru_cache(maxsize=None)
def safe_components(row: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for speed in row:
        intervals.extend(speed_danger_intervals(speed))
    danger = union_intervals(intervals)
    components: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for left, right in danger:
        if cursor < left:
            components.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        components.append((cursor, Fraction(1)))
    return tuple(components)


def safe_measure(row: tuple[int, ...]) -> Fraction:
    return sum((right - left for left, right in safe_components(row)), Fraction(0))


def dist_to_integer(x: Fraction) -> Fraction:
    y = x % 1
    return min(y, Fraction(1) - y)


def min_distance(row: tuple[int, ...], t: Fraction) -> Fraction:
    return min(dist_to_integer(speed * t) for speed in row)


def boundary_units(row: tuple[int, ...]) -> tuple[int, ...]:
    out: list[int] = []
    for a in (1, 3, 5, 9, 11, 13):
        if min_distance(row, Fraction(a, N)) >= THRESHOLD:
            out.append(a)
    return tuple(out)


def max_adjacent_ratio(row: tuple[int, ...]) -> Fraction:
    vals = sorted(row)
    return max(Fraction(b, a) for a, b in zip(vals, vals[1:]))


def ratio_bucket(row: tuple[int, ...]) -> str:
    r = max_adjacent_ratio(row)
    if r <= 2:
        return "<=2"
    if r <= 3:
        return "(2,3]"
    return ">3"


@dataclass(frozen=True)
class Row:
    label: str
    speeds: tuple[int, ...]
    family: str

    @property
    def mu(self) -> Fraction:
        return safe_measure(self.speeds)

    @property
    def status(self) -> str:
        return "boundary" if self.mu == 0 else "open"

    @property
    def components(self) -> int:
        return len(safe_components(self.speeds))

    @property
    def largest_component(self) -> Fraction:
        comps = safe_components(self.speeds)
        if not comps:
            return Fraction(0)
        return max(right - left for left, right in comps)


def make_rows() -> list[Row]:
    rows = [
        Row("AP13", tuple(range(1, 14)), "named"),
        Row("GW_12_to_24", tuple(list(range(1, 12)) + [13, 24]), "named"),
        Row("K33_12_to_36", tuple(list(range(1, 12)) + [13, 36]), "named"),
        Row("loose_12_to_26", tuple(list(range(1, 12)) + [13, 26]), "named"),
        Row("loose_12_to_96", tuple(list(range(1, 12)) + [13, 96]), "named"),
        Row("petal_10_to_20", tuple(list(range(1, 10)) + [11, 12, 13, 20]), "named"),
        Row("petal_13_to_26", tuple(list(range(1, 13)) + [26]), "named"),
        Row("cover_12_to_84", tuple(list(range(1, 12)) + [13, 84]), "named"),
    ]
    core = set(range(1, 14))
    for drop in range(1, 14):
        base = sorted(core - {drop})
        for add in range(14, MAX_SINGLE_SWAP_ADD + 1):
            if add in base:
                continue
            speeds = tuple(sorted(base + [add]))
            rows.append(Row(f"swap_{drop}_to_{add}", speeds, f"single_swap_drop_{drop}"))

    dedup: dict[tuple[int, ...], Row] = {}
    for row in rows:
        if gcd(*row.speeds) != 1:
            continue
        dedup.setdefault(row.speeds, row)
    return list(dedup.values())


def key_mfc_counts(row: Row) -> tuple[int, int, int]:
    c = Counter(language_letter(v) for v in row.speeds)
    return c["M"], c["F"], c["C"]


def key_mfc_word(row: Row) -> str:
    return "".join(language_letter(v) for v in row.speeds)


def key_terminal_word(row: Row) -> tuple[str, ...]:
    return tuple(terminal_state(v) for v in row.speeds)


def key_residue_multiset(row: Row) -> tuple[int, ...]:
    return tuple(sorted(v % N for v in row.speeds))


def key_residue_mfc_pairs(row: Row) -> tuple[tuple[int, str], ...]:
    return tuple(sorted((v % N, language_letter(v)) for v in row.speeds))


def key_residue_terminal_pairs(row: Row) -> tuple[tuple[int, str], ...]:
    return tuple(sorted((v % N, terminal_state(v)) for v in row.speeds))


def key_residue_mfc_bitlen_pairs(row: Row) -> tuple[tuple[int, str, int], ...]:
    return tuple(sorted((v % N, language_letter(v), bit_length(v)) for v in row.speeds))


def key_perfect_power_word(row: Row) -> tuple[str, ...]:
    return tuple(perfect_power_tag(v) for v in row.speeds)


def key_gap_bucket(row: Row) -> str:
    return ratio_bucket(row.speeds)


def key_exact(row: Row) -> tuple[int, ...]:
    return row.speeds


QUOTIENTS = (
    ("mfc_counts", key_mfc_counts),
    ("mfc_word", key_mfc_word),
    ("terminal_word", key_terminal_word),
    ("residue_multiset_mod14", key_residue_multiset),
    ("residue_mfc_pairs", key_residue_mfc_pairs),
    ("residue_terminal_pairs", key_residue_terminal_pairs),
    ("residue_mfc_bitlen_pairs", key_residue_mfc_bitlen_pairs),
    ("perfect_power_word", key_perfect_power_word),
    ("gap_ratio_bucket", key_gap_bucket),
    ("exact_speed_tuple", key_exact),
)


@dataclass
class FiberReport:
    name: str
    fibers: int
    mixed_status_fibers: int
    max_measure_spread: Fraction
    max_component_spread: int
    example_labels: tuple[str, ...]
    example_measures: tuple[Fraction, ...]
    example_statuses: tuple[str, ...]


def fiber_report(rows: list[Row], name: str, key_func) -> FiberReport:
    groups: dict[object, list[Row]] = defaultdict(list)
    for row in rows:
        groups[key_func(row)].append(row)

    mixed = 0
    max_spread = Fraction(-1)
    max_comp_spread = -1
    best: list[Row] = []
    for group in groups.values():
        statuses = {row.status for row in group}
        measures = [row.mu for row in group]
        components = [row.components for row in group]
        spread = max(measures) - min(measures)
        comp_spread = max(components) - min(components)
        if len(statuses) > 1:
            mixed += 1
        if (len(statuses) > 1, spread, comp_spread, len(group)) > (
            bool(best and len({r.status for r in best}) > 1),
            max_spread,
            max_comp_spread,
            len(best),
        ):
            max_spread = spread
            max_comp_spread = comp_spread
            best = group

    ordered = sorted(best, key=lambda r: (r.status, r.mu, r.label))
    sample = ordered[:3] + ordered[-3:] if len(ordered) > 6 else ordered
    return FiberReport(
        name=name,
        fibers=len(groups),
        mixed_status_fibers=mixed,
        max_measure_spread=max_spread,
        max_component_spread=max_comp_spread,
        example_labels=tuple(row.label for row in sample),
        example_measures=tuple(row.mu for row in sample),
        example_statuses=tuple(row.status for row in sample),
    )


def named_row_table(rows: list[Row]) -> list[Row]:
    wanted = {
        "AP13",
        "GW_12_to_24",
        "K33_12_to_36",
        "loose_12_to_26",
        "loose_12_to_96",
        "petal_10_to_20",
        "petal_13_to_26",
        "cover_12_to_84",
    }
    return [row for row in rows if row.label in wanted]


def collision_family(rows: list[Row], target: Row, key_func, limit: int = 12) -> list[Row]:
    key = key_func(target)
    hits = [row for row in rows if key_func(row) == key]
    hits.sort(key=lambda row: (row.mu, max(row.speeds), row.label))
    return hits[:limit]


@dataclass(frozen=True)
class Carrier:
    name: str
    fiber_purity: int
    finite_state: int
    keeps_magnitude: int
    keeps_residue_endpoint: int
    route_compatibility: int
    anti_scalar_guard: int
    proof_cost: int

    def vector(self) -> tuple[int, ...]:
        # proof_cost is inverted below so larger is better.
        return (
            self.fiber_purity,
            self.finite_state,
            self.keeps_magnitude,
            self.keeps_residue_endpoint,
            self.route_compatibility,
            self.anti_scalar_guard,
            6 - self.proof_cost,
        )


CARRIERS = (
    Carrier("exact_labelled_packet", 5, 2, 5, 5, 5, 5, 5),
    Carrier("farey_magnitude_height", 4, 2, 5, 3, 5, 5, 4),
    Carrier("residue_mfc_bitlen_pairs", 3, 4, 3, 5, 3, 4, 3),
    Carrier("residue_terminal_pairs", 2, 5, 1, 5, 2, 4, 2),
    Carrier("residue_mfc_pairs", 2, 5, 1, 5, 2, 4, 2),
    Carrier("terminal_word", 1, 5, 0, 1, 1, 3, 1),
    Carrier("mfc_word", 1, 5, 0, 1, 1, 3, 1),
    Carrier("perfect_power_word", 1, 2, 1, 1, 2, 2, 1),
    Carrier("gap_ratio_bucket", 1, 1, 2, 1, 2, 2, 1),
    Carrier("raw_counts", 0, 3, 0, 0, 0, 1, 1),
)

TIE_PATH = tuple(c.name for c in CARRIERS)


def orient(a: Carrier, b: Carrier) -> str:
    av = a.vector()
    bv = b.vector()
    awins = sum(x > y for x, y in zip(av, bv))
    bwins = sum(x < y for x, y in zip(av, bv))
    if awins != bwins:
        return a.name if awins > bwins else b.name
    if av != bv:
        return a.name if av > bv else b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def edges() -> dict[tuple[str, str], str]:
    out: dict[tuple[str, str], str] = {}
    for a, b in combinations(CARRIERS, 2):
        out[(a.name, b.name)] = orient(a, b)
    return out


def outdegree_scores(edge_map: dict[tuple[str, str], str]) -> dict[str, int]:
    scores = {c.name: 0 for c in CARRIERS}
    for winner in edge_map.values():
        scores[winner] += 1
    return scores


def directed_3cycles(edge_map: dict[tuple[str, str], str]) -> list[tuple[str, str, str]]:
    names = [c.name for c in CARRIERS]

    def wins(x: str, y: str) -> bool:
        key = (x, y) if (x, y) in edge_map else (y, x)
        return edge_map[key] == x

    cycles = []
    for a, b, c in combinations(names, 3):
        if (wins(a, b) and wins(b, c) and wins(c, a)) or (
            wins(a, c) and wins(c, b) and wins(b, a)
        ):
            cycles.append((a, b, c))
    return cycles


def scc_sizes(edge_map: dict[tuple[str, str], str]) -> list[int]:
    names = [c.name for c in CARRIERS]
    adj = {name: [] for name in names}
    for a, b in combinations(names, 2):
        winner = edge_map[(a, b)]
        loser = b if winner == a else a
        adj[winner].append(loser)

    def reach(start: str) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        return seen

    remaining = set(names)
    sizes = []
    reaches = {name: reach(name) for name in names}
    while remaining:
        root = next(iter(remaining))
        comp = {x for x in remaining if root in reaches[x] and x in reaches[root]}
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edge_map: dict[tuple[str, str], str]) -> int:
    names = tuple(c.name for c in CARRIERS)

    def wins(x: str, y: str) -> bool:
        key = (x, y) if (x, y) in edge_map else (y, x)
        return edge_map[key] == x

    dp: dict[tuple[frozenset[str], str], int] = {}
    for name in names:
        dp[(frozenset([name]), name)] = 1
    for size in range(2, len(names) + 1):
        next_dp: dict[tuple[frozenset[str], str], int] = {}
        for subset_end, count in dp.items():
            subset, end = subset_end
            if len(subset) != size - 1:
                continue
            for name in names:
                if name in subset or not wins(end, name):
                    continue
                new_subset = frozenset((*subset, name))
                next_dp[(new_subset, name)] = next_dp.get((new_subset, name), 0) + count
        dp.update(next_dp)
    full = frozenset(names)
    return sum(count for (subset, _), count in dp.items() if subset == full)


def main() -> None:
    rows = make_rows()
    print("=== LRC14 automaton quotient fiber-mixing scout S175 ===")
    print(f"threshold=1/14 max_single_swap_add={MAX_SINGLE_SWAP_ADD}")
    print(f"row_bank={len(rows)} distinct primitive rows")
    status_counts = Counter(row.status for row in rows)
    print(f"status_counts={dict(status_counts)}")

    positives = [row for row in rows if row.mu > 0]
    min_positive = min(positives, key=lambda row: row.mu)
    print(
        "min_positive_safe_measure="
        f"{fmt(min_positive.mu)} at {min_positive.label} "
        f"components={min_positive.components} largest={fmt(min_positive.largest_component)}"
    )

    print("\n## Named Rows")
    for row in named_row_table(rows):
        print(
            f"{row.label:18s} status={row.status:8s} mu={fmt(row.mu):>12s} "
            f"components={row.components:3d} largest={fmt(row.largest_component):>12s} "
            f"boundary_units={boundary_units(row.speeds)} "
            f"mfc={key_mfc_word(row)} residues={key_residue_multiset(row)}"
        )

    print("\n## Quotient Fiber-Mixing Reports")
    reports = [fiber_report(rows, name, key_func) for name, key_func in QUOTIENTS]
    for rep in reports:
        ex = ", ".join(
            f"{label}:{status}:{fmt(mu)}"
            for label, status, mu in zip(rep.example_labels, rep.example_statuses, rep.example_measures)
        )
        print(
            f"{rep.name:28s} fibers={rep.fibers:4d} "
            f"mixed_boundary_open={rep.mixed_status_fibers:4d} "
            f"max_mu_spread={fmt(rep.max_measure_spread):>12s} "
            f"max_component_spread={rep.max_component_spread:3d} "
            f"example=[{ex}]"
        )

    row_by_label = {row.label: row for row in rows}
    print("\n## Collision Fibers")
    for target_label, key_name, key_func in [
        ("AP13", "residue_mfc_pairs", key_residue_mfc_pairs),
        ("GW_12_to_24", "residue_mfc_pairs", key_residue_mfc_pairs),
        ("AP13", "residue_terminal_pairs", key_residue_terminal_pairs),
    ]:
        target = row_by_label[target_label]
        print(f"{target_label} under {key_name}:")
        for row in collision_family(rows, target, key_func, limit=14):
            tail = max(row.speeds)
            print(
                f"  {row.label:16s} tail={tail:3d} status={row.status:8s} "
                f"mu={fmt(row.mu):>12s} components={row.components:3d} "
                f"largest={fmt(row.largest_component):>12s}"
            )

    print("\n## Automaton Closure Stress Inside Single-Swap Bank")
    unit_excess_hits = []
    for p in range(1, 128):
        q = 14 * p - 1
        tags = []
        if is_fibbinary(q):
            tags.append("F")
        if is_moser(q):
            tags.append("M")
        if tags:
            row = Row(f"unit_excess_q{q}", tuple(list(range(1, 12)) + [13, q]), "unit_excess")
            unit_excess_hits.append((p, q, "".join(tags), row.mu, row.status))
    print(f"unit_excess_automaton_hits_p<128={len(unit_excess_hits)}")
    for p, q, tags, mu, status in unit_excess_hits[:14]:
        print(f"  p={p:3d} q={q:5d} tags={tags:2s} status={status:8s} mu={fmt(mu)}")

    print("\n## Tournament Analysis")
    print(
        "vertices_are=quotient/proof carriers, not runners; "
        "observable=fiber purity, finite-state checkability, magnitude retention, "
        "residue/endpoint retention, route compatibility, anti-scalar guard, proof cost"
    )
    print("tie_hamiltonian_path=" + " > ".join(TIE_PATH))
    edge_map = edges()
    scores = outdegree_scores(edge_map)
    hist = Counter(scores.values())
    cycles = directed_3cycles(edge_map)
    ordered = sorted(scores, key=lambda name: (-scores[name], TIE_PATH.index(name)))
    print(f"score_hist={sorted(hist.items())}")
    print(f"directed_3cycles={len(cycles)}")
    print(f"scc_sizes={scc_sizes(edge_map)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(edge_map)}")
    print("score_order=" + " > ".join(ordered))

    print("\n## Proof Readout")
    print(
        "1. The AP residue+automaton fiber is mixed: AP13 is boundary-only, "
        "while residue- and automaton-indistinguishable tails such as 12->26 "
        "and 12->96 are strictly open.  Thus mod-14 residue plus Moser/"
        "fibbinary terminal state cannot prove tightness."
    )
    print(
        "2. The GW one-dipole residue+automaton fiber is also mixed by later "
        "tails.  The Jacobsthal/doubling site 12 is visible only after the "
        "magnitude/Farey height coordinate is retained."
    )
    print(
        "3. Adding crude bit-length reduces some collisions but is not a proof "
        "invariant; it is a proxy for the real coordinate: exact speed/Farey "
        "scale and safe-component endpoint ownership."
    )
    print(
        "4. A theorem-safe packet quotient for HYP-2963 should therefore carry "
        "automaton state as a side channel, but discharge or retain the "
        "magnitude cocycle measuring how far a tail has moved inside a fixed "
        "residue-language fiber."
    )


if __name__ == "__main__":
    main()
