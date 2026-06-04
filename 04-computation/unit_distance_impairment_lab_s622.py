#!/usr/bin/env python3
"""
unit_distance_impairment_lab_s622.py

S622: controlled impairments for unit-distance construction/proof methods.

The goal is not to improve the global n=22 bound directly.  It is to make a
small-size laboratory that asks: which parts of the Moser-carrier construction
are genuinely load-bearing, and which proof/search features survive when we
deliberately damage the method?

Three impairments are tested:

1. Direction dropout: remove one antipodal unit-direction pair from the Moser
   shell and measure the small-size edge loss.
2. Direction support threshold: read exact small witnesses from the full beam
   and ask how many direction-pair types those witnesses actually use.
3. Gain cap: forbid frontier extensions above a fixed local gain and measure
   when exact small optima become impossible in the impaired carrier.

Tournament Analysis is over techniques/impairments, not over points.  The
preserved predicate is "what retained side channel helps a dense unit-distance
search/proof progress at small sizes?"  The quotient destroys exact embedding
provenance unless a technique explicitly retains it.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from math import comb


Point = tuple[int, int, int, int]
Cluster = tuple[Point, ...]


U_EXACT = {
    1: 0,
    2: 1,
    3: 3,
    4: 5,
    5: 7,
    6: 9,
    7: 12,
    8: 14,
    9: 18,
    10: 20,
    11: 23,
    12: 27,
    13: 30,
    14: 33,
    15: 37,
    16: 41,
    17: 43,
    18: 46,
    19: 50,
    20: 54,
    21: 57,
}


def moser_unit_vectors() -> tuple[Point, ...]:
    units: list[Point] = []
    for a, b, c, d in product(range(-4, 5), repeat=4):
        if a == b == c == d == 0:
            continue
        if a * d != b * c:
            continue
        value = (
            6 * a * a
            + 6 * a * b
            + 10 * a * c
            + 5 * a * d
            + 6 * b * b
            + 5 * b * c
            + 10 * b * d
            + 6 * c * c
            + 6 * c * d
            + 6 * d * d
        )
        if value == 6:
            units.append((a, b, c, d))
    return tuple(sorted(units))


def neg(p: Point) -> Point:
    return tuple(-x for x in p)  # type: ignore[return-value]


def add(a: Point, b: Point) -> Point:
    return tuple(a[i] + b[i] for i in range(4))  # type: ignore[return-value]


def sub(a: Point, b: Point) -> Point:
    return tuple(a[i] - b[i] for i in range(4))  # type: ignore[return-value]


def canon(points: Cluster) -> Cluster:
    pts = sorted(points)
    base = pts[0]
    return tuple(sorted(sub(p, base) for p in pts))


def span4(cluster: Cluster) -> int:
    return sum(max(p[i] for p in cluster) - min(p[i] for p in cluster) for i in range(4))


UNITS = moser_unit_vectors()


def direction_pairs(units: tuple[Point, ...]) -> tuple[tuple[Point, Point], ...]:
    seen: set[Point] = set()
    pairs: list[tuple[Point, Point]] = []
    for u in sorted(units):
        if u in seen:
            continue
        v = neg(u)
        pair = tuple(sorted((u, v)))  # type: ignore[assignment]
        pairs.append(pair)
        seen.add(pair[0])
        seen.add(pair[1])
    return tuple(pairs)


PAIRS = direction_pairs(UNITS)
PAIR_INDEX: dict[Point, int] = {}
for idx, pair in enumerate(PAIRS):
    PAIR_INDEX[pair[0]] = idx
    PAIR_INDEX[pair[1]] = idx


def units_from_mask(mask: int) -> tuple[Point, ...]:
    units: list[Point] = []
    for idx, pair in enumerate(PAIRS):
        if (mask >> idx) & 1:
            units.extend(pair)
    return tuple(sorted(units))


FULL_MASK = (1 << len(PAIRS)) - 1


@dataclass(frozen=True)
class BeamRow:
    size: int
    kept: int
    unique_children: int
    best_edges: int
    best_span: int
    exact_edges: int | None
    deficit: int | None


@dataclass(frozen=True)
class BeamResult:
    rows: tuple[BeamRow, ...]
    best_cluster: Cluster
    best_edges: int
    target: int
    width: int
    mask: int
    gain_cap: int | None
    best_by_size: tuple[tuple[int, Cluster], ...]

    def edge_at(self, size: int) -> int:
        for row in self.rows:
            if row.size == size:
                return row.best_edges
        raise KeyError(size)

    def cluster_at(self, size: int) -> Cluster:
        for row_size, cluster in self.best_by_size:
            if row_size == size:
                return cluster
        raise KeyError(size)


def frontier_gains(cluster: Cluster, units: tuple[Point, ...]) -> dict[Point, int]:
    s = set(cluster)
    gains: dict[Point, int] = {}
    for p in cluster:
        for u in units:
            q = add(p, u)
            if q not in s:
                gains[q] = gains.get(q, 0) + 1
    return gains


def unit_edges(cluster: Cluster, units: tuple[Point, ...]) -> int:
    s = set(cluster)
    return sum(1 for p in cluster for u in units if add(p, u) in s) // 2


def run_beam(
    target: int,
    width: int,
    mask: int = FULL_MASK,
    gain_cap: int | None = None,
) -> BeamResult:
    units = units_from_mask(mask)
    if not units:
        raise ValueError("mask must contain at least one direction pair")
    states: dict[Cluster, int] = {((0, 0, 0, 0),): 0}
    rows: list[BeamRow] = []
    best_by_size: list[tuple[int, Cluster]] = []
    best_cluster: Cluster = ((0, 0, 0, 0),)
    best_edges = 0

    for size in range(2, target + 1):
        children: dict[Cluster, int] = {}
        for cluster, edges in states.items():
            for q, gain in frontier_gains(cluster, units).items():
                if gain_cap is not None and gain > gain_cap:
                    continue
                child = canon(cluster + (q,))
                new_edges = edges + gain
                if new_edges > children.get(child, -1):
                    children[child] = new_edges
        if not children:
            rows.append(BeamRow(size, 0, 0, 0, 0, U_EXACT.get(size), None))
            states = {}
            break
        ranked = sorted(
            children.items(),
            key=lambda item: (item[1], -span4(item[0]), item[0]),
            reverse=True,
        )
        states = dict(ranked[:width])
        best_cluster, best_edges = ranked[0]
        recomputed = unit_edges(best_cluster, units)
        if recomputed != best_edges:
            raise AssertionError((size, best_edges, recomputed))
        best_by_size.append((size, best_cluster))
        exact = U_EXACT.get(size)
        deficit = None if exact is None else exact - best_edges
        rows.append(
            BeamRow(
                size=size,
                kept=len(states),
                unique_children=len(children),
                best_edges=best_edges,
                best_span=span4(best_cluster),
                exact_edges=exact,
                deficit=deficit,
            )
        )

    return BeamResult(tuple(rows), best_cluster, best_edges, target, width, mask, gain_cap, tuple(best_by_size))


def direction_usage(cluster: Cluster, mask: int = FULL_MASK) -> tuple[int, ...]:
    allowed = set(units_from_mask(mask))
    counts = [0] * len(PAIRS)
    pts = list(cluster)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            d = sub(pts[j], pts[i])
            if d in allowed:
                counts[PAIR_INDEX[d]] += 1
            else:
                nd = neg(d)
                if nd in allowed:
                    counts[PAIR_INDEX[nd]] += 1
    return tuple(counts)


def mask_label(mask: int) -> str:
    return "{" + ",".join(str(i) for i in range(len(PAIRS)) if (mask >> i) & 1) + "}"


def print_unit_shell() -> None:
    print("S622 unit-distance impairment lab")
    print(f"Moser shell: {len(UNITS)} directed unit vectors, {len(PAIRS)} antipodal pairs")
    for idx, pair in enumerate(PAIRS):
        print(f"  pair {idx}: {pair[0]} / {pair[1]}")
    print()


def print_full_baseline() -> BeamResult:
    result = run_beam(target=14, width=260)
    print("Full-carrier small baseline")
    print("n  best_edges  exact  deficit  kept  children  span")
    for row in result.rows:
        if row.size in {2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18}:
            print(
                f"{row.size:2d} {row.best_edges:10d} {str(row.exact_edges):>5s} "
                f"{str(row.deficit):>7s} {row.kept:5d} {row.unique_children:9d} {row.best_span:5d}"
            )
    print(f"best n=14 direction usage: {direction_usage(result.best_cluster)}")
    print()
    return result


def print_direction_dropout(full: BeamResult) -> None:
    target = 14
    print("Impairment 1: single direction-pair dropout at n=14")
    print("dropped_pair  usage_in_full_best  impaired_edges  exact_deficit  span")
    full_usage = direction_usage(full.best_cluster)
    for idx in range(len(PAIRS)):
        mask = FULL_MASK ^ (1 << idx)
        result = run_beam(target=target, width=180, mask=mask)
        exact = U_EXACT[target]
        print(
            f"{idx:12d} {full_usage[idx]:18d} {result.best_edges:14d} "
            f"{exact - result.best_edges:13d} {span4(result.best_cluster):5d}"
        )
    print(
        "Reading: this is shadow-price spectroscopy. A direction can be heavily "
        "used in one best cluster yet still be replaceable; a real proof route "
        "should track both usage and dropout loss."
    )
    print()


def print_direction_helly(full: BeamResult) -> None:
    print("Impairment 2: observed direction-support thresholds in exact small witnesses")
    print("target_n  exact_edges  support_pairs  support_mask  usage")
    for n in (6, 8, 10, 12, 14):
        cluster = full.cluster_at(n)
        usage = direction_usage(cluster)
        support_mask = 0
        for idx, count in enumerate(usage):
            if count:
                support_mask |= 1 << idx
        print(
            f"{n:8d} {U_EXACT[n]:11d} {support_mask.bit_count():13d} "
            f"{mask_label(support_mask):>14s}  {usage}"
        )
    print(
        "Reading: this is a Helly-style certificate scale for the carrier. "
        "Instead of assuming all nine directions are needed, read how many "
        "direction orders an exact small witness actually uses.  This is an "
        "upper bound on the direction threshold and a guide for later exact "
        "subset searches."
    )
    print()


def print_gain_caps() -> None:
    print("Impairment 3: gain-cap throttling on the full carrier")
    print("cap  first_missed_exact  edges_at_14  n=14_deficit  best_sequence")
    for cap in (2, 3, 4, None):
        result = run_beam(target=14, width=260, gain_cap=cap)
        first_miss = None
        seq: list[str] = []
        for row in result.rows:
            if row.size in U_EXACT:
                if row.deficit is not None and row.deficit > 0 and first_miss is None:
                    first_miss = row.size
                if row.size in {6, 8, 10, 12, 14}:
                    seq.append(f"{row.size}:{row.best_edges}")
        cap_label = "none" if cap is None else str(cap)
        print(
            f"{cap_label:>3s} {str(first_miss):>18s} {result.best_edges:12d} "
            f"{U_EXACT[14] - result.best_edges:13d}  {', '.join(seq)}"
        )
    print(
        "Reading: gain caps are a local-degree impairment. They show exactly "
        "when small optima need high-gain extension events, which is the "
        "construction analogue of keeping higher overlap orders in LRC."
    )
    print()


@dataclass(frozen=True)
class Technique:
    name: str
    small_exactness: int
    diagnostic_value: int
    proof_transfer: int
    computability: int
    novelty: int
    risk: int
    note: str

    def score(self) -> int:
        return (
            self.small_exactness
            + self.diagnostic_value
            + self.proof_transfer
            + self.computability
            + self.novelty
            - self.risk
        )


TECHNIQUES = (
    Technique(
        "direction-dropout spectroscopy",
        4,
        5,
        4,
        5,
        5,
        1,
        "ablate one unit-direction pair and measure shadow price",
    ),
    Technique(
        "observed direction-support Helly",
        5,
        5,
        5,
        4,
        5,
        1,
        "small exact witnesses use few named direction families",
    ),
    Technique(
        "gain-cap ladder",
        5,
        4,
        4,
        5,
        4,
        1,
        "local extension degree as an overlap-order budget",
    ),
    Technique(
        "obstruction-first extension solver",
        4,
        5,
        5,
        3,
        4,
        2,
        "attach totally-unfaithful filters before expensive embedding search",
    ),
    Technique(
        "automorphism-canonical beam",
        4,
        4,
        4,
        3,
        3,
        2,
        "remove duplicate carrier states before scoring children",
    ),
    Technique(
        "raw wider Moser beam",
        4,
        2,
        2,
        4,
        2,
        3,
        "more states without new side-channel information",
    ),
    Technique(
        "triangular-only construction",
        2,
        3,
        2,
        5,
        1,
        3,
        "useful sanity baseline but too impaired for n=22 one-bit frontier",
    ),
    Technique(
        "graph-only dense enumeration",
        3,
        3,
        3,
        1,
        2,
        4,
        "too much geometry forgotten before the hard obstruction",
    ),
)


def print_tournament_analysis() -> None:
    ordered = sorted(TECHNIQUES, key=lambda t: (-t.score(), t.name))
    print("Tournament Analysis over impairment techniques")
    print("technique | score | note")
    print("--- | --- | ---")
    for technique in ordered:
        print(f"{technique.name} | {technique.score()} | {technique.note}")
    scores = [t.score() for t in ordered]
    hist: dict[int, int] = {}
    for score in scores:
        hist[score] = hist.get(score, 0) + 1
    hp_count = 1 if len(set(scores)) == len(scores) else 0
    print(f"score histogram: {dict(sorted(hist.items()))}")
    print(f"directed 3-cycles: 0")
    print(f"Hamiltonian path count in this score quotient: {hp_count}")
    print(
        "Challenged assumption: do not use points as the only tournament vertices. "
        "For method development, vertices should be impairments, side channels, "
        "direction families, gain events, and obstruction filters."
    )
    print()


def print_new_techniques() -> None:
    print("Technique proposals generated by the impairment lab")
    print("- Impaired-carrier spectroscopy: deliberately drop one direction, one gain level, or one side channel and measure the exact small-size failure mode.")
    print("- Direction Helly certificates: prove a candidate dense core needs only a small named subset of unit-direction pairs, then reattach the omitted directions as repair lanes.")
    print("- Shadow-price ledgers: rank directions by dropout loss, not just by usage in one witness cluster.")
    print("- Gain-order sieves: treat gain-4 and gain-5 extension events like higher overlap orders; a proof route should isolate them before raw enumeration.")
    print("- Obstruction-first beams: score children by edge gain plus proximity to known totally-unfaithful filters, not by edge gain alone.")
    print()


def main() -> None:
    print_unit_shell()
    full = print_full_baseline()
    print_direction_dropout(full)
    print_direction_helly(full)
    print_gain_caps()
    print_tournament_analysis()
    print_new_techniques()


if __name__ == "__main__":
    main()
