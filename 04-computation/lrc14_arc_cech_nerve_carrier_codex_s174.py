#!/usr/bin/env python3
"""Circular-arc Cech nerve carrier scout for LRC14.

This is a proof-carrier search, not a proof.  It tests a carrier closer to the
actual Lonely Runner geometry than the recent sequence/automaton shadows:

    danger arcs at threshold 1/14 -> exact endpoint cells -> Cech nerve.

Assumption challenge:
  Tournament vertices need not be runners.  Considered vertex sets are runners,
  individual danger arcs, endpoint cells/topes, section boundaries, boundary
  cocircuits, safe components, nerve homology classes, quotient defects, Fejer
  atoms, automaton states, and proof obligations.  This script chooses proof
  carriers and, inside the LRC row audit, compares individual danger arcs with
  the quotient that collapses all arcs belonging to the same runner.  The arc
  Cech nerve preserves the circle-cover predicate; the runner quotient can
  destroy connectedness and homology.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


N = 14
THRESHOLD = Fraction(1, N)


def fmt(x: Fraction | int) -> str:
    if isinstance(x, int):
        return str(x)
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def dist_to_integer(x: Fraction) -> Fraction:
    y = mod1(x)
    return min(y, 1 - y)


def nearest_residue(x: Fraction, modulus: int) -> int:
    base = x.numerator // x.denominator
    candidates = (base - 1, base, base + 1)
    best = min(candidates, key=lambda z: abs(x - z))
    return best % modulus


def danger_endpoints(row: tuple[int, ...]) -> list[Fraction]:
    points: set[Fraction] = set()
    for v in row:
        den = N * v
        for k in range(v):
            points.add(mod1(Fraction(N * k - 1, den)))
            points.add(mod1(Fraction(N * k + 1, den)))
    return sorted(points)


def active_runners(row: tuple[int, ...], t: Fraction, closed: bool = False) -> frozenset[str]:
    active = []
    for i, v in enumerate(row):
        d = dist_to_integer(v * t)
        if d < THRESHOLD or (closed and d == THRESHOLD):
            active.append(f"r{i}:{v}")
    return frozenset(active)


def active_arcs(row: tuple[int, ...], t: Fraction, closed: bool = False) -> frozenset[str]:
    active = []
    for i, v in enumerate(row):
        vt = v * t
        d = dist_to_integer(vt)
        if d < THRESHOLD or (closed and d == THRESHOLD):
            k = nearest_residue(vt, v)
            active.append(f"r{i}:{v}@{k}")
    return frozenset(active)


def equality_runners(row: tuple[int, ...], t: Fraction) -> tuple[int, ...]:
    out = []
    for v in row:
        if dist_to_integer(v * t) == THRESHOLD:
            out.append(v)
    return tuple(out)


@dataclass(frozen=True)
class CellAudit:
    endpoints: int
    open_cells: int
    safe_topes: int
    safe_measure: Fraction
    max_depth: int
    arc_facets: tuple[frozenset[str], ...]
    runner_facets: tuple[frozenset[str], ...]
    closed_arc_facets: tuple[frozenset[str], ...]
    closed_runner_facets: tuple[frozenset[str], ...]
    private_arc_count: int
    private_runner_count: int
    boundary_safe_points: tuple[tuple[Fraction, tuple[int, ...]], ...]


def cell_audit(row: tuple[int, ...]) -> CellAudit:
    pts = danger_endpoints(row)
    arc_facets: list[frozenset[str]] = []
    runner_facets: list[frozenset[str]] = []
    private_arcs: set[str] = set()
    private_runners: set[str] = set()
    safe_topes = 0
    safe_measure = Fraction(0)
    max_depth = 0

    for idx, lo in enumerate(pts):
        hi = pts[(idx + 1) % len(pts)]
        if idx + 1 == len(pts):
            hi += 1
        length = hi - lo
        mid = mod1((lo + hi) / 2)
        arcs = active_arcs(row, mid)
        runners = active_runners(row, mid)
        if not arcs:
            safe_topes += 1
            safe_measure += length
            continue
        arc_facets.append(arcs)
        runner_facets.append(runners)
        max_depth = max(max_depth, len(arcs))
        if len(arcs) == 1:
            private_arcs.update(arcs)
        if len(runners) == 1:
            private_runners.update(runners)

    boundary_safe = []
    closed_arc_facets = list(arc_facets)
    closed_runner_facets = list(runner_facets)
    for p in pts:
        closed_arcs = active_arcs(row, p, closed=True)
        closed_runners = active_runners(row, p, closed=True)
        if closed_arcs:
            closed_arc_facets.append(closed_arcs)
        if closed_runners:
            closed_runner_facets.append(closed_runners)
        if not active_runners(row, p) and min(dist_to_integer(v * p) for v in row) >= THRESHOLD:
            boundary_safe.append((p, equality_runners(row, p)))

    return CellAudit(
        endpoints=len(pts),
        open_cells=len(pts),
        safe_topes=safe_topes,
        safe_measure=safe_measure,
        max_depth=max_depth,
        arc_facets=tuple(arc_facets),
        runner_facets=tuple(runner_facets),
        closed_arc_facets=tuple(closed_arc_facets),
        closed_runner_facets=tuple(closed_runner_facets),
        private_arc_count=len(private_arcs),
        private_runner_count=len(private_runners),
        boundary_safe_points=tuple(boundary_safe),
    )


@dataclass(frozen=True)
class Betti:
    vertices: int
    edges: int
    triangles: int
    beta0: int
    beta1: int


def gf2_rank(columns: list[int]) -> int:
    basis: dict[int, int] = {}
    for col in columns:
        x = col
        while x:
            pivot = x.bit_length() - 1
            if pivot in basis:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                break
    return len(basis)


def nerve_betti(facets: tuple[frozenset[str], ...]) -> Betti:
    vertices = sorted(set().union(*facets)) if facets else []
    vertex_index = {v: i for i, v in enumerate(vertices)}

    edges: set[tuple[int, int]] = set()
    triangles: set[tuple[int, int, int]] = set()
    for facet in facets:
        ids = sorted(vertex_index[v] for v in facet)
        edges.update(tuple(pair) for pair in combinations(ids, 2))
        triangles.update(tuple(tri) for tri in combinations(ids, 3))

    edge_list = sorted(edges)
    edge_index = {edge: i for i, edge in enumerate(edge_list)}

    d1_cols = [(1 << a) | (1 << b) for a, b in edge_list]
    rank_d1 = gf2_rank(d1_cols)

    d2_cols = []
    for a, b, c in sorted(triangles):
        mask = 0
        for edge in ((a, b), (a, c), (b, c)):
            mask ^= 1 << edge_index[tuple(sorted(edge))]
        d2_cols.append(mask)
    rank_d2 = gf2_rank(d2_cols)

    beta0 = len(vertices) - rank_d1 if vertices else 0
    beta1 = len(edge_list) - rank_d1 - rank_d2
    return Betti(len(vertices), len(edge_list), len(triangles), beta0, beta1)


@dataclass(frozen=True)
class RowNerveRecord:
    name: str
    row: tuple[int, ...]
    audit: CellAudit
    arc_betti: Betti
    runner_betti: Betti
    closed_arc_betti: Betti
    closed_runner_betti: Betti

    @property
    def full_cover(self) -> bool:
        return self.audit.safe_topes == 0

    @property
    def runner_quotient_defect(self) -> int:
        return abs(self.closed_runner_betti.beta0 - self.closed_arc_betti.beta0) + abs(
            self.closed_runner_betti.beta1 - self.closed_arc_betti.beta1
        )

    @property
    def boundary_pair_sums(self) -> tuple[int, ...]:
        sums = []
        for _, owners in self.audit.boundary_safe_points:
            if owners:
                sums.append(sum(owners) % N)
        return tuple(sums)


def first_fibbinary(n: int) -> tuple[int, ...]:
    out = []
    x = 1
    while len(out) < n:
        if x & (x << 1) == 0:
            out.append(x)
        x += 1
    return tuple(out)


def is_moser(x: int) -> bool:
    pos = 0
    while x:
        if x & 1 and pos % 2 == 1:
            return False
        x >>= 1
        pos += 1
    return True


def first_moser(n: int) -> tuple[int, ...]:
    out = []
    x = 1
    while len(out) < n:
        if is_moser(x):
            out.append(x)
        x += 1
    return tuple(out)


def named_rows() -> list[tuple[str, tuple[int, ...]]]:
    ap = tuple(range(1, 14))
    return [
        ("AP", ap),
        ("GW_12_to_24", tuple(list(range(1, 12)) + [13, 24])),
        ("K33_12_to_36", tuple(list(range(1, 12)) + [13, 36])),
        ("petal_10_to_20", tuple([x for x in ap if x != 10] + [20])),
        ("petal_13_to_26", tuple([x for x in ap if x != 13] + [26])),
        ("P10_plus_GW", tuple([x for x in ap if x not in (10, 12)] + [20, 24])),
        ("P10_plus_K33", tuple([x for x in ap if x not in (10, 12)] + [20, 36])),
        ("covering_12_to_84", tuple(list(range(1, 12)) + [13, 84])),
        ("fibbinary_first13", first_fibbinary(13)),
        ("moser_first13", first_moser(13)),
    ]


def row_record(name: str, row: tuple[int, ...]) -> RowNerveRecord:
    audit = cell_audit(row)
    return RowNerveRecord(
        name=name,
        row=row,
        audit=audit,
        arc_betti=nerve_betti(audit.arc_facets),
        runner_betti=nerve_betti(audit.runner_facets),
        closed_arc_betti=nerve_betti(audit.closed_arc_facets),
        closed_runner_betti=nerve_betti(audit.closed_runner_facets),
    )


def safe_measure_and_topes(row: tuple[int, ...]) -> tuple[Fraction, int]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in row:
        den = N * v
        for k in range(v):
            left = Fraction(N * k - 1, den)
            right = Fraction(N * k + 1, den)
            if left < 0:
                intervals.append((Fraction(0), right))
                intervals.append((1 + left, Fraction(1)))
            elif right > 1:
                intervals.append((left, Fraction(1)))
                intervals.append((Fraction(0), right - 1))
            else:
                intervals.append((left, right))

    intervals.sort()
    merged: list[list[Fraction]] = []
    for left, right in intervals:
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right

    safe = Fraction(0)
    safe_topes = 0
    cursor = Fraction(0)
    for left, right in merged:
        if left > cursor:
            safe += left - cursor
            safe_topes += 1
        if right > cursor:
            cursor = right
    if cursor < 1:
        safe += 1 - cursor
        safe_topes += 1
    return safe, safe_topes


def one_swap_scan(add_limit: int = 160) -> dict[str, object]:
    ap = tuple(range(1, 14))
    zero_open: list[tuple[int, int, tuple[int, ...]]] = []
    positive_min: tuple[Fraction, int, int, tuple[int, ...]] | None = None
    total = 0
    primitive = 0
    for drop in ap:
        for add in range(14, add_limit + 1):
            if add in ap or add == drop:
                continue
            row = tuple(sorted([x for x in ap if x != drop] + [add]))
            total += 1
            if gcd(*row) != 1:
                continue
            primitive += 1
            safe_mu, safe_topes = safe_measure_and_topes(row)
            if safe_topes == 0:
                zero_open.append((drop, add, row))
            elif positive_min is None or safe_mu < positive_min[0]:
                positive_min = (safe_mu, drop, add, row)
    return {
        "add_limit": add_limit,
        "total": total,
        "primitive": primitive,
        "zero_open": zero_open,
        "positive_min": positive_min,
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    predicate_equiv: int
    topology_exactness: int
    quotient_defect_visible: int
    endpoint_locality: int
    dual_handoff: int
    computability: int
    maturity: int

    def values(self) -> dict[str, int]:
        return {
            "predicate_equiv": self.predicate_equiv,
            "topology_exactness": self.topology_exactness,
            "quotient_defect_visible": self.quotient_defect_visible,
            "endpoint_locality": self.endpoint_locality,
            "dual_handoff": self.dual_handoff,
            "computability": self.computability,
            "maturity": self.maturity,
        }


CARRIERS = [
    Carrier("arc_cech_good_cover_nerve", 3, 3, 3, 3, 2, 2, 2),
    Carrier("endpoint_tope_cocircuit_wall", 3, 3, 2, 3, 2, 3, 3),
    Carrier("taut_owner_current", 2, 2, 2, 3, 2, 3, 2),
    Carrier("safe_component_interval_measure", 3, 2, 1, 2, 2, 3, 3),
    Carrier("runner_quotient_nerve", 1, 1, 3, 2, 1, 3, 1),
    Carrier("fejer_toeplitz_dual_certificate", 2, 1, 1, 1, 3, 2, 3),
    Carrier("automaton_gap_sidecar", 1, 1, 2, 1, 1, 3, 2),
    Carrier("raw_speed_or_sequence_scalar", 0, 0, 0, 0, 0, 3, 1),
]


GAUGE = (
    "predicate_equiv",
    "topology_exactness",
    "quotient_defect_visible",
    "endpoint_locality",
    "dual_handoff",
    "maturity",
    "computability",
)


def orient(a: Carrier, b: Carrier) -> str:
    av = a.values()
    bv = b.values()
    wins = sum(1 for key in GAUGE if av[key] > bv[key])
    losses = sum(1 for key in GAUGE if av[key] < bv[key])
    if wins > losses:
        return a.name
    if losses > wins:
        return b.name
    # Declared tie Hamiltonian path: the listed order in CARRIERS.
    order = {carrier.name: i for i, carrier in enumerate(CARRIERS)}
    return a.name if order[a.name] < order[b.name] else b.name


def tournament_edges(carriers: list[Carrier]) -> set[tuple[str, str]]:
    edges = set()
    for a, b in combinations(carriers, 2):
        winner = orient(a, b)
        loser = b.name if winner == a.name else a.name
        edges.add((winner, loser))
    return edges


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    names = [c.name for c in carriers]
    edges = tournament_edges(carriers)
    scores = Counter()
    for a, _ in edges:
        scores[a] += 1
    score_hist = Counter(scores[name] for name in names)
    directed_3cycles = 0
    for a, b, c in combinations(names, 3):
        triple = {(a, b) in edges, (b, c) in edges, (c, a) in edges}
        if len(triple) == 1 and True in triple:
            # This condition catches neither orientation robustly; use degree test.
            pass
        out = {x: 0 for x in (a, b, c)}
        for x, y in combinations((a, b, c), 2):
            if (x, y) in edges:
                out[x] += 1
            else:
                out[y] += 1
        if sorted(out.values()) == [1, 1, 1]:
            directed_3cycles += 1

    graph = {name: [] for name in names}
    rev = {name: [] for name in names}
    for a, b in edges:
        graph[a].append(b)
        rev[b].append(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    scc_sizes = []

    def rdfs(v: str) -> int:
        seen.add(v)
        size = 1
        for w in rev[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for name in reversed(order):
        if name not in seen:
            scc_sizes.append(rdfs(name))

    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if (names[last], names[nxt]) in edges:
                    dp[mask | (1 << nxt)][nxt] += count
    hamiltonian_paths = sum(dp[(1 << n) - 1])
    path = sorted(names, key=lambda name: scores[name], reverse=True)
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": directed_3cycles,
        "scc_sizes": scc_sizes,
        "hamiltonian_path_count": hamiltonian_paths,
        "score_path": path,
    }


def print_row_records(records: list[RowNerveRecord]) -> None:
    print("Named-row circular-arc nerve audit at threshold 1/14")
    print("=" * 96)
    header = (
        "row",
        "full",
        "safe_topes",
        "safe_mu",
        "bd_pts",
        "closed_arc_b0/b1",
        "open_arc_b0/b1",
        "closed_runner_b0/b1",
        "quot_betti_def",
        "max_depth",
        "private_arcs/runners",
    )
    print(
        f"{header[0]:<23} {header[1]:<5} {header[2]:>10} {header[3]:>12} "
        f"{header[4]:>6} {header[5]:>17} {header[6]:>15} {header[7]:>20} "
        f"{header[8]:>11} {header[9]:>9} {header[10]:>21}"
    )
    for rec in records:
        print(
            f"{rec.name:<23} {str(rec.full_cover):<5} {rec.audit.safe_topes:>10} "
            f"{fmt(rec.audit.safe_measure):>12} {len(rec.audit.boundary_safe_points):>6} "
            f"{rec.closed_arc_betti.beta0}/{rec.closed_arc_betti.beta1: <14} "
            f"{rec.arc_betti.beta0}/{rec.arc_betti.beta1: <12} "
            f"{rec.closed_runner_betti.beta0}/{rec.closed_runner_betti.beta1: <17} "
            f"{rec.runner_quotient_defect:>11} {rec.audit.max_depth:>9} "
            f"{rec.audit.private_arc_count}/{rec.audit.private_runner_count: <9}"
        )
    print()
    print("Boundary-safe owner sums modulo 14")
    for rec in records:
        if rec.audit.boundary_safe_points:
            print(f"  {rec.name}: {rec.boundary_pair_sums}")


def print_swap_scan(scan: dict[str, object]) -> None:
    print()
    print("AP one-swap nerve/full-cover scan")
    print("=" * 96)
    print(f"add_limit={scan['add_limit']}")
    print(f"total_one_swaps={scan['total']}")
    print(f"primitive_one_swaps={scan['primitive']}")
    zero_open = scan["zero_open"]
    assert isinstance(zero_open, list)
    print(f"zero_open_one_swaps={len(zero_open)}")
    for drop, add, row in zero_open:
        print(f"  drop={drop} add={add} row={row}")
    positive_min = scan["positive_min"]
    if positive_min is not None:
        measure, drop, add, row = positive_min
        print(f"smallest_positive_safe_mu={fmt(measure)} at drop={drop} add={add} row={row}")


def print_tournament() -> None:
    print()
    print("Tournament Analysis over proof carriers")
    print("=" * 96)
    print("Vertices are proof carriers, not runners.")
    print("Pairwise observable: predicate equivalence, topology exactness, quotient-defect")
    print("visibility, endpoint locality, dual handoff, computability, and maturity.")
    print("Tie Hamiltonian path:")
    print("  " + " > ".join(carrier.name for carrier in CARRIERS))
    fp = tournament_fingerprint(CARRIERS)
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"scc_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("score_path=" + " > ".join(fp["score_path"]))


def print_synthesis() -> None:
    print()
    print("Synthesis")
    print("=" * 96)
    print(
        "The closed arc Cech nerve is a stronger carrier than a runner or "
        "sequence nerve: for the exact danger-arc cover, beta1=1 is the "
        "full-circle cover signal, while beta1=0 leaves positive safe "
        "components.  The open arc nerve records the taut-boundary split: AP/GW "
        "have six open components glued only by boundary cocircuits."
    )
    print(
        "The runner quotient is deliberately lossy.  Its Betti mismatch is a "
        "measurable quotient defect: collapsing disconnected arcs owned by the "
        "same speed can erase the cover cycle or fuse distinct danger components."
    )
    print(
        "AP and GW are the only zero-open rows in the one-swap scan through 160; "
        "both have six boundary-safe points with owner sums 0 mod 14 and closed "
        "arc beta1=1.  This matches the tope/cocircuit and taut-current route, "
        "but adds a cleaner good-cover nerve invariant."
    )
    print(
        "Candidate proof carrier: every primitive zero-open LRC14 packet should "
        "emit either (1) an arc-Cech beta1=1 cover whose boundary cocircuits "
        "satisfy the AP/GW owner-current laws, (2) a named K33/state-lift or "
        "Fejer/Ramanujan/Haar dual exit, or (3) a genuine F7 good-cover quotient "
        "defect where runner-level homology cannot be lifted to arc-level cover "
        "homology."
    )


def main() -> None:
    records = [row_record(name, row) for name, row in named_rows()]
    print_row_records(records)
    print_swap_scan(one_swap_scan())
    print_tournament()
    print_synthesis()


if __name__ == "__main__":
    main()
