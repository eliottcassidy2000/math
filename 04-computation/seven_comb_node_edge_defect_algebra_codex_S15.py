#!/usr/bin/env python3
"""Exact node/edge defect audit for the six THM-856 seven-comb pilots.

For the prefix-safe set

    E_[5] = {t in R/Z : ||qt|| > 1/13 for q=1,...,5},

the seven remaining danger events are A_i=E_[5] cap D_(x_i).  This script
checks, using only ``fractions.Fraction`` interval arithmetic, the coloured
Hunter decomposition

    s_i   = mu(A_i) - 2e/13,
    h_ij  = [Q(a+b)-Q(b-a)]/(169ab),
    eta_ij= mu(A_i cap A_j) - e*mu(D_a cap D_b),
    c_ij  = e*h_ij + eta_ij,

where e=mu(E_[5]), g=gcd(x_i,x_j), (a,b)=(x_i/g,x_j/g), and
Q(r)=r(13-r) for r reduced modulo 13.  It then verifies exactly that

    direct Hunter lower bound = 11e/169 - sum_i s_i + MST(c).

The Kruskal threshold-rank word is the sequence of graphic-matroid rank
increments as the distinct c-levels are added in descending order.  It is a
strictly smaller sufficient statistic for the MST value than the entire
edge ordering, but it still retains incidence at every contributing level.

Tournament Analysis uses proof-carrier vertices A_i, not physical runners.
For a pair with p_ij=mu(A_i cap A_j), compare the two conditional supports

    omega(i,j) = p_ij/mu(A_i) - p_ij/mu(A_j).

Orient i -> j when omega(i,j)>0.  Equal values (including zero-overlap pairs)
are switched by the vertex-colour gauge (mu(A_i),i).  The declared tie
Hamiltonian path is the exact increasing order of that gauge.  This quotient preserves conditional-support order.  It
destroys the magnitudes and incidence of c_ij and therefore cannot replace
the coloured graph or certify Hunter positivity by itself.  That information
loss is an intentional assumption challenge required by AGENTS.md.

Scope: these are six finite exact packet audits.  Positive Hunter values in
four packets and negative values in two packets do not prove a uniform
radius-seven theorem and do not classify all projective-ratio packets.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd


BASE_PAIR = F(4, 169)
PREFIX = (1, 2, 3, 4, 5)
RESIDUES = (6, 7, 8, 9, 10, 11, 12)
PACKETS = (
    ("random-1", (292, 150, 359, 74, 88, 479, 116)),
    ("random-2", (331, 514, 73, 451, 205, 63, 103)),
    ("random-3", (383, 371, 86, 230, 101, 492, 389)),
    ("random-4", (71, 501, 125, 217, 517, 76, 506)),
    ("consecutive-high", (499, 500, 501, 502, 503, 504, 505)),
    ("consecutive-low", (32, 33, 34, 35, 36, 37, 38)),
)

Interval = tuple[F, F]
Intervals = tuple[Interval, ...]
Edge = tuple[int, int]


def measure(intervals: Intervals) -> F:
    return sum((hi - lo for lo, hi in intervals), F(0))


def meet(left: Intervals, right: Intervals) -> Intervals:
    """Intersection of two sorted, internally disjoint interval unions."""
    out: list[Interval] = []
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


def subtract(left: Intervals, right: Intervals) -> Intervals:
    """Measure-exact left-minus-right for sorted interval unions."""
    out: list[Interval] = []
    j = 0
    for lo, hi in left:
        cursor = lo
        while j < len(right) and right[j][1] <= lo:
            j += 1
        k = j
        while k < len(right) and right[k][0] < hi:
            cut_lo, cut_hi = right[k]
            if cursor < cut_lo:
                out.append((cursor, min(cut_lo, hi)))
            cursor = max(cursor, cut_hi)
            if cursor >= hi:
                break
            k += 1
        if cursor < hi:
            out.append((cursor, hi))
    return tuple((lo, hi) for lo, hi in out if lo < hi)


def danger(speed: int) -> Intervals:
    """D_speed={t: ||speed*t||<=1/13}, split at the circle cut."""
    assert speed > 0
    radius = F(1, 13 * speed)
    intervals: list[Interval] = [(F(0), radius)]
    for k in range(1, speed):
        centre = F(k, speed)
        intervals.append((centre - radius, centre + radius))
    intervals.append((F(1) - radius, F(1)))
    return tuple(intervals)


def safe_set(speeds: tuple[int, ...]) -> Intervals:
    current: Intervals = ((F(0), F(1)),)
    for speed in speeds:
        one_safe = tuple(
            (
                F(13 * k + 1, 13 * speed),
                F(13 * (k + 1) - 1, 13 * speed),
            )
            for k in range(speed)
        )
        current = meet(current, one_safe)
    return current


def q13(value: int) -> int:
    residue = value % 13
    return residue * (13 - residue)


def projective_pair(x: int, y: int) -> tuple[int, int, int]:
    common = gcd(x, y)
    a, b = x // common, y // common
    if a > b:
        a, b = b, a
    assert gcd(a, b) == 1
    return common, a, b


def projective_defect(a: int, b: int) -> F:
    assert 1 <= a <= b and gcd(a, b) == 1
    return F(q13(a + b) - q13(b - a), 169 * a * b)


class DSU:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.components = size

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> bool:
        left_root, right_root = self.find(left), self.find(right)
        if left_root == right_root:
            return False
        self.parent[right_root] = left_root
        self.components -= 1
        return True


def descending_edges(weights: dict[Edge, F]) -> list[Edge]:
    return sorted(
        weights,
        key=lambda edge: (weights[edge], -edge[0], -edge[1]),
        reverse=True,
    )


def kruskal_max_tree(weights: dict[Edge, F], order: int) -> tuple[F, tuple[Edge, ...]]:
    dsu = DSU(order)
    chosen: list[Edge] = []
    total = F(0)
    for edge in descending_edges(weights):
        if dsu.union(*edge):
            chosen.append(edge)
            total += weights[edge]
            if len(chosen) == order - 1:
                break
    assert len(chosen) == order - 1 and dsu.components == 1
    return total, tuple(chosen)


def prim_max_tree(weights: dict[Edge, F], order: int) -> F:
    """Independent maximum-tree implementation used as an exact oracle."""
    inside = {0}
    total = F(0)
    while len(inside) < order:
        candidates: list[tuple[F, int, int]] = []
        for left in inside:
            for right in range(order):
                if right not in inside:
                    edge = (min(left, right), max(left, right))
                    candidates.append((weights[edge], -left, -right))
        weight, neg_left, neg_right = max(candidates)
        del neg_left
        total += weight
        inside.add(-neg_right)
    return total


@dataclass(frozen=True)
class ThresholdRow:
    level: F
    edges_at_level: tuple[Edge, ...]
    components: int
    rank_increment: int


def threshold_rank_word(
    weights: dict[Edge, F], order: int
) -> tuple[tuple[ThresholdRow, ...], F]:
    """Graphic rank increments for descending distinct weight thresholds."""
    levels = sorted(set(weights.values()), reverse=True)
    dsu = DSU(order)
    previous_rank = 0
    rows: list[ThresholdRow] = []
    rank_integral = F(0)
    for level in levels:
        level_edges = tuple(sorted(edge for edge, value in weights.items() if value == level))
        for edge in level_edges:
            dsu.union(*edge)
        rank = order - dsu.components
        increment = rank - previous_rank
        assert increment >= 0
        rows.append(ThresholdRow(level, level_edges, dsu.components, increment))
        rank_integral += level * increment
        previous_rank = rank
    assert previous_rank == order - 1
    return tuple(rows), rank_integral


def directed_three_cycles(adjacency: tuple[tuple[bool, ...], ...]) -> int:
    count = 0
    for i, j, k in combinations(range(len(adjacency)), 3):
        if (
            adjacency[i][j] and adjacency[j][k] and adjacency[k][i]
        ) or (
            adjacency[j][i] and adjacency[i][k] and adjacency[k][j]
        ):
            count += 1
    return count


def scc_sizes(adjacency: tuple[tuple[bool, ...], ...]) -> tuple[int, ...]:
    """Tarjan SCC fingerprint for the small exact tournament."""
    order = len(adjacency)
    index = 0
    indices = [-1] * order
    lowlink = [0] * order
    stack: list[int] = []
    on_stack = [False] * order
    sizes: list[int] = []

    def visit(vertex: int) -> None:
        nonlocal index
        indices[vertex] = lowlink[vertex] = index
        index += 1
        stack.append(vertex)
        on_stack[vertex] = True
        for target in range(order):
            if not adjacency[vertex][target]:
                continue
            if indices[target] == -1:
                visit(target)
                lowlink[vertex] = min(lowlink[vertex], lowlink[target])
            elif on_stack[target]:
                lowlink[vertex] = min(lowlink[vertex], indices[target])
        if lowlink[vertex] == indices[vertex]:
            size = 0
            while True:
                target = stack.pop()
                on_stack[target] = False
                size += 1
                if target == vertex:
                    break
            sizes.append(size)

    for vertex in range(order):
        if indices[vertex] == -1:
            visit(vertex)
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count(adjacency: tuple[tuple[bool, ...], ...]) -> int:
    """Subset DP count of directed Hamiltonian paths."""
    order = len(adjacency)
    dp = [[0] * order for _ in range(1 << order)]
    for vertex in range(order):
        dp[1 << vertex][vertex] = 1
    for mask in range(1 << order):
        for last in range(order):
            ways = dp[mask][last]
            if not ways:
                continue
            for target in range(order):
                if not mask & (1 << target) and adjacency[last][target]:
                    dp[mask | (1 << target)][target] += ways
    return sum(dp[-1])


@dataclass(frozen=True)
class TournamentFingerprint:
    adjacency: tuple[tuple[bool, ...], ...]
    tie_count: int
    tie_path: tuple[int, ...]
    score_histogram: tuple[tuple[int, int], ...]
    three_cycles: int
    sccs: tuple[int, ...]
    canonical_flips: tuple[Edge, ...]
    hamiltonian_paths: int


def proof_carrier_tournament(
    masses: tuple[F, ...], pair_masses: dict[Edge, F]
) -> TournamentFingerprint:
    order = len(masses)
    adjacency = [[False] * order for _ in range(order)]
    tie_count = 0
    for i, j in combinations(range(order), 2):
        pair = pair_masses[(i, j)]
        assert pair >= 0 and masses[i] > 0 and masses[j] > 0
        observable = pair / masses[i] - pair / masses[j]
        if observable > 0:
            adjacency[i][j] = True
        elif observable < 0:
            adjacency[j][i] = True
        else:
            tie_count += 1
            if (masses[i], i) < (masses[j], j):
                adjacency[i][j] = True
            else:
                adjacency[j][i] = True
    frozen = tuple(tuple(row) for row in adjacency)
    tie_path = tuple(sorted(range(order), key=lambda vertex: (masses[vertex], vertex)))
    for left, right in zip(tie_path, tie_path[1:]):
        assert frozen[left][right]
    scores = tuple(sum(row) for row in frozen)
    score_histogram = tuple(sorted(Counter(scores).items()))
    canonical_flips = tuple(
        (i, j) for i, j in combinations(range(order), 2) if frozen[j][i]
    )
    return TournamentFingerprint(
        frozen,
        tie_count,
        tie_path,
        score_histogram,
        directed_three_cycles(frozen),
        scc_sizes(frozen),
        canonical_flips,
        hamiltonian_path_count(frozen),
    )


@dataclass(frozen=True)
class EdgeColour:
    common_scale: int
    a: int
    b: int
    h: F
    eta: F
    credit: F
    restricted_mass: F


@dataclass(frozen=True)
class PacketAudit:
    name: str
    speeds: tuple[int, ...]
    masses: tuple[F, ...]
    anomalies: tuple[F, ...]
    edge_colours: dict[Edge, EdgeColour]
    hunter: F
    decomposition: F
    uncovered: F
    mst_credit: F
    mst_edges: tuple[Edge, ...]
    thresholds: tuple[ThresholdRow, ...]
    tournament: TournamentFingerprint


def audit_packet(name: str, speeds: tuple[int, ...], core: Intervals) -> PacketAudit:
    order = len(speeds)
    assert order == 7
    expected_residues = (5, 6, 7, 8, 9, 10, 11) if name == "consecutive-high" else RESIDUES
    # The S312 high-consecutive stress test is deliberately shifted one
    # residue from the scale-one 6,...,12 chart; retain it exactly as recorded
    # instead of silently relabelling the historical finite pilot.
    assert tuple(speed % 13 for speed in speeds) == expected_residues
    e = measure(core)
    events = tuple(meet(core, danger(speed)) for speed in speeds)
    masses = tuple(measure(event) for event in events)
    anomalies = tuple(mass - F(2, 13) * e for mass in masses)

    edge_colours: dict[Edge, EdgeColour] = {}
    pair_masses: dict[Edge, F] = {}
    credits: dict[Edge, F] = {}
    for i, j in combinations(range(order), 2):
        edge = (i, j)
        restricted = measure(meet(events[i], danger(speeds[j])))
        assert restricted == measure(meet(events[j], danger(speeds[i])))
        common, a, b = projective_pair(speeds[i], speeds[j])
        h = projective_defect(a, b)
        global_pair = measure(meet(danger(a), danger(b)))
        assert global_pair == BASE_PAIR + h
        eta = restricted - e * global_pair
        assert abs(eta) <= F(2 * len(core), common)
        credit = e * h + eta
        assert credit == restricted - e * BASE_PAIR
        colour = EdgeColour(common, a, b, h, eta, credit, restricted)
        edge_colours[edge] = colour
        pair_masses[edge] = restricted
        credits[edge] = credit

    pair_tree_prim = prim_max_tree(pair_masses, order)
    pair_tree_kruskal, _ = kruskal_max_tree(pair_masses, order)
    assert pair_tree_prim == pair_tree_kruskal
    hunter = e - sum(masses, F(0)) + pair_tree_prim

    mst_credit, mst_edges = kruskal_max_tree(credits, order)
    assert prim_max_tree(credits, order) == mst_credit
    thresholds, threshold_integral = threshold_rank_word(credits, order)
    assert threshold_integral == mst_credit
    assert pair_tree_prim == F(4 * (order - 1), 169) * e + mst_credit

    decomposition = F(11, 169) * e - sum(anomalies, F(0)) + mst_credit
    assert hunter == decomposition

    uncovered_intervals = core
    for speed in speeds:
        uncovered_intervals = subtract(uncovered_intervals, danger(speed))
    uncovered = measure(uncovered_intervals)
    assert uncovered >= hunter

    tournament = proof_carrier_tournament(masses, pair_masses)
    # This conditional-support orientation is a total mass order.  Freeze the
    # loss explicitly: the tournament is transitive even when Hunter fails.
    assert tournament.three_cycles == 0
    assert tournament.sccs == (1,) * order
    assert tournament.score_histogram == tuple((score, 1) for score in range(order))
    assert tournament.hamiltonian_paths == 1

    return PacketAudit(
        name,
        speeds,
        masses,
        anomalies,
        edge_colours,
        hunter,
        decomposition,
        uncovered,
        mst_credit,
        mst_edges,
        thresholds,
        tournament,
    )


def edge_name(edge: Edge) -> str:
    return f"A{edge[0] + 1}A{edge[1] + 1}"


def print_packet(audit: PacketAudit, reference: TournamentFingerprint | None) -> None:
    print(f"\nPACKET {audit.name}")
    print(f"speeds={audit.speeds}")
    print("VERTEX COLOURS (proof events, not runners)")
    for index, (speed, mass, anomaly) in enumerate(
        zip(audit.speeds, audit.masses, audit.anomalies), start=1
    ):
        print(
            f"  A{index} residue={speed % 13:2d} x={speed:3d} "
            f"mass={mass} s={anomaly}"
        )

    print("EDGE COLOURS (reduced ratio, scale, h, eta, c, restricted mass)")
    mst_set = set(audit.mst_edges)
    for edge in sorted(audit.edge_colours):
        colour = audit.edge_colours[edge]
        selected = "*" if edge in mst_set else "."
        print(
            f" {selected}{edge_name(edge)} ratio={colour.a}:{colour.b} "
            f"g={colour.common_scale} h={colour.h} eta={colour.eta} "
            f"c={colour.credit} pair={colour.restricted_mass}"
        )

    rank_word = tuple(row.rank_increment for row in audit.thresholds)
    print(f"mst_c={audit.mst_credit} mst_edges={tuple(edge_name(e) for e in audit.mst_edges)}")
    print(f"kruskal_threshold_rank_word={rank_word}")
    print("KRUSKAL LEVELS (descending c; incidence retained)")
    for row in audit.thresholds:
        print(
            f"  lambda={row.level} edges={tuple(edge_name(e) for e in row.edges_at_level)} "
            f"components={row.components} delta_rank={row.rank_increment}"
        )

    print(
        f"hunter_direct={audit.hunter} (~{float(audit.hunter):+.8f}) "
        f"decomposition={audit.decomposition} actual_uncovered={audit.uncovered} "
        f"coercive={audit.hunter > 0}"
    )

    tournament = audit.tournament
    reference_flips: tuple[Edge, ...]
    if reference is None:
        reference_flips = ()
    else:
        reference_flips = tuple(
            (i, j)
            for i, j in combinations(range(7), 2)
            if tournament.adjacency[i][j] != reference.adjacency[i][j]
        )
    print("TOURNAMENT ANALYSIS")
    print("  vertices=A_i=E_[5] cap D_xi (proof carriers; not runner identities)")
    print("  observable=pair/mass_i-pair/mass_j; tie gauge=(mass_i,i)")
    print(
        f"  tie_path={tuple(vertex + 1 for vertex in tournament.tie_path)} "
        f"ties={tournament.tie_count} score_hist={tournament.score_histogram}"
    )
    print(
        f"  directed_3cycles={tournament.three_cycles} scc_sizes={tournament.sccs} "
        f"hamiltonian_paths={tournament.hamiltonian_paths}"
    )
    print(
        f"  flips_from_slot_path={tuple(edge_name(e) for e in tournament.canonical_flips)} "
        f"flips_from_packet_1={tuple(edge_name(e) for e in reference_flips)}"
    )


def main() -> None:
    core = safe_set(PREFIX)
    e = measure(core)
    assert len(core) == 10
    assert e == F(7, 15)
    print("SEVEN-COMB NODE/EDGE DEFECT ALGEBRA: EXACT FRACTION AUDIT")
    print(f"prefix={PREFIX} E_components={len(core)} e={e}")
    print("vertex=s_i=mu(E cap D_xi)-2e/13")
    print("edge=h(projective ratio), eta(endpoint pullback), c=e*h+eta")
    print("identity=Hunter=11e/169-sum(s_i)+MST(c)")
    print("kruskal_word=graphic-rank increments at descending distinct c levels")

    audits = tuple(audit_packet(name, speeds, core) for name, speeds in PACKETS)
    signs = tuple(audit.hunter > 0 for audit in audits)
    assert signs == (True, True, True, True, False, False)
    expected_hunter = (
        F(6350550702737, 192614324316700),
        F(13322563785667, 393758910940395),
        F(3207706171151467, 88285303438984470),
        F(15342920255059, 404839718122125),
        -F(3248280870311, 5180297684976375),
        -F(48097487, 2871108240),
    )
    expected_uncovered = (
        F(584564587283, 4377598279925),
        F(10891281944419, 78751782188079),
        F(4137098311279637, 29428434479661490),
        F(131331423591653, 944626008951625),
        F(1190249283463042, 21757250276900775),
        F(364233, 9202270),
    )
    assert tuple(audit.hunter for audit in audits) == expected_hunter
    assert tuple(audit.uncovered for audit in audits) == expected_uncovered

    reference = audits[0].tournament
    for index, audit in enumerate(audits):
        print_packet(audit, None if index == 0 else reference)

    print("\nCROSS-PACKET VERDICT")
    print(f"hunter_signs={signs} expected=(True,True,True,True,False,False)")
    print("all_interval_masses=exact Fraction; direct/decomposition identities=PASS")
    print("Prim=Kruskal and Kruskal threshold-rank integral=MST(c): PASS")
    print("Tournament quotient is transitive in all six packets and is NOT sufficient for MST(c)")
    print("HONEST SCOPE: six finite pilots only; uniform radius-seven closure remains open")


if __name__ == "__main__":
    main()
