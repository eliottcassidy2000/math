#!/usr/bin/env python3
"""Permutohedral geometry audit for an n=14 LRC proof attempt.

codex-2026-06-01 S526

The geometric lift:

    points x_i(t) = v_i t mod 1, together with observer x_0=0

determine a chamber of the affine braid arrangement whenever no two points
collide.  The circular order is a permutohedral chamber; collision walls are
adjacent-swap facets.  THM-384 says the observer is lonely exactly when the
two gaps adjacent to vertex 0 are both at least 1/14.  Thus n=14 asks whether
the one-parameter rational sweep must hit a long-long source facet.

This script audits that idea on the initial row, the known seven-ladder and
gate-ladder near-misses, and a small bounded family.  The point is not to claim
a proof; it is to find the exact permutohedral obstruction a proof must beat.

Tournament Analysis declaration:
    vertices: selected n=14 speed systems
    pairwise observable: permutohedral proof-strength profile
        (open source measure, wall source count, negative raw handoff defects,
         negative blocker SCC size, negative chamber count, best margin)
    switch/gauge: lexicographic comparison of the profile vector
    tie Hamiltonian path: listed selected-row order

Stored output:
    05-knowledge/results/lrc_n14_permutohedron_s526.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd


N = 14
THR = Fraction(1, N)
ZERO = Fraction(0)
ONE = Fraction(1)


def frac(value: Fraction) -> Fraction:
    return value - (value.numerator // value.denominator)


def dist0(value: Fraction) -> Fraction:
    f = frac(value)
    return min(f, ONE - f)


def fmt_frac(value: Fraction | None) -> str:
    if value is None:
        return "-"
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fmt_float(value: Fraction | float | None, places: int = 6) -> str:
    if value is None:
        return "-"
    return f"{float(value):.{places}f}"


def primitive(speeds: tuple[int, ...]) -> bool:
    return reduce(gcd, speeds) == 1


def positions(speeds: tuple[int, ...], time: Fraction) -> list[Fraction]:
    return [ZERO] + [frac(Fraction(speed) * time) for speed in speeds]


def endpoint_walls(speeds: tuple[int, ...]) -> set[Fraction]:
    walls = {ZERO, ONE}
    for speed in speeds:
        for k in range(speed):
            walls.add(Fraction(k * N + 1, N * speed))
            walls.add(Fraction(k * N + N - 1, N * speed))
    return {w for w in walls if ZERO <= w <= ONE}


def collision_walls(speeds: tuple[int, ...]) -> set[Fraction]:
    all_speeds = (0,) + speeds
    walls = {ZERO, ONE}
    for i, j in combinations(range(len(all_speeds)), 2):
        d = abs(all_speeds[i] - all_speeds[j])
        if d == 0:
            continue
        for k in range(d + 1):
            walls.add(Fraction(k, d))
    return {w for w in walls if ZERO <= w <= ONE}


def all_walls(speeds: tuple[int, ...]) -> list[Fraction]:
    return sorted(endpoint_walls(speeds) | collision_walls(speeds))


def circular_order(pos: list[Fraction]) -> tuple[int, ...]:
    return tuple(sorted(range(len(pos)), key=lambda idx: (pos[idx], idx)))


def order_distance(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    """Kendall distance between two orders; simultaneous walls may exceed 1."""
    rank = {v: i for i, v in enumerate(b)}
    seq = [rank[v] for v in a]
    inv = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                inv += 1
    return inv


@dataclass(frozen=True)
class State:
    time: Fraction
    order: tuple[int, ...]
    left_gap: Fraction
    right_gap: Fraction
    fiber: str
    left_neighbor: int
    right_neighbor: int
    blocker_sides: tuple[str, ...]
    blocker_classes: tuple[int, ...]
    qsafe: int
    qclasses: int

    def source(self) -> bool:
        return self.fiber == "LL"


def state_at(speeds: tuple[int, ...], time: Fraction) -> State:
    pos = positions(speeds, time)
    order = circular_order(pos)
    obs = order.index(0)
    left = order[(obs - 1) % len(order)]
    right = order[(obs + 1) % len(order)]
    left_gap = frac(pos[0] - pos[left])
    right_gap = frac(pos[right] - pos[0])
    fiber = ("L" if left_gap >= THR else "S") + (
        "L" if right_gap >= THR else "S"
    )

    blocker_sides: list[str] = []
    blocker_classes: list[int] = []
    if left_gap < THR and left != 0:
        blocker_sides.append(f"L:{speeds[left - 1]}")
        blocker_classes.append(speeds[left - 1] % 7)
    if right_gap < THR and right != 0:
        blocker_sides.append(f"R:{speeds[right - 1]}")
        blocker_classes.append(speeds[right - 1] % 7)

    classes: defaultdict[int, list[int]] = defaultdict(list)
    for speed in speeds:
        classes[ speed % 7 ].append(speed)
    qsafe = 0
    for members in classes.values():
        if all(dist0(Fraction(speed) * time) >= THR for speed in members):
            qsafe += 1

    return State(
        time=time,
        order=order,
        left_gap=left_gap,
        right_gap=right_gap,
        fiber=fiber,
        left_neighbor=left,
        right_neighbor=right,
        blocker_sides=tuple(blocker_sides),
        blocker_classes=tuple(sorted(set(blocker_classes))),
        qsafe=qsafe,
        qclasses=len(classes),
    )


@dataclass
class RowProfile:
    label: str
    speeds: tuple[int, ...]
    chambers: int
    distinct_orders: int
    source_measure: Fraction
    open_source_cells: int
    wall_source_count: int
    best_margin: Fraction
    max_qsafe_open: int
    max_qsafe_wall: int
    crt_class_count: int
    collision_events: int
    total_adjacent_swaps: int
    observer_neighbor_swaps: int
    blocker_handoffs: int
    raw_handoff_defects: int
    single_owner_handoffs: int
    blocker_scc_size: int
    blocker_edge_count: int
    blocker_hist: Counter[tuple[int, ...]]
    sample_source_times: list[Fraction]

    def strength(self) -> tuple:
        return (
            self.source_measure,
            self.wall_source_count,
            -self.raw_handoff_defects,
            -self.blocker_scc_size,
            -self.chambers,
            self.best_margin,
            Fraction(self.max_qsafe_open, self.crt_class_count),
            Fraction(self.max_qsafe_wall, self.crt_class_count),
        )


def scc_largest(edges: set[tuple[tuple[int, ...], tuple[int, ...]]]) -> int:
    nodes: set[tuple[int, ...]] = set()
    graph: defaultdict[tuple[int, ...], set[tuple[int, ...]]] = defaultdict(set)
    rgraph: defaultdict[tuple[int, ...], set[tuple[int, ...]]] = defaultdict(set)
    for u, v in edges:
        nodes.add(u)
        nodes.add(v)
        graph[u].add(v)
        rgraph[v].add(u)
    if not nodes:
        return 0

    def reach(start: tuple[int, ...], g: dict) -> set[tuple[int, ...]]:
        seen = {start}
        q = deque([start])
        while q:
            u = q.popleft()
            for v in g.get(u, ()):
                if v not in seen:
                    seen.add(v)
                    q.append(v)
        return seen

    largest = 1
    remaining = set(nodes)
    while remaining:
        node = next(iter(remaining))
        comp = reach(node, graph) & reach(node, rgraph)
        largest = max(largest, len(comp))
        remaining -= comp
    return largest


def profile_row(label: str, speeds: tuple[int, ...]) -> RowProfile:
    assert len(speeds) == 13, (label, speeds)
    assert primitive(speeds), (label, speeds)
    walls = all_walls(speeds)
    wall_set = set(walls)
    cells: list[tuple[Fraction, Fraction, State]] = []
    source_measure = ZERO
    open_source_cells = 0
    distinct_orders: set[tuple[int, ...]] = set()
    max_qsafe_open = 0
    crt_class_count = len({speed % 7 for speed in speeds})
    best_margin: Fraction | None = None

    for left, right in zip(walls, walls[1:]):
        if left >= right:
            continue
        mid = (left + right) / 2
        st = state_at(speeds, mid)
        cells.append((left, right, st))
        distinct_orders.add(st.order)
        max_qsafe_open = max(max_qsafe_open, st.qsafe)
        assert st.qclasses == crt_class_count
        margin = min(st.left_gap, st.right_gap) - THR
        best_margin = margin if best_margin is None else max(best_margin, margin)
        if st.source():
            source_measure += right - left
            open_source_cells += 1

    wall_source_count = 0
    max_qsafe_wall = 0
    sample_source_times: list[Fraction] = []
    for wall in walls:
        st = state_at(speeds, wall)
        assert st.qclasses == crt_class_count
        max_qsafe_wall = max(max_qsafe_wall, st.qsafe)
        margin = min(st.left_gap, st.right_gap) - THR
        best_margin = margin if best_margin is None else max(best_margin, margin)
        if st.source():
            wall_source_count += 1
            if len(sample_source_times) < 6:
                sample_source_times.append(wall)

    collision_set = collision_walls(speeds)
    collision_events = 0
    total_adjacent_swaps = 0
    observer_neighbor_swaps = 0
    blocker_handoffs = 0
    raw_handoff_defects = 0
    single_owner_handoffs = 0
    blocker_edges: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    blocker_hist: Counter[tuple[int, ...]] = Counter()

    for _, _, st in cells:
        if not st.source():
            blocker_hist[st.blocker_classes] += 1

    for idx in range(len(cells) - 1):
        _, _, prev = cells[idx]
        _, _, nxt = cells[idx + 1]
        wall = cells[idx][1]
        if wall in collision_set and prev.order != nxt.order:
            collision_events += 1
            dist = order_distance(prev.order, nxt.order)
            total_adjacent_swaps += dist
            if (
                prev.left_neighbor != nxt.left_neighbor
                or prev.right_neighbor != nxt.right_neighbor
            ):
                observer_neighbor_swaps += 1

        if prev.source() or nxt.source():
            continue
        if prev.blocker_classes != nxt.blocker_classes:
            blocker_handoffs += 1
            if prev.blocker_classes and nxt.blocker_classes:
                blocker_edges.add((prev.blocker_classes, nxt.blocker_classes))
            if len(prev.blocker_classes) == 1 and len(nxt.blocker_classes) == 1:
                single_owner_handoffs += 1
            wall_state = state_at(speeds, wall)
            if not wall_state.source():
                raw_handoff_defects += 1

    return RowProfile(
        label=label,
        speeds=speeds,
        chambers=len(cells),
        distinct_orders=len(distinct_orders),
        source_measure=source_measure,
        open_source_cells=open_source_cells,
        wall_source_count=wall_source_count,
        best_margin=best_margin if best_margin is not None else -THR,
        max_qsafe_open=max_qsafe_open,
        max_qsafe_wall=max_qsafe_wall,
        crt_class_count=crt_class_count,
        collision_events=collision_events,
        total_adjacent_swaps=total_adjacent_swaps,
        observer_neighbor_swaps=observer_neighbor_swaps,
        blocker_handoffs=blocker_handoffs,
        raw_handoff_defects=raw_handoff_defects,
        single_owner_handoffs=single_owner_handoffs,
        blocker_scc_size=scc_largest(blocker_edges),
        blocker_edge_count=len(blocker_edges),
        blocker_hist=blocker_hist,
        sample_source_times=sample_source_times,
    )


def selected_speed_sets() -> list[tuple[str, tuple[int, ...]]]:
    seven_ladder = (
        1,
        7,
        14,
        21,
        28,
        35,
        49,
        56,
        63,
        70,
        77,
        84,
        91,
    )
    gate_ladder = (1,) + tuple(14 * q for q in range(1, 14) if q != 6)
    double_gate = (1,) + tuple(28 * q for q in range(1, 14) if q != 6)
    return [
        ("initial_1_to_13", tuple(range(1, 14))),
        ("seven_ladder_drop_42", seven_ladder),
        ("gate_ladder_s380", gate_ladder),
        ("double_gate_ladder", double_gate),
        ("primes_row", (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)),
        ("hard_sporadic_s524", (1, 3, 4, 5, 9, 11, 13, 15, 17, 19, 23, 25, 29)),
    ]


def tournament_profile(rows: list[RowProfile]) -> tuple[Counter[int], int, list[int], int]:
    m = len(rows)
    adj = [[0] * m for _ in range(m)]
    for i, j in combinations(range(m), 2):
        if rows[i].strength() >= rows[j].strength():
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    score_hist = Counter(sum(row) for row in adj)
    c3 = 0
    for a, b, c in combinations(range(m), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            c3 += 1

    # SCC sizes by reachability.
    def reach(start: int) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            u = q.popleft()
            for v in range(m):
                if adj[u][v] and v not in seen:
                    seen.add(v)
                    q.append(v)
        return seen

    remaining = set(range(m))
    sccs: list[int] = []
    while remaining:
        node = next(iter(remaining))
        fwd = reach(node)
        bwd = {v for v in range(m) if node in reach(v)}
        comp = fwd & bwd
        sccs.append(len(comp))
        remaining -= comp

    full = (1 << m) - 1
    dp = [[0] * m for _ in range(1 << m)]
    for i in range(m):
        dp[1 << i][i] = 1
    for mask in range(1 << m):
        for last in range(m):
            cur = dp[mask][last]
            if not cur:
                continue
            for nxt in range(m):
                if not ((mask >> nxt) & 1) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += cur
    hp = sum(dp[full])
    return score_hist, c3, sorted(sccs, reverse=True), hp


def bounded_family(max_speed: int = 16) -> dict[str, int | Fraction]:
    systems = []
    for combo in combinations(range(1, max_speed + 1), 13):
        if primitive(combo):
            systems.append(combo)

    no_source = 0
    wall_only = 0
    max_raw_defects = 0
    max_blocker_scc = 0
    hardest: tuple[str, Fraction, int, int] | None = None
    for idx, speeds in enumerate(systems):
        prof = profile_row(f"bounded_{idx}", speeds)
        if prof.source_measure == 0 and prof.wall_source_count == 0:
            no_source += 1
        if prof.source_measure == 0 and prof.wall_source_count > 0:
            wall_only += 1
        max_raw_defects = max(max_raw_defects, prof.raw_handoff_defects)
        max_blocker_scc = max(max_blocker_scc, prof.blocker_scc_size)
        key = (prof.source_measure, -prof.raw_handoff_defects)
        if hardest is None or key < (hardest[1], -hardest[2]):
            hardest = (str(speeds), prof.source_measure, prof.raw_handoff_defects, prof.blocker_scc_size)

    return {
        "systems": len(systems),
        "no_source": no_source,
        "wall_only": wall_only,
        "max_raw_defects": max_raw_defects,
        "max_blocker_scc": max_blocker_scc,
        "hardest": hardest[0] if hardest else "-",
        "hardest_source_measure": hardest[1] if hardest else ZERO,
        "hardest_raw_defects": hardest[2] if hardest else 0,
        "hardest_blocker_scc": hardest[3] if hardest else 0,
    }


def print_row(prof: RowProfile) -> None:
    print(f"{prof.label}:")
    print(f"  speeds={prof.speeds}")
    print(
        "  chambers={ch} distinct_orders={orders} collisions={coll} "
        "adjacent_swaps={swaps} obs_neighbor_swaps={obs}".format(
            ch=prof.chambers,
            orders=prof.distinct_orders,
            coll=prof.collision_events,
            swaps=prof.total_adjacent_swaps,
            obs=prof.observer_neighbor_swaps,
        )
    )
    print(
        "  source_measure={measure} ({measure_f}) open_cells={open_cells} "
        "wall_sources={walls} best_margin={margin} ({margin_f})".format(
            measure=fmt_frac(prof.source_measure),
            measure_f=fmt_float(prof.source_measure),
            open_cells=prof.open_source_cells,
            walls=prof.wall_source_count,
            margin=fmt_frac(prof.best_margin),
            margin_f=fmt_float(prof.best_margin),
        )
    )
    print(
        "  CRT qsafe max: open={open_q}/{classes} wall={wall_q}/{classes}".format(
            open_q=prof.max_qsafe_open,
            wall_q=prof.max_qsafe_wall,
            classes=prof.crt_class_count,
        )
    )
    print(
        "  blocker handoffs={h} single_owner={single} raw_handoff_defects={bad} "
        "blocker_edges={edges} blocker_SCC_max={scc}".format(
            h=prof.blocker_handoffs,
            single=prof.single_owner_handoffs,
            bad=prof.raw_handoff_defects,
            edges=prof.blocker_edge_count,
            scc=prof.blocker_scc_size,
        )
    )
    top = prof.blocker_hist.most_common(6)
    print(f"  top blocker class states: {top}")
    if prof.sample_source_times:
        print(
            "  sample wall/source times: "
            + ", ".join(fmt_frac(t) for t in prof.sample_source_times)
        )
    print()


def main() -> None:
    print("LRC n=14 permutohedron audit -- codex-2026-06-01 S526")
    print()
    print("GEOMETRIC MODEL")
    print("=" * 72)
    print("The circular order of {0, v_i t} is a chamber of the affine braid fan.")
    print("Collision times are adjacent-swap facets of the permutohedron.")
    print("THM-384 makes the LRC target the long-long pair of gaps adjacent to")
    print("the observer.  A counterexample would be a closed one-parameter sweep")
    print("through these chambers avoiding every long-long source facet.")
    print()

    rows = [profile_row(label, speeds) for label, speeds in selected_speed_sets()]

    print("PART A: selected n=14 rows")
    print("=" * 72)
    for prof in rows:
        print_row(prof)

    print("PART B: bounded primitive family, max_speed=16")
    print("=" * 72)
    bounded = bounded_family(16)
    print(f"systems={bounded['systems']}")
    print(f"no_source={bounded['no_source']} wall_only={bounded['wall_only']}")
    print(
        "max_raw_handoff_defects={bad} max_blocker_SCC={scc}".format(
            bad=bounded["max_raw_defects"],
            scc=bounded["max_blocker_scc"],
        )
    )
    print(
        "hardest_by_source_measure={hardest} source_measure={measure} "
        "raw_defects={bad} blocker_SCC={scc}".format(
            hardest=bounded["hardest"],
            measure=fmt_frac(bounded["hardest_source_measure"]),
            bad=bounded["hardest_raw_defects"],
            scc=bounded["hardest_blocker_scc"],
        )
    )
    print()

    print("PART C: Tournament Analysis over selected rows")
    print("=" * 72)
    score_hist, c3, sccs, hp = tournament_profile(rows)
    print(f"score_hist={tuple(sorted(score_hist.items()))}")
    print(f"c3={c3}")
    print(f"SCCs={tuple(sccs)}")
    print(f"Hamiltonian_paths={hp}")
    print("ranking by proof-strength profile:")
    for prof in sorted(rows, key=lambda row: row.strength(), reverse=True):
        print(
            "  {label}: strength={strength}".format(
                label=prof.label,
                strength=(
                    fmt_frac(prof.source_measure),
                    prof.wall_source_count,
                    -prof.raw_handoff_defects,
                    -prof.blocker_scc_size,
                    -prof.chambers,
                    fmt_frac(prof.best_margin),
                    f"{prof.max_qsafe_open}/{prof.crt_class_count}",
                    f"{prof.max_qsafe_wall}/{prof.crt_class_count}",
                ),
            )
        )
    print()

    print("SYNTHESIS")
    print("=" * 72)
    print("1. The permutohedron lift is the right geometry: each row is a")
    print("   one-parameter path through adjacent-swap chambers, and the LRC")
    print("   target is exactly the long-long facet at the observer.")
    print("2. The tempting raw lemma 'every blocker-owner handoff crosses source'")
    print("   is false in the hard rows: raw_handoff_defects are nonzero and")
    print("   blocker owner states can have nontrivial SCCs.  So the free graph")
    print("   obstruction from HYP-1997 survives in permutohedral language.")
    print("3. The useful refinement is labelled debt: hard handoffs that avoid")
    print("   source are precisely where endpoint debt is exported to deeper")
    print("   quotient layers (seven-ladder -> gate-ladder -> double-gate).")
    print("4. A plausible n=14 proof should combine this chamber ledger with a")
    print("   Farkas/Hall endpoint-debt certificate: no closed blocker-owner")
    print("   circulation in the permutohedral fan should be realizable after all")
    print("   private endpoint leaves are peeled.")
    print("5. This session therefore does not prove n=14.  It narrows the proof")
    print("   target to a labelled permutohedral handoff theorem, rather than a")
    print("   scalar gap argument or a pure metagraph feedback-vertex argument.")


if __name__ == "__main__":
    main()
