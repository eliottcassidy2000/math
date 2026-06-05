#!/usr/bin/env python3
"""Unit-distance Moser fixed-quantum carrier scout.

codex-2026-06-04-S648

S646 turned the perfect-number fixed-point idea into an LRC gcd-shell clock:
the pair-sum clock C=27 is fixed by a divisor/Pillai carrier.  This script
ports the same habit to the unit-distance side.

In the THM-408 Moser spine ladder, adding one full slab eventually has a fixed
edge-channel increment

    Delta E = 27 = 8 spine edges + 19 bulk edges.

The fixed object is therefore not a scalar point count.  It is a direction
pair ledger for the unit shell, split into the retained unit-spine section and
the hidden bulk carrier.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

import unit_distance_impairment_lab_s622 as lab  # noqa: E402
import unit_distance_spine_ladder_s628 as s628  # noqa: E402


Point = tuple[int, int, int, int]


def direction_reps() -> tuple[Point, ...]:
    reps: list[Point] = []
    for u in sorted(lab.UNITS):
        rep = min(u, tuple(-x for x in u))
        if rep not in reps:
            reps.append(rep)
    return tuple(reps)


REPS = direction_reps()
REP_OF: dict[Point, Point] = {}
for unit in lab.UNITS:
    REP_OF[unit] = min(unit, tuple(-x for x in unit))


def edge_rep(p: Point, q: Point) -> Point | None:
    d = lab.sub(q, p)
    if d in REP_OF:
        return REP_OF[d]
    d = lab.sub(p, q)
    if d in REP_OF:
        return REP_OF[d]
    return None


def edge_decomposition(path: tuple[Point, ...]) -> tuple[Counter[Point], Counter[Point], Counter[Point]]:
    """Return total, spine, and bulk direction-pair counters."""
    path_edges = {
        frozenset((path[i], path[i + 1]))
        for i in range(len(path) - 1)
    }
    total: Counter[Point] = Counter()
    spine: Counter[Point] = Counter()
    bulk: Counter[Point] = Counter()
    for p, q in combinations(tuple(sorted(path)), 2):
        rep = edge_rep(p, q)
        if rep is None:
            continue
        total[rep] += 1
        if frozenset((p, q)) in path_edges:
            spine[rep] += 1
        else:
            bulk[rep] += 1
    return total, spine, bulk


def vector(counter: Counter[Point]) -> tuple[int, ...]:
    return tuple(counter[rep] for rep in REPS)


def vector_sub(a: Counter[Point], b: Counter[Point]) -> tuple[int, ...]:
    return tuple(a[rep] - b[rep] for rep in REPS)


def vector_sum(v: tuple[int, ...]) -> int:
    return sum(v)


def degree_hist(path: tuple[Point, ...]) -> Counter[int]:
    points = tuple(sorted(path))
    point_set = set(points)
    out: Counter[int] = Counter()
    for p in points:
        deg = sum(1 for u in lab.UNITS if lab.add(p, u) in point_set)
        out[deg] += 1
    return out


def deletion_hist(path: tuple[Point, ...]) -> tuple[Counter[int], Counter[int]]:
    points = tuple(sorted(path))
    edge_count = s628.unit_edges(points)
    losses: Counter[int] = Counter()
    remaining: Counter[int] = Counter()
    for p in points:
        child = tuple(q for q in points if q != p)
        loss = edge_count - s628.unit_edges(child)
        losses[loss] += 1
        remaining[edge_count - loss] += 1
    return losses, remaining


def classify_direction(total: int, spine: int, bulk: int) -> str:
    if total == 0:
        return "unused"
    if spine and bulk:
        return "mixed"
    if spine:
        return "pure_spine"
    return "pure_bulk"


def family(name: str, builder, max_m: int = 5) -> list[tuple[int, tuple[Counter[Point], Counter[Point], Counter[Point]]]]:
    rows = []
    print(f"{name} family")
    print("m | vertices | edges | spine | bulk | delta_total | delta_spine | delta_bulk")
    print("--- | --- | --- | --- | --- | --- | --- | ---")
    previous: tuple[Counter[Point], Counter[Point], Counter[Point]] | None = None
    for m in range(max_m + 1):
        path = builder(m)
        total, spine, bulk = edge_decomposition(path)
        edges = sum(total.values())
        spine_edges = sum(spine.values())
        bulk_edges = sum(bulk.values())
        if previous is None:
            delta_total = delta_spine = delta_bulk = "-"
        else:
            delta_total = str(vector_sub(total, previous[0]))
            delta_spine = str(vector_sub(spine, previous[1]))
            delta_bulk = str(vector_sub(bulk, previous[2]))
        print(
            f"{m} | {len(path)} | {edges} | {spine_edges} | {bulk_edges} | "
            f"{delta_total} | {delta_spine} | {delta_bulk}"
        )
        rows.append((m, (total, spine, bulk)))
        previous = (total, spine, bulk)
    print()
    return rows


def stable_quantum(rows: list[tuple[int, tuple[Counter[Point], Counter[Point], Counter[Point]]]]) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    deltas = []
    for i in range(2, len(rows)):
        prev = rows[i - 1][1]
        cur = rows[i][1]
        deltas.append(
            (
                vector_sub(cur[0], prev[0]),
                vector_sub(cur[1], prev[1]),
                vector_sub(cur[2], prev[2]),
            )
        )
    first = deltas[0]
    if any(delta != first for delta in deltas):
        raise AssertionError(deltas)
    return first


def print_quantum_summary(plus_rows, minus_rows) -> None:
    plus_q = stable_quantum(plus_rows)
    minus_q = stable_quantum(minus_rows)
    if plus_q != minus_q:
        raise AssertionError((plus_q, minus_q))

    total_q, spine_q, bulk_q = plus_q
    print("=" * 72)
    print("Fixed Moser layer quantum")
    print("=" * 72)
    print(f"direction pairs = {len(REPS)}")
    print(f"total_quantum = {total_q}, sum={vector_sum(total_q)}")
    print(f"spine_quantum = {spine_q}, sum={vector_sum(spine_q)}")
    print(f"bulk_quantum  = {bulk_q}, sum={vector_sum(bulk_q)}")
    print(f"identity: {vector_sum(total_q)} = {vector_sum(spine_q)} + {vector_sum(bulk_q)}")
    print(f"also: {vector_sum(total_q)} = 3 * {len(REPS)} direction pairs")
    uniform = tuple(3 for _ in REPS)
    defect = tuple(total_q[i] - uniform[i] for i in range(len(REPS)))
    print(f"uniform_3_defect = {defect}, L1={sum(abs(x) for x in defect)}")
    print()
    print("Direction pair representatives:")
    for idx, rep in enumerate(REPS):
        print(f"  d{idx}: {rep}")
    print()


def print_named_instance(name: str, path: tuple[Point, ...]) -> None:
    total, spine, bulk = edge_decomposition(path)
    edges = sum(total.values())
    spine_edges = sum(spine.values())
    bulk_edges = sum(bulk.values())
    losses, remaining = deletion_hist(path)
    print("=" * 72)
    print(f"{name}")
    print("=" * 72)
    print(f"vertices={len(path)} edges={edges} spine_edges={spine_edges} bulk_edges={bulk_edges}")
    print(f"degree_hist={dict(sorted(degree_hist(path).items()))}")
    print(f"deletion_loss_hist={dict(sorted(losses.items()))}")
    print(f"remaining_edge_hist_after_one_deletion={dict(sorted(remaining.items()))}")
    print()
    print("direction | rep | total | spine | bulk | class")
    print("--- | --- | --- | --- | --- | ---")
    pure_bulk = 0
    pure_spine = 0
    mixed = 0
    for idx, rep in enumerate(REPS):
        t = total[rep]
        s = spine[rep]
        b = bulk[rep]
        cls = classify_direction(t, s, b)
        if cls == "pure_bulk":
            pure_bulk += t
        elif cls == "pure_spine":
            pure_spine += t
        elif cls == "mixed":
            mixed += t
        print(f"d{idx} | {rep} | {t} | {s} | {b} | {cls}")
    print()
    print(f"pure_bulk_edges={pure_bulk}")
    print(f"pure_spine_edges={pure_spine}")
    print(f"mixed_edges={mixed}")
    print(f"pure_bulk_minus_spine_length={pure_bulk - spine_edges}")
    print()


@dataclass(frozen=True)
class Route:
    name: str
    preserves_ud_predicate: int
    carrier_transfer: int
    n22_use: int
    computable: int
    scalar_forgetfulness: int


ROUTES = [
    Route("fixed_27_edge_quantum", 5, 5, 5, 5, 1),
    Route("spine_bulk_direction_split", 5, 5, 5, 4, 1),
    Route("cap_endpoint_repair_channel", 5, 4, 5, 4, 1),
    Route("pure_bulk_direction_jackknife", 4, 5, 4, 5, 1),
    Route("degree_deletion_core_ledger", 4, 3, 5, 5, 1),
    Route("traceable_section_word", 4, 4, 4, 5, 2),
    Route("triangular_lattice_baseline", 3, 3, 2, 4, 3),
    Route("raw_edge_count_only", 1, 1, 2, 5, 5),
]


def route_score(r: Route) -> int:
    return (
        3 * r.preserves_ud_predicate
        + 2 * r.carrier_transfer
        + 2 * r.n22_use
        + r.computable
        - 2 * r.scalar_forgetfulness
    )


def beats(a: Route, b: Route) -> bool:
    sa = route_score(a)
    sb = route_score(b)
    if sa != sb:
        return sa > sb
    return a.name < b.name


def tournament_adj(routes: tuple[Route, ...]) -> list[list[int]]:
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and beats(routes[i], routes[j]):
                adj[i][j] = 1
    return adj


def directed_3cycles(adj: list[list[int]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            total += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            total += 1
    return total


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            count = dp[mask][v]
            if not count:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += count
    return sum(dp[-1])


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, reverse: bool = False) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u in range(n):
                edge = adj[u][v] if reverse else adj[v][u]
                if edge and u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen

    remaining = set(range(n))
    sizes = []
    while remaining:
        root = next(iter(remaining))
        comp = reach(root) & reach(root, reverse=True)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def print_tournament_analysis() -> None:
    print("=" * 72)
    print("Tournament Analysis")
    print("=" * 72)
    adj = tournament_adj(tuple(ROUTES))
    scores = sorted(sum(row) for row in adj)
    print("Pairwise observable:")
    print("  Which retained unit-distance carrier best preserves dense-core and")
    print("  possible n=22 proof obligations under quotienting.")
    print("Switch/gauge:")
    print("  fixed-carrier transfer from S646, n=22 obstruction value,")
    print("  computability, and scalar-forgetfulness penalty.")
    print("Challenged assumption:")
    print("  Vertices need not be points or unit edges.  Considered points,")
    print("  unit edges, direction pairs, slabs, caps, ears, deletion cores,")
    print("  obstruction labels, and proof obligations.  This scout uses")
    print("  carrier lenses as vertices; it preserves section/bulk proof data")
    print("  while destroying raw point order.")
    order = sorted(ROUTES, key=lambda r: (-route_score(r), r.name))
    print("Tie Hamiltonian path:")
    print("  " + " > ".join(r.name for r in order))
    print(f"score_hist={dict(Counter(scores))}")
    print(f"directed_3_cycles={directed_3cycles(adj)}")
    print(f"scc_sizes={scc_sizes(adj)}")
    print(f"hamiltonian_paths={hamiltonian_paths(adj)}")
    print("majority_score_order=" + " > ".join(r.name for r in order))
    print()


def main() -> None:
    print("=" * 72)
    print("Unit-distance Moser fixed-quantum carrier S648")
    print("=" * 72)
    print("S646 fixed the LRC C=27 gcd-shell clock.  The unit-distance")
    print("analogue is the THM-408 Moser add-one-slab operator.")
    print()

    plus_rows = family("P_m^+", s628.path_plus)
    minus_rows = family("P_m^-", s628.path_minus)
    print_quantum_summary(plus_rows, minus_rows)

    print_named_instance("P_2^- exact n=21 Moser row", s628.path_minus(2))
    print_named_instance("P_2^+ n=22 60-edge Moser lane", s628.path_plus(2))

    print("=" * 72)
    print("n=22 proof-facing interpretation")
    print("=" * 72)
    print("S614's graph-only deletion reduction says a 61-edge 22-point row")
    print("cannot have a degree <=3 vertex, because deleting it would leave")
    print("more than u(21)=57 edges.  P_2^+ has one degree-3 cap vertex.")
    print("Therefore a one-edge improvement in this carrier must repair that")
    print("cap endpoint channel, not merely add anonymous bulk.")
    print()
    print("The n=21 and n=22 m=2 rows also have a striking side-channel")
    print("balance: pure-bulk direction edges equal the unit-spine length.")
    print("This equality is not a general m identity; it singles out the")
    print("frontier m=2 rows as especially balanced section/bulk carriers.")
    print()

    print_tournament_analysis()

    print("=" * 72)
    print("Research takeaways")
    print("=" * 72)
    print("1. The unit-distance analogue of S646 is a fixed Moser layer")
    print("   quantum: after the cap transient, every added slab contributes")
    print("   the same direction vector with total 27 edges.")
    print("2. The 27 quantum splits as 8 spine edges plus 19 hidden bulk")
    print("   edges.  A point-order tournament can forget this split; a proof")
    print("   carrier should retain it.")
    print("3. At n=21 and n=22, pure-bulk direction mass equals the spine")
    print("   length.  That marks the current frontier as a balanced section/")
    print("   bulk fixed pocket, not just a high edge count.")
    print("4. The next no-leak target is cap-endpoint conservativity: a")
    print("   hypothetical 61-edge n=22 row must repair the degree-3 cap")
    print("   endpoint or leave the fixed Moser carrier.")


if __name__ == "__main__":
    main()
