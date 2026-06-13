#!/usr/bin/env python3
"""S669: equinumerosity versus outer-extension usability.

Prompt:

    outer extension usability theorem, outer extension embedding theorem;
    Friedman finite trees/homeomorphic embeddings; open intervals and Euclidean
    spaces are equinumerous; tournaments are equinumerous with even graphs.

The finite toy keeps cardinality deliberately fixed.  A path with 27 vertices,
a 3x9 grid, and a 3x3x3 grid all have the same number of points, just as
R, R^2, and R^3 all have continuum cardinality.  The target predicate is the
dimension label.  Cardinality alone leaks; deletion/link/growth/embedding
profiles repair the leak.

Tournament Analysis:
  Vertices are proof-transfer channels, not points, runners, arcs, or graphs.
  Pairwise observable is majority over (cardinal fidelity, extension purity,
  dimension retention, tournament fiber retention, Friedman embedding fit,
  actionability).  Ties follow the listed Hamiltonian priority order.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations, product


Graph = dict[int, set[int]]


def path_graph(n: int) -> Graph:
    g = {i: set() for i in range(n)}
    for i in range(n - 1):
        g[i].add(i + 1)
        g[i + 1].add(i)
    return g


def grid_graph(shape: tuple[int, ...]) -> Graph:
    coords = list(product(*[range(s) for s in shape]))
    index = {c: i for i, c in enumerate(coords)}
    g = {i: set() for i in range(len(coords))}
    for c in coords:
        i = index[c]
        for axis, limit in enumerate(shape):
            if c[axis] + 1 < limit:
                d = list(c)
                d[axis] += 1
                j = index[tuple(d)]
                g[i].add(j)
                g[j].add(i)
    return g


def component_count_after_deleting(g: Graph, removed: int | None = None) -> int:
    live = set(g)
    if removed is not None:
        live.discard(removed)
    if not live:
        return 0
    components = 0
    while live:
        components += 1
        start = live.pop()
        q = deque([start])
        while q:
            v = q.popleft()
            for w in g[v]:
                if w in live:
                    live.remove(w)
                    q.append(w)
    return components


def degree_hist(g: Graph) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(len(nbrs) for nbrs in g.values()).items()))


def deletion_hist(g: Graph) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(component_count_after_deleting(g, v) for v in g).items()))


def ball_count(g: Graph, root: int, radius: int) -> int:
    seen = {root}
    q = deque([(root, 0)])
    while q:
        v, d = q.popleft()
        if d == radius:
            continue
        for w in g[v]:
            if w not in seen:
                seen.add(w)
                q.append((w, d + 1))
    return len(seen)


def growth_profile(g: Graph) -> tuple[tuple[tuple[int, int, int], int], ...]:
    profiles = Counter((ball_count(g, v, 0), ball_count(g, v, 1), ball_count(g, v, 2)) for v in g)
    return tuple(sorted(profiles.items()))


def has_four_cycle(g: Graph) -> bool:
    vertices = list(g)
    for a, b in combinations(vertices, 2):
        if len(g[a] & g[b]) >= 2:
            return True
    return False


def small_embedding_profile(g: Graph) -> tuple[object, ...]:
    max_degree = max(len(nbrs) for nbrs in g.values())
    max_ball2 = max(ball_count(g, v, 2) for v in g)
    return (
        ("P3", max_degree >= 2),
        ("branch3", max_degree >= 3),
        ("branch5", max_degree >= 5),
        ("C4", has_four_cycle(g)),
        ("puncture_splits", any(k > 1 for k, _ in deletion_hist(g))),
        ("max_ball2", max_ball2),
    )


def dimension_signature(g: Graph) -> tuple[object, ...]:
    """Boundary-stable finite shadow of local dimension for this toy."""
    prof = dict(small_embedding_profile(g))
    return (
        ("branch3", prof["branch3"]),
        ("branch5", prof["branch5"]),
        ("C4", prof["C4"]),
        ("puncture_splits", prof["puncture_splits"]),
    )


@dataclass(frozen=True)
class Carrier:
    name: str
    dimension: int
    graph: Graph
    extension_name: str
    extension: Graph


def carriers() -> list[Carrier]:
    return [
        Carrier("line_P27", 1, path_graph(27), "line_P29", path_graph(29)),
        Carrier("plane_3x9", 2, grid_graph((3, 9)), "plane_5x9", grid_graph((5, 9))),
        Carrier("space_3x3x3", 3, grid_graph((3, 3, 3)), "space_5x3x3", grid_graph((5, 3, 3))),
    ]


def audit(carriers_: list[Carrier], key_name: str, key_fn) -> tuple[int, int, int, dict[object, set[int]]]:
    groups: dict[object, set[int]] = defaultdict(set)
    sizes: Counter[object] = Counter()
    for c in carriers_:
        key = key_fn(c)
        groups[key].add(c.dimension)
        sizes[key] += 1
    mixed = sum(1 for dims in groups.values() if len(dims) > 1)
    max_bucket = max(sizes.values(), default=0)
    return len(groups), mixed, max_bucket, groups


def channel_rows(carriers_: list[Carrier]) -> list[tuple[str, int, int, int]]:
    channels = [
        ("cardinality", lambda c: len(c.graph)),
        ("cardinality+deletion", lambda c: (len(c.graph), deletion_hist(c.graph))),
        ("cardinality+degree", lambda c: (len(c.graph), degree_hist(c.graph))),
        ("growth_profile", lambda c: growth_profile(c.graph)),
        ("small_embedding_profile", lambda c: small_embedding_profile(c.graph)),
    ]
    out = []
    for name, fn in channels:
        groups, mixed, max_bucket, _ = audit(carriers_, name, fn)
        out.append((name, groups, mixed, max_bucket))
    return out


S617_COUNTS = [
    (3, 4, 2, 2, 2, 2),
    (4, 11, 4, 4, 3, 7),
    (5, 34, 12, 12, 7, 22),
    (6, 156, 56, 56, 16, 100),
    (7, 1044, 456, None, None, 588),
    (8, 12346, 6880, None, None, 5466),
    (9, 274668, 191536, None, None, 83132),
    (10, 12005168, 9733056, None, None, 2272112),
    (11, 1018997864, 903753248, None, None, 115244616),
    (12, 165091172592, 154108311168, None, None, 10982861424),
]


S617_SPLITS = {
    5: {9: [0, 1]},
    6: {17: [0, 1], 23: [0, 1], 33: [0, 1], 37: [0, 1], 45: [0, 1]},
}


@dataclass(frozen=True)
class Lane:
    name: str
    vector: tuple[int, int, int, int, int, int]


def tournament_edges(lanes: list[Lane]) -> dict[tuple[int, int], int]:
    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(len(lanes)), 2):
        vi = lanes[i].vector
        vj = lanes[j].vector
        wins_i = sum(a > b for a, b in zip(vi, vj))
        wins_j = sum(b > a for a, b in zip(vi, vj))
        if wins_i > wins_j:
            winner = i
        elif wins_j > wins_i:
            winner = j
        else:
            winner = min(i, j)
        edges[(i, j)] = winner
    return edges


def score_hist(lanes: list[Lane], edges: dict[tuple[int, int], int]) -> dict[int, int]:
    scores = [0] * len(lanes)
    for (i, j), winner in edges.items():
        scores[winner] += 1
    return dict(sorted(Counter(scores).items()))


def directed_3cycles(lanes: list[Lane], edges: dict[tuple[int, int], int]) -> int:
    def beats(a: int, b: int) -> bool:
        i, j = sorted((a, b))
        return edges[(i, j)] == a

    total = 0
    for a, b, c in combinations(range(len(lanes)), 3):
        if beats(a, b) and beats(b, c) and beats(c, a):
            total += 1
        if beats(a, c) and beats(c, b) and beats(b, a):
            total += 1
    return total


def scc_sizes(lanes: list[Lane], edges: dict[tuple[int, int], int]) -> list[int]:
    n = len(lanes)
    adj = [[] for _ in range(n)]
    radj = [[] for _ in range(n)]
    for (i, j), winner in edges.items():
        loser = j if winner == i else i
        adj[winner].append(loser)
        radj[loser].append(winner)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: int) -> int:
        seen.add(v)
        total = 1
        for w in radj[v]:
            if w not in seen:
                total += rdfs(w)
        return total

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(lanes: list[Lane], edges: dict[tuple[int, int], int]) -> int:
    n = len(lanes)

    def beats(a: int, b: int) -> bool:
        i, j = sorted((a, b))
        return edges[(i, j)] == a

    dp: dict[tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp.get((mask, last), 0)
            if not cur:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(last, nxt):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + cur
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def print_tournament() -> None:
    lanes = [
        Lane("address_plus_embedding", (4, 5, 5, 4, 5, 5)),
        Lane("small_embedding_profile", (4, 5, 5, 3, 5, 5)),
        Lane("local_growth_profile", (3, 5, 5, 1, 3, 4)),
        Lane("(H,beta1,packet)_fiber", (4, 4, 1, 5, 3, 4)),
        Lane("deletion_profile", (3, 3, 4, 1, 2, 4)),
        Lane("Royle_even_cardinal_theorem", (5, 1, 0, 1, 0, 3)),
        Lane("bijection_without_profile", (5, 1, 0, 1, 0, 2)),
        Lane("raw_cardinality", (5, 0, 0, 0, 0, 1)),
    ]
    edges = tournament_edges(lanes)
    scores = [0] * len(lanes)
    for (_, _), winner in edges.items():
        scores[winner] += 1

    print("F. Tournament Analysis over transfer channels")
    print("  vertices=proof-transfer channels")
    print("  observable=(cardinal fidelity, extension purity, dimension retention, tournament fiber retention, Friedman embedding fit, actionability)")
    print("  switch=majority; tie Hamiltonian path=listed priority order")
    print(f"  score_hist={score_hist(lanes, edges)}")
    print(f"  directed_3cycles={directed_3cycles(lanes, edges)}")
    print(f"  scc_sizes={scc_sizes(lanes, edges)}")
    print(f"  hamiltonian_paths={hamiltonian_path_count(lanes, edges)}")
    print("  top_order:")
    for score, lane in sorted(zip(scores, lanes), key=lambda item: (-item[0], item[1].name)):
        print(f"    {lane.name:34s} score={score} vector={lane.vector}")


def main() -> None:
    cs = carriers()
    print("=" * 78)
    print("S669 equinumerosity versus outer-extension usability")
    print("=" * 78)
    print()
    print("Source anchors")
    print("  Friedman official pages: homeomorphic embeddings, internal finite tree embeddings, finite functions / large cardinals, invariant maximality.")
    print("  Royle et al.: tournaments and even graphs are equinumerous.")
    print("  Topology anchor: cardinality does not preserve dimension/homeomorphism; deletion/link profiles are finite shadows.")
    print()

    print("A. Equal-cardinality Euclidean toy")
    for c in cs:
        print(
            f"  {c.name:14s} dim={c.dimension} |V|={len(c.graph):2d} "
            f"degree={degree_hist(c.graph)} deletion={deletion_hist(c.graph)} "
            f"profile={small_embedding_profile(c.graph)}"
        )
    print()
    print("  Interpretation: all three carriers have the same finite cardinal shadow |V|=27.")
    print("  The 1D carrier is puncture-disconnected; 2D/3D are not.  Degree and growth profiles separate all three.")
    print()

    print("B. Usability audit for the dimension predicate")
    print(f"{'channel':30s} {'groups':>6s} {'mixed':>6s} {'max_bucket':>10s}")
    for name, groups, mixed, max_bucket in channel_rows(cs):
        print(f"{name:30s} {groups:6d} {mixed:6d} {max_bucket:10d}")
    print("  Reading: cardinality alone has one mixed bucket containing dimensions 1,2,3.")
    print("  Deletion repairs the 1D-vs-higher split but still mixes 2D/3D; local growth or embedding profiles repair the toy completely.")
    print()

    print("C. Outer-extension check")
    print("  extension operation=add an outer layer while preserving the internal carrier type")
    for c in cs:
        before = small_embedding_profile(c.graph)
        after = small_embedding_profile(c.extension)
        signature_before = dimension_signature(c.graph)
        signature_after = dimension_signature(c.extension)
        print(
            f"  {c.name:14s} -> {c.extension_name:14s} "
            f"|V| {len(c.graph):2d}->{len(c.extension):2d} dim={c.dimension} "
            f"full_profile_changed={before != after} "
            f"dimension_signature_changed={signature_before != signature_after}"
        )
    print("  Reading: extension changes raw size and sometimes boundary/growth counts.")
    print("  The dimension signature remains stable, which is the usability coordinate the cardinal shadow lacked.")
    print()

    print("D. Royle/even-graph count shadow, imported from S617")
    print(f"{'n':>4s} {'graphs':>12s} {'tournaments':>12s} {'degree_even':>12s} {'naive_even':>12s}")
    for n, graphs, t_burnside, _t_enum, degree_even, naive_even in S617_COUNTS:
        degree = "-" if degree_even is None else str(degree_even)
        naive = "-" if naive_even is None else str(naive_even)
        print(f"{n:4d} {graphs:12d} {t_burnside:12d} {degree:>12s} {naive:>12s}")
    print("  S617 finding: degree-even and naive Burnside-even are not the Royle-even theorem.")
    print("  Even when a true equinumerous class exists, HYP-2187 still asks for a fiber functor preserving H/beta1/packet data.")
    print("  H splits by beta1 in the stored audit:")
    for n, splits in sorted(S617_SPLITS.items()):
        print(f"    n={n}: {splits}")
    print()

    print("E. S668 Friedman-tree compatibility")
    print("  S668 finite-tree target=contains rooted homeomorphic color chain 0 -> 1 -> 2.")
    print("  Coarse channels leaked: size_height mixed=6, color_hist=10, frontier=46, raw outer_address=55.")
    print("  Embedding channels repaired: small_embedding_profile=0, full_embedding_profile=0, address_plus_small_embedding=0 mixed fibers.")
    print("  S669 reading: cardinal/count equality is exactly the wrong level for outer-extension theorem transfer.")
    print()

    print_tournament()
    print()
    print("G. Repo theorem shape")
    print("  Equinumerosity-to-usability repair:")
    print("    A count-equivalent quotient can transfer a predicate only after attaching")
    print("    the smallest deletion/link/embedding/fiber profile that is pure under")
    print("    the allowed outer extensions.")
    print("  LRC14 use:")
    print("    visible Res_27/cardinal shadows are allowed to guide search, but proof")
    print("    transfer must live in owner-private deletion, carry cocycle, endpoint")
    print("    protector, and bounded proof-obligation embedding profiles.")


if __name__ == "__main__":
    main()
