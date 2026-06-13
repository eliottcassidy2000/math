#!/usr/bin/env python3
"""
equinum_equidecomp_bridge_s617.py

S617 scout: separate equinumerosity from equidecomposability.

The old repo trail has two nearby but different ideas:

* Equinumerosity: two quotient worlds have the same number of objects.
  Royle-even graphs and tournaments are the advertised cardinal shadow, while
  degree-even cycle-space graphs and the naive even-order Burnside complement
  are tempting but wrong substitutes.

* Equidecomposability: two tournaments are equivalent only after retaining
  enough scissors-style invariants.  Here we audit small tournament classes
  by H, by (H,beta1), and by a small odd-cycle packet polynomial.

Tournament Analysis is run on quotient lenses rather than runners or arcs.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from math import comb, factorial, gcd, lcm


PAIR_CACHE: dict[int, list[tuple[int, int]]] = {}
PERM_CACHE: dict[int, list[tuple[int, ...]]] = {}
TOURNAMENT_CLASS_CACHE: dict[int, list["TournamentClass"]] = {}


def pairs(n: int) -> list[tuple[int, int]]:
    if n not in PAIR_CACHE:
        PAIR_CACHE[n] = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return PAIR_CACHE[n]


def perms(n: int) -> list[tuple[int, ...]]:
    if n not in PERM_CACHE:
        PERM_CACHE[n] = list(permutations(range(n)))
    return PERM_CACHE[n]


def tournament_adj(mask: int, n: int) -> list[list[int]]:
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs(n)):
        if (mask >> k) & 1:
            a[i][j] = 1
        else:
            a[j][i] = 1
    return a


def graph_adj(mask: int, n: int) -> list[list[int]]:
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs(n)):
        if (mask >> k) & 1:
            a[i][j] = 1
            a[j][i] = 1
    return a


def canonical_form(a: list[list[int]], n: int) -> tuple[int, ...]:
    best: tuple[int, ...] | None = None
    for p in perms(n):
        form = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    assert best is not None
    return best


@dataclass(frozen=True)
class TournamentClass:
    canon: tuple[int, ...]
    adj: tuple[tuple[int, ...], ...]
    labeled_size: int


def tournament_classes(n: int) -> list[TournamentClass]:
    if n in TOURNAMENT_CLASS_CACHE:
        return TOURNAMENT_CLASS_CACHE[n]

    by_canon: dict[tuple[int, ...], list[object]] = {}
    for mask in range(1 << comb(n, 2)):
        a = tournament_adj(mask, n)
        c = canonical_form(a, n)
        if c not in by_canon:
            by_canon[c] = [tuple(tuple(row) for row in a), 0]
        by_canon[c][1] = int(by_canon[c][1]) + 1

    out = [
        TournamentClass(canon=c, adj=adj, labeled_size=int(size))
        for c, (adj, size) in by_canon.items()
    ]
    out.sort(key=lambda d: d.canon)
    TOURNAMENT_CLASS_CACHE[n] = out
    return out


def degree_even_graph_class_count(n: int) -> int:
    classes: set[tuple[int, ...]] = set()
    for mask in range(1 << comb(n, 2)):
        deg = [0] * n
        for k, (i, j) in enumerate(pairs(n)):
            if (mask >> k) & 1:
                deg[i] += 1
                deg[j] += 1
        if any(d % 2 for d in deg):
            continue
        classes.add(canonical_form(graph_adj(mask, n), n))
    return len(classes)


def int_partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if max_part is None:
        max_part = n
    if n == 0:
        return [()]
    out: list[tuple[int, ...]] = []
    for first in range(min(n, max_part), 0, -1):
        for rest in int_partitions(n - first, first):
            out.append((first,) + rest)
    return out


def conjugacy_class_size(n: int, cycle_type: tuple[int, ...]) -> int:
    counts = Counter(cycle_type)
    denom = 1
    for length, multiplicity in counts.items():
        denom *= (length**multiplicity) * factorial(multiplicity)
    return factorial(n) // denom


def edge_orbits(cycle_type: tuple[int, ...]) -> int:
    total = 0
    for i, a in enumerate(cycle_type):
        total += a // 2
        for b in cycle_type[i + 1 :]:
            total += gcd(a, b)
    return total


def perm_order(cycle_type: tuple[int, ...]) -> int:
    order = 1
    for part in cycle_type:
        order = lcm(order, part)
    return order


def burnside_counts(n: int) -> tuple[int, int, int]:
    """Return graph classes, tournament classes, naive even-order complement."""
    graph_sum = 0
    odd_order_sum = 0
    even_order_sum = 0
    for ct in int_partitions(n):
        class_size = conjugacy_class_size(n, ct)
        fixed_graphs = 2 ** edge_orbits(ct)
        contribution = class_size * fixed_graphs
        graph_sum += contribution
        if perm_order(ct) % 2:
            odd_order_sum += contribution
        else:
            even_order_sum += contribution
    nf = factorial(n)
    return graph_sum // nf, odd_order_sum // nf, even_order_sum // nf


def hamiltonian_paths(a: tuple[tuple[int, ...], ...]) -> int:
    n = len(a)
    total = 0
    for p in perms(n):
        if all(a[p[i]][p[i + 1]] for i in range(n - 1)):
            total += 1
    return total


def modular_rank(matrix: list[list[int]], prime: int) -> int:
    if not matrix or not matrix[0]:
        return 0
    a = [[x % prime for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot = None
        for r in range(rank, rows):
            if a[r][col] % prime:
                pivot = r
                break
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], prime - 2, prime)
        a[rank] = [(x * inv) % prime for x in a[rank]]
        for r in range(rows):
            if r == rank or not a[r][col]:
                continue
            factor = a[r][col]
            a[r] = [(a[r][c] - factor * a[rank][c]) % prime for c in range(cols)]
        rank += 1
        if rank == rows:
            break
    return rank


def rank_over_q(matrix: list[list[int]]) -> int:
    ranks = [modular_rank(matrix, p) for p in (1_000_000_007, 1_000_000_009)]
    if ranks[0] != ranks[1]:
        raise RuntimeError(f"rank changed across primes: {ranks}")
    return ranks[0]


def beta1(a: tuple[tuple[int, ...], ...]) -> int:
    n = len(a)
    edges = [(i, j) for i in range(n) for j in range(n) if i != j and a[i][j]]
    edge_index = {edge: i for i, edge in enumerate(edges)}
    ne = len(edges)

    d1 = [[0] * ne for _ in range(n)]
    for col, (u, v) in enumerate(edges):
        d1[v][col] += 1
        d1[u][col] -= 1

    paths2: list[tuple[int, int, int]] = []
    for u in range(n):
        for v in range(n):
            if u == v or not a[u][v]:
                continue
            for w in range(n):
                if w == u or w == v or not a[v][w]:
                    continue
                if a[u][w]:
                    paths2.append((u, v, w))

    d2 = [[0] * len(paths2) for _ in range(ne)]
    for col, (u, v, w) in enumerate(paths2):
        d2[edge_index[(v, w)]][col] += 1
        d2[edge_index[(u, w)]][col] -= 1
        d2[edge_index[(u, v)]][col] += 1

    return ne - rank_over_q(d1) - rank_over_q(d2)


def c3_independence_poly(a: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    n = len(a)
    c3s: list[frozenset[int]] = []
    for i, j, k in combinations(range(n), 3):
        forward = a[i][j] and a[j][k] and a[k][i]
        backward = a[i][k] and a[k][j] and a[j][i]
        if forward or backward:
            c3s.append(frozenset((i, j, k)))

    coeff = [1, len(c3s)]
    disjoint_pairs = 0
    for x, y in combinations(range(len(c3s)), 2):
        if not (c3s[x] & c3s[y]):
            disjoint_pairs += 1
    if disjoint_pairs:
        coeff.append(disjoint_pairs)
    return tuple(coeff)


@dataclass(frozen=True)
class ClassInvariants:
    h: int
    b1: int
    c3_poly: tuple[int, ...]
    labeled_size: int


def class_invariants(n: int) -> list[ClassInvariants]:
    rows: list[ClassInvariants] = []
    for cls in tournament_classes(n):
        rows.append(
            ClassInvariants(
                h=hamiltonian_paths(cls.adj),
                b1=beta1(cls.adj),
                c3_poly=c3_independence_poly(cls.adj),
                labeled_size=cls.labeled_size,
            )
        )
    return rows


def summarize_equicounts() -> None:
    print("Equinumerosity audit: three noninterchangeable quotients")
    print("  T_burnside is the A000568 tournament count from odd-order cycle types.")
    print("  degree-even means Euler/cycle-space graph classes, not Royle-even graphs.")
    print("  naive-even is the even-order Burnside complement; it is also not Royle-even.")
    print()
    print(
        f"  {'n':>2} {'graphs':>10} {'T_burnside':>12} "
        f"{'T_enum':>8} {'degree-even':>12} {'naive-even':>12}"
    )
    for n in range(3, 13):
        graphs, tourn, naive_even = burnside_counts(n)
        t_enum = len(tournament_classes(n)) if n <= 6 else None
        deg_even = degree_even_graph_class_count(n) if n <= 6 else None
        print(
            f"  {n:>2} {graphs:>10} {tourn:>12} "
            f"{str(t_enum) if t_enum is not None else '-':>8} "
            f"{str(deg_even) if deg_even is not None else '-':>12} "
            f"{naive_even:>12}"
        )
    print()
    print("  Findings:")
    print("    degree-even classes agree with tournaments only at n=3 in this range.")
    print("    naive-even equals graphs - tournaments, not tournaments, once n>=4.")
    print("    Therefore equinumerosity needs the Royle-even intrinsic graph property,")
    print("    not the cycle-space projection and not the Burnside parity split.")
    print()


def summarize_equidecomposability() -> None:
    print("Equidecomposability audit: retained tournament invariants")
    print("  H is the volume shadow. beta1 is the first Dehn/scissors obstruction.")
    print("  c3-poly is the disjoint directed-3-cycle packet polynomial 1+c1*x+c2*x^2.")
    print()
    print(
        f"  {'n':>2} {'iso':>5} {'H':>5} {'(H,b1)':>8} "
        f"{'(H,b1,c3)':>12} {'split H-values by beta1'}"
    )
    for n in range(3, 7):
        inv = class_invariants(n)
        h_values = {r.h for r in inv}
        hb_values = {(r.h, r.b1) for r in inv}
        hbc_values = {(r.h, r.b1, r.c3_poly) for r in inv}
        split: dict[int, list[int]] = {}
        for h in sorted(h_values):
            beta_values = sorted({r.b1 for r in inv if r.h == h})
            if len(beta_values) > 1:
                split[h] = beta_values
        print(
            f"  {n:>2} {len(inv):>5} {len(h_values):>5} {len(hb_values):>8} "
            f"{len(hbc_values):>12} {split}"
        )
    print()

    for n in (5, 6):
        inv = class_invariants(n)
        print(f"  n={n} split ledger")
        by_h: dict[int, list[ClassInvariants]] = defaultdict(list)
        for row in inv:
            by_h[row.h].append(row)
        for h in sorted(by_h):
            beta_counts = Counter(row.b1 for row in by_h[h])
            if len(beta_counts) > 1:
                labeled = {
                    b: sum(row.labeled_size for row in by_h[h] if row.b1 == b)
                    for b in sorted(beta_counts)
                }
                print(
                    f"    H={h}: beta1 iso-counts={dict(sorted(beta_counts.items()))}, "
                    f"labeled-mass={labeled}"
                )
        print()

    print("  Findings:")
    print("    At n=5, H=9 already splits into beta1=0 and beta1=1 scissors classes.")
    print("    At n=6, five H-values split by beta1: 17, 23, 33, 37, 45.")
    print("    Adding the c3 packet polynomial refines 24 (H,beta1) classes to 29.")
    print("    Thus equal H is equinumerosity-like volume data, not equidecomposition.")
    print()


@dataclass(frozen=True)
class Lens:
    name: str
    cardinal_evidence: int
    constructive_map: int
    preserves_h: int
    scissors_refinement: int
    cross_domain: int
    low_risk: int


LENSES = (
    Lens("Royle-even cardinal theorem", 5, 1, 1, 1, 5, 2),
    Lens("Degree-even cycle projection", 2, 5, 2, 2, 4, 3),
    Lens("Naive Burnside parity split", 1, 4, 0, 0, 2, 1),
    Lens("H-volume quotient", 4, 5, 5, 1, 4, 4),
    Lens("(H,beta1) scissors quotient", 4, 4, 5, 5, 5, 4),
    Lens("(H,beta1,c3) packet quotient", 4, 4, 5, 5, 4, 3),
    Lens("Full tournament isomorphism", 5, 5, 5, 5, 2, 5),
)


def lens_votes(a: Lens, b: Lens, include_risk: bool = True) -> tuple[int, int]:
    criteria = [
        (a.cardinal_evidence > b.cardinal_evidence, b.cardinal_evidence > a.cardinal_evidence),
        (a.constructive_map > b.constructive_map, b.constructive_map > a.constructive_map),
        (a.preserves_h > b.preserves_h, b.preserves_h > a.preserves_h),
        (a.scissors_refinement > b.scissors_refinement, b.scissors_refinement > a.scissors_refinement),
        (a.cross_domain > b.cross_domain, b.cross_domain > a.cross_domain),
    ]
    if include_risk:
        criteria.append((a.low_risk > b.low_risk, b.low_risk > a.low_risk))
    av = sum(1 for x, y in criteria if x and not y)
    bv = sum(1 for x, y in criteria if y and not x)
    return av, bv


def lens_tournament(include_risk: bool = True) -> list[list[int]]:
    n = len(LENSES)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        iv, jv = lens_votes(LENSES[i], LENSES[j], include_risk=include_risk)
        if iv > jv or (iv == jv and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def score_hist(adj: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for row in adj:
        score = sum(row)
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            total += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, forward: bool) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in range(n):
                edge = adj[v][w] if forward else adj[w][v]
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    left = set(range(n))
    sizes: list[int] = []
    while left:
        v = next(iter(left))
        comp = reach(v, True) & reach(v, False)
        sizes.append(len(comp))
        left -= comp
    return sorted(sizes, reverse=True)


def hp_count(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            val = dp[mask][v]
            if not val:
                continue
            for w in range(n):
                if not (mask >> w) & 1 and adj[v][w]:
                    dp[mask | (1 << w)][w] += val
    return sum(dp[-1])


def summarize_tournament_analysis() -> None:
    print("Tournament Analysis on quotient lenses")
    print("  Vertices are quotient choices, not runners, arcs, or graph vertices.")
    print("  Pairwise observable: majority over cardinal evidence, constructivity,")
    print("  H preservation, scissors refinement, cross-domain reach, and risk.")
    adj = lens_tournament(include_risk=True)
    no_risk = lens_tournament(include_risk=False)
    flips = sum(
        1
        for i, j in combinations(range(len(LENSES)), 2)
        if adj[i][j] != no_risk[i][j]
    )
    scores = [sum(row) for row in adj]
    for lens, score in sorted(zip(LENSES, scores), key=lambda item: (-item[1], item[0].name)):
        print(
            f"  score={score} lens={lens.name}; "
            f"features=(card={lens.cardinal_evidence}, map={lens.constructive_map}, "
            f"H={lens.preserves_h}, scissors={lens.scissors_refinement}, "
            f"cross={lens.cross_domain}, lowrisk={lens.low_risk})"
        )
    print(f"  score histogram: {score_hist(adj)}")
    print(f"  directed 3-cycles: {directed_triangles(adj)}")
    print(f"  SCC sizes: {scc_sizes(adj)}")
    print(f"  Hamiltonian path count: {hp_count(adj)}")
    print(f"  edge flips when risk is removed: {flips}")
    print("  Reading: full isomorphism is safest but too fine; (H,beta1) is the")
    print("    best working scissors quotient; Royle-even is cardinal signal only")
    print("    until a constructive predicate-preserving map is supplied.")
    print()


def assumption_challenge() -> None:
    print("Assumption challenge")
    print("  Do not make tournament vertices be runners or arcs by reflex.")
    print("  Alternative vertex sets considered:")
    print("    Burnside cycle types, tournament iso classes, degree-even graph classes,")
    print("    Royle-even graph classes, H-fibers, beta1 chains, disjoint c3 packets,")
    print("    OCF proof obligations, and quotient lenses.")
    print("  Preserved predicate:")
    print("    Equinumerosity preserves only a count unless the map also preserves H,")
    print("    beta1, packet data, or the target scissors class.")
    print("  Destroyed information:")
    print("    Degree-even projection forgets cut/score data; H forgets beta1;")
    print("    (H,beta1) still forgets some c3 packet structure at n=6.")
    print("  Challenged assumption:")
    print("    A000568-sized quotients are not automatically equivalent proof objects.")
    print()


def main() -> None:
    print("==== S617 equinumerosity / equidecomposability bridge ====")
    summarize_equicounts()
    summarize_equidecomposability()
    summarize_tournament_analysis()
    assumption_challenge()
    print("Hypothesis HYP-2187:")
    print("  Equinumerosity is a cardinal shadow; equidecomposability is a retained")
    print("  invariant quotient.  The bridge must be a fiber functor that records which")
    print("  side channels survive: score/cut data, beta1, and odd-cycle packets.")


if __name__ == "__main__":
    main()
