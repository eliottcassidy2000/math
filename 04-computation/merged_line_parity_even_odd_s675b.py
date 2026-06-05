#!/usr/bin/env python3
"""S675b: black-even / blue-odd parity in the merged tournament line lift.

The user conjectured that the black-line part of the merged tournament
metagraph is an even graph, and asked whether the blue-line part is always an
odd graph.  The repository has two nearby conventions:

1. The corrected tiling-complement convention: a blue/black line pairs a
   fixed-base tiling with its all-free-tile complement.  This line lift carries
   multiplicity and self-loops.
2. The older class-level convention: a single-flip merged metagraph edge is
   blue when its endpoints have the same SC/non-SC type and black otherwise.

This script audits both through n=7.  The parity theorem lives in (1), not in
the simple graph projection (2).
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations, permutations
from math import comb


A000568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
KNOWN_MERGED_EDGES = {3: 1, 4: 3, 5: 21, 6: 143, 7: 2123}


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


@lru_cache(maxsize=None)
def pair_list(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(n), 2))


@lru_cache(maxsize=None)
def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {pair: i for i, pair in enumerate(pair_list(n))}


def bit_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return pair_index(n)[(i, j)]


def has_arc(mask: int, n: int, i: int, j: int) -> bool:
    if i == j:
        return False
    if i < j:
        return ((mask >> bit_index(n, i, j)) & 1) == 1
    return ((mask >> bit_index(n, j, i)) & 1) == 0


@lru_cache(maxsize=None)
def outdegrees(mask: int, n: int) -> tuple[int, ...]:
    scores = [0] * n
    for i, j in pair_list(n):
        if has_arc(mask, n, i, j):
            scores[i] += 1
        else:
            scores[j] += 1
    return tuple(scores)


@lru_cache(maxsize=None)
def grouped_perms(scores: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    groups: dict[int, list[int]] = defaultdict(list)
    for v, score in enumerate(scores):
        groups[score].append(v)
    blocks = [tuple(groups[s]) for s in sorted(groups)]

    def rec(k: int):
        if k == len(blocks):
            yield ()
            return
        for p in permutations(blocks[k]):
            for tail in rec(k + 1):
                yield p + tail

    return tuple(rec(0))


def encode_under_perm(mask: int, n: int, perm: tuple[int, ...]) -> int:
    """New label a is old vertex perm[a]."""
    out = 0
    for a, b in pair_list(n):
        if has_arc(mask, n, perm[a], perm[b]):
            out |= 1 << bit_index(n, a, b)
    return out


@lru_cache(maxsize=None)
def canon_tournament(mask: int, n: int) -> int:
    best = None
    for perm in grouped_perms(outdegrees(mask, n)):
        code = encode_under_perm(mask, n, perm)
        if best is None or code < best:
            best = code
    assert best is not None
    return best


@lru_cache(maxsize=None)
def free_pairs(n: int) -> tuple[tuple[int, int], ...]:
    """Free arcs after fixing the Hamiltonian path 0->1->...->n-1."""
    return tuple((i, j) for i, j in pair_list(n) if j != i + 1)


def tournament_from_free(n: int, free_mask: int) -> int:
    mask = 0
    for i in range(n - 1):
        mask |= 1 << bit_index(n, i, i + 1)
    for k, (i, j) in enumerate(free_pairs(n)):
        if (free_mask >> k) & 1:
            mask |= 1 << bit_index(n, i, j)
    return mask


@lru_cache(maxsize=None)
def reflection_map(n: int) -> tuple[int, ...]:
    index = {pair: i for i, pair in enumerate(free_pairs(n))}
    reflected = []
    for i, j in free_pairs(n):
        ri, rj = n - 1 - j, n - 1 - i
        if ri > rj:
            ri, rj = rj, ri
        reflected.append(index[(ri, rj)])
    return tuple(reflected)


def is_grid_symmetric(n: int, free_mask: int) -> bool:
    refl = reflection_map(n)
    for k, rk in enumerate(refl):
        if ((free_mask >> k) & 1) != ((free_mask >> rk) & 1):
            return False
    return True


def even_graph_edges_from_free(n: int, free_mask: int) -> frozenset[tuple[int, int]]:
    """A labelled linear bijection F_2^C(n-1,2) -> even graphs on n vertices."""
    nonroot = tuple(combinations(range(1, n), 2))
    edges: set[tuple[int, int]] = set()
    degree = [0] * n
    for k, (i, j) in enumerate(nonroot):
        if (free_mask >> k) & 1:
            edges.add((i, j))
            degree[i] ^= 1
            degree[j] ^= 1
    for v in range(1, n):
        if degree[v]:
            edges.add((0, v))
            degree[0] ^= 1
            degree[v] ^= 1
    assert all(d == 0 for d in degree)
    return frozenset(edges)


def xor_edges(a: frozenset[tuple[int, int]], b: frozenset[tuple[int, int]]):
    return frozenset(set(a) ^ set(b))


def add_multigraph_degree(degree: list[int], a: int, b: int) -> None:
    if a == b:
        degree[a] += 2
    else:
        degree[a] += 1
        degree[b] += 1


def build_fixed_base_quotient(n: int) -> dict[str, object]:
    m = len(free_pairs(n))
    total = 1 << m
    class_id: dict[int, int] = {}
    class_of_free: dict[int, int] = {}
    class_rep_canon: dict[int, int] = {}
    members: dict[int, list[int]] = defaultdict(list)

    for free_mask in range(total):
        canon = canon_tournament(tournament_from_free(n, free_mask), n)
        if canon not in class_id:
            cid = len(class_id)
            class_id[canon] = cid
            class_rep_canon[cid] = canon
        cid = class_id[canon]
        class_of_free[free_mask] = cid
        members[cid].append(free_mask)

    all_arc_mask = (1 << comb(n, 2)) - 1
    complement_of: dict[int, int] = {}
    is_sc: dict[int, bool] = {}
    for cid, canon in class_rep_canon.items():
        comp_canon = canon_tournament(canon ^ all_arc_mask, n)
        complement_of[cid] = class_id[comp_canon]
        is_sc[cid] = complement_of[cid] == cid

    merged: dict[int, int] = {}
    mid = 0
    for cid in range(len(class_rep_canon)):
        if cid in merged:
            continue
        merged[cid] = mid
        merged[complement_of[cid]] = mid
        mid += 1

    merged_sc = [False] * mid
    for cid, node in merged.items():
        merged_sc[node] = is_sc[cid]

    return {
        "n": n,
        "m": m,
        "total": total,
        "class_of_free": class_of_free,
        "members": members,
        "class_count": len(class_rep_canon),
        "merged": merged,
        "merged_count": mid,
        "merged_sc": tuple(merged_sc),
    }


def audit_line_lift(data: dict[str, object]) -> dict[str, object]:
    n = int(data["n"])
    m = int(data["m"])
    total = int(data["total"])
    full_free_mask = total - 1
    class_of_free = data["class_of_free"]
    merged = data["merged"]
    merged_count = int(data["merged_count"])

    blue_degree = [0] * merged_count
    black_degree = [0] * merged_count
    blue_edges = Counter()
    black_edges = Counter()
    blue_lines = black_lines = 0
    grid_symmetric_tilings = 0

    for free_mask in range(total):
        comp = free_mask ^ full_free_mask
        if free_mask > comp:
            continue
        a = merged[class_of_free[free_mask]]
        b = merged[class_of_free[comp]]
        edge = (a, b) if a <= b else (b, a)
        if is_grid_symmetric(n, free_mask):
            grid_symmetric_tilings += 2
            blue_lines += 1
            blue_edges[edge] += 1
            add_multigraph_degree(blue_degree, a, b)
        else:
            black_lines += 1
            black_edges[edge] += 1
            add_multigraph_degree(black_degree, a, b)

    active_blue = [d for d in blue_degree if d > 0]
    active_black = [d for d in black_degree if d > 0]
    return {
        "grid_symmetric_tilings": grid_symmetric_tilings,
        "blue_lines": blue_lines,
        "black_lines": black_lines,
        "blue_active_nodes": len(active_blue),
        "black_active_nodes": len(active_black),
        "blue_all_active_odd": all(d % 2 == 1 for d in active_blue),
        "black_all_even": all(d % 2 == 0 for d in black_degree),
        "blue_degree_parity": Counter(d % 2 for d in active_blue),
        "black_degree_parity": Counter(d % 2 for d in black_degree),
        "blue_loop_edges": sum(v for (a, b), v in blue_edges.items() if a == b),
        "black_loop_edges": sum(v for (a, b), v in black_edges.items() if a == b),
    }


def audit_simple_sc_type_edges(data: dict[str, object]) -> dict[str, object]:
    n = int(data["n"])
    m = int(data["m"])
    total = int(data["total"])
    class_of_free = data["class_of_free"]
    merged = data["merged"]
    merged_sc = data["merged_sc"]
    merged_count = int(data["merged_count"])

    unmerged_edges: set[tuple[int, int]] = set()
    for free_mask in range(total):
        a = class_of_free[free_mask]
        for k in range(m):
            b = class_of_free[free_mask ^ (1 << k)]
            if a != b:
                unmerged_edges.add((a, b) if a < b else (b, a))

    merged_edges: set[tuple[int, int]] = set()
    collapsed = 0
    for a, b in unmerged_edges:
        ma, mb = merged[a], merged[b]
        if ma == mb:
            collapsed += 1
        else:
            merged_edges.add((ma, mb) if ma < mb else (mb, ma))

    blue_degree = [0] * merged_count
    black_degree = [0] * merged_count
    blue_edges = black_edges = 0
    for a, b in merged_edges:
        if merged_sc[a] == merged_sc[b]:
            blue_edges += 1
            blue_degree[a] += 1
            blue_degree[b] += 1
        else:
            black_edges += 1
            black_degree[a] += 1
            black_degree[b] += 1

    black_odd = [i for i, d in enumerate(black_degree) if d % 2]
    blue_even_active = [i for i, d in enumerate(blue_degree) if d > 0 and d % 2 == 0]
    return {
        "unmerged_edges": len(unmerged_edges),
        "merged_edges": len(merged_edges),
        "collapsed": collapsed,
        "blue_edges": blue_edges,
        "black_edges": black_edges,
        "black_all_even": not black_odd,
        "active_blue_all_odd": not blue_even_active,
        "black_odd_witnesses": black_odd[:8],
        "blue_even_active_witnesses": blue_even_active[:8],
    }


def audit_even_graph_carrier(n: int) -> dict[str, object]:
    m = len(free_pairs(n))
    full = (1 << m) - 1
    zero_graph = even_graph_edges_from_free(n, 0)
    one_graph = even_graph_edges_from_free(n, full)
    toggler = xor_edges(zero_graph, one_graph)
    degrees = [0] * n
    for i, j in toggler:
        degrees[i] += 1
        degrees[j] += 1
    sample_ok = True
    for x in range(min(1 << m, 256)):
        lhs = xor_edges(even_graph_edges_from_free(n, x), even_graph_edges_from_free(n, x ^ full))
        if lhs != toggler:
            sample_ok = False
            break
    return {
        "even_graph_dimension": m,
        "complement_toggler_edges": len(toggler),
        "complement_toggler_degree_sequence": tuple(sorted(degrees)),
        "translation_constant_on_sample": sample_ok,
    }


def c3_count(adj: list[list[int]]) -> int:
    n = len(adj)
    total = 0
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)
    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rev = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order: list[int] = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in graph[v]:
            if not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    seen = [False] * n
    sizes = []

    def rdfs(v: int) -> int:
        seen[v] = True
        size = 1
        for w in rev[v]:
            if not seen[w]:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if not seen[v]:
            sizes.append(rdfs(v))
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for used in range(1 << n):
        for v in range(n):
            cur = dp[used][v]
            if not cur:
                continue
            for w in range(n):
                if ((used >> w) & 1) == 0 and adj[v][w]:
                    dp[used | (1 << w)][w] += cur
    return sum(dp[-1])


def tournament_analysis() -> None:
    lanes = [
        ("tiling_line_lift", (5, 5, 5, 4, 4)),
        ("even_graph_cycle_space", (5, 4, 5, 5, 3)),
        ("line_multiplicity_address", (5, 5, 4, 4, 4)),
        ("simple_sc_type_edges", (1, 2, 2, 2, 2)),
        ("raw_node_cardinality", (0, 1, 1, 1, 1)),
        ("lrc_owner_filter_transfer", (4, 4, 4, 3, 5)),
    ]
    n = len(lanes)
    adj = [[0] * n for _ in range(n)]
    for i, (_, si) in enumerate(lanes):
        for j, (_, sj) in enumerate(lanes):
            if i == j:
                continue
            if si > sj:
                adj[i][j] = 1
    print("\nTournament Analysis over parity carriers")
    print("  Pairwise observable: lexicographic preservation of parity, multiplicity,")
    print("  quotient descent, cycle-space meaning, and LRC transfer value.")
    scores = [sum(row) for row in adj]
    order = sorted(range(n), key=lambda i: scores[i], reverse=True)
    print(f"  score_hist={dict(Counter(scores))}")
    print(f"  directed_3cycles={c3_count(adj)}")
    print(f"  scc_sizes={scc_sizes(adj)}")
    print(f"  hamiltonian_paths={hamiltonian_paths(adj)}")
    print("  Hamiltonian order:")
    for i in order:
        print(f"    {scores[i]}  {lanes[i][0]}")


def main() -> None:
    banner("Merged line parity audit, S675b")
    print("HYP-2251 / T748")
    print("Question: black-line portion even? blue-line portion odd?")
    print("Answer preview: yes for the complement-line lift; no for the simple")
    print("SC/NS-colored merged edge projection.")

    rows = []
    for n in range(3, 8):
        data = build_fixed_base_quotient(n)
        line = audit_line_lift(data)
        simple = audit_simple_sc_type_edges(data)
        even = audit_even_graph_carrier(n)
        rows.append((n, data, line, simple, even))

    banner("Exact table")
    header = (
        "n  m  cls  merged  E  blueL blackL  blueOdd  blackEven  "
        "simpleBlackEven simpleBlueOdd"
    )
    print(header)
    print("-" * len(header))
    for n, data, line, simple, _even in rows:
        cls_ok = "ok" if data["class_count"] == A000568[n] else "BAD"
        edge_ok = "ok" if simple["merged_edges"] == KNOWN_MERGED_EDGES[n] else "BAD"
        print(
            f"{n:1d} {data['m']:2d} {data['class_count']:4d}{cls_ok:>4} "
            f"{data['merged_count']:6d} {simple['merged_edges']:4d}{edge_ok:>4} "
            f"{line['blue_lines']:6d} {line['black_lines']:6d} "
            f"{str(line['blue_all_active_odd']):>8} {str(line['black_all_even']):>10} "
            f"{str(simple['black_all_even']):>15} {str(simple['active_blue_all_odd']):>13}"
        )

    banner("Witnesses against the simple projected version")
    for n, _data, _line, simple, _even in rows:
        if simple["black_all_even"] and simple["active_blue_all_odd"]:
            continue
        print(
            f"n={n}: simple SC/NS edge split has "
            f"black_odd_witnesses={simple['black_odd_witnesses']} and "
            f"blue_even_active_witnesses={simple['blue_even_active_witnesses']}"
        )

    banner("Line-lift parity details")
    for n, data, line, _simple, even in rows:
        print(
            f"n={n}: grid-symmetric tilings={line['grid_symmetric_tilings']}, "
            f"blue active nodes={line['blue_active_nodes']}, "
            f"black active nodes={line['black_active_nodes']}, "
            f"blue loops={line['blue_loop_edges']}, black loops={line['black_loop_edges']}"
        )
        print(
            f"     blue parity={dict(line['blue_degree_parity'])}, "
            f"black parity={dict(line['black_degree_parity'])}, "
            f"SC merged nodes={sum(data['merged_sc'])}, "
            f"NS merged nodes={data['merged_count'] - sum(data['merged_sc'])}"
        )
        print(
            f"     even-graph carrier dim={even['even_graph_dimension']}, "
            f"complement toggler edges={even['complement_toggler_edges']}, "
            f"degree seq={even['complement_toggler_degree_sequence']}, "
            f"constant translation sample={even['translation_constant_on_sample']}"
        )

    banner("Proof sketch")
    print("Let Q_m be the fixed-base-path tiling cube, and let C(x)=x+1 be")
    print("all-free-tile complement.  The blue/black line graph has one line")
    print("for each pair {x,C(x)}.")
    print()
    print("Blue side: grid-reflection-fixed tilings occur only over SC merged")
    print("nodes.  In an SC class the non-fixed tilings pair under reflection,")
    print("while the total fixed-base tiling count is odd, so the blue endpoint")
    print("count in each active node is odd.  With loops counted twice, the")
    print("blue multigraph has odd degree at every active node.")
    print()
    print("Black side: in an SC node, black = odd total - odd blue, hence even.")
    print("In an NS merged pair, black = odd class A + odd class B, also even.")
    print("Thus the black line graph is Eulerian: a genuine even graph.")
    print()
    print("The simple merged metagraph forgets the line multiplicity/address.")
    print("That projection is exactly where the parity theorem fails.")

    banner("Assumption challenge")
    print("Considered vertex sets: merged tournament nodes, tiling-complement")
    print("pairs, single-flip edges, even-graph cycle coordinates, and proof")
    print("obligation addresses.  The quotient predicate preserved by the line")
    print("lift is degree parity of complement-pair endpoints.  The simple")
    print("merged graph preserves adjacency but destroys multiplicity and loop")
    print("half-edge data.  Challenged assumption: metagraph edges themselves")
    print("need not be the right vertices for parity; the lines are.")

    tournament_analysis()

    banner("Synthesis")
    print("Even numbers/graphs are kernel objects: boundary zero.  Odd line")
    print("graphs are affine boundary objects: boundary one on the SC support.")
    print("So the useful identity is not just 'tournaments are equinumerous with")
    print("even graphs'; it is 'the black complement-line lift is already an")
    print("even graph, while the blue lift is the odd coset whose boundary tells")
    print("which SC nodes kept a grid-symmetric witness.'")


if __name__ == "__main__":
    main()
