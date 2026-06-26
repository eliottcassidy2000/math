#!/usr/bin/env python3
"""S211: A000568, k-depth perspectives, and controlled-forgetting carriers.

The prompt asks where the old count coincidence stops:

    A000568(n) = total node perspectives on all (n-1)-tournament classes.

Writing U(n) for unlabelled tournaments and P(m) for rooted/node-perspective
classes on m vertices, this asks when U(m+1) = P(m) fails.  The answer is
m=5, i.e. shifted n=6:

    P(5)=48 but U(6)=56.

This script treats that failure as the first controlled-forgetting alarm.  It
computes small exact counts and several finite "perspective carrier" lifts:
node k-depth refinement, directed-edge perspectives, directed-cycle
perspectives, transitive-clique insertion perspectives, and a qualitative
Tournament Analysis over proof carriers.

Methodological guardrail: because tournaments are complete, graph distance is
not the useful notion of depth.  The k-depth node signature is an iterative
role-coloring: a vertex remembers the multiset of depth-(k-1) perspectives on
its out-neighbors and in-neighbors.  This is a controlled-forgetting ladder,
not a replacement for exact rooted isomorphism.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations

from tournament_rigidity_cascade_s589 import (
    canonical,
    classes,
    directed_three_cycles,
    edge,
    ham_paths_adj,
    relabel,
    rooted_canonical,
    scc_sizes,
    vertex_orbits,
)


KNOWN_U = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
}


def u_count(n: int) -> int:
    """A000568/unlabelled tournament count."""
    if n in KNOWN_U:
        return KNOWN_U[n]
    return len(classes(n))


def rooted_count(n: int) -> int:
    """Number of exact node perspectives on all n-tournament classes."""
    return sum(len(vertex_orbits(rep, n)) for rep in classes(n))


def source_root_count(n: int) -> int:
    """Number of source-rooted n-tournament perspectives."""
    total = 0
    for rep in classes(n):
        for orbit in vertex_orbits(rep, n):
            root = orbit[0]
            if all(edge(rep, n, root, v) for v in range(n) if v != root):
                total += 1
    return total


def node_depth_colors(n: int, max_depth: int) -> list[dict[tuple[int, int], int]]:
    """Global k-depth node colors over all n-tournament class representatives."""
    reps = classes(n)
    colors: list[dict[tuple[int, int], int]] = []
    raw0: dict[tuple[int, int], tuple[int]] = {}
    for rep in reps:
        for v in range(n):
            outdeg = sum(edge(rep, n, v, u) for u in range(n) if u != v)
            raw0[(rep, v)] = (outdeg,)
    colors.append(compress(raw0))

    for depth in range(1, max_depth + 1):
        prev = colors[-1]
        raw: dict[tuple[int, int], tuple[object, ...]] = {}
        for rep in reps:
            for v in range(n):
                out_colors = sorted(prev[(rep, u)] for u in range(n) if u != v and edge(rep, n, v, u))
                in_colors = sorted(prev[(rep, u)] for u in range(n) if u != v and edge(rep, n, u, v))
                raw[(rep, v)] = (prev[(rep, v)], tuple(out_colors), tuple(in_colors))
        colors.append(compress(raw))
        if len(set(colors[-1].values())) == len(set(colors[-2].values())):
            # Continue to max_depth anyway so the printed table is rectangular.
            pass
    return colors


def compress(raw: dict[tuple[int, int], tuple[object, ...]]) -> dict[tuple[int, int], int]:
    mapping = {value: i for i, value in enumerate(sorted(set(raw.values())))}
    return {key: mapping[value] for key, value in raw.items()}


def exact_rooted_labels(n: int) -> dict[tuple[int, int], int]:
    raw = {(rep, v): rooted_canonical(rep, n, v) for rep in classes(n) for v in range(n)}
    mapping = {value: i for i, value in enumerate(sorted(set(raw.values())))}
    return {key: mapping[value] for key, value in raw.items()}


def edge_rooted_canonical(mask: int, n: int, tail: int, tip: int) -> int:
    """Canonical tournament with ordered directed edge tail -> tip fixed first."""
    others = tuple(v for v in range(n) if v not in (tail, tip))
    return min(relabel(mask, n, (tail, tip) + p) for p in permutations(others))


def oriented_edges(mask: int, n: int) -> list[tuple[int, int]]:
    out = []
    for a, b in combinations(range(n), 2):
        out.append((a, b) if edge(mask, n, a, b) else (b, a))
    return out


def exact_edge_perspective_count(n: int) -> int:
    values = set()
    for rep in classes(n):
        for tail, tip in oriented_edges(rep, n):
            values.add(edge_rooted_canonical(rep, n, tail, tip))
    return len(values)


def edge_depth_signature(mask: int, n: int, tail: int, tip: int, node_color: dict[tuple[int, int], int]) -> tuple[object, ...]:
    groups: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for x in range(n):
        if x in (tail, tip):
            continue
        sector = (int(edge(mask, n, tail, x)), int(edge(mask, n, tip, x)))
        groups[(sector[0], sector[1], node_color[(mask, x)])].append(x)

    group_keys = tuple(sorted((key, len(value)) for key, value in groups.items()))
    cross = []
    sorted_keys = sorted(groups)
    for i, ka in enumerate(sorted_keys):
        for kb in sorted_keys[i + 1 :]:
            wins = sum(1 for a in groups[ka] for b in groups[kb] if edge(mask, n, a, b))
            cross.append((ka, kb, wins, len(groups[ka]) * len(groups[kb])))
    return (group_keys, tuple(cross))


def edge_depth_counts(n: int, node_colors: list[dict[tuple[int, int], int]]) -> list[int]:
    counts = []
    for color in node_colors:
        sigs = set()
        for rep in classes(n):
            for tail, tip in oriented_edges(rep, n):
                sigs.add(edge_depth_signature(rep, n, tail, tip, color))
        counts.append(len(sigs))
    return counts


def is_directed_cycle(mask: int, n: int, cyc: tuple[int, int, int]) -> bool:
    a, b, c = cyc
    return edge(mask, n, a, b) and edge(mask, n, b, c) and edge(mask, n, c, a)


def cycle_rooted_canonical(mask: int, n: int, cyc: tuple[int, int, int]) -> int:
    orders = (cyc, cyc[1:] + cyc[:1], cyc[2:] + cyc[:2])
    best = None
    for order in orders:
        others = tuple(v for v in range(n) if v not in order)
        for p in permutations(others):
            value = relabel(mask, n, order + p)
            best = value if best is None else min(best, value)
    assert best is not None
    return best


def exact_cycle_perspective_count(n: int) -> int:
    values = set()
    for rep in classes(n):
        for a, b, c in permutations(range(n), 3):
            cyc = (a, b, c)
            if is_directed_cycle(rep, n, cyc):
                values.add(cycle_rooted_canonical(rep, n, cyc))
    return len(values)


def transitive_chain(triple: tuple[int, int, int], mask: int, n: int) -> tuple[int, int, int] | None:
    scores = []
    for v in triple:
        scores.append((sum(edge(mask, n, v, u) for u in triple if u != v), v))
    if sorted(score for score, _ in scores) != [0, 1, 2]:
        return None
    ordered = tuple(v for _, v in sorted(scores, reverse=True))
    return ordered


def chain_rooted_canonical(mask: int, n: int, chain: tuple[int, int, int]) -> int:
    others = tuple(v for v in range(n) if v not in chain)
    return min(relabel(mask, n, chain + p) for p in permutations(others))


def exact_transitive_clique_perspective_count(n: int) -> int:
    values = set()
    for rep in classes(n):
        for triple in combinations(range(n), 3):
            chain = transitive_chain(triple, rep, n)
            if chain is not None:
                values.add(chain_rooted_canonical(rep, n, chain))
    return len(values)


def carrier_tournament() -> tuple[list[tuple[str, dict[str, int]]], list[list[int]]]:
    """Qualitative Tournament Analysis over proposed LRC proof carriers."""
    carriers = [
        (
            "raw_A000568_class",
            dict(source=0, incident=0, pair=0, cycle=0, insertion=0, owner=0, automaton=0, cost=0),
        ),
        (
            "exact_rooted_node",
            dict(source=1, incident=0, pair=0, cycle=0, insertion=0, owner=0, automaton=0, cost=1),
        ),
        (
            "k_depth_node_ladder",
            dict(source=1, incident=0, pair=0, cycle=0, insertion=0, owner=0, automaton=1, cost=2),
        ),
        (
            "directed_edge_perspective",
            dict(source=1, incident=1, pair=1, cycle=0, insertion=0, owner=0, automaton=1, cost=3),
        ),
        (
            "directed_cycle_conflict",
            dict(source=0, incident=1, pair=1, cycle=1, insertion=0, owner=0, automaton=1, cost=3),
        ),
        (
            "transitive_clique_insertion",
            dict(source=1, incident=1, pair=1, cycle=0, insertion=1, owner=0, automaton=1, cost=3),
        ),
        (
            "edge_cycle_incidence_conflict",
            dict(source=1, incident=1, pair=1, cycle=1, insertion=1, owner=0, automaton=1, cost=4),
        ),
        (
            "endpoint_owner_packet_sheaf",
            dict(source=1, incident=1, pair=1, cycle=1, insertion=1, owner=1, automaton=1, cost=5),
        ),
    ]

    def score(data: dict[str, int]) -> tuple[int, int]:
        retained = (
            3 * data["source"]
            + 3 * data["incident"]
            + 2 * data["pair"]
            + 2 * data["cycle"]
            + 2 * data["insertion"]
            + 4 * data["owner"]
            + data["automaton"]
        )
        return retained, -data["cost"]

    n = len(carriers)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        si = score(carriers[i][1])
        sj = score(carriers[j][1])
        if si > sj or (si == sj and carriers[i][0] < carriers[j][0]):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return carriers, adj


def print_count_tables(max_n: int, max_depth: int) -> None:
    print("1. SHIFTED A000568 / NODE-PERSPECTIVE COUNT")
    print("   m is the tournament size; shifted_n=m+1 is the user's A000568 index.")
    print("   The first failure is shifted_n=6, equivalently m=5.")
    print()
    print("   m  shifted_n  U(m)  P(m)  U(m+1)  U(m+1)-P(m)  source_roots  U(m-1)")
    for m in range(1, max_n + 1):
        u_prev = u_count(m - 1) if m > 1 else 1
        print(
            f"   {m:1d}     {m+1:2d}      {u_count(m):4d}  {rooted_count(m):4d}  "
            f"{u_count(m+1):7d}  {u_count(m+1)-rooted_count(m):13d}  "
            f"{source_root_count(m):12d}  {u_prev:6d}"
        )

    print()
    print("2. k-DEPTH NODE LADDER")
    print("   depth 0 is score.  depth k+1 sees in/out multisets of depth-k node types.")
    for m in range(3, max_n + 1):
        colors = node_depth_colors(m, max_depth)
        exact = len(set(exact_rooted_labels(m).values()))
        counts = [len(set(c.values())) for c in colors]
        print(
            f"   m={m}: depth_counts={counts} exact_rooted={exact} "
            f"U(m+1)-exact={u_count(m+1)-exact}"
        )

    print()
    print("3. EDGE / CYCLE / CLIQUE PERSPECTIVE LIFTS")
    print("   Edge depth uses sector colors around the directed edge tail -> tip.")
    for m in range(3, max_n + 1):
        colors = node_depth_colors(m, max_depth)
        edge_counts = edge_depth_counts(m, colors)
        exact_edge = exact_edge_perspective_count(m)
        exact_cycle = exact_cycle_perspective_count(m)
        exact_chain = exact_transitive_clique_perspective_count(m)
        print(
            f"   m={m}: edge_depth_counts={edge_counts} exact_edge={exact_edge} "
            f"directed_cycle={exact_cycle} transitive_chain={exact_chain} "
            f"target_U(m+1)={u_count(m+1)}"
        )


def print_carrier_index() -> None:
    carriers, adj = carrier_tournament()
    scores = [sum(row) for row in adj]
    print()
    print("4. CREATIVE PROOF-CARRIER INDEX")
    lines = {
        "k_depth_node_ladder": "Controlled local memory.  Good for diagnosing the break; insufficient by itself because exact rooted nodes still cap at P(5)=48.",
        "directed_edge_perspective": "Tip/tail dual carrier.  Other vertices sit in four sectors (both beaten, tail-only, tip-only, both beat endpoints); cross-sector orientation is the missing coupling payload.",
        "directed_cycle_conflict": "Frustration carrier.  A directed 3-cycle remembers chirality and catches the kind of orientation debt that the eight-class gap suggests.",
        "transitive_clique_insertion": "Insertion-cut carrier.  A chain/rooted clique classifies external vertices by cut position, matching the add-a-new-observer problem.",
        "edge_cycle_incidence_conflict": "Bipartite conflict carrier between directed edges and cyclic witnesses; useful when local edge sectors and cycle chirality must agree.",
        "endpoint_owner_packet_sheaf": "LRC-facing lift: source predicate plus incident threshold word plus endpoint-owner/gap-pressure labels.",
        "raw_A000568_class": "Useful cache only.  It forgets observer and incident word.",
        "exact_rooted_node": "The old perspective object.  It preserves observer choice but forgets the incident word/cocycle needed at the first failure.",
    }
    for name, _ in carriers:
        print(f"   - {name}: {lines[name]}")

    print()
    print("5. TOURNAMENT ANALYSIS OVER CARRIERS")
    print("   Vertices are proof carriers, not runners or arcs.")
    print("   Observable: retained source/incident/pair/cycle/insertion/owner/automaton payload minus proof cost.")
    print("   Switch: orient toward the carrier retaining more missing observer-coupling payload; name order breaks exact ties.")
    print(f"   vertices={[name for name, _ in carriers]}")
    print(f"   score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"   directed_3_cycles={directed_three_cycles(adj)}")
    print(f"   scc_sizes={scc_sizes(adj)}")
    print(f"   hamiltonian_paths={ham_paths_adj(adj)}")
    order = sorted(range(len(carriers)), key=lambda i: -scores[i])
    print("   one_hamiltonian_path=" + " -> ".join(carriers[i][0] for i in order))


def print_reading() -> None:
    print()
    print("READING")
    print(
        "  The old equality is a shallow coincidence between observer-coupled m-objects "
        "and observer-blind (m+1)-objects.  It fails at shifted n=6 because a rooted "
        "five-tournament has forgotten an incident-word/cross-coupling coordinate."
    )
    print(
        "  The source slice remains exact: source-rooted m-tournaments delete to "
        "U(m-1), so source deletion is the controlled-forgetting map that LRC can "
        "actually use.  Non-source perspectives need sidecars."
    )
    print(
        "  Edge perspectives are the first natural lift because an observer added to "
        "a rooted node is really an incident word.  Directed edges expose the four "
        "tip/tail sectors and their cross-sector orientations.  Cycle and clique "
        "perspectives then split chirality and insertion-order debt."
    )
    print(
        "  LRC translation: do not use raw A000568 or exact rooted nodes as final "
        "proof states.  Use them as cache layers below source threshold arcs, "
        "endpoint-owner packets, gap-pressure fibers, and proof-obligation automata."
    )
    print()
    print("ASSUMPTION CHALLENGE")
    print(
        "  Considered vertices: tournament nodes, rooted perspectives, directed edges, "
        "directed cycles, transitive clique insertion cuts, edge-cycle conflicts, "
        "endpoint owners, gap-pressure fibers, automaton states, and proof obligations."
    )
    print(
        "  Preserved predicate in the selected quotient: whether a carrier keeps the "
        "observer-source/incident-coupling data needed for an LRC safe-box hit."
    )
    print(
        "  Destroyed information: exact labelled runner identities, full extension "
        "rows, and scalar ordering.  These losses are acceptable only when sidecars "
        "retain/reconstruct/annihilate the forgotten coordinate or emit named debt."
    )
    print(
        "  Challenged assumption: tournament vertices should be only runners or raw "
        "arcs.  For the LRC proof they should often be proof carriers, observer "
        "sidecars, wall-crossing events, endpoint-owner packets, or conflict cells."
    )


def main() -> None:
    print("=" * 80)
    print("S211: A000568 k-depth perspective ladder and LRC carriers")
    print("=" * 80)
    print_count_tables(max_n=5, max_depth=4)
    print_carrier_index()
    print_reading()


if __name__ == "__main__":
    main()
