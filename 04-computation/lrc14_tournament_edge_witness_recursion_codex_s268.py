#!/usr/bin/env python3
"""HYP-3124 tournament edge-witness recursion scout.

This is an executable synthesis for the prompt "tournament edges as witness,
think recursively on both tip and tail."  It is not an LRC14 proof.

Tournament Analysis declaration:
  vertices: edge-witness proof reframes, not runners, roots, or scalar edges;
  pairwise observable: retained proof payload across named axes;
  switch/gauge: majority of axis scores, with a fixed tie Hamiltonian path;
  tie path: coordinate repair and two-ended witness packets before local
            sector decks, one-sided recursions, and raw edge counts.

The exact census layer uses actual oriented edges.  For every directed edge
tail -> tip in every labelled tournament through n=5, it compares:
  1. the outside four-sector word around the edge;
  2. one-sided endpoint deletion children;
  3. the paired tail-deletion and tip-deletion children;
  4. one recursive layer inside both endpoint-deletion children.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations
from typing import Dict, Iterable, List, Sequence, Tuple


AXES = (
    "lrc_predicate_retention",
    "tail_tip_symmetry",
    "recursive_closure",
    "sector_resolution",
    "observer_gluing",
    "missing_input_repair",
    "domain_wall_alignment",
    "formal_readiness",
    "scalar_guardrail",
)


@dataclass(frozen=True)
class Reframe:
    name: str
    axes: Dict[str, int]
    preserves: str
    destroys: str
    next_use: str
    locality_rank: int


def pair_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return i * (2 * n - i - 1) // 2 + (j - i - 1)


def edge(mask: int, n: int, i: int, j: int) -> bool:
    """Return True exactly when i -> j."""
    if i == j:
        raise ValueError("loops are not tournament edges")
    if i < j:
        return bool((mask >> pair_index(n, i, j)) & 1)
    return not edge(mask, n, j, i)


def directed_edges(mask: int, n: int) -> Iterable[Tuple[int, int]]:
    for i, j in combinations(range(n), 2):
        if edge(mask, n, i, j):
            yield i, j
        else:
            yield j, i


def transform(mask: int, n: int, perm: Tuple[int, ...]) -> int:
    """Relabel old vertex perm[new_vertex] as new_vertex."""
    out = 0
    for i, j in combinations(range(n), 2):
        if edge(mask, n, perm[i], perm[j]):
            out |= 1 << pair_index(n, i, j)
    return out


@lru_cache(maxsize=None)
def canonical(mask: int, n: int) -> int:
    return min(transform(mask, n, perm) for perm in permutations(range(n)))


def induced(mask: int, n: int, keep: Tuple[int, ...]) -> int:
    out = 0
    m = len(keep)
    for i, j in combinations(range(m), 2):
        if edge(mask, n, keep[i], keep[j]):
            out |= 1 << pair_index(m, i, j)
    return out


def delete_vertex(mask: int, n: int, v: int) -> Tuple[int, int]:
    keep = tuple(i for i in range(n) if i != v)
    return induced(mask, n, keep), n - 1


def outdegree(mask: int, n: int, v: int) -> int:
    return sum(1 for w in range(n) if w != v and edge(mask, n, v, w))


@lru_cache(maxsize=None)
def score_sequence(mask: int, n: int) -> Tuple[int, ...]:
    cmask = canonical(mask, n)
    return tuple(sorted(outdegree(cmask, n, v) for v in range(n)))


@lru_cache(maxsize=None)
def triad_profile(mask: int, n: int) -> Tuple[int, int]:
    cmask = canonical(mask, n)
    transitive = 0
    cyclic = 0
    for a, b, c in combinations(range(n), 3):
        local = [
            sum(1 for y in (a, b, c) if y != x and edge(cmask, n, x, y))
            for x in (a, b, c)
        ]
        if sorted(local) == [0, 1, 2]:
            transitive += 1
        else:
            cyclic += 1
    return transitive, cyclic


@lru_cache(maxsize=None)
def tournament_sig(mask: int, n: int) -> Tuple[object, ...]:
    cmask = canonical(mask, n)
    return (n, score_sequence(cmask, n), triad_profile(cmask, n), cmask)


def sector_word(mask: int, n: int, tail: int, tip: int) -> Tuple[int, int, int, int]:
    """Four outside sectors around a directed edge tail -> tip.

    Order:
      0. tail -> w and tip -> w
      1. tail -> w and w -> tip
      2. w -> tail and tip -> w
      3. w -> tail and w -> tip
    """
    sectors = [0, 0, 0, 0]
    for w in range(n):
        if w in (tail, tip):
            continue
        tail_to_w = edge(mask, n, tail, w)
        tip_to_w = edge(mask, n, tip, w)
        if tail_to_w and tip_to_w:
            sectors[0] += 1
        elif tail_to_w and not tip_to_w:
            sectors[1] += 1
        elif not tail_to_w and tip_to_w:
            sectors[2] += 1
        else:
            sectors[3] += 1
    return tuple(sectors)


def one_sided_child(mask: int, n: int, deleted: int) -> Tuple[object, ...]:
    child, child_n = delete_vertex(mask, n, deleted)
    return tournament_sig(child, child_n)


@lru_cache(maxsize=None)
def edge_signature(mask: int, n: int, tail: int, tip: int, depth: int) -> Tuple[object, ...]:
    if not edge(mask, n, tail, tip):
        raise ValueError("edge_signature expects tail -> tip")
    base = (
        n,
        sector_word(mask, n, tail, tip),
        outdegree(mask, n, tail),
        outdegree(mask, n, tip),
        one_sided_child(mask, n, tail),
        one_sided_child(mask, n, tip),
    )
    if depth <= 0 or n <= 3:
        return base

    tail_deleted, m = delete_vertex(mask, n, tail)
    tip_deleted, _ = delete_vertex(mask, n, tip)
    tail_child_deck = tuple(
        sorted(edge_signature(tail_deleted, m, a, b, depth - 1) for a, b in directed_edges(tail_deleted, m))
    )
    tip_child_deck = tuple(
        sorted(edge_signature(tip_deleted, m, a, b, depth - 1) for a, b in directed_edges(tip_deleted, m))
    )
    return base + (tail_child_deck, tip_child_deck)


def edge_census(n: int, depth: int = 1) -> Dict[str, object]:
    total_tournaments = 1 << (n * (n - 1) // 2)
    sector_counter: Counter[Tuple[int, int, int, int]] = Counter()
    tail_child_counter: Counter[Tuple[object, ...]] = Counter()
    tip_child_counter: Counter[Tuple[object, ...]] = Counter()
    paired_child_counter: Counter[Tuple[object, ...]] = Counter()
    depth_counter: Counter[Tuple[object, ...]] = Counter()
    sector_to_tail_child: Dict[Tuple[int, int, int, int], set] = defaultdict(set)
    sector_to_tip_child: Dict[Tuple[int, int, int, int], set] = defaultdict(set)
    sector_to_pair: Dict[Tuple[int, int, int, int], set] = defaultdict(set)
    pair_to_depth: Dict[Tuple[object, ...], set] = defaultdict(set)

    for mask in range(total_tournaments):
        for tail, tip in directed_edges(mask, n):
            sector = sector_word(mask, n, tail, tip)
            tail_child = one_sided_child(mask, n, tail)
            tip_child = one_sided_child(mask, n, tip)
            child_pair = (sector, tail_child, tip_child)
            depth_sig = edge_signature(mask, n, tail, tip, depth)

            sector_counter[sector] += 1
            tail_child_counter[(sector, tail_child)] += 1
            tip_child_counter[(sector, tip_child)] += 1
            paired_child_counter[child_pair] += 1
            depth_counter[depth_sig] += 1
            sector_to_tail_child[sector].add(tail_child)
            sector_to_tip_child[sector].add(tip_child)
            sector_to_pair[sector].add(child_pair)
            pair_to_depth[child_pair].add(depth_sig)

    top_sector_splits = sorted(
        (
            len(sector_to_pair[sector]),
            len(sector_to_tail_child[sector]),
            len(sector_to_tip_child[sector]),
            sector_counter[sector],
            sector,
        )
        for sector in sector_counter
    )[-5:]
    top_sector_splits.reverse()

    return {
        "n": n,
        "total_tournaments": total_tournaments,
        "directed_edge_instances": total_tournaments * n * (n - 1) // 2,
        "unique_sector_words": len(sector_counter),
        "unique_sector_plus_tail_child": len(tail_child_counter),
        "unique_sector_plus_tip_child": len(tip_child_counter),
        "unique_sector_plus_child_pair": len(paired_child_counter),
        "unique_depth1_signatures": len(depth_counter),
        "sector_groups_split_by_tail_child": sum(1 for vals in sector_to_tail_child.values() if len(vals) > 1),
        "sector_groups_split_by_tip_child": sum(1 for vals in sector_to_tip_child.values() if len(vals) > 1),
        "sector_groups_split_by_child_pair": sum(1 for vals in sector_to_pair.values() if len(vals) > 1),
        "child_pair_groups_split_by_depth1": sum(1 for vals in pair_to_depth.values() if len(vals) > 1),
        "max_child_pairs_inside_one_sector": max(len(vals) for vals in sector_to_pair.values()),
        "max_depth1_sigs_inside_one_child_pair": max(len(vals) for vals in pair_to_depth.values()),
        "top_sector_splits": top_sector_splits,
    }


def reframe(name: str, scores: Sequence[int], preserves: str, destroys: str, next_use: str, locality_rank: int) -> Reframe:
    if len(scores) != len(AXES):
        raise ValueError(name)
    return Reframe(name, dict(zip(AXES, scores)), preserves, destroys, next_use, locality_rank)


REFRAMES: List[Reframe] = [
    reframe(
        "coordinate_resurrection_edge_sheaf",
        (5, 5, 5, 4, 5, 5, 4, 5, 5),
        "destroyed-coordinate vector, repair cover, adjoint section, and live endpoint section",
        "nothing essential if every lost edge coordinate is routed to repair or named debt",
        "make HYP-3118 repair covers accept edge-tail/tip coordinates as live sections",
        8,
    ),
    reframe(
        "edge_witness_two_ended_packet",
        (5, 5, 5, 5, 4, 4, 4, 4, 5),
        "tail packet, tip packet, four-sector deck, both endpoint-deletion children, repair/debt",
        "raw labels and scalar edge counts",
        "promote `edge_witness(tail->tip)` as a packet row field in HYP-2963/HYP-3098/HYP-3112",
        7,
    ),
    reframe(
        "paired_tail_tip_deletion_recursion",
        (4, 5, 5, 3, 3, 3, 3, 3, 5),
        "two recursive child obligations after deleting the tail and after deleting the tip",
        "outside cross-sector orientation if the child pair is stored alone",
        "use as the recursive kernel beneath the full edge witness packet",
        6,
    ),
    reframe(
        "proof_circuit_edge_gate",
        (4, 4, 3, 3, 3, 4, 3, 5, 5),
        "edge witness as a proof-circuit gate with missing-input vector and terminal exit",
        "geometric sector detail unless included as gate input",
        "add edge-witness inputs to HYP-3116/HYP-3117 circuit-certificate vectors",
        9,
    ),
    reframe(
        "cross_sector_orientation_word",
        (4, 3, 2, 5, 5, 3, 3, 4, 5),
        "HYP-3050/HYP-3054/HYP-3106 cross-sector orientation around tail -> tip",
        "recursive endpoint child obligations",
        "retain as the local observer-gluing coordinate inside the two-ended packet",
        3,
    ),
    reframe(
        "domain_wall_edge_classifier",
        (4, 3, 3, 4, 3, 3, 5, 3, 5),
        "HYP-3115 one-swap / Ising domain-wall edge labels",
        "tail-tip recursion if a wall is treated as one scalar cut",
        "classify domain walls by both endpoint-deletion children before discharging them",
        4,
    ),
    reframe(
        "ear_payload_edge_mass",
        (4, 3, 3, 3, 3, 2, 4, 3, 5),
        "HYP-3112 one-runner ear payload attached to an edge boundary",
        "missing finite-address, observer-gluing, and repair coordinates",
        "feed edge-witness sectors into root-motion/ear-payload ledgers",
        5,
    ),
    reframe(
        "outside_four_sector_deck",
        (3, 4, 1, 4, 4, 2, 2, 3, 5),
        "the four outside sectors around tail -> tip",
        "both recursive children and most repair obligations",
        "use only as the local observable, not as the whole witness",
        2,
    ),
    reframe(
        "tail_deleted_one_sided_recursion",
        (3, 2, 3, 2, 2, 2, 2, 2, 4),
        "the child after deleting the tail endpoint",
        "tip-side reciprocity and the paired compatibility condition",
        "negative control for one-ended recursion",
        1,
    ),
    reframe(
        "tip_deleted_one_sided_recursion",
        (3, 2, 3, 2, 2, 2, 2, 2, 4),
        "the child after deleting the tip endpoint",
        "tail-side reciprocity and the paired compatibility condition",
        "negative control for one-ended recursion",
        1,
    ),
    reframe(
        "raw_edge_count_scalar",
        (1, 1, 0, 0, 0, 0, 1, 1, 0),
        "only the number of oriented edges or boundary edges",
        "tail/tip roles, sector word, recursive children, packet route, and repair debt",
        "keep as a scalar warning, never as a proof carrier",
        0,
    ),
]


TIE_PATH = [r.name for r in REFRAMES]


def compare(a: Reframe, b: Reframe) -> str:
    a_votes = 0
    b_votes = 0
    for axis in AXES:
        if a.axes[axis] > b.axes[axis]:
            a_votes += 1
        elif b.axes[axis] > a.axes[axis]:
            b_votes += 1
    if a_votes > b_votes:
        return a.name
    if b_votes > a_votes:
        return b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_edges(vertices: Sequence[Reframe]) -> Dict[str, set]:
    adj: Dict[str, set] = {v.name: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        winner = compare(a, b)
        loser = b.name if winner == a.name else a.name
        adj[winner].add(loser)
    return adj


def score_histogram(adj: Dict[str, set]) -> Dict[int, int]:
    return dict(sorted(Counter(len(v) for v in adj.values()).items()))


def directed_3cycles(adj: Dict[str, set]) -> int:
    count = 0
    names = sorted(adj)
    for a, b, c in combinations(names, 3):
        ab = b in adj[a]
        bc = c in adj[b]
        ca = a in adj[c]
        ba = a in adj[b]
        cb = b in adj[c]
        ac = c in adj[a]
        if (ab and bc and ca) or (ba and cb and ac):
            count += 1
    return count


def scc_sizes(adj: Dict[str, set]) -> List[int]:
    index = 0
    stack: List[str] = []
    on_stack = set()
    indices: Dict[str, int] = {}
    low: Dict[str, int] = {}
    sizes: List[int] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in adj[v]:
            if w not in indices:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            size = 0
            while True:
                w = stack.pop()
                on_stack.remove(w)
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for v in adj:
        if v not in indices:
            strongconnect(v)
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: Dict[str, set]) -> int:
    names = list(adj)
    idx = {name: i for i, name in enumerate(names)}
    n = len(names)
    dp: Dict[Tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            paths = dp.get((mask, last), 0)
            if not paths:
                continue
            for nxt_name in adj[names[last]]:
                nxt = idx[nxt_name]
                if mask & (1 << nxt):
                    continue
                dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + paths
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def one_hamiltonian_path(adj: Dict[str, set]) -> List[str]:
    names = list(adj)

    def search(path: List[str], unused: set) -> List[str] | None:
        if not unused:
            return path
        last = path[-1]
        for nxt in sorted(unused, key=lambda name: -len(adj[name])):
            if nxt in adj[last]:
                found = search(path + [nxt], unused - {nxt})
                if found:
                    return found
        return None

    for start in sorted(names, key=lambda name: -len(adj[name])):
        found = search([start], set(names) - {start})
        if found:
            return found
    return []


def locality_gauge_flips(adj: Dict[str, set], vertices: Sequence[Reframe]) -> int:
    by_name = {v.name: v for v in vertices}
    flips = 0
    for a, b in combinations(sorted(by_name), 2):
        if by_name[a].locality_rank == by_name[b].locality_rank:
            continue
        locality_winner = a if by_name[a].locality_rank > by_name[b].locality_rank else b
        if (locality_winner == a and b in adj[a]) or (locality_winner == b and a in adj[b]):
            continue
        flips += 1
    return flips


def tournament_report() -> List[str]:
    adj = tournament_edges(REFRAMES)
    path = one_hamiltonian_path(adj)
    lines = [
        "3. TOURNAMENT ANALYSIS OVER EDGE-WITNESS REFRAMES",
        "vertices=edge-witness proof reframes, not runners, roots, spins, or scalar edges",
        "pairwise_observable=majority of retained proof payload over axes "
        + ",".join(AXES),
        "switch=orient A->B when A wins more axes; ties follow the declared Hamiltonian path",
        "tie_hamiltonian_path=" + " -> ".join(TIE_PATH),
        f"score_hist={score_histogram(adj)}",
        f"directed_3cycles={directed_3cycles(adj)}",
        f"scc_sizes={scc_sizes(adj)}",
        f"hamiltonian_path_count={hamiltonian_path_count(adj)}",
        f"edge_flips_against_locality_first_gauge={locality_gauge_flips(adj, REFRAMES)}",
        "selected_hamiltonian_path=" + " -> ".join(path),
        "",
        "Top reframe fields:",
    ]
    for name in path[:5]:
        ref = next(r for r in REFRAMES if r.name == name)
        lines.append(f"- {ref.name}: preserves={ref.preserves}; next={ref.next_use}")
    return lines


def census_report() -> List[str]:
    lines = [
        "2. EXACT SMALL-TOURNAMENT EDGE CENSUS",
        "sector order=(tail->w & tip->w, tail->w & w->tip, w->tail & tip->w, w->tail & w->tip)",
        (
            "n  labelled_T  directed_edges  sectors  sector+tail_child  "
            "sector+tip_child  sector+child_pair  depth1_sigs  "
            "sector_groups_split_by_child_pair"
        ),
    ]
    summaries = []
    for n in range(2, 6):
        summary = edge_census(n)
        summaries.append(summary)
        lines.append(
            f"{n}  {summary['total_tournaments']:10d}  "
            f"{summary['directed_edge_instances']:14d}  "
            f"{summary['unique_sector_words']:7d}  "
            f"{summary['unique_sector_plus_tail_child']:17d}  "
            f"{summary['unique_sector_plus_tip_child']:16d}  "
            f"{summary['unique_sector_plus_child_pair']:17d}  "
            f"{summary['unique_depth1_signatures']:11d}  "
            f"{summary['sector_groups_split_by_child_pair']:34d}"
        )
    lines.append("")
    lines.append("n=5 largest sector refinements (child_pairs, tail_children, tip_children, instances, sector):")
    for child_pairs, tail_children, tip_children, instances, sector in summaries[-1]["top_sector_splits"]:
        lines.append(
            f"- {sector}: child_pairs={child_pairs}, tail_children={tail_children}, "
            f"tip_children={tip_children}, instances={instances}"
        )
    lines.append("")
    lines.append(
        "Readout: the four-sector word is a useful local observable, but every "
        "nontrivial n=5 sector group splits after both endpoint-deletion "
        "children are kept.  One-sided recursion sees only one projection; the "
        "paired child object is the first natural witness carrier."
    )
    return lines


def packet_schema_report() -> List[str]:
    return [
        "1. EDGE-WITNESS PACKET SCHEMA",
        "namespace=HYP-3124/T1198/LTI-259/LTT-157; provisional S268 HYP-3119 was renumbered after rebase",
        "edge_witness(tail -> tip) = (",
        "  endpoint_role_word,",
        "  outside_four_sector_deck,",
        "  tail_deletion_child_signature,",
        "  tip_deletion_child_signature,",
        "  recursive_tail_child_edge_deck,",
        "  recursive_tip_child_edge_deck,",
        "  observer_gluing_payload_orbit,",
        "  missing_input_vector,",
        "  coordinate_resurrection_sidecar_or_named_debt,",
        "  terminal_exit",
        ")",
        "",
        "Preserved LRC predicate: a local edge quotient is legal only if the "
        "LRC14 packet route remains constant on the quotient fiber, or the "
        "lost tail/tip coordinate is reconstructed, dual-annihilated, "
        "recursed into a smaller child, or named as residual debt.",
        "Destroyed information if scalarized: endpoint role, four-sector "
        "orientation, which side owns the obstruction, both endpoint-deletion "
        "children, and the proof-circuit/coordinate-repair inputs.",
        "Challenged assumption: tournament vertices need not be runners or "
        "tournament classes.  Here the candidate vertices are directed edges, "
        "endpoint-deletion children, outside sectors, proof obligations, and "
        "repair sidecars.",
        "",
    ]


def proof_route_report() -> List[str]:
    return [
        "4. PROOF-ROUTE CONSEQUENCES",
        "- HYP-3050/HYP-3054/HYP-3106 supply the local observable: the edge tail/tip sector word and cross-sector orientation.",
        "- HYP-3112 supplies an analytic payload: edge boundaries can carry one-runner ear mass and root-motion debt.",
        "- HYP-3115 supplies a finite wall payload: domain-wall and one-swap edges should be classified by their two endpoint-deletion children before using wall counts.",
        "- HYP-3116 supplies the circuit guardrail: an edge shortcut must name its missing-input vector rather than advertise small local complexity.",
        "- HYP-3118 supplies the repair rule: any quotient that loses endpoint role, child recursion, or sector orientation must add a resurrection sidecar or name residual debt.",
        "",
        "Candidate invariant:",
        "  edge_witness_certificate = four_sector_deck + paired_endpoint_deletion_recursion + repair_sidecar",
        "",
        "Next executable test: attach this packet to HYP-3115 one-swap/domain-wall edges and ask which walls become legal observer-gluing discharges, which recurse to smaller tail/tip children, and which remain named HYP-2963/HYP-3098 debt.",
    ]


def main() -> None:
    print("HYP-3124 / codex-2026-06-27-S268")
    print("Tournament edge witness recursion scout; executable synthesis, not a proof.")
    print()
    for block in (packet_schema_report(), census_report(), tournament_report(), proof_route_report()):
        for line in block:
            print(line)
        print()


if __name__ == "__main__":
    main()
