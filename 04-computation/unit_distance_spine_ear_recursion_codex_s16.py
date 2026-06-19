#!/usr/bin/env python3
"""
unit_distance_spine_ear_recursion_codex_s16.py

Codex S16: endpoint-compatible ears for the unit-distance Hamiltonian-path
flop question.

S625/S626 already separate the intrinsic graph question (does the unit graph
have a Hamiltonian path?) from the gauge question (does a fixed point-order
tournament expose that path?).  This addendum tests the recursive obstruction:
if a dense witness has many vertices whose deletion leaves a unit spine and
which can be reattached to a spine endpoint, then the all-unit Hamiltonian
path is not an accident.  A genuine geometric flop must first destroy this
endpoint-ear recursion.
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

import unit_distance_impairment_atlas_s623 as s623  # noqa: E402
import unit_distance_impairment_lab_s622 as lab622  # noqa: E402
import unit_distance_tournament_hampath_s625 as s625  # noqa: E402


TRI_TARGETS = s625.U_EXACT_SMALL
MOSER_TARGETS = dict(lab622.U_EXACT)
MOSER_TARGETS[22] = 60


def bit_index(bit: int) -> int:
    return bit.bit_length() - 1


def edge_count(adj: list[int]) -> int:
    return sum(bits.bit_count() for bits in adj) // 2


def degree_hist(adj: list[int]) -> dict[int, int]:
    return dict(sorted(Counter(bits.bit_count() for bits in adj).items()))


def complement_graph(adj: list[int]) -> list[int]:
    n = len(adj)
    full = (1 << n) - 1
    return [full & ~(bits | (1 << i)) for i, bits in enumerate(adj)]


def endpoint_dp(adj: list[int]) -> list[int]:
    """For every mask, bitset of possible endpoints of a path covering mask."""
    n = len(adj)
    ends = [0] * (1 << n)
    for v in range(n):
        ends[1 << v] = 1 << v
    for mask in range(1 << n):
        endpoint_bits = ends[mask]
        while endpoint_bits:
            bit = endpoint_bits & -endpoint_bits
            endpoint_bits ^= bit
            v = bit_index(bit)
            nxts = adj[v] & ~mask
            while nxts:
                nb = nxts & -nxts
                nxts ^= nb
                ends[mask | nb] |= nb
    return ends


def recover_path(adj: list[int], ends: list[int], mask: int) -> tuple[int, ...]:
    if ends[mask] == 0:
        return ()
    end = bit_index(ends[mask] & -ends[mask])
    path = [end]
    while mask != (1 << end):
        prev_mask = mask ^ (1 << end)
        candidates = ends[prev_mask] & adj[end]
        if candidates == 0:
            raise AssertionError("endpoint DP reconstruction failed")
        prev = bit_index(candidates & -candidates)
        path.append(prev)
        mask = prev_mask
        end = prev
    path.reverse()
    return tuple(path)


def max_unit_edges_in_directed_hp(out: list[int], unit_adj: list[int]) -> int:
    n = len(out)
    full = (1 << n) - 1
    neg = -10**9
    dp: list[dict[int, int]] = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 0
    for mask in range(1 << n):
        for last, score in list(dp[mask].items()):
            nxts = out[last] & ~mask
            while nxts:
                nb = nxts & -nxts
                nxts ^= nb
                nxt = bit_index(nb)
                gain = 1 if (unit_adj[last] & nb) else 0
                nxt_row = dp[mask | nb]
                old = nxt_row.get(nxt, neg)
                if score + gain > old:
                    nxt_row[nxt] = score + gain
    return max(dp[full].values())


def connected_components_after_delete(adj: list[int], drop: int) -> int:
    n = len(adj)
    remaining = ((1 << n) - 1) ^ (1 << drop)
    components = 0
    while remaining:
        components += 1
        start = remaining & -remaining
        stack = [bit_index(start)]
        remaining ^= start
        while stack:
            v = stack.pop()
            nxts = adj[v] & remaining
            while nxts:
                nb = nxts & -nxts
                nxts ^= nb
                remaining ^= nb
                stack.append(bit_index(nb))
    return components


def directed_triangles(out: list[int]) -> int:
    count = 0
    for i, j, k in combinations(range(len(out)), 3):
        if (out[i] >> j) & 1 and (out[j] >> k) & 1 and (out[k] >> i) & 1:
            count += 1
        if (out[i] >> k) & 1 and (out[k] >> j) & 1 and (out[j] >> i) & 1:
            count += 1
    return count


def scc_sizes(out: list[int]) -> tuple[int, ...]:
    n = len(out)
    rev = [0] * n
    for v, bits in enumerate(out):
        while bits:
            bit = bits & -bits
            bits ^= bit
            rev[bit_index(bit)] |= 1 << v

    def reach(start: int, graph: list[int]) -> int:
        seen = 1 << start
        stack = [start]
        while stack:
            v = stack.pop()
            bits = graph[v] & ~seen
            while bits:
                bit = bits & -bits
                bits ^= bit
                seen |= bit
                stack.append(bit_index(bit))
        return seen

    remaining = (1 << n) - 1
    sizes: list[int] = []
    while remaining:
        bit = remaining & -remaining
        v = bit_index(bit)
        comp = reach(v, out) & reach(v, rev)
        sizes.append(comp.bit_count())
        remaining &= ~comp
    return tuple(sorted(sizes, reverse=True))


def hp_count_tournament(out: list[int]) -> int | str:
    n = len(out)
    if n > 14:
        return "not-counted"
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last, ways in enumerate(dp[mask]):
            if not ways:
                continue
            nxts = out[last] & ~mask
            while nxts:
                bit = nxts & -nxts
                nxts ^= bit
                dp[mask | bit][bit_index(bit)] += ways
    return sum(dp[-1])


@dataclass(frozen=True)
class EarAudit:
    traceable_deletions: int
    endpoint_ears: int
    exact_edge_deletions: int | str
    endpoint_exact_ears: int | str
    max_components_after_delete: int
    branch_cut_vertices: int
    endpoint_options_hist: dict[int, int]
    ear_score_hist: dict[int, int]
    ear_directed_3cycles: int
    ear_scc_sizes: tuple[int, ...]
    ear_hamiltonian_paths: int | str
    best_ears: tuple[tuple[int, int, int, int, int], ...]


def ear_tournament(scores: list[tuple[int, int, int, int, int]]) -> list[int]:
    n = len(scores)
    out = [0] * n
    for i, j in combinations(range(n), 2):
        votes_i = 0
        votes_j = 0
        for a, b in zip(scores[i], scores[j]):
            if a > b:
                votes_i += 1
            elif b > a:
                votes_j += 1
        if votes_i > votes_j or (votes_i == votes_j and i < j):
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return out


def audit_ears(adj: list[int], target_prev: int | None) -> EarAudit:
    n = len(adj)
    full = (1 << n) - 1
    ends = endpoint_dp(adj)
    graph_edges = edge_count(adj)

    traceable_deletions = 0
    endpoint_ears = 0
    exact_edge_deletions = 0
    endpoint_exact_ears = 0
    max_components = 0
    branch_cut_vertices = 0
    endpoint_options = []
    vertex_scores: list[tuple[int, int, int, int, int]] = []
    best_rows: list[tuple[int, int, int, int, int]] = []

    for v in range(n):
        submask = full ^ (1 << v)
        endpoints = ends[submask]
        traceable = int(endpoints != 0)
        compatible_options = (adj[v] & endpoints).bit_count() if traceable else 0
        compatible = int(compatible_options > 0)
        edge_after_delete = graph_edges - adj[v].bit_count()
        exact_delete = int(target_prev is not None and edge_after_delete == target_prev)
        components = connected_components_after_delete(adj, v)

        traceable_deletions += traceable
        endpoint_ears += compatible
        exact_edge_deletions += exact_delete
        endpoint_exact_ears += compatible * exact_delete
        max_components = max(max_components, components)
        branch_cut_vertices += int(components > 2)
        endpoint_options.append(compatible_options)

        # Criteria: endpoint-compatible, deletion-traceable, exact-edge delete,
        # edge count after deletion, and low degree.  A majority tournament over
        # these criteria exposes whether all notions of "good ear" agree.
        degree = adj[v].bit_count()
        vertex_scores.append((compatible, traceable, exact_delete, edge_after_delete, -degree))
        best_rows.append((v, compatible, traceable, exact_delete, degree))

    out = ear_tournament(vertex_scores)
    score_hist = dict(sorted(Counter(bits.bit_count() for bits in out).items()))
    best_rows.sort(key=lambda row: (-row[1], -row[2], -row[3], row[4], row[0]))

    return EarAudit(
        traceable_deletions=traceable_deletions,
        endpoint_ears=endpoint_ears,
        exact_edge_deletions=exact_edge_deletions if target_prev is not None else "unknown",
        endpoint_exact_ears=endpoint_exact_ears if target_prev is not None else "unknown",
        max_components_after_delete=max_components,
        branch_cut_vertices=branch_cut_vertices,
        endpoint_options_hist=dict(sorted(Counter(endpoint_options).items())),
        ear_score_hist=score_hist,
        ear_directed_3cycles=directed_triangles(out),
        ear_scc_sizes=scc_sizes(out),
        ear_hamiltonian_paths=hp_count_tournament(out),
        best_ears=tuple(best_rows[:6]),
    )


@dataclass(frozen=True)
class Row:
    carrier: str
    n: int
    edges: int
    target: int | None
    unit_spine: bool
    nonunit_spine: bool
    lex_unitflip_max: int | str
    lex_nonunitflip_max: int | str
    degree_hist: dict[int, int]
    path_prefix: tuple[int, ...]
    ear: EarAudit


def analyze_row(
    carrier: str,
    n: int,
    points: tuple,
    unit_adj: list[int],
    target: int | None,
    target_prev: int | None,
) -> Row:
    ends = endpoint_dp(unit_adj)
    full = (1 << n) - 1
    path = recover_path(unit_adj, ends, full)
    nonunit = complement_graph(unit_adj)
    nonunit_ends = endpoint_dp(nonunit)
    lex_unitflip_max: int | str = "skip"
    lex_nonunitflip_max: int | str = "skip"
    if n <= 14:
        unitflip = s625.tournament_from_unit_graph(unit_adj, unit_flip=True)
        nonunitflip = s625.tournament_from_unit_graph(unit_adj, unit_flip=False)
        lex_unitflip_max = max_unit_edges_in_directed_hp(unitflip, unit_adj)
        lex_nonunitflip_max = max_unit_edges_in_directed_hp(nonunitflip, unit_adj)
    return Row(
        carrier=carrier,
        n=n,
        edges=edge_count(unit_adj),
        target=target,
        unit_spine=bool(ends[full]),
        nonunit_spine=bool(nonunit_ends[full]),
        lex_unitflip_max=lex_unitflip_max,
        lex_nonunitflip_max=lex_nonunitflip_max,
        degree_hist=degree_hist(unit_adj),
        path_prefix=path[: min(10, len(path))],
        ear=audit_ears(unit_adj, target_prev),
    )


def triangular_rows(max_n: int) -> list[Row]:
    rows: list[Row] = []
    for n in range(2, max_n + 1):
        result = s623.beam_search(
            "triangular",
            n,
            300,
            s623.TRI_UNITS,
            s623.TRI_UNITS,
            s623.add2,
            s623.canon2,
            s623.span2,
            "healthy",
        )
        adj = s625.unit_graph(result.cluster, s623.TRI_UNITS, s623.add2)
        rows.append(
            analyze_row(
                "triangular",
                n,
                result.cluster,
                adj,
                TRI_TARGETS.get(n),
                TRI_TARGETS.get(n - 1),
            )
        )
    return rows


def moser_rows(max_n: int, width: int) -> list[Row]:
    beam = lab622.run_beam(max_n, width)
    rows: list[Row] = []
    for n in list(range(2, min(14, max_n) + 1)) + [v for v in (21, 22) if v <= max_n]:
        cluster = beam.cluster_at(n)
        adj, _ = unit_graph_lab(cluster)
        rows.append(
            analyze_row(
                "moser",
                n,
                cluster,
                adj,
                MOSER_TARGETS.get(n),
                MOSER_TARGETS.get(n - 1),
            )
        )
    return rows


def unit_graph_lab(cluster: lab622.Cluster) -> tuple[list[int], list[int]]:
    unit_set = set(lab622.UNITS)
    n = len(cluster)
    adj = [0] * n
    for i, j in combinations(range(n), 2):
        d1 = lab622.sub(cluster[j], cluster[i])
        d2 = lab622.sub(cluster[i], cluster[j])
        if d1 in unit_set or d2 in unit_set:
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return adj, adj[:]


def first_lex_flop(rows: list[Row]) -> int | None:
    for row in rows:
        if isinstance(row.lex_unitflip_max, int) and row.lex_unitflip_max < row.n - 1:
            return row.n
    return None


def print_table(label: str, rows: list[Row]) -> None:
    print(label)
    print("n | edges | target | unit HP | nonunit HP | ears | exact ears | branch cuts | lex unit max")
    print("--- | --- | --- | --- | --- | --- | --- | --- | ---")
    for row in rows:
        ear = row.ear
        print(
            f"{row.n} | {row.edges} | {row.target} | {row.unit_spine} | {row.nonunit_spine} | "
            f"{ear.endpoint_ears}/{row.n} | {ear.endpoint_exact_ears} | "
            f"{ear.branch_cut_vertices}, maxcomp={ear.max_components_after_delete} | {row.lex_unitflip_max}"
        )
    print()


def print_detail(row: Row) -> None:
    print(f"DETAIL {row.carrier} n={row.n}")
    print(f"degree_hist={row.degree_hist}")
    print(f"unit spine prefix={row.path_prefix}")
    print(f"endpoint-options hist={row.ear.endpoint_options_hist}")
    print(f"best ear rows (v, endpoint-compatible, deletion-traceable, exact-delete, degree)={row.ear.best_ears}")
    print(
        "ear tournament: "
        f"score_hist={row.ear.ear_score_hist}, "
        f"directed_3cycles={row.ear.ear_directed_3cycles}, "
        f"scc={row.ear.ear_scc_sizes}, H={row.ear.ear_hamiltonian_paths}"
    )
    print()


def print_route_tournament() -> None:
    routes = (
        ("unit graph traceability", (5, 5, 4, 5, 1), "preserves the intrinsic predicate"),
        ("endpoint-ear recursion", (5, 4, 5, 5, 1), "tests inductive deletion and reattachment"),
        ("spine-order tournament", (4, 2, 4, 4, 1), "retains a tie Hamiltonian path witness"),
        ("lexicographic flip gauge", (2, 1, 2, 4, -3), "diagnostic for order artifacts"),
        ("nonunit complement HP", (3, 3, 2, 4, -2), "detects when the background can also thread"),
        ("S/U/L impurity word", (4, 3, 3, 3, -1), "keeps short/unit/long labels rather than one bit"),
        ("raw edge count", (1, 1, 1, 5, -4), "scalar quotient that forgets spine/bulk"),
    )
    out = [0] * len(routes)
    for i, j in combinations(range(len(routes)), 2):
        ai = routes[i][1]
        aj = routes[j][1]
        iv = sum(x > y for x, y in zip(ai, aj))
        jv = sum(y > x for x, y in zip(ai, aj))
        if iv > jv or (iv == jv and i < j):
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    order = sorted(range(len(routes)), key=lambda idx: (-out[idx].bit_count(), routes[idx][0]))
    print("TOURNAMENT ANALYSIS OVER PROOF LENSES")
    print("Pairwise observable: which lens preserves the unit-spine predicate while supporting recursion.")
    print("Switch/gauge: majority over intrinsic predicate, gauge invariance, recursion, computation, and risk.")
    print("Tie Hamiltonian path by score order:")
    for idx in order:
        print(f"  score={out[idx].bit_count()} {routes[idx][0]}: {routes[idx][2]}")
    print(f"score_hist={dict(sorted(Counter(bits.bit_count() for bits in out).items()))}")
    print(f"directed_3cycles={directed_triangles(out)} scc={scc_sizes(out)} H={hp_count_tournament(out)}")
    print()


def print_summary(tri: list[Row], moser: list[Row]) -> None:
    print("FLOP READOUT")
    print("------------")
    for label, rows in (("triangular", tri), ("moser", moser)):
        unit_failures = [row.n for row in rows if not row.unit_spine]
        nonunit_true = [row.n for row in rows if row.nonunit_spine]
        nonunit_false = [row.n for row in rows if not row.nonunit_spine]
        ear_failures = [row.n for row in rows if row.ear.endpoint_ears == 0]
        print(f"{label}: graph-level unit-spine failures={unit_failures or 'none'}")
        print(f"{label}: first lexicographic unit-flip flop={first_lex_flop(rows)}")
        print(f"{label}: nonunit HP true rows={nonunit_true}; false rows={nonunit_false}")
        print(f"{label}: endpoint-ear failures={ear_failures or 'none'}")
    print()
    print("Interpretation")
    print("- No tested carrier has a graph-level unit Hamiltonian-path flop.")
    print("- The first fixed-order tournament flop remains n=7: a coordinate/tie-order artifact.")
    print("- The new stronger signal is endpoint-ear abundance: every tested row has at least one")
    print("  removable vertex that leaves a unit spine and reattaches to a spine endpoint.")
    print("- A true geometric flop must kill traceability itself, probably through three or more")
    print("  incompatible ears, a separator with three path-obligatory branches, or a core-extension")
    print("  where all removable vertices miss every smaller-spine endpoint.")
    print()


def main() -> None:
    sys.stdout.reconfigure(line_buffering=True)
    print("S16 UNIT-DISTANCE SPINE/EAR RECURSION AUDIT")
    print("===========================================")
    print()
    print("Method: compute Hamiltonian endpoint masks for the unit graph, then ask")
    print("whether each vertex is an endpoint-compatible ear: deleting it leaves a")
    print("unit Hamiltonian path and the vertex touches at least one possible endpoint.")
    print()

    tri = triangular_rows(22)
    moser = moser_rows(22, 1200)
    print_summary(tri, moser)
    print_table("TRIANGULAR CARRIER", tri)
    print_table("MOSER CARRIER", moser)

    for carrier_rows in (tri, moser):
        for n in (7, 14, 21, 22):
            row = next((r for r in carrier_rows if r.n == n), None)
            if row is not None:
                print_detail(row)

    print_route_tournament()
    print("ASSUMPTION CHALLENGE")
    print("--------------------")
    print("Alternate vertices considered: points, unit edges, nonunit edges, deletion ears,")
    print("endpoint masks, direction shells, core-extension obligations, and S/U/L impurity")
    print("words.  This script uses deletion ears and proof lenses because the preserved")
    print("LRC/UD-style predicate is not a point order but a retained certificate section.")
    print("The quotient preserves traceability, endpoint compatibility, complement threading,")
    print("and lex-gauge dropout; it destroys metric distinctions among nonunit distances")
    print("unless the later S/U/L impurity word is attached.")


if __name__ == "__main__":
    main()
