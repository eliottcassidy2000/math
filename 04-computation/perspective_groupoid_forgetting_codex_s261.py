#!/usr/bin/env python3
"""S261: perspective groupoid scout for controlled forgetting.

This script extends the A000568 perspective ladder into a typed functor
ledger.  The main point is negative and useful: k-depth node views recover
rooted memory, but the first shifted A000568 defect needs an
observer-extension/cut sidecar.  Edge sectors, cycle/chirality, clique
insertion cuts, component-bound packets, first-obstruction cocycles, and
miss-count PGF roots are different functors because they survive different
next operations.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import combinations
from math import atan2, prod

import numpy as np

import a000568_edge_perspective_extension_codex_s213 as edge_lift
import perspective_depth_ladder_codex_s214 as depth
import tournament_rigidity_cascade_s589 as tour


A000568 = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
}


def classes_by_n(n: int) -> tuple[int, ...]:
    if n == 6:
        # Generated as 5-vertex parent + incident word, avoiding full n=6
        # labelled enumeration while still producing all 56 classes.
        return edge_lift.generated_classes6()
    return tour.classes(n)


def opposite(mask: int, n: int) -> int:
    out = 0
    for i, j in combinations(range(n), 2):
        out = tour.set_edge(out, n, i, j, not tour.edge(mask, n, i, j))
    return out


def exact_tables(nmax: int = 6, depth_max: int = 4) -> list[str]:
    lines: list[str] = []
    lines.append("1. SHIFTED A000568 / ROOTED-PERSPECTIVE LADDER")
    lines.append("m  U(m)  P_node(m)  U(m+1)  defect  node_depths")
    for n in range(1, nmax + 1):
        reps = classes_by_n(n)
        p_node = 0
        node_depths = [0 for _ in range(depth_max)]
        for rep in reps:
            auts = depth.automorphisms(rep, n)
            p_node += depth.vertex_orbit_count(rep, n, auts)
            for d in range(1, depth_max + 1):
                node_depths[d - 1] += depth.node_depth_count(rep, n, d)
        target = A000568.get(n + 1)
        defect = None if target is None else target - p_node
        lines.append(
            f"{n:1d}  {len(reps):4d}  {p_node:9d}  {target!s:7s}  "
            f"{defect!s:6s}  {node_depths}"
        )
    lines.append("")
    lines.append("Readout: the prompt's comparison is P_node(m) versus U(m+1).")
    lines.append("The first failure is m=5: P_node(5)=48, U(6)=56, defect=8.")
    lines.append("At m=6 the gap has widened: P_node(6)=296 versus U(7)=456.")
    return lines


def non_node_tables(nmax: int = 6) -> list[str]:
    lines: list[str] = []
    lines.append("2. EXACT NON-NODE PERSPECTIVE CARRIERS")
    lines.append(
        "m  arc_orbits  triple_orbits  transitive_triples  cyclic_triples  "
        "conflict_pair_orbits"
    )
    for n in range(1, nmax + 1):
        arc = triple = trans = cyc = conflict = 0
        for rep in classes_by_n(n):
            auts = depth.automorphisms(rep, n)
            arc += depth.arc_orbit_count(rep, n, auts)
            triple += depth.triple_orbit_count(rep, n, auts)
            trans += depth.triple_orbit_count(rep, n, auts, "transitive")
            cyc += depth.triple_orbit_count(rep, n, auts, "cyclic")
            conflict += depth.conflict_pair_orbit_count(rep, n, auts)
        lines.append(f"{n:1d}  {arc:10d}  {triple:13d}  {trans:18d}  {cyc:14d}  {conflict:19d}")
    lines.append("")
    lines.append(
        "Conflict/Omega-style carriers first become nonzero at m=6; that is "
        "exactly where the shifted rooted-node gap is no longer a one-step toy."
    )
    return lines


def sector_and_converse_audit() -> list[str]:
    classes6 = classes_by_n(6)
    lines: list[str] = []
    lines.append("3. EDGE-SECTOR AND DIHEDRAL/CONVERSE LOSS AT U(6)")
    for mode in ("size", "internal", "cross", "full"):
        decks = {edge_lift.class_sector_deck(rep, 6, mode) for rep in classes6}
        sigs = {
            edge_lift.sector_signature(rep, 6, a, b, mode)
            for rep in classes6
            for a in range(6)
            for b in range(6)
            if a != b
        }
        lines.append(
            f"{mode:8s}: individual_sigs={len(sigs):3d}; "
            f"class_decks={len(decks):2d}/56"
        )

    internal_decks: dict[tuple[tuple[object, ...], ...], list[int]] = defaultdict(list)
    for rep in classes6:
        internal_decks[edge_lift.class_sector_deck(rep, 6, "internal")].append(rep)
    collisions = [reps for reps in internal_decks.values() if len(reps) > 1]
    sc_count = sum(1 for rep in classes6 if tour.canonical(opposite(rep, 6), 6) == rep)
    lines.append(f"self_converse_classes_U6={sc_count}")
    for reps in collisions:
        a, b = sorted(reps)
        lines.append(
            "only_size_internal_collision="
            f"{a},{b}; converse_pair={tour.canonical(opposite(a, 6), 6) == b}"
        )
    lines.append(
        "Interpretation: sector size/internal quotients forget a dihedral/"
        "converse chirality coordinate; cross-sector orientation restores it."
    )
    return lines


def sector_of(x: Fraction) -> int:
    return int((x % 1) * 7)


def miss_distribution(speeds: tuple[int, ...]) -> list[Fraction]:
    speeds = tuple(sorted(set(speeds)))
    breaks = {Fraction(0), Fraction(1)}
    for e in speeds:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            breaks.add(Fraction(m, 7 * e))
    pts = sorted(breaks)
    q = [Fraction(0) for _ in range(7)]
    for a, b in zip(pts, pts[1:]):
        if b <= a:
            continue
        t = (a + b) / 2
        covered = {sector_of(Fraction(e) * t) for e in speeds}
        missing = 7 - len(covered)
        if 0 <= missing <= 6:
            q[missing] += b - a
    return q


def pgf_roots(q: list[Fraction]) -> tuple[int, list[complex]]:
    coeffs = [float(q[i]) for i in range(6, -1, -1)]
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-14:
        coeffs = coeffs[1:]
    roots = list(np.roots(coeffs))
    real_count = sum(1 for z in roots if abs(z.imag) < 1e-7)
    return real_count, sorted(roots, key=lambda z: (abs(z), z.real, z.imag))


def pgf_sidecar_audit() -> list[str]:
    named = {
        "consec_8": tuple(range(8)),
        "even_AP": tuple(2 * i for i in range(8)),
        "top_cluster": (0, 7, 8, 9, 10, 11, 12, 13),
    }
    lines = ["4. MISS-COUNT PGF ROOT SIDECAR (S66 INTEGRATION)"]
    for name, speeds in named.items():
        q = miss_distribution(speeds)
        real_count, roots = pgf_roots(q)
        extreme = float(q[0] + q[6])
        ly = float(10 * q[0] + q[3] + 10 * q[6])
        root_summary = [
            (
                round(float(abs(z)), 4),
                round(float(atan2(z.imag, z.real) * 180 / np.pi), 1),
            )
            for z in roots
        ]
        lines.append(
            f"{name:11s}: #real={real_count}/{len(roots)}  "
            f"extreme_mass={extreme:.4f}  L_yK8={ly:.4f}  roots(|z|,arg)={root_summary}"
        )
    q = miss_distribution(tuple(range(8)))
    _, roots = pgf_roots(q)
    lines.append(f"consec_8_root_mod_product={prod(abs(z) for z in roots):.4f}; q0/q6={float(q[0] / q[6]):.4f}")
    lines.append(
        "Interpretation: PGF root structure is another sidecar.  A scalar "
        "miss-count moment is unsafe unless root-realness / zero-confinement "
        "is retained or proven irrelevant."
    )
    return lines


def functor_ledger() -> list[tuple[str, dict[str, object]]]:
    return [
        (
            "raw_A000568_class",
            dict(root="none", action="S_n", next="add observer", sidecar="none", score=(0, 0, 0, 0, 0, 0, 0, 0, 0), cost=1),
        ),
        (
            "node_k_depth_cache",
            dict(root="vertex color", action="Aut(T) on vertices", next="recover rooted view", sidecar="depth index", score=(1, 1, 0, 0, 0, 0, 0, 0, 0), cost=2),
        ),
        (
            "exact_rooted_node",
            dict(root="marked vertex", action="Aut(T)_root", next="add observer", sidecar="incident word", score=(2, 1, 0, 0, 0, 0, 0, 0, 0), cost=2),
        ),
        (
            "ordered_pair_extension",
            dict(root="old root + new observer", action="Aut(T) on ordered pair", next="unroot/delete", sidecar="endpoint role", score=(3, 3, 1, 0, 0, 1, 0, 0, 1), cost=4),
        ),
        (
            "directed_edge_sector",
            dict(root="tail -> tip", action="edge stabilizer", next="sector quotient", sidecar="cross-sector orientation", score=(2, 3, 3, 1, 0, 1, 0, 1, 1), cost=4),
        ),
        (
            "cycle_chirality",
            dict(root="directed 3-cycle", action="C3/D3 quotient", next="Omega conflict", sidecar="chirality word", score=(1, 1, 1, 3, 0, 1, 0, 2, 0), cost=3),
        ),
        (
            "transitive_clique_insertion",
            dict(root="transitive triple", action="chain stabilizer", next="insert observer", sidecar="cut position", score=(1, 2, 1, 1, 0, 1, 0, 0, 0), cost=3),
        ),
        (
            "conflict_omega",
            dict(root="cycle conflict pair", action="Omega automorphism", next="H/Omega certificate", sidecar="conflict incidence", score=(1, 1, 1, 3, 0, 2, 0, 0, 0), cost=4),
        ),
        (
            "dihedral_reflection_quotient",
            dict(root="reflection orbit", action="D_n / converse", next="time reversal", sidecar="orientation owner", score=(1, 1, 1, 2, 0, 1, 0, 3, 0), cost=3),
        ),
        (
            "normal_fan_component_packet",
            dict(root="safe component", action="normal-fan chamber", next="component bound", sidecar="Cech/barcode word", score=(1, 1, 0, 1, 3, 1, 0, 1, 2), cost=5),
        ),
        (
            "first_obstruction_cocycle",
            dict(root="visible fiber pair", action="payload orbit", next="chart gluing", sidecar="obstruction syndrome", score=(1, 2, 1, 2, 1, 3, 0, 1, 2), cost=5),
        ),
        (
            "miss_count_pgf_roots",
            dict(root="miss-count distribution", action="root multiset", next="moment bound", sidecar="real-root count", score=(0, 1, 0, 0, 1, 1, 3, 0, 0), cost=3),
        ),
        (
            "endpoint_owner_packet",
            dict(root="observer endpoint owner", action="packet automorphism", next="formal certificate", sidecar="full owner ledger", score=(3, 3, 3, 2, 3, 3, 2, 2, 3), cost=7),
        ),
    ]


def carrier_tournament() -> tuple[list[tuple[str, dict[str, object]]], list[list[int]]]:
    carriers = functor_ledger()

    def wins(a: dict[str, object], b: dict[str, object]) -> bool:
        sa = a["score"]
        sb = b["score"]
        assert isinstance(sa, tuple) and isinstance(sb, tuple)
        majority = sum(1 for x, y in zip(sa, sb) if x > y) - sum(1 for x, y in zip(sa, sb) if y > x)
        if majority:
            return majority > 0
        ta = sum(sa) - int(a["cost"])
        tb = sum(sb) - int(b["cost"])
        return ta > tb

    n = len(carriers)
    adj = [[0 for _ in range(n)] for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if wins(carriers[i][1], carriers[j][1]):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return carriers, adj


def directed_three_cycles(adj: list[list[int]]) -> int:
    n = len(adj)
    total = 0
    for a, b, c in combinations(range(n), 3):
        total += int(adj[a][b] and adj[b][c] and adj[c][a])
        total += int(adj[a][c] and adj[c][b] and adj[b][a])
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, graph: list[list[int]]) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for u, bit in enumerate(graph[v]):
                if bit and u not in seen:
                    seen.add(u)
                    q.append(u)
        return seen

    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    out = []
    while remaining:
        v = min(remaining)
        comp = reach(v, adj) & reach(v, rev)
        out.append(len(comp))
        remaining -= comp
    return sorted(out, reverse=True)


def hamiltonian_path_count(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            if not dp[seen][v]:
                continue
            for u in range(n):
                if seen & (1 << u):
                    continue
                if adj[v][u]:
                    dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[-1])


def tournament_analysis_lines() -> list[str]:
    carriers, adj = carrier_tournament()
    out_scores = [sum(row) for row in adj]
    order = sorted(range(len(carriers)), key=lambda i: (-out_scores[i], carriers[i][0]))
    lines = ["5. TOURNAMENT ANALYSIS OVER PERSPECTIVE FUNCTORS"]
    lines.append("Pairwise observable axes: root, extension, edge, cycle, topology, obstruction, analytic roots, dihedral, owner.")
    lines.append("Switch/gauge: orient by majority of retained axes; total retention minus cost breaks ties.")
    lines.append(f"vertices={[name for name, _ in carriers]}")
    lines.append(f"score_hist={dict(sorted(Counter(out_scores).items()))}")
    lines.append(f"directed_3_cycles={directed_three_cycles(adj)}")
    lines.append(f"scc_sizes={scc_sizes(adj)}")
    lines.append(f"hamiltonian_path_count={hamiltonian_path_count(adj)}")
    lines.append("one_priority_path=" + " -> ".join(carriers[i][0] for i in order))
    lines.append("")
    lines.append("Functor summary:")
    for name, data in carriers:
        lines.append(
            f"  {name}: root={data['root']}; action={data['action']}; "
            f"next={data['next']}; sidecar={data['sidecar']}"
        )
    return lines


def main() -> None:
    print("Perspective groupoid controlled-forgetting scout - codex S261")
    print("=" * 78)
    print()
    blocks = [
        exact_tables(),
        non_node_tables(),
        sector_and_converse_audit(),
        pgf_sidecar_audit(),
        tournament_analysis_lines(),
    ]
    for block in blocks:
        print("\n".join(block))
        print()
    print("SYNTHESIS")
    print(
        "A k-depth node perspective is a cache functor: it recovers rooted "
        "memory and then stops.  Edge, cycle, clique, conflict, component, "
        "obstruction, dihedral, owner, and PGF-root perspectives are not "
        "deeper versions of that cache; they are sidecars for different next "
        "operations."
    )
    print(
        "For LRC14, the rule is to choose the functor whose root/action pair "
        "survives the immediate operation: add observer, delete/unroot, reflect, "
        "bound components, glue charts, push H/Omega certificates, or close a "
        "moment/PGF-root inequality."
    )
    print()
    print("ASSUMPTION CHALLENGE")
    print(
        "Considered vertices: runners, tournament nodes, directed edges, "
        "ordered pairs, cycles, cliques, conflict pairs, dihedral orbits, "
        "safe components, obstruction cochains, miss-count roots, endpoint "
        "owners, matrix sidecar columns, and proof obligations."
    )
    print(
        "Chosen vertices: perspective functors.  Preserved predicate: the "
        "quotient's legality for the next LRC operation.  Destroyed data: "
        "labels, endpoint roles, cross-sector orientation, chirality, topology, "
        "obstruction, and analytic root structure unless the named sidecar "
        "carries it."
    )


if __name__ == "__main__":
    main()
