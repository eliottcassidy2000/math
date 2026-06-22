#!/usr/bin/env python3
"""
Finite signed packet graph for the HYP-2632 repeated-residue kernel.

HYP-2632 found the d=9 two-large support-six kernel:

  4+2 packets:   (1,1,1,1,a,a) have negative weights.
  4+1+1 packets: (1,1,1,1,a,b) have positive or zero weights.

This script forgets the floating coimage coefficient unit U and keeps the exact
integer packet weights.  The useful new observation is a Kirchhoff-style
balance: for each residue vertex a in {0,2,3,4,5,6}, the negative repeated
packet at a is exactly cancelled by the positive edge packets incident to a.

This is the first concrete version of HYP-2883's signed packet graph.  It is a
cycle-space address for the support-six LRC residual, not a runner graph.
"""
from __future__ import annotations

import itertools
from collections import Counter, defaultdict, deque

MOD = 7
VERTICES = (0, 2, 3, 4, 5, 6)


def chi7(x: int) -> int:
    x %= MOD
    if x == 0:
        return 0
    return 1 if x in {1, 2, 4} else -1


def q_selector(a: int, b: int) -> int:
    return (a * b * (1 + 3 * ((a + b) % MOD)) - 1) % MOD


def loop_weight(a: int) -> int:
    """Weight of packet (1,1,1,1,a,a), measured in HYP-2632 units U."""
    if a == 0:
        return -4
    return (-43 - 7 * chi7(a)) // 2


def edge_weight(a: int, b: int) -> int:
    """Weight of packet (1,1,1,1,a,b), measured in HYP-2632 units U."""
    if (a + b) % MOD == 2:
        return 0
    return 8 if chi7(q_selector(a, b)) == 1 else 1


def edge_rows() -> list[tuple[int, int, int, str]]:
    rows = []
    for a, b in itertools.combinations(VERTICES, 2):
        w = edge_weight(a, b)
        if w == 0:
            kind = "zero"
        elif w == 8:
            kind = "high"
        else:
            kind = "low"
        rows.append((a, b, w, kind))
    return rows


def simple_edges(kind: str) -> set[tuple[int, int]]:
    return {
        (a, b)
        for a, b, _w, k in edge_rows()
        if k == kind
    }


def degree_hist(edges: set[tuple[int, int]]) -> dict[int, int]:
    deg = Counter({v: 0 for v in VERTICES})
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    return dict(sorted(Counter(deg.values()).items()))


def triangle_count(edges: set[tuple[int, int]]) -> int:
    edge_set = {tuple(sorted(e)) for e in edges}
    count = 0
    for tri in itertools.combinations(VERTICES, 3):
        if all(tuple(sorted(e)) in edge_set for e in itertools.combinations(tri, 2)):
            count += 1
    return count


def c4_count(edges: set[tuple[int, int]]) -> int:
    edge_set = {tuple(sorted(e)) for e in edges}
    count = 0
    for quad in itertools.combinations(VERTICES, 4):
        # Three cyclic orders modulo reversal for four labelled vertices.
        a, b, c, d = quad
        cycles = [
            ((a, b), (b, c), (c, d), (d, a)),
            ((a, b), (b, d), (d, c), (c, a)),
            ((a, c), (c, b), (b, d), (d, a)),
        ]
        for cyc in cycles:
            if all(tuple(sorted(e)) in edge_set for e in cyc):
                count += 1
    return count


def components(edges: set[tuple[int, int]]) -> list[list[int]]:
    adj = {v: set() for v in VERTICES}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen: set[int] = set()
    comps: list[list[int]] = []
    for root in VERTICES:
        if root in seen:
            continue
        q = deque([root])
        seen.add(root)
        comp = []
        while q:
            v = q.popleft()
            comp.append(v)
            for u in sorted(adj[v]):
                if u not in seen:
                    seen.add(u)
                    q.append(u)
        comps.append(comp)
    return comps


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def packet_table() -> None:
    section("SIGNED PACKET TABLE")
    print("Vertex packets are 4+2 repeated residues; edge packets are 4+1+1 residues.")
    print(f"{'v':>2} {'chi7(v)':>7} {'loop':>6} {'incident':>8} {'balance':>8}")
    rows = edge_rows()
    failures = []
    for v in VERTICES:
        incident = sum(w for a, b, w, _k in rows if v in (a, b))
        bal = loop_weight(v) + incident
        if bal:
            failures.append((v, bal))
        print(f"{v:>2} {chi7(v):>7} {loop_weight(v):>6} {incident:>8} {bal:>8}")
    print(f"balance failures: {failures}")

    print("\nEdges:")
    print(f"{'edge':>8} {'a+b':>5} {'Q':>3} {'chi7(Q)':>8} {'weight':>6} {'kind':>6}")
    for a, b, w, kind in rows:
        q = q_selector(a, b)
        print(
            f"{str((a,b)):>8} {(a+b)%MOD:>5} {q:>3} {chi7(q):>8} {w:>6} {kind:>6}"
        )


def graph_fingerprints() -> None:
    section("GRAPH FINGERPRINTS")
    rows = edge_rows()
    by_kind = defaultdict(list)
    for a, b, _w, kind in rows:
        by_kind[kind].append((a, b))

    loop_sum = sum(loop_weight(v) for v in VERTICES)
    edge_sum = sum(w for _a, _b, w, _kind in rows)
    print(f"loop sum: {loop_sum}")
    print(f"edge sum: {edge_sum}")
    print(f"loop + edge: {loop_sum + edge_sum}")
    print(f"loop + 2*edge incidence balance: {loop_sum + 2 * edge_sum}")
    print(
        "Readout: scalar edge-counting leaves net -54 U, but incidence-counting "
        "gives exact zero.  The support-six packet graph is locally balanced."
    )

    for kind in ("zero", "low", "high"):
        edges = set(by_kind[kind])
        print(f"\n{kind} edges: {sorted(edges)}")
        print(f"  degree_hist: {degree_hist(edges)}")
        print(f"  components: {components(edges)}")
        print(f"  triangles: {triangle_count(edges)}")
        print(f"  c4_cycles: {c4_count(edges)}")

    nonzero_edges = {(a, b) for a, b, w, _kind in rows if w != 0}
    print("\nnonzero graph:")
    print(f"  degree_hist: {degree_hist(nonzero_edges)}")
    print(f"  components: {components(nonzero_edges)}")
    print(f"  triangles: {triangle_count(nonzero_edges)}")
    print(f"  c4_cycles: {c4_count(nonzero_edges)}")


def tournament_analysis_note() -> None:
    section("TOURNAMENT ANALYSIS NOTE")
    print("Pairwise observable:")
    print("  W(a,b) = S_9(1,1,1,1,a,b)/U in {0,1,8}.")
    print("Switch/gauge:")
    print("  first split by affine zero lane a+b=2 mod 7; off that lane split by chi7(Q).")
    print("Tie Hamiltonian path:")
    print("  additive_fourier_kernel > affine_zero_matching > Legendre_Q_edges >")
    print("  local_balance_identity > reciprocal_lift > raw_runner_labels")
    print(
        "Why this is not forced into an orientation tournament: W(a,b) is a "
        "symmetric packet coupling.  Orienting it would destroy the affine "
        "zero matching, which is the point of the audit.  The correct object is "
        "a signed graph plus a proof-order tournament on quotient stages."
    )
    print("Challenged assumption:")
    print(
        "  The repeated-residue tail should be bounded by absolute packet mass. "
        "  The exact local balance says the finite kernel is already a "
        "signed-current graph; the analytic tail should preserve this incidence "
        "before applying triangle inequalities."
    )


def main() -> None:
    section("LRC14 REPEATED PACKET GRAPH - CODEX S101")
    packet_table()
    graph_fingerprints()
    tournament_analysis_note()
    section("S101 READING")
    print(
        "The HYP-2632 finite kernel is a perfectly balanced signed graph: every "
        "negative 4+2 vertex packet is exactly cancelled by the positive "
        "4+1+1 edge packets incident to that residue."
    )
    print(
        "This gives HYP-2883 a concrete proof object.  The remaining analytic "
        "task is to lift this local finite balance through the reciprocal "
        "hyperplane sums after finite low-height wall deletion."
    )


if __name__ == "__main__":
    main()
