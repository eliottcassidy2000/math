#!/usr/bin/env python3
"""
lrc_endpoint_pressure_formal_s503.py

codex-2026-06-01 S503

Finite audits for the LRC endpoint-pressure formalization.

The mathematical core is tiny but important:

  finite directed graph + indegree >= 1 at every active owner => directed cycle.

For LRC, a nonempty endpoint protection core gives one incoming owner edge for
each active owner once strict protection incidences are recorded as
protector-owner -> endpoint-owner pressure edges.  Thus a pressure DAG can only
coexist with a counterexample if the endpoint core is empty or some protection
incidence is not represented in the pressure lift.
"""

from __future__ import annotations

from itertools import product


def has_cycle(n: int, edges: set[tuple[int, int]]) -> bool:
    visiting = [False] * n
    seen = [False] * n
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)

    def dfs(v: int) -> bool:
        if visiting[v]:
            return True
        if seen[v]:
            return False
        visiting[v] = True
        for w in adj[v]:
            if dfs(w):
                return True
        visiting[v] = False
        seen[v] = True
        return False

    return any(dfs(v) for v in range(n))


def indegrees(n: int, edges: set[tuple[int, int]]) -> list[int]:
    deg = [0] * n
    for _, b in edges:
        deg[b] += 1
    return deg


def all_digraphs_without_loops(n: int):
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    for bits in product([0, 1], repeat=len(pairs)):
        yield {edge for edge, bit in zip(pairs, bits) if bit}


def audit_min_indegree_cycle(max_n: int = 4) -> None:
    print("FINITE OWNER-GRAPH LEMMA")
    print("Every loopless finite digraph with indegree >= 1 at each vertex has a directed cycle.")
    print(f"{'n':>3} {'digraphs':>10} {'min-in>=1':>10} {'acyclic among them':>18}")
    for n in range(1, max_n + 1):
        total = good = bad = 0
        for edges in all_digraphs_without_loops(n):
            total += 1
            if min(indegrees(n, edges), default=0) >= 1:
                good += 1
                if not has_cycle(n, edges):
                    bad += 1
        print(f"{n:3d} {total:10d} {good:10d} {bad:18d}")
    print()


def audit_protector_selectors(max_n: int = 8) -> None:
    print("PROTECTOR SELECTOR MODEL")
    print("Each active owner chooses one strict protector owner different from itself.")
    print("This functional subgraph is the minimal shadow of an endpoint core.")
    print(f"{'owners':>6} {'selectors':>12} {'acyclic selectors':>18}")
    for n in range(2, max_n + 1):
        choices = []
        for owner in range(n):
            choices.append([p for p in range(n) if p != owner])
        total = bad = 0
        for selector in product(*choices):
            total += 1
            edges = {(selector[owner], owner) for owner in range(n)}
            if not has_cycle(n, edges):
                bad += 1
        print(f"{n:6d} {total:12d} {bad:18d}")
    print()


def print_lrc_trilemma() -> None:
    print("LRC FORMAL CERTIFICATE TRILEMMA")
    print("For a speed set V at denominator n, a no-lonely open cover must pass:")
    print("  1. denominator sieve completeness (THM-366/369),")
    print("  2. nonempty terminal endpoint protection core (THM-357/359),")
    print("  3. a cyclic owner-pressure graph if core protection is pressure-realized (THM-379/380).")
    print()
    print("Therefore a proof search can stop on any of these certificates:")
    print("  missing small denominator -> rational lonely witness;")
    print("  empty endpoint core       -> no all-protected open cover;")
    print("  pressure DAG + realized core incidences -> no counterexample core.")
    print()
    print("The remaining disproof-shaped object is narrow:")
    print("  sieve-complete + nonempty endpoint core + nontrivial labelled pressure SCC.")


def main() -> None:
    print("LRC endpoint-pressure formalization (codex-2026-06-01 S503)")
    print("=" * 76)
    audit_min_indegree_cycle()
    audit_protector_selectors()
    print_lrc_trilemma()


if __name__ == "__main__":
    main()
