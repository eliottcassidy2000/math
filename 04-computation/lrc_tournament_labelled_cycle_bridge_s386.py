#!/usr/bin/env python3
"""
lrc_tournament_labelled_cycle_bridge_s386.py

codex-2026-05-31 S386

Creative bridge between THM-365's labelled LRC endpoint cycles and the
tournament machinery already living in the repository.

The central observation is exact rather than metaphorical:

    Along a Hamiltonian path of a tournament, every backward arc v_j -> v_i
    protects the whole cut interval i < k <= j.

So the good-cut/SCC theorem is a tournament cut-protection system on a line,
while the Lonely Runner endpoint program is a labelled endpoint-protection
system on a circle.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class TournamentRow:
    label: str
    n: int
    protected_cuts: tuple[int, ...]
    bad_cuts: tuple[int, ...]
    backward_arcs: tuple[tuple[int, int], ...]
    components: int
    protection_cycle: tuple[int, ...]


def circulant_tournament(n: int, steps: set[int]) -> set[tuple[int, int]]:
    arcs: set[tuple[int, int]] = set()
    for i in range(n):
        for j in range(i + 1, n):
            if (j - i) % n in steps:
                arcs.add((i, j))
            else:
                arcs.add((j, i))
    return arcs


def transitive_tournament(n: int) -> set[tuple[int, int]]:
    return {(i, j) for i in range(n) for j in range(i + 1, n)}


def one_core_tail() -> set[tuple[int, int]]:
    """A 3-cycle followed by two singleton condensation components."""

    arcs = {(0, 1), (1, 2), (2, 0), (3, 4)}
    for i in range(3):
        for j in (3, 4):
            arcs.add((i, j))
    return arcs


def protects_cut(arc: tuple[int, int], cut: int) -> bool:
    tail, head = arc
    return head < cut <= tail


def backward_arcs(n: int, arcs: set[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((tail, head) for tail, head in arcs if tail > head))


def cut_protectors(
    n: int, arcs: set[tuple[int, int]]
) -> dict[int, tuple[tuple[int, int], ...]]:
    backs = backward_arcs(n, arcs)
    return {
        cut: tuple(arc for arc in backs if protects_cut(arc, cut))
        for cut in range(1, n)
    }


def protection_graph(
    n: int, arcs: set[tuple[int, int]]
) -> dict[int, set[int]]:
    graph = {cut: set() for cut in range(1, n)}
    for cut, protectors in cut_protectors(n, arcs).items():
        for tail, head in protectors:
            # The protected interval of cuts is [head+1, tail].
            # Its boundary cuts are the places where the protector starts/stops
            # being visible inside the Hamiltonian path cut line.
            for endpoint_cut in (head + 1, tail):
                if 1 <= endpoint_cut <= n - 1:
                    graph[cut].add(endpoint_cut)
    return graph


def find_cycle(graph: dict[int, set[int]]) -> tuple[int, ...]:
    seen: set[int] = set()
    stack: list[int] = []
    on_stack: set[int] = set()

    def visit(node: int) -> tuple[int, ...] | None:
        seen.add(node)
        stack.append(node)
        on_stack.add(node)
        for nxt in sorted(graph[node]):
            if nxt not in seen:
                found = visit(nxt)
                if found:
                    return found
            elif nxt in on_stack:
                i = stack.index(nxt)
                return tuple(stack[i:] + [nxt])
        stack.pop()
        on_stack.remove(node)
        return None

    for node in sorted(graph):
        if node not in seen:
            found = visit(node)
            if found:
                return found
    return tuple()


def scc_count(n: int, arcs: set[tuple[int, int]]) -> int:
    adj = {i: [] for i in range(n)}
    radj = {i: [] for i in range(n)}
    for a, b in arcs:
        adj[a].append(b)
        radj[b].append(a)

    order: list[int] = []
    seen: set[int] = set()

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
    components = 0

    def rdfs(v: int) -> None:
        seen.add(v)
        for w in radj[v]:
            if w not in seen:
                rdfs(w)

    for v in reversed(order):
        if v not in seen:
            components += 1
            rdfs(v)
    return components


def tournament_row(label: str, n: int, arcs: set[tuple[int, int]]) -> TournamentRow:
    protectors = cut_protectors(n, arcs)
    protected = tuple(cut for cut, ps in protectors.items() if ps)
    bad = tuple(cut for cut, ps in protectors.items() if not ps)
    return TournamentRow(
        label=label,
        n=n,
        protected_cuts=protected,
        bad_cuts=bad,
        backward_arcs=backward_arcs(n, arcs),
        components=scc_count(n, arcs),
        protection_cycle=find_cycle(protection_graph(n, arcs)),
    )


def print_schema() -> None:
    print("LABELLED PROTECTION SCHEMA")
    print("=" * 88)
    rows = [
        (
            "LRC endpoint program",
            "endpoint e",
            "forbidden interval that strictly contains e",
            "speed labels (u,p,m,eps,a) and strict slack",
            "nonempty core => endpoint cycle",
        ),
        (
            "Tournament good cuts",
            "Hamiltonian-path cut k",
            "backward arc v_j -> v_i crossing the cut",
            "root label e_j-e_i and interval i<k<=j",
            "all cuts protected inside each SCC",
        ),
        (
            "Endpoint transfer",
            "parent quotient row",
            "odd child column / private pivot",
            "incidence matrix entry over F2",
            "rank core, private-pivot or peel order",
        ),
        (
            "OCF/Omega",
            "odd-cycle packet",
            "vertex-disjoint compatible packet",
            "root-relation support and fugacity weight",
            "hard-core independent packet family",
        ),
    ]
    print("system                 atom                  protector")
    print("                       label                 core/cycle meaning")
    print("-" * 88)
    for system, atom, protector, label, meaning in rows:
        print(f"{system:<22} {atom:<22} {protector}")
        print(f"{'':<22} {label:<44} {meaning}")
    print()


def print_tournament_cut_audit() -> None:
    print("TOURNAMENT CUT-PROTECTION AUDIT")
    print("=" * 88)
    examples = [
        tournament_row("transitive T4", 4, transitive_tournament(4)),
        tournament_row("directed C3", 3, {(0, 1), (1, 2), (2, 0)}),
        tournament_row("3-core plus tail", 5, one_core_tail()),
        tournament_row("circulant C5", 5, circulant_tournament(5, {1, 2})),
    ]
    print(
        "label              n scc goodCuts badCuts backward_arcs "
        "cut_cycle"
    )
    print("-" * 88)
    for row in examples:
        cycle = "none" if not row.protection_cycle else "->".join(
            str(x) for x in row.protection_cycle
        )
        print(
            f"{row.label:<18} {row.n:>1} {row.components:>3} "
            f"{len(row.protected_cuts):>8} {str(row.bad_cuts):<10} "
            f"{row.backward_arcs!s:<28} {cycle}"
        )
    print()
    print(
        "THM-354 says goodCuts = n - scc.  Read through this protection lens, "
        "bad cuts are unprotected endpoints and strong components are the "
        "surviving protection cores."
    )
    print()


def print_deep_analogies() -> None:
    print("DEEP ANALOGIES")
    print("=" * 88)
    analogies = [
        (
            "Good-cut/SCC",
            "LRC endpoints are circular cuts; protector intervals are backward arcs on a circle.",
            "Try a condensation theorem for endpoint cores: quotient layers should form a transitive leak order.",
        ),
        (
            "Endpoint-transfer rank",
            "The labelled LRC cycle is an incidence object, not ordinary graph adjacency.",
            "Use private endpoints, triangular minors, or leaf-peelable hypergraphs instead of overlap density.",
        ),
        (
            "Root-sign representation",
            "Tournament cycles are zero sums of type-A roots; LRC protection cycles are near-zero relations among characters t -> vt.",
            "Look for a root/character packet module whose boundary is the strict slack vector.",
        ),
        (
            "OCF hard-core packets",
            "Omega counts compatible odd-cycle packets; LRC disproof would be a compatible endpoint-protector packet.",
            "Build an LRC Omega graph where vertices are labelled protection arrows and conflicts are inconsistent endpoints/speeds.",
        ),
        (
            "Projection residue",
            "Both theories fail when a scalar shadow forgets incidence: score forgets cycle packets, speed lists forget endpoint labels.",
            "TDA features should keep residue rank, peel depth, private labels, and slack, not just H or gap width.",
        ),
        (
            "Single-core anomalies",
            "H=63 complete-Omega single-core and LRC all-protected core are opposite concentrated-core regimes.",
            "Search for endpoint single-core signatures with weighted language gaps, like r_core(signature).",
        ),
    ]
    for name, analogy, next_move in analogies:
        print(f"[{name}]")
        print(f"  analogy:   {analogy}")
        print(f"  next move: {next_move}")
    print()


def print_new_questions() -> None:
    print("NEW QUESTIONS")
    print("=" * 88)
    questions = [
        "Can THM-354 be rewritten as an endpoint-protection theorem, then generalized from path cuts to circular LRC endpoints?",
        "Is there an LRC analogue of condensation order: a quotient-layer partial order that every labelled protection cycle must violate?",
        "Can labelled endpoint cycles be triangularized by private endpoint witnesses, just like THM-356 private child columns?",
        "What is the LRC Omega graph: labelled protection arrows as vertices, conflicts from incompatible speed/endpoint labels, fugacity measuring repair cost?",
        "Does the n=14 seven-ladder fail because its putative cycle lives in the wrong projected graph, exactly like endpoint-collision triples that are not parent-metagraph triangles?",
        "Can cycle slack be represented as a character-valued boundary map, with positive leak analogous to nonzero rank or SCC defect?",
    ]
    for i, question in enumerate(questions, 1):
        print(f"{i}. {question}")
    print()


def main() -> None:
    print("LRC labelled cycles vs tournament machinery (codex-2026-05-31 S386)")
    print()
    print_schema()
    print_tournament_cut_audit()
    print_deep_analogies()
    print_new_questions()


if __name__ == "__main__":
    main()
