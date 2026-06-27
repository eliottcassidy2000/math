#!/usr/bin/env python3
"""Root-curve, relation-lattice, and ear-extremality maps for HYP-3113.

This is a synthesis scout, not a proof.  It converts the prompt cues
Savitch/Bravais/Lee-Yang/ears into two proof-obligation tournaments.

Map A: root-curve extremality.  Vertices are signals from a single value
through the full PGF/root locus.  The observable asks which signal retains
more of the LRC predicate and the extremal mechanism.

Map B: memory-lattice-ear certificates.  Vertices are proof carriers from raw
runner sets through midpoint recursion, relation-lattice shape, and ear
decomposition certificates.

The output records Tournament Analysis fingerprints: score histograms, SCCs,
directed 3-cycles, Hamiltonian path counts, and a declared tie path.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class Vertex:
    name: str
    scores: dict[str, int]
    preserves: str
    destroys: str
    next_signal: str


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            row = adj[last]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if row[nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def sccs(adj: list[list[bool]]) -> list[list[int]]:
    n = len(adj)
    index = 0
    stack: list[int] = []
    on_stack = [False] * n
    idx = [-1] * n
    low = [0] * n
    out: list[list[int]] = []

    def visit(v: int) -> None:
        nonlocal index
        idx[v] = low[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        for w, has_edge in enumerate(adj[v]):
            if not has_edge:
                continue
            if idx[w] == -1:
                visit(w)
                low[v] = min(low[v], low[w])
            elif on_stack[w]:
                low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp: list[int] = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            out.append(sorted(comp))

    for v in range(n):
        if idx[v] == -1:
            visit(v)
    return out


def directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    n = len(adj)
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                if (adj[a][b] and adj[b][c] and adj[c][a]) or (
                    adj[a][c] and adj[c][b] and adj[b][a]
                ):
                    total += 1
    return total


def orient(vertices: list[Vertex], axes: list[str]) -> tuple[list[list[bool]], dict[str, int], dict[str, tuple[int, int]]]:
    n = len(vertices)
    tie_rank = {v.name: i for i, v in enumerate(vertices)}
    adj = [[False] * n for _ in range(n)]
    edge_votes: dict[str, tuple[int, int]] = {}
    for i, a in enumerate(vertices):
        for j in range(i + 1, n):
            b = vertices[j]
            votes_a = 0
            votes_b = 0
            for axis in axes:
                av = a.scores.get(axis, 0)
                bv = b.scores.get(axis, 0)
                if av > bv:
                    votes_a += 1
                elif bv > av:
                    votes_b += 1
            if votes_a > votes_b or (votes_a == votes_b and tie_rank[a.name] < tie_rank[b.name]):
                adj[i][j] = True
                edge_votes[f"{a.name} -> {b.name}"] = (votes_a, votes_b)
            else:
                adj[j][i] = True
                edge_votes[f"{b.name} -> {a.name}"] = (votes_b, votes_a)
    scores = {vertices[i].name: sum(adj[i]) for i in range(n)}
    return adj, scores, edge_votes


def print_map(title: str, vertices: list[Vertex], axes: list[str]) -> None:
    adj, scores, edge_votes = orient(vertices, axes)
    score_hist = Counter(scores.values())
    comps = [[vertices[i].name for i in comp] for comp in sccs(adj)]
    path_count = hamiltonian_paths(adj)
    cycles = directed_3cycles(adj)
    path = sorted(vertices, key=lambda v: (-scores[v.name], v.name))

    print("=" * 100)
    print(title)
    print("=" * 100)
    print("axes=" + ", ".join(axes))
    print("tie_path=" + " > ".join(v.name for v in vertices))
    print("score_hist=" + str(dict(sorted(score_hist.items()))))
    print("scores=" + str(dict(sorted(scores.items(), key=lambda kv: (-kv[1], kv[0])))))
    print(f"directed_3cycles={cycles}")
    print("SCCs=" + str(comps))
    print(f"Hamiltonian_path_count={path_count}")
    print("selected_retention_path=" + " > ".join(v.name for v in path))
    print()
    print("vertex ledger:")
    for v in vertices:
        print(f"- {v.name}")
        print(f"  preserves: {v.preserves}")
        print(f"  destroys: {v.destroys}")
        print(f"  next_signal: {v.next_signal}")
    print()
    print("sample edge votes:")
    for edge, vote in list(edge_votes.items())[:14]:
        print(f"  {edge} by {vote[0]}:{vote[1]}")
    print()


def root_curve_vertices() -> list[Vertex]:
    return [
        Vertex(
            "root_curve_zero_locus",
            dict(predicate=3, curve=3, zeros=3, cumulants=2, extremality=3, transfer=2, computable=2, proof_exit=3),
            "whole fugacity curve and zero geometry",
            "direct interval topology unless packet sidecars are attached",
            "Lee-Yang zero-free region and discriminant boundary",
        ),
        Vertex(
            "Lee_Yang_zero_free_region",
            dict(predicate=3, curve=2, zeros=3, cumulants=1, extremality=3, transfer=2, computable=1, proof_exit=3),
            "coverage extremality as absence of dangerous roots near the LRC evaluation",
            "coefficient anatomy away from the certified region",
            "interval-certified zero exclusion for bounded banks",
        ),
        Vertex(
            "PGF_discriminant_break",
            dict(predicate=2, curve=2, zeros=3, cumulants=1, extremality=3, transfer=2, computable=3, proof_exit=2),
            "root-collision transition between extremal regimes",
            "why a collision happens unless paired with sidecars",
            "track k=8/k=9 break as discriminant sign crossing",
        ),
        Vertex(
            "miss_count_PGF_root_stratum",
            dict(predicate=2, curve=2, zeros=3, cumulants=1, extremality=2, transfer=2, computable=3, proof_exit=2),
            "sector correlation and real-root count stratum",
            "full coefficient source without q_t ledger",
            "join #real roots with exchange traps and gK8 bins",
        ),
        Vertex(
            "phi4_quartic_cumulant_stabilizer",
            dict(predicate=2, curve=2, zeros=2, cumulants=3, extremality=2, transfer=1, computable=2, proof_exit=2),
            "fourth-moment stabilization and tail penalty",
            "exact packet family unless tied to row labels",
            "measure quartic cumulant against boundary_mass_charge",
        ),
        Vertex(
            "tournament_Iomega_root_spectrum",
            dict(predicate=1, curve=3, zeros=3, cumulants=2, extremality=2, transfer=3, computable=2, proof_exit=2),
            "parallel tournament hard-core gas root spectrum",
            "LRC sector PGF semantics",
            "compare zero arcs of G_N and I(Omega,x)",
        ),
        Vertex(
            "full_PGF_coefficients",
            dict(predicate=2, curve=2, zeros=1, cumulants=2, extremality=2, transfer=1, computable=3, proof_exit=1),
            "all q_t bins and moments",
            "root geometry and phase-boundary signal",
            "compute roots, log-concavity, and nearest-zero distance",
        ),
        Vertex(
            "fugacity_rank_curve",
            dict(predicate=2, curve=3, zeros=1, cumulants=1, extremality=2, transfer=2, computable=2, proof_exit=1),
            "rank changes as z moves from LRC to tournament-like fugacity",
            "root mechanism unless zeros are attached",
            "rank-crossing table for representative packet families",
        ),
        Vertex(
            "single_value_p0",
            dict(predicate=2, curve=0, zeros=0, cumulants=0, extremality=2, transfer=0, computable=3, proof_exit=0),
            "coverage value at z=0",
            "all curve, root, and collision data",
            "use only with root-stratum or coefficient sidecar",
        ),
        Vertex(
            "raw_scalar_rank",
            dict(predicate=1, curve=0, zeros=0, cumulants=0, extremality=1, transfer=0, computable=3, proof_exit=0),
            "cheap ordering",
            "mechanism, sidecar legality, and proof exits",
            "demote to scout unless quotient is certified",
        ),
    ]


def memory_lattice_ear_vertices() -> list[Vertex]:
    return [
        Vertex(
            "packet_sheaf_legal_exit",
            dict(predicate=3, compression=3, lattice=2, ears=2, parity=2, quotient=3, computable=2, proof_exit=3),
            "the LRC predicate with named reconstruction, annihilation, descent, or debt",
            "nothing essential if all sidecars are populated",
            "use as the target interface for root/lattice/ear signals",
        ),
        Vertex(
            "first_obstruction_cocycle",
            dict(predicate=3, compression=2, lattice=2, ears=2, parity=2, quotient=3, computable=2, proof_exit=3),
            "the first hidden payload over a visible quotient fiber",
            "direct value if the obstruction is only diagnostic",
            "generate or kill each root/lattice/ear syndrome",
        ),
        Vertex(
            "Savitch_midpoint_certificate",
            dict(predicate=2, compression=3, lattice=1, ears=2, parity=1, quotient=3, computable=2, proof_exit=2),
            "recursive reachability/compression with named midpoint sidecars",
            "metric detail inside each subproblem",
            "bounded-depth witness recursion for observer charts",
        ),
        Vertex(
            "Bravais_shape_tensor",
            dict(predicate=2, compression=2, lattice=3, ears=1, parity=1, quotient=3, computable=2, proof_exit=2),
            "relation-lattice shape beyond shortest vector or covolume",
            "ear/order certificates unless attached",
            "record Gram shape and automorphism chamber of relation lattice",
        ),
        Vertex(
            "successive_minima_covolume",
            dict(predicate=2, compression=2, lattice=3, ears=0, parity=0, quotient=2, computable=3, proof_exit=2),
            "Minkowski-style volume pressure with rank/covolume/minima",
            "root curve and parity ear legality",
            "tie HYP-3062 fields to extremal rows",
        ),
        Vertex(
            "relation_lattice_basis",
            dict(predicate=2, compression=2, lattice=3, ears=0, parity=1, quotient=2, computable=3, proof_exit=1),
            "integer relations, low-height walls, signed tails",
            "shape class and proof-certificate order",
            "add Bravais-shape and wall-deletion fields",
        ),
        Vertex(
            "strong_ear_spine",
            dict(predicate=2, compression=2, lattice=0, ears=3, parity=1, quotient=2, computable=3, proof_exit=2),
            "strong-connectivity certificates by labelled ears",
            "lattice height and root curve unless coupled",
            "use ears as certificates for proof-obligation digraphs",
        ),
        Vertex(
            "odd_ear_parity_spine",
            dict(predicate=2, compression=2, lattice=0, ears=3, parity=3, quotient=2, computable=2, proof_exit=2),
            "factor-critical style odd parity certificate",
            "metric relation-lattice geometry",
            "test parity of residual proof-obligation ears",
        ),
        Vertex(
            "nested_ear_series_parallel_spine",
            dict(predicate=2, compression=2, lattice=0, ears=3, parity=2, quotient=2, computable=2, proof_exit=2),
            "nested reducibility and no-crossing proof grammar",
            "non-nested obstruction payload",
            "separate nested-refinement O2 from K33 cross-handoff O3",
        ),
        Vertex(
            "raw_runner_vertices",
            dict(predicate=1, compression=0, lattice=0, ears=0, parity=0, quotient=0, computable=3, proof_exit=0),
            "labelled speeds before quotient",
            "proof obligation structure when used alone",
            "only use after a preserved predicate is declared",
        ),
    ]


def print_alignment() -> None:
    print("=" * 100)
    print("TWO-MAP ALIGNMENT TABLE")
    print("=" * 100)
    rows: Iterable[tuple[str, str, str]] = [
        (
            "single_value_p0 -> root_curve_zero_locus",
            "raw_runner_vertices -> packet_sheaf_legal_exit",
            "the single value is not the object; the object is the curve plus the legal-forgetting ledger",
        ),
        (
            "PGF_discriminant_break",
            "strong_ear_spine / nested_ear_series_parallel_spine",
            "a root collision should name the structural ear where the proof graph changes",
        ),
        (
            "Lee_Yang_zero_free_region",
            "Savitch_midpoint_certificate",
            "zero exclusion can be certified recursively if each midpoint keeps the missing sidecar",
        ),
        (
            "phi4_quartic_cumulant_stabilizer",
            "Bravais_shape_tensor / successive_minima_covolume",
            "quartic stabilization is a fourth-moment face of the relation-lattice shape, not only a scalar tail",
        ),
        (
            "tournament_Iomega_root_spectrum",
            "odd_ear_parity_spine",
            "tournament root asymmetry and odd-ear parity both measure how cycle-packing avoids factorization",
        ),
        (
            "miss_count_PGF_root_stratum",
            "first_obstruction_cocycle",
            "near-real roots are first-obstruction alarms: attach the payload before promoting the quotient",
        ),
    ]
    for left, right, readout in rows:
        print(f"- {left}")
        print(f"  paired with: {right}")
        print(f"  readout: {readout}")
    print()
    print("new signals:")
    signals = [
        "PGF_zero_locus_signature",
        "nearest_zero_to_LRC_evaluation",
        "Lee_Yang_confinement_margin",
        "PGF_discriminant_break_index",
        "quartic_cumulant_stabilizer",
        "root_curve_vs_Iomega_arc_distance",
        "Savitch_midpoint_sidecar_depth",
        "Bravais_relation_shape_class",
        "successive_minima_anisotropy",
        "ear_certificate_type",
        "odd_ear_parity_debt",
        "nested_ear_crossing_defect",
        "root_lattice_ear_resonance_portfolio",
    ]
    for signal in signals:
        print(f"  - {signal}")
    print()
    print("assumption challenge:")
    print("  considered vertices: runners, gaps, sections, wall events, residues, cover arcs,")
    print("  Fourier modes, PGF roots, relation lattices, Bravais cells, ears, obstruction")
    print("  syndromes, and proof obligations.")
    print("  selected vertices: signal families / proof carriers.  The quotient preserves")
    print("  the named proof predicate only after the destroyed coordinate is recorded.")


def main() -> None:
    print_map(
        "MAP A: ROOT-CURVE EXTREMALITY TOURNAMENT",
        root_curve_vertices(),
        ["predicate", "curve", "zeros", "cumulants", "extremality", "transfer", "computable", "proof_exit"],
    )
    print_map(
        "MAP B: MEMORY-LATTICE-EAR CERTIFICATE TOURNAMENT",
        memory_lattice_ear_vertices(),
        ["predicate", "compression", "lattice", "ears", "parity", "quotient", "computable", "proof_exit"],
    )
    print_alignment()


if __name__ == "__main__":
    main()
