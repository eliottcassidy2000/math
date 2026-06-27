#!/usr/bin/env python3
"""S216: diagonal-layer transport and half-tiling quotient laws.

The prompt asks what happens when a tournament grows from n to n+1 in the
fixed-path tiling model.  In that model the new vertex appends one diagonal
word of length n.  Between consecutive geometric layers of lengths k and k+1
we can put the complete bipartite position carrier K_{k,k+1}, with k(k+1)
position-lines.

This script keeps two objects separate:

* the geometric position carrier K_{k,k+1}, whose vertices are tile positions;
* the binary word layers {0,1}^k and {0,1}^{k+1}, whose labels are the actual
  orientations.

That separation exposes the algebraic efficiency: K_{k,k+1} has k(k+1) lines,
but its binary pair-label matrix is an outer product of two diagonal words.
The line counts obey rank-one laws and the aligned diagonal plus the link bit
updates the triangles that use the two newest vertices.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import comb

from tournament_rigidity_cascade_s589 import (
    automorphisms,
    canonical,
    classes,
    delete_vertex,
    directed_three_cycles,
    edge,
    ham_paths_adj,
    scc_sizes,
    set_edge,
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


def all_words(length: int) -> tuple[tuple[int, ...], ...]:
    return tuple(product((0, 1), repeat=length))


EXTENSION_CACHE: dict[int, dict[str, object]] = {}


def lift_mask(mask: int, n: int) -> int:
    """Copy an n-vertex mask into the first n vertices of an (n+1)-mask."""
    out = 0
    for i, j in combinations(range(n), 2):
        if edge(mask, n, i, j):
            out = set_edge(out, n + 1, i, j, True)
    return out


def extend_by_diagonal_word(mask: int, n: int, word: tuple[int, ...]) -> int:
    """Append vertex n.  word[i]=1 means old vertex i beats the new vertex."""
    if len(word) != n:
        raise ValueError("word length must match parent size")
    out = lift_mask(mask, n)
    for i, bit in enumerate(word):
        out = set_edge(out, n + 1, i, n, bool(bit))
    return out


def diagonal_word(mask: int, n: int, vertex: int) -> tuple[int, ...]:
    """Return the fixed-path layer word for vertex, against all earlier vertices."""
    return tuple(int(edge(mask, n, i, vertex)) for i in range(vertex))


def act_on_word(word: tuple[int, ...], aut: tuple[int, ...]) -> tuple[int, ...]:
    """Action induced by the parent automorphism tuple used by relabel()."""
    return tuple(word[aut[i]] for i in range(len(word)))


def word_orbits(parent: int, n: int) -> list[tuple[tuple[int, ...], set[tuple[int, ...]]]]:
    auts = automorphisms(parent, n)
    unseen = set(all_words(n))
    out: list[tuple[tuple[int, ...], set[tuple[int, ...]]]] = []
    while unseen:
        word = min(unseen)
        orbit = {act_on_word(word, aut) for aut in auts}
        rep = min(orbit)
        out.append((rep, orbit))
        unseen -= orbit
    return out


@lru_cache(maxsize=None)
def u_count(n: int) -> int:
    if n in KNOWN_U:
        return KNOWN_U[n]
    return len(classes(n))


def rooted_count_from_reps(reps: tuple[int, ...], n: int) -> int:
    return sum(len(vertex_orbits(rep, n)) for rep in reps)


def rooted_count(n: int) -> int:
    return rooted_count_from_reps(classes(n), n)


def class_reps_fast(n: int) -> tuple[int, ...]:
    if n <= 5:
        return classes(n)
    if n == 6:
        return tuple(extension_census(5)["child_reps"])  # type: ignore[index]
    return classes(n)


def self_converse_count_from_reps(reps: tuple[int, ...], n: int) -> int:
    return sum(1 for rep in reps if canonical(converse(rep, n), n) == rep)


def deletion_profile(mask: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(canonical(delete_vertex(mask, n, v), n - 1) for v in range(n)))


def converse(mask: int, n: int) -> int:
    return mask ^ ((1 << (n * (n - 1) // 2)) - 1)


def self_converse_count(n: int) -> int:
    return sum(1 for rep in classes(n) if canonical(converse(rep, n), n) == rep)


def half_bits(n: int) -> int:
    return ((n - 1) ** 2) // 4


def full_bits(n: int) -> int:
    return n * (n - 1) // 2


def extension_census(parent_n: int) -> dict[str, object]:
    if parent_n in EXTENSION_CACHE:
        return EXTENSION_CACHE[parent_n]

    parents = classes(parent_n)
    raw_child_fibers: Counter[int] = Counter()
    orbit_child_fibers: Counter[int] = Counter()
    orbit_size_hist: Counter[int] = Counter()
    parent_orbit_hist: Counter[int] = Counter()
    source_children: set[int] = set()
    sink_children: set[int] = set()

    for parent in parents:
        aut_size = len(automorphisms(parent, parent_n))
        parent_orbit_hist[aut_size] += 1

        for word in all_words(parent_n):
            child = canonical(extend_by_diagonal_word(parent, parent_n, word), parent_n + 1)
            raw_child_fibers[child] += 1

        for rep_word, orbit in word_orbits(parent, parent_n):
            orbit_size_hist[len(orbit)] += 1
            child = canonical(extend_by_diagonal_word(parent, parent_n, rep_word), parent_n + 1)
            orbit_child_fibers[child] += 1

        source_children.add(canonical(extend_by_diagonal_word(parent, parent_n, (0,) * parent_n), parent_n + 1))
        sink_children.add(canonical(extend_by_diagonal_word(parent, parent_n, (1,) * parent_n), parent_n + 1))

    child_reps = tuple(sorted(orbit_child_fibers))
    child_n = parent_n + 1
    rooted = rooted_count_from_reps(child_reps, child_n)
    vertex_orbit_match = all(
        orbit_child_fibers[rep] == len(vertex_orbits(rep, child_n)) for rep in child_reps
    )
    unique_parent_hist = Counter(len(set(deletion_profile(rep, child_n))) for rep in child_reps)
    repeated_parent_hist = Counter(
        max(Counter(deletion_profile(rep, child_n)).values()) for rep in child_reps
    )

    out = {
        "child_reps": child_reps,
        "parent_n": parent_n,
        "child_n": child_n,
        "parent_classes": len(parents),
        "raw_extensions": len(parents) * (1 << parent_n),
        "child_classes": len(child_reps),
        "layer_orbits": sum(orbit_child_fibers.values()),
        "rooted_count": rooted,
        "vertex_orbit_match": vertex_orbit_match,
        "raw_sink_hist": dict(sorted(Counter(raw_child_fibers.values()).items())),
        "orbit_sink_hist": dict(sorted(Counter(orbit_child_fibers.values()).items())),
        "orbit_size_hist": dict(sorted(orbit_size_hist.items())),
        "parent_aut_hist": dict(sorted(parent_orbit_hist.items())),
        "unique_parent_hist": dict(sorted(unique_parent_hist.items())),
        "repeated_parent_hist": dict(sorted(repeated_parent_hist.items())),
        "source_children": len(source_children),
        "sink_children": len(sink_children),
    }
    EXTENSION_CACHE[parent_n] = out
    return out


def two_layer_profiles(k: int) -> dict[str, object]:
    line_profiles: set[tuple[int, int, int, int]] = set()
    aligned_profiles: set[tuple[int, tuple[tuple[tuple[int, int], int], ...]]] = set()
    score_cycle_profiles: set[tuple[int, int, int, int]] = set()
    cycle_counts: set[int] = set()

    for w in all_words(k):
        for u in all_words(k + 1):
            a = sum(w)
            b = sum(u)
            # Counts of K_{k,k+1} line labels (w_i, u_j), ordered 00,01,10,11.
            line_profiles.add(((k - a) * (k + 1 - b), (k - a) * b, a * (k + 1 - b), a * b))

            link = u[-1]
            aligned = Counter((w[i], u[i]) for i in range(k))
            aligned_key = tuple(sorted(aligned.items()))
            aligned_profiles.add((link, aligned_key))

            two_newest_cycles = sum(1 for i in range(k) if w[i] == link and u[i] != link)
            cycle_counts.add(two_newest_cycles)
            score_cycle_profiles.add((a, b, link, two_newest_cycles))

    return {
        "k": k,
        "raw_two_words": 1 << (2 * k + 1),
        "position_lines": k * (k + 1),
        "stored_bits_for_words": 2 * k + 1,
        "line_profiles": len(line_profiles),
        "line_profile_formula": (k + 1) * (k + 2),
        "aligned_profiles": len(aligned_profiles),
        "aligned_profile_formula": 2 * comb(k + 3, 3),
        "score_cycle_profiles": len(score_cycle_profiles),
        "cycle_count_values": len(cycle_counts),
    }


def carrier_tournament() -> tuple[list[tuple[str, dict[str, int]]], list[list[int]]]:
    carriers = [
        ("raw_labelled_diagonal_word", dict(iso=0, aut=0, deletion=0, line=0, cycle=0, half=0, owner=0, automaton=0, cost=0)),
        ("binary_half_tiling_shadow", dict(iso=0, aut=0, deletion=0, line=0, cycle=0, half=1, owner=0, automaton=0, cost=1)),
        ("K_position_line_profile", dict(iso=0, aut=0, deletion=0, line=1, cycle=0, half=0, owner=0, automaton=0, cost=1)),
        ("parent_aut_layer_orbit", dict(iso=1, aut=1, deletion=0, line=0, cycle=0, half=0, owner=0, automaton=1, cost=2)),
        ("aligned_triangle_flow", dict(iso=0, aut=0, deletion=0, line=1, cycle=1, half=0, owner=0, automaton=1, cost=2)),
        ("deletion_parent_fiber", dict(iso=1, aut=1, deletion=1, line=0, cycle=0, half=0, owner=0, automaton=1, cost=3)),
        ("diagonal_transport_sidecar", dict(iso=1, aut=1, deletion=1, line=1, cycle=1, half=1, owner=0, automaton=1, cost=4)),
        ("endpoint_owner_packet_sheaf", dict(iso=1, aut=1, deletion=1, line=1, cycle=1, half=0, owner=1, automaton=1, cost=5)),
        ("proof_obligation_automaton", dict(iso=1, aut=1, deletion=1, line=1, cycle=1, half=1, owner=1, automaton=1, cost=6)),
    ]

    def score(data: dict[str, int]) -> tuple[int, int]:
        retained = (
            3 * data["iso"]
            + 3 * data["aut"]
            + 3 * data["deletion"]
            + 2 * data["line"]
            + 2 * data["cycle"]
            + data["half"]
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


def one_hamiltonian_path(adj: list[list[int]], names: list[str]) -> list[str]:
    n = len(names)
    path: list[int] = []
    used = [False] * n

    def dfs(v: int) -> bool:
        path.append(v)
        used[v] = True
        if len(path) == n:
            return True
        for u in range(n):
            if not used[u] and adj[v][u] and dfs(u):
                return True
        used[v] = False
        path.pop()
        return False

    for start in range(n):
        if dfs(start):
            return [names[i] for i in path]
    return []


def print_extension_census(max_parent: int = 5) -> None:
    print("1. DIAGONAL EXTENSION CENSUS")
    print("   A parent n-class plus a binary diagonal word of length n builds a child (n+1)-class.")
    print("   Word orbits under Aut(parent) are exactly rooted child classes; child sinks are unrooted classes.")
    print()
    print(
        "   n -> n+1  U(n)  raw=U(n)2^n  layer_orbits  rooted  U(n+1)  "
        "raw/rooted  rooted/U  source  sink  orbit=vertex_orbits"
    )
    for n in range(1, max_parent + 1):
        c = extension_census(n)
        raw = int(c["raw_extensions"])
        rooted = int(c["rooted_count"])
        child = int(c["child_classes"])
        print(
            f"   {n} -> {n+1:<2d}  {c['parent_classes']:4d}  {raw:12d}  "
            f"{c['layer_orbits']:12d}  {rooted:6d}  {child:6d}  "
            f"{raw / rooted:10.3f}  {rooted / child:8.3f}  "
            f"{c['source_children']:6d}  {c['sink_children']:4d}  {c['vertex_orbit_match']}"
        )

    print()
    print("   DUPLICATION HISTOGRAMS AT THE FIRST SHIFTED FAILURE SURFACE (5 -> 6)")
    c = extension_census(5)
    print(f"   parent automorphism sizes among 5-classes: {c['parent_aut_hist']}")
    print(f"   diagonal word orbit sizes under parent automorphisms: {c['orbit_size_hist']}")
    print(f"   rooted extension orbits per 6-class: {c['orbit_sink_hist']}")
    print(f"   raw labelled extensions per 6-class: {c['raw_sink_hist']}")
    print(f"   unique deletion-parent classes per 6-class: {c['unique_parent_hist']}")
    print(f"   max repeated deletion-parent multiplicity per 6-class: {c['repeated_parent_hist']}")


def print_two_layer_laws(max_k: int = 6) -> None:
    print()
    print("2. TWO-LAYER POSITION-LINE FLOW")
    print("   Consecutive geometric layers have k and k+1 tile positions.")
    print("   The complete carrier K_{k,k+1} has k(k+1) position-lines, while the binary labels are two words with 2k+1 bits.")
    print()
    print(
        "   k  K_edges  raw_word_pairs  stored_bits  line_profiles  aligned_profiles  "
        "score_cycle_profiles  cycle_values"
    )
    for k in range(1, max_k + 1):
        p = two_layer_profiles(k)
        print(
            f"   {k:1d}  {p['position_lines']:7d}  {p['raw_two_words']:14d}  "
            f"{p['stored_bits_for_words']:11d}  {p['line_profiles']:13d}  "
            f"{p['aligned_profiles']:16d}  {p['score_cycle_profiles']:20d}  "
            f"{p['cycle_count_values']:12d}"
        )

    print()
    print("   Exact laws:")
    print("   - If a=sum(w) and b=sum(u), the K-line label counts are")
    print("     N00=(k-a)(k+1-b), N01=(k-a)b, N10=a(k+1-b), N11=ab.")
    print("   - Therefore N00*N11=N01*N10: every line-count matrix has rank one.")
    print("   - The triangles using the two newest vertices are read only from aligned pairs")
    print("     (w_i,u_i) and the link bit ell=u_k:")
    print("     c3_newest_pair = #{i<k : w_i=ell and u_i!=ell}.")
    print("   - Score flow is affine: old vertex i gains w_i+u_i, vertex k gains")
    print("     k-sum(w)+u_k, and the newest vertex has outdegree k+1-sum(u).")


def print_half_tiling_checks(max_n: int = 7) -> None:
    print()
    print("3. HALF-TILING / DIAGONAL-SYMMETRY CHECK")
    print("   The half-tiling is the fixed-path converse fold.  It is a branch/fundamental-domain input, not generally an isomorphism-class count.")
    print()
    print(
        "   n  U(n)  rooted_P(n)  2^half(n-1)  2^half(n)  rho_orbits(n)  self_converse_classes"
    )
    for n in range(3, max_n + 1):
        reps = class_reps_fast(n) if n <= 6 else ()
        rooted = rooted_count_from_reps(reps, n) if n <= 6 else -1
        sc = self_converse_count_from_reps(reps, n) if n <= 6 else -1
        print(
            f"   {n:1d}  {u_count(n):4d}  {rooted:11d}  "
            f"{1 << half_bits(n - 1):12d}  {1 << half_bits(n):9d}  "
            f"{((1 << full_bits(n)) + (1 << half_bits(n))) // 2:13d}  {sc:21d}"
        )
    print()
    print("   Readout:")
    print("   - U(n) is not the count of binary half-tilings of size n-1 or n.")
    print("   - Symmetry along the fold detects the self-converse branch locus; even that class count is not U(n-1) in general.")
    print("   - The exact growth object is the diagonal extension orbit DAG: parent class + word orbit -> rooted child -> unrooted child sink.")


def print_carrier_analysis() -> None:
    print()
    print("4. TOURNAMENT ANALYSIS OVER DIAGONAL-TRANSPORT CARRIERS")
    carriers, adj = carrier_tournament()
    names = [name for name, _ in carriers]
    scores = [sum(row) for row in adj]
    print("   Vertices are proof carriers, not runners.")
    print("   Observable: retained iso/aut/deletion/line/cycle/half/owner/automaton payload minus proof cost.")
    print("   Switch: orient toward the carrier retaining more recoverable diagonal-transport payload; name order breaks exact ties.")
    print(f"   vertices={names}")
    print(f"   score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"   directed_3_cycles={directed_three_cycles(adj)}")
    print(f"   scc_sizes={scc_sizes(adj)}")
    print(f"   hamiltonian_paths={ham_paths_adj(adj)}")
    print(f"   one_hamiltonian_path={' -> '.join(one_hamiltonian_path(adj, names))}")


def print_reading() -> None:
    print()
    print("READING")
    print("  The user's k^2+k lines are best interpreted as K_{k,k+1} on geometric tile positions, not as the complete graph between all binary words.")
    print("  Once that is fixed, the redundancy is sharp: k(k+1) line labels are generated by only 2k+1 bits, and their count matrix satisfies a rank-one determinant law.")
    print("  The recursive tournament DP should therefore carry parent automorphism word-orbits, aligned triangle-flow data, and deletion-parent fibers, then merge into unrooted class sinks.")
    print("  At 5 -> 6, raw labelled diagonal extensions are 384, parent-aut word orbits are 296, and the unrooted sinks are 56.  The 296 rooted orbits match the sum of vertex-orbits over the 56 child classes.")
    print("  Source and sink words are exact deletion slices: all-zero and all-one diagonal words give U(n) source/sink-rooted children.  Arbitrary words need the sidecar.")
    print("  The half-tiling idea is still valuable, but as a converse-fold branch locus and complement-pair quotient.  It does not by itself enumerate tournament isomorphism classes.")

    print()
    print("ASSUMPTION CHALLENGE")
    print("  Considered vertices: runners, diagonal tile positions, K position-lines, binary layer words, parent-aut word orbits, deletion fibers, fixed circle sections, section boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits, endpoint owners, and proof obligations.")
    print("  Selected quotient: diagonal transport sidecars preserve observer-source, incident-word, two-newest-triangle, deletion-parent, and endpoint-owner payloads needed for LRC safe-box routing.")
    print("  Destroyed information: full labelled runner identity, full binary word order beyond the chosen sidecar, and raw scalar half-tiling counts.")
    print("  Challenged assumption: isomorphism classes should equal half-tiling counts.  The safer statement is that class growth factors through diagonal word-orbits and deletion fibers, while half-tilings mark the converse-symmetric branch of that transport.")


def main() -> None:
    print("=" * 80)
    print("S216: TOURNAMENT DIAGONAL-LAYER TRANSPORT")
    print("=" * 80)
    print_extension_census()
    print_two_layer_laws()
    print_half_tiling_checks()
    print_carrier_analysis()
    print_reading()


if __name__ == "__main__":
    main()
