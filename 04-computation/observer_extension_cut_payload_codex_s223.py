#!/usr/bin/env python3
"""S223: observer-extension/cut payload synthesis.

This is a follow-on to S211--S218.  The finite question is the user's
``48 + 12 = 56`` observation around the first A000568/rooted-perspective
failure.  The arithmetic as written is false, so this script audits the nearby
objects exactly:

    P(5)=48, U(6)=56, defect=8, and 56 = 48 + 12 - 4.

The goal is not to promote that last identity as a theorem by itself.  The goal
is to identify which ``12`` layers are real payload carriers and which ``4``
is the previous parent-layer correction when observer-extension data is
blended with a rooted cache.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations

import a000568_edge_perspective_extension_codex_s213 as s213
import tournament_rigidity_cascade_s589 as tour


def converse(mask: int, n: int) -> int:
    return mask ^ ((1 << (n * (n - 1) // 2)) - 1)


def is_self_converse(mask: int, n: int) -> bool:
    return tour.canonical(converse(mask, n), n) == tour.canonical(mask, n)


def rooted_count(n: int) -> int:
    return sum(len(tour.vertex_orbits(rep, n)) for rep in tour.classes(n))


def rooted_count_from_reps(reps: tuple[int, ...], n: int) -> int:
    return sum(len(tour.vertex_orbits(rep, n)) for rep in reps)


def extend_with_word(parent: int, n: int, word: int) -> int:
    """Add a new vertex n.  Bit i=1 means the new vertex beats old i."""
    out = 0
    for i, j in combinations(range(n), 2):
        out = tour.set_edge(out, n + 1, i, j, tour.edge(parent, n, i, j))
    for i in range(n):
        out = tour.set_edge(out, n + 1, n, i, bool((word >> i) & 1))
    return out


def act_word(word: int, n: int, aut: tuple[int, ...]) -> int:
    """Aut action compatible with tour.relabel(old_for_new)."""
    out = 0
    for i in range(n):
        if (word >> aut[i]) & 1:
            out |= 1 << i
    return out


def word_orbit_reps(parent: int, n: int) -> tuple[int, ...]:
    auts = tour.automorphisms(parent, n)
    unseen = set(range(1 << n))
    reps = []
    while unseen:
        word = min(unseen)
        orb = {act_word(word, n, aut) for aut in auts}
        reps.append(min(orb))
        unseen -= orb
    return tuple(reps)


def deletion_profile(mask: int, n: int) -> tuple[int, ...]:
    return tuple(
        sorted(tour.canonical(tour.delete_vertex(mask, n, v), n - 1) for v in range(n))
    )


def directed_three_cycles(mask: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        total += int(
            (tour.edge(mask, n, a, b) and tour.edge(mask, n, b, c) and tour.edge(mask, n, c, a))
            or (
                tour.edge(mask, n, a, c)
                and tour.edge(mask, n, c, b)
                and tour.edge(mask, n, b, a)
            )
        )
    return total


def local_rectangle_redundancy(k: int) -> int:
    """Cycle-space rank of K_{k,k+1}; S217's local rectangle count."""
    return k * (k - 1)


def full_hourglass_redundancy(n: int) -> int:
    """S217 global adjacent-layer redundancy."""
    # 2*C(n-1,3)+C(n-2,2); written without math.comb for old Python parity.
    c1 = (n - 1) * (n - 2) * (n - 3) // 6
    c2 = (n - 2) * (n - 3) // 2
    return 2 * c1 + c2


def frac_text(num: int, den: int) -> str:
    f = Fraction(num, den)
    return f"{f.numerator}/{f.denominator}"


def extension_surface(parent_n: int = 5) -> dict[str, object]:
    parents = tour.classes(parent_n)
    child_n = parent_n + 1
    child_reps = tuple(
        sorted(
            {
                tour.canonical(extend_with_word(parent, parent_n, word), child_n)
                for parent in parents
                for word in range(1 << parent_n)
            }
        )
    )
    source_word = (1 << parent_n) - 1
    sink_word = 0

    label = {rep: i + 1 for i, rep in enumerate(parents)}
    child_label = {rep: i + 1 for i, rep in enumerate(child_reps)}

    source_children: set[int] = set()
    sink_children: set[int] = set()
    targets_by_parent: dict[int, set[int]] = defaultdict(set)
    orbit_targets_by_parent: dict[int, set[int]] = defaultdict(set)
    word_orbit_hist: Counter[int] = Counter()
    parent_aut_hist: Counter[int] = Counter()
    parent_rows = []
    self_converse_children_by_parent: dict[int, set[int]] = defaultdict(set)

    for parent in parents:
        aut_size = len(tour.automorphisms(parent, parent_n))
        parent_aut_hist[aut_size] += 1
        wreps = word_orbit_reps(parent, parent_n)
        word_orbit_hist[len(wreps)] += 1

        for word in range(1 << parent_n):
            child = tour.canonical(extend_with_word(parent, parent_n, word), child_n)
            targets_by_parent[parent].add(child)
            if is_self_converse(child, child_n):
                self_converse_children_by_parent[parent].add(child)

        for word in wreps:
            child = tour.canonical(extend_with_word(parent, parent_n, word), child_n)
            orbit_targets_by_parent[parent].add(child)

        source = tour.canonical(extend_with_word(parent, parent_n, source_word), child_n)
        sink = tour.canonical(extend_with_word(parent, parent_n, sink_word), child_n)
        source_children.add(source)
        sink_children.add(sink)
        parent_rows.append(
            {
                "parent_id": label[parent],
                "mask": parent,
                "aut": aut_size,
                "word_orbits": len(wreps),
                "raw_targets": len(targets_by_parent[parent]),
                "orbit_targets": len(orbit_targets_by_parent[parent]),
                "source_child": child_label[source],
                "sink_child": child_label[sink],
                "source_self_converse": is_self_converse(source, child_n),
                "sink_self_converse": is_self_converse(sink, child_n),
                "self_converse_targets": len(self_converse_children_by_parent[parent]),
            }
        )

    sc6 = {rep for rep in child_reps if is_self_converse(rep, child_n)}
    deletion_decks = {rep: deletion_profile(rep, child_n) for rep in child_reps}
    sc_parent_union = {p for rep in sc6 for p in set(deletion_decks[rep])}
    sc_deck_hist = Counter(len(set(deletion_decks[rep])) for rep in sc6)
    all_deck_hist = Counter(len(set(deletion_decks[rep])) for rep in child_reps)
    deck_to_children: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for rep, deck in deletion_decks.items():
        deck_to_children[deck].append(rep)

    return {
        "parents": parents,
        "children": child_reps,
        "label": label,
        "child_label": child_label,
        "source_children": source_children,
        "sink_children": sink_children,
        "source_sink_overlap": source_children & sink_children,
        "word_orbit_hist": dict(sorted(word_orbit_hist.items())),
        "parent_aut_hist": dict(sorted(parent_aut_hist.items())),
        "parent_rows": parent_rows,
        "targets_by_parent_hist": dict(
            sorted(Counter(len(v) for v in targets_by_parent.values()).items())
        ),
        "orbit_targets_by_parent_hist": dict(
            sorted(Counter(len(v) for v in orbit_targets_by_parent.values()).items())
        ),
        "self_converse_6": sc6,
        "sc_source_overlap": sc6 & source_children,
        "sc_sink_overlap": sc6 & sink_children,
        "sc_parent_union": sc_parent_union,
        "sc_deck_hist": dict(sorted(sc_deck_hist.items())),
        "all_deck_hist": dict(sorted(all_deck_hist.items())),
        "unique_decks": len(deck_to_children),
        "deck_collision_hist": dict(sorted(Counter(len(v) for v in deck_to_children.values()).items())),
        "deck_collisions": [tuple(v) for v in deck_to_children.values() if len(v) > 1],
    }


APPLICATIONS = [
    (
        "A000568/rooted perspective",
        "rooted node type",
        "new observer incident word and cross-sector chirality",
        "ordered-pair sector deck / edge tail-tip sector word",
        "exact finite first-failure audit",
    ),
    (
        "LRC threshold packet",
        "safe/open status near a boundary phase",
        "endpoint owner, period deck, route label",
        "labelled packet sheaf with hidden-coordinate ledger",
        "HYP-3039/HYP-3043 controlled-forgetting protocol",
    ),
    (
        "Residual capacitor cuts",
        "mixed route pair survives coarse boundary type",
        "which generator/cut class blocks it",
        "finite min-cut payload plus owner strip",
        "pair-good decoy repair and capacitor first-cut IDs",
    ),
    (
        "Haar drop/add square",
        "row/column margin shadow",
        "mixed zeta T00-T01-T10+T11",
        "rectangle product residue",
        "exact zeta sign or dual annihilation",
    ),
    (
        "S217 fixed-path diagonal flow",
        "line-count or half-tiling shadow",
        "rectangle/hourglass cycle residue",
        "GF(2) cut-space basis on K_{k,k+1}",
        "sidecar admissibility; not a scalar count",
    ),
    (
        "Matrix sidecar observability",
        "rank/spectrum/scalar invariant",
        "hidden coordinate separating route/status fibers",
        "observability matrix column",
        "separate, reconstruct, annihilate, descend, or name debt",
    ),
    (
        "Moser/fibbinary automata",
        "regular-language membership",
        "scale, parity-position, carry boundary",
        "automaton state plus valuation/topology sidecar",
        "automatic telemetry before LRC handoff",
    ),
    (
        "Perfect-number/divisor lane",
        "abundancy or Euler-product scalar",
        "prime-power lane and divisor-lattice fiber",
        "divisor-lattice/SNF payload",
        "squarefree blindness repair",
    ),
    (
        "Fixed-path tournament atlas",
        "Hamiltonian-path presentation",
        "fiber H(T)/Aut(T) and path-reversal branch data",
        "presentation-fiber payload",
        "unrooted class quotient only after fiber audit",
    ),
    (
        "Proof-obligation automaton",
        "coarse theorem stage",
        "which obligation remains after a quotient",
        "named residual state plus legal exit",
        "agent coordination and next lemma target",
    ),
]


def application_tournament() -> tuple[list[tuple[str, dict[str, int]]], list[list[int]]]:
    carriers = [
        ("raw_scalar_or_count", dict(payload=0, exact=0, lrc=0, automaton=0, cost=0)),
        ("rooted_node_cache", dict(payload=1, exact=2, lrc=1, automaton=0, cost=1)),
        ("incident_word_orbit", dict(payload=2, exact=2, lrc=1, automaton=1, cost=2)),
        ("ordered_pair_sector_deck", dict(payload=3, exact=3, lrc=2, automaton=1, cost=3)),
        ("rectangle_hourglass_residue", dict(payload=3, exact=2, lrc=3, automaton=1, cost=3)),
        ("deletion_fiber_payload", dict(payload=3, exact=3, lrc=2, automaton=1, cost=3)),
        ("endpoint_owner_packet", dict(payload=4, exact=2, lrc=4, automaton=1, cost=4)),
        ("sidecar_observability_matrix", dict(payload=4, exact=3, lrc=4, automaton=1, cost=5)),
        ("proof_obligation_automaton", dict(payload=5, exact=2, lrc=5, automaton=3, cost=6)),
    ]

    def score(data: dict[str, int]) -> tuple[int, int]:
        retained = (
            3 * data["payload"]
            + 2 * data["exact"]
            + 3 * data["lrc"]
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


def reps_for_n(surface: dict[str, object], n: int) -> tuple[int, ...]:
    if n == 6:
        return surface["children"]  # type: ignore[return-value]
    return tour.classes(n)


def rooted_count_fast(surface: dict[str, object], n: int) -> int:
    reps = reps_for_n(surface, n)
    return rooted_count_from_reps(reps, n)


def print_count_tables(surface: dict[str, object]) -> None:
    print("1. COUNT ORIGINS AND THE 12-LAYER")
    print("   U(n)=unrooted tournament classes, P(n)=rooted/node perspectives, SC(n)=self-converse classes.")
    print("   n  U(n)  P(n)  SC(n)  P(n-1)->U(n)  local_K_{n-1,n}_rect")
    for n in range(1, 7):
        reps = reps_for_n(surface, n)
        u = len(reps)
        p = rooted_count_from_reps(reps, n)
        sc = sum(1 for rep in reps if is_self_converse(rep, n))
        shifted = "-" if n == 1 else str(u - rooted_count_fast(surface, n - 1))
        rect = "-" if n == 1 else str(local_rectangle_redundancy(n - 1))
        print(f"   {n:1d}  {u:4d}  {p:4d}  {sc:5d}  {shifted:12s}  {rect:20s}")

    print()
    print("   Exact 12 carriers near the first defect:")
    print(f"   P(4) rooted perspectives                         = {rooted_count_fast(surface, 4)}")
    print(f"   U(5) unrooted classes                            = {len(reps_for_n(surface, 5))}")
    print(
        "   source-word slice in 5->6 extension              = "
        f"{len(surface['source_children'])}"
    )
    print(
        "   sink-word slice in 5->6 extension                = "
        f"{len(surface['sink_children'])}"
    )
    print(
        "   self-converse classes at n=6                     = "
        f"{len(surface['self_converse_6'])}"
    )
    print(f"   K_{{4,5}} rectangle cycle-space redundancy        = {local_rectangle_redundancy(4)}")


def print_arithmetic(surface: dict[str, object]) -> None:
    p5 = rooted_count_fast(surface, 5)
    u6 = len(reps_for_n(surface, 6))
    u5 = len(reps_for_n(surface, 5))
    u4 = len(reps_for_n(surface, 4))
    print()
    print("2. ARITHMETIC CORRECTION")
    print(f"   48 + 12 is {p5 + u5}, not {u6}.")
    print(f"   The first defect is U(6)-P(5) = {u6}-{p5} = {u6-p5}.")
    print(f"   The exact splice is U(6) = P(5) + U(5) - U(4) = {p5}+{u5}-{u4} = {u6}.")
    print("   This identity is a useful alarm, not yet a theorem: the 4 is the previous")
    print("   parent/source layer that must be quotient-corrected when a source/sink")
    print("   slice is merged into the arbitrary rooted-perspective cache.")
    print(
        "   Fractions of U(6): "
        f"P(5)/U(6)={frac_text(p5, u6)}, "
        f"defect/U(6)={frac_text(u6-p5, u6)}, "
        f"U(5)/U(6)={frac_text(u5, u6)}, "
        f"U(4)/U(6)={frac_text(u4, u6)}."
    )
    print(
        "   Hence the correction reads 1 = 6/7 + 3/14 - 1/14; "
        "the net missing mass is 1/7."
    )


def print_extension_fibers(surface: dict[str, object]) -> None:
    print()
    print("3. INCIDENT WORDS, SOURCE/SINK SLICES, AND DELETION FIBERS")
    print(f"   parent Aut-size histogram at U(5): {surface['parent_aut_hist']}")
    print(f"   Aut(parent)-word-orbit count histogram: {surface['word_orbit_hist']}")
    print(f"   raw target classes per U(5) parent: {surface['targets_by_parent_hist']}")
    print(f"   orbit-target classes per U(5) parent: {surface['orbit_targets_by_parent_hist']}")
    print(
        "   source/sink child sets: "
        f"|source|={len(surface['source_children'])}, |sink|={len(surface['sink_children'])}, "
        f"overlap={len(surface['source_sink_overlap'])}"
    )
    print(
        "   SC(6) intersections: "
        f"with source={len(surface['sc_source_overlap'])}, "
        f"with sink={len(surface['sc_sink_overlap'])}, "
        f"with source_or_sink={len(surface['self_converse_6'] & (surface['source_children'] | surface['sink_children']))}"
    )
    print(
        "   deletion parent diversity over all U(6) classes: "
        f"{surface['all_deck_hist']}"
    )
    print(
        "   deletion parent diversity over SC(6) classes: "
        f"{surface['sc_deck_hist']}; parent classes touched={len(surface['sc_parent_union'])}/12"
    )
    print(
        "   deletion deck collision histogram: "
        f"{surface['deck_collision_hist']} (unique decks={surface['unique_decks']}/56)"
    )

    print()
    print("   Parent rows (id aut word_orbits raw_targets orbit_targets source sink SC-targets):")
    for row in surface["parent_rows"]:  # type: ignore[index]
        print(
            f"   {row['parent_id']:2d} {row['aut']:3d} {row['word_orbits']:11d}"
            f" {row['raw_targets']:11d} {row['orbit_targets']:13d}"
            f" {row['source_child']:6d} {row['sink_child']:4d}"
            f" {row['self_converse_targets']:10d}"
        )


def print_sector_and_s217(surface: dict[str, object]) -> None:
    print()
    print("4. ORDERED-PAIR/EDGE SECTORS AND S217 RESIDUES")
    classes6 = reps_for_n(surface, 6)
    for mode in ("size", "internal", "cross", "full"):
        decks = {s213.class_sector_deck(rep, 6, mode) for rep in classes6}
        print(f"   ordered-pair sector deck mode={mode:8s}: separates {len(decks)}/56 classes")
    sc_c3 = Counter(directed_three_cycles(rep, 6) for rep in surface["self_converse_6"])  # type: ignore[arg-type]
    all_c3 = Counter(directed_three_cycles(rep, 6) for rep in classes6)
    print(f"   C3 histogram over all U(6): {dict(sorted(all_c3.items()))}")
    print(f"   C3 histogram over SC(6):   {dict(sorted(sc_c3.items()))}")
    print(
        "   S217 local line law: K_{k,k+1} has k(k+1) lines, rank 2k, "
        "rectangle residue rank k(k-1)."
    )
    print(
        "   For k=4 this residue rank is 12; for n=6 global rectangle/hourglass "
        f"redundancy is {full_hourglass_redundancy(6)}."
    )
    print(
        "   Reading: ordered-pair cross sectors and S217 rectangle/hourglass cycles "
        "are the same kind of cut payload: a mixed coordinate killed by a scalar "
        "shadow but recoverable as a sidecar."
    )


def print_payload_abstraction() -> None:
    print()
    print("5. OBSERVER-EXTENSION / CUT PAYLOAD ABSTRACT")
    print("   Data form:")
    print("     base quotient Q, observer/cut word sigma, stabilizer Aut(Q),")
    print("     sidecar C(sigma), sink map Phi(Q,sigma), and legality exit.")
    print("   Controlled-forgetting recurrence:")
    print("     retained_{r+1} = (retained_r, C_r) / Aut(retained_r)")
    print("     debt_r = kernel(C_r -> LRC predicate)")
    print("     allow forgetting only if debt_r is reconstructed, annihilated,")
    print("     descended to a smaller family, AP/GW equality, or named residual debt.")
    print()
    print("   Applications across this repo:")
    print("   problem                         retained shadow                    hidden payload")
    for name, retained, hidden, sidecar, exit_rule in APPLICATIONS:
        print(f"   {name:30s} {retained:34s} {hidden}")
        print(f"     sidecar={sidecar}; exit={exit_rule}")


def print_tournament_analysis() -> None:
    carriers, adj = application_tournament()
    out_scores = [sum(row) for row in adj]
    order = sorted(range(len(carriers)), key=lambda i: -out_scores[i])
    print()
    print("6. TOURNAMENT ANALYSIS OVER PAYLOAD CARRIERS")
    print("   Vertices are proof carriers, not runners or tournament nodes.")
    print("   Observable: retained cut payload, exact finite checkability, LRC predicate")
    print("   relevance, and automaton/handoff value; proof cost breaks ties.")
    print(f"   vertices={[name for name, _ in carriers]}")
    print(f"   score_hist={dict(sorted(Counter(out_scores).items()))}")
    print(f"   directed_3_cycles={tour.directed_three_cycles(adj)}")
    print(f"   scc_sizes={tour.scc_sizes(adj)}")
    print(f"   hamiltonian_paths={tour.ham_paths_adj(adj)}")
    print("   one_hamiltonian_path=" + " -> ".join(carriers[i][0] for i in order))


def print_reading() -> None:
    print()
    print("READING")
    print(
        "  The user's 12 observation is real, but not as 48+12=56.  The exact "
        "nearby arithmetic is P(5)=48, U(6)=56, defect=8, with the splice "
        "56=48+12-4 and fractions 1=6/7+3/14-1/14."
    )
    print(
        "  The recurring 12 layers are different facets of the same payload theme: "
        "P(4), U(5), source/sink extension slices in 5->6, SC(6), and the k=4 "
        "rectangle-cycle rank in S217.  Each is a controlled-forgetting boundary "
        "where a quotient is legal only after the cut coordinate is named."
    )
    print(
        "  The proposed carrier is observer-extension/cut payload: a rooted or "
        "coarse object plus the incident word, sector/cross-sector coupling, "
        "deletion fiber, rectangle/hourglass residue, endpoint owner, or proof "
        "obligation that was destroyed by the quotient."
    )
    print(
        "  LRC use: do not treat A000568, a half-tiling count, a matrix scalar, "
        "or an automaton state as the proof object unless its observer/cut payload "
        "has been retained, reconstructed, annihilated, descended, or routed to "
        "named debt."
    )
    print()
    print("ASSUMPTION CHALLENGE")
    print(
        "  Considered vertices: runners, rooted nodes, directed edges, ordered "
        "pairs, incident words, deletion fibers, source/sink slices, self-converse "
        "branch loci, rectangle/hourglass cycles, endpoint owners, automaton "
        "states, divisor lanes, matrix columns, and proof obligations."
    )
    print(
        "  Selected quotient: proof carriers with explicit observer/cut payload.  "
        "Preserved predicate: LRC safe-box route/status information and tournament "
        "extension recoverability.  Destroyed information: labelled runner identity, "
        "raw word order, full half-tiling presentation, and scalar-only counts."
    )


def main() -> None:
    surface = extension_surface(5)
    print("=" * 80)
    print("S223: observer-extension/cut payload and the 12-layer audit")
    print("=" * 80)
    print_count_tables(surface)
    print_arithmetic(surface)
    print_extension_fibers(surface)
    print_sector_and_s217(surface)
    print_payload_abstraction()
    print_tournament_analysis()
    print_reading()


if __name__ == "__main__":
    main()
