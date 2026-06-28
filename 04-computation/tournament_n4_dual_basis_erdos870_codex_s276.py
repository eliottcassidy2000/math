#!/usr/bin/env python3
"""n=4 tournament class bases, Erdős-870 subbasis analogy, and LRC packet tests.

This scout is deliberately finite and exact.  It verifies the two n=4
tournament modelling schemes from the prompt, then reads them as quotient
guardrails for the LRC14 labelled-packet program.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, product


VERTICES = tuple(range(4))
EDGES = tuple(combinations(VERTICES, 2))
CLASS_BY_SCORE = {
    (0, 1, 2, 3): "T",
    (1, 1, 2, 2): "S",
    (0, 2, 2, 2): "+",
    (1, 1, 1, 3): "-",
}


def edge_name(edge: tuple[int, int]) -> str:
    return f"{edge[0]}{edge[1]}"


def scores(bits: tuple[int, ...]) -> tuple[int, ...]:
    """Outdegree vector.  bit=0 means i->j for edge (i,j), bit=1 means j->i."""
    out = [0, 0, 0, 0]
    for bit, (i, j) in zip(bits, EDGES):
        if bit == 0:
            out[i] += 1
        else:
            out[j] += 1
    return tuple(out)


def klass(bits: tuple[int, ...]) -> str:
    return CLASS_BY_SCORE[tuple(sorted(scores(bits)))]


def flip(base: tuple[int, ...], indices: tuple[int, ...]) -> tuple[int, ...]:
    bits = list(base)
    for idx in indices:
        bits[idx] = 1 - bits[idx]
    return tuple(bits)


def format_table(headers: tuple[str, ...], rows: list[tuple[str, list[str]]]) -> str:
    widths = [max(len("*"), *(len(label) for label, _ in rows))]
    widths.extend(max(len(headers[i]), *(len(row[i]) for _, row in rows)) for i in range(len(headers)))
    first = "*".ljust(widths[0]) + " | " + " | ".join(headers[i].ljust(widths[i + 1]) for i in range(len(headers)))
    rule = "-+-".join("-" * w for w in widths)
    body = []
    for label, row in rows:
        body.append(label.ljust(widths[0]) + " | " + " | ".join(row[i].ljust(widths[i + 1]) for i in range(len(row))))
    return first + "\n" + rule + "\n" + "\n".join(body)


def all_class_counts() -> Counter[str]:
    counts: Counter[str] = Counter()
    for bits in product([0, 1], repeat=len(EDGES)):
        counts[klass(tuple(bits))] += 1
    return counts


def scheme_a() -> dict[str, object]:
    base = (0, 0, 0, 0, 0, 0)
    # Fixed Hamiltonian path: 0->1, 1->2, 2->3.
    # The remaining chords are named to match the prompt's one-flip classes.
    generators = {
        "a": (EDGES.index((0, 2)),),
        "b": (EDGES.index((1, 3)),),
        "c": (EDGES.index((0, 3)),),
    }
    labels = ("E", "a", "b", "c")
    states = {"E": ()} | generators
    rows: list[tuple[str, list[str]]] = []
    for left in labels:
        row = []
        for right in labels:
            idxs = tuple(sorted(set(states[left]).symmetric_difference(states[right])))
            row.append(klass(flip(base, idxs)))
        rows.append((left, row))

    full_cube: dict[str, tuple[int, ...]] = {}
    for mask in range(8):
        active = []
        name = []
        for bit, gen in enumerate(("a", "b", "c")):
            if mask & (1 << bit):
                active.extend(generators[gen])
                name.append(gen)
        full_cube["".join(name) or "E"] = tuple(sorted(active))
    fiber = Counter(klass(flip(base, idxs)) for idxs in full_cube.values())
    by_state = {
        name: (idxs, scores(flip(base, idxs)), klass(flip(base, idxs)))
        for name, idxs in full_cube.items()
    }
    return {
        "base": base,
        "generators": generators,
        "table": rows,
        "full_cube_fiber": fiber,
        "by_state": by_state,
    }


@dataclass(frozen=True)
class SchemeBWitness:
    free_edges: tuple[int, int]
    fixed_edges: tuple[int, ...]
    fixed_bits: tuple[int, ...]
    default_free_bits: tuple[int, int]
    classes: dict[str, str]
    state_scores: dict[str, tuple[int, ...]]


def make_bits_for_scheme_b(w: SchemeBWitness, flip_x: bool, flip_y: bool) -> tuple[int, ...]:
    bits: list[int | None] = [None] * len(EDGES)
    for bit, idx in zip(w.fixed_bits, w.fixed_edges):
        bits[idx] = bit
    free_bits = list(w.default_free_bits)
    if flip_x:
        free_bits[0] = 1 - free_bits[0]
    if flip_y:
        free_bits[1] = 1 - free_bits[1]
    for bit, idx in zip(free_bits, w.free_edges):
        bits[idx] = bit
    return tuple(int(bit) for bit in bits)


def find_scheme_b_witnesses() -> list[SchemeBWitness]:
    witnesses: list[SchemeBWitness] = []
    wanted = {"E": "T", "x": "+", "y": "-", "xy": "S"}
    flags = {"E": (False, False), "x": (True, False), "y": (False, True), "xy": (True, True)}

    for free_edges in combinations(range(len(EDGES)), 2):
        fixed_edges = tuple(idx for idx in range(len(EDGES)) if idx not in free_edges)
        for fixed_bits in product([0, 1], repeat=4):
            partial = [0, 0, 0, 0]
            for bit, idx in zip(fixed_bits, fixed_edges):
                i, j = EDGES[idx]
                if bit == 0:
                    partial[i] += 1
                else:
                    partial[j] += 1
            if tuple(sorted(partial)) != (0, 1, 1, 2):
                continue
            for default_free_bits in product([0, 1], repeat=2):
                probe = SchemeBWitness(
                    free_edges=tuple(free_edges),
                    fixed_edges=fixed_edges,
                    fixed_bits=tuple(fixed_bits),
                    default_free_bits=tuple(default_free_bits),
                    classes={},
                    state_scores={},
                )
                classes = {name: klass(make_bits_for_scheme_b(probe, *flag)) for name, flag in flags.items()}
                if classes == wanted:
                    state_scores = {name: scores(make_bits_for_scheme_b(probe, *flag)) for name, flag in flags.items()}
                    witnesses.append(
                        SchemeBWitness(
                            free_edges=tuple(free_edges),
                            fixed_edges=fixed_edges,
                            fixed_bits=tuple(fixed_bits),
                            default_free_bits=tuple(default_free_bits),
                            classes=classes,
                            state_scores=state_scores,
                        )
                    )
    return witnesses


def scheme_b_table(w: SchemeBWitness) -> list[tuple[str, list[str]]]:
    labels = ("E", "x", "y")
    state_flags = {"E": (False, False), "x": (True, False), "y": (False, True)}
    rows: list[tuple[str, list[str]]] = []
    for left in labels:
        row = []
        for right in labels:
            fx = state_flags[left][0] ^ state_flags[right][0]
            fy = state_flags[left][1] ^ state_flags[right][1]
            row.append(klass(make_bits_for_scheme_b(w, fx, fy)))
        rows.append((left, row))
    return rows


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves_class: int
    injective_payload: int
    lower_order_exclusion: int
    sidecar_repair: int
    lrc_transfer: int
    proof_readiness: int
    novelty: int

    def score(self) -> int:
        return (
            3 * self.preserves_class
            + 3 * self.injective_payload
            + 3 * self.lower_order_exclusion
            + 2 * self.sidecar_repair
            + 3 * self.lrc_transfer
            + 2 * self.proof_readiness
            + self.novelty
        )


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    edge_flips = 0
    directed_3cycles = 0
    hpaths = 0

    def beats(a: Carrier, b: Carrier) -> bool:
        if a.score() != b.score():
            return a.score() > b.score()
        return a.name < b.name

    for a, b in combinations(carriers, 2):
        if not beats(a, b):
            edge_flips += 1

    for a, b, c in combinations(carriers, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            directed_3cycles += 1

    ordered = sorted(carriers, key=lambda item: (-item.score(), item.name))
    hpaths = 1
    score_hist = Counter(carrier.score() for carrier in carriers)
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": directed_3cycles,
        "hamiltonian_path_count": hpaths,
        "selected_path": [carrier.name for carrier in ordered],
        "edge_flips_against_name_tie_path": edge_flips,
    }


def carrier_tournament() -> dict[str, object]:
    carriers = [
        Carrier("lrc_packet_subbasis_theorem", 5, 5, 5, 4, 5, 4, 5),
        Carrier("two_bit_matching_basis_C4_filler", 5, 5, 5, 3, 4, 5, 4),
        Carrier("lower_order_exclusion_audit", 5, 4, 5, 4, 5, 4, 3),
        Carrier("erdos870_random_core_unique_h2", 4, 5, 5, 3, 4, 4, 4),
        Carrier("erdos870_power_two_filler_D", 4, 4, 5, 4, 4, 4, 3),
        Carrier("HYP3141_tip_tail_edge_packet", 5, 4, 4, 5, 5, 4, 3),
        Carrier("HYP3142_U4_sidecar_exit", 4, 4, 4, 5, 5, 5, 3),
        Carrier("A000568_global_class_shadow", 4, 3, 3, 4, 4, 4, 4),
        Carrier("schemeA_three_chord_tiling_cube", 5, 2, 1, 5, 3, 5, 5),
        Carrier("schemeA_S_collision_fiber", 4, 1, 0, 5, 3, 4, 4),
        Carrier("erdos870_h3_cluster_repair", 3, 3, 4, 5, 3, 3, 5),
        Carrier("raw_score_sequence_table", 3, 2, 1, 1, 2, 5, 1),
    ]
    return tournament_fingerprint(carriers)


def main() -> None:
    print("HYP-3143 n=4 dual tournament bases and Erdős-870 subbasis lens")
    print("source=codex-2026-06-28-S276")
    print("prompt_url=https://github.com/davidturturean/erdos-870")
    print()

    print("0. GLOBAL N=4 TOURNAMENT CLASS COUNTS")
    counts = all_class_counts()
    print(f"edge_order={tuple(edge_name(edge) for edge in EDGES)}")
    print(f"class_counts={dict(sorted(counts.items()))}")
    print("score_class_key={(0,1,2,3):T, (1,1,2,2):S, (0,2,2,2):+, (1,1,1,3):-}")
    print()

    print("1. SCHEME A: HAMILTONIAN-PATH TILING MODEL WITH THREE FREE CHORDS")
    a = scheme_a()
    print("fixed_path=(0->1, 1->2, 2->3)")
    print("free_chords={a:02, b:13, c:03}")
    print(format_table(("E", "a", "b", "c"), a["table"]))  # type: ignore[arg-type]
    print(f"full_3bit_class_fiber={dict(sorted(a['full_cube_fiber'].items()))}")
    for name, (idxs, out, cls) in a["by_state"].items():  # type: ignore[union-attr]
        chord_word = tuple(edge_name(EDGES[idx]) for idx in idxs)
        print(f"  state={name:3s} flipped={chord_word!s:18s} scores={out} class={cls}")
    print("readout=faithful_to_the_tiling_chords_but_not_a_unique_class_basis")
    print("lower_order_exclusion_for_S=FAILS: S occurs at Hamming weights 1,2,3")
    print()

    print("2. SCHEME B: PARTIAL SCORE 0,1,1,2 WITH TWO FREE ARCS")
    witnesses = find_scheme_b_witnesses()
    preferred = next(w for w in witnesses if w.free_edges == (EDGES.index((0, 2)), EDGES.index((1, 3))) and w.fixed_bits == (0, 0, 0, 0))
    free_pair_counter = Counter(tuple(edge_name(EDGES[idx]) for idx in w.free_edges) for w in witnesses)
    disjoint_count = sum(1 for w in witnesses if set(EDGES[w.free_edges[0]]).isdisjoint(EDGES[w.free_edges[1]]))
    print(f"matching_witness_count={len(witnesses)}")
    print(f"all_witnesses_have_free_edges_disjoint={disjoint_count == len(witnesses)}")
    print(f"free_pair_orbit_counts={dict(free_pair_counter)}")
    print("preferred_fixed_edges=(01,03,12,23) all forward; preferred_free_edges={x:02,y:13}")
    print(f"partial_score_sorted={(0, 1, 1, 2)}")
    print(format_table(("E", "x", "y"), scheme_b_table(preferred)))
    for name in ("E", "x", "y", "xy"):
        print(f"  state={name:2s} scores={preferred.state_scores[name]} class={preferred.classes[name]}")
    print("readout=minimal_two_bit_class_basis_after_a_fixed_C4_filler")
    print("lower_order_exclusion_for_S=HOLDS in squarefree flip-order: S occurs only as x*y")
    print()

    print("3. ERDŐS-870 SUBBASIS TRANSFER")
    print("- The erdos-870 proof strategy is representation economy: build unique order-h representations and exclude lower-order representations.")
    print("- Scheme B is the n=4 tournament analogue of an order-2 subbasis: fixed C4 filler + two matching free bits represents T,+,-,S exactly once.")
    print("- Scheme A is the useful but unsafe tiling quotient: it keeps all three chord directions, but S has a one-generator representative and several higher-order representatives.")
    print("- Therefore the LRC packet theorem should ask for a quotient basis only after checking exact-order representation, not merely after checking class coverage.")
    print()

    print("4. LRC14 PACKET THEOREM CANDIDATE")
    print("packet_subbasis(P) = fixed_filler_sidecar + free_obstruction_basis + lower_order_exclusion + coordinate_resurrection_or_named_debt")
    print("proposed_rule=an LRC quotient is legal only if each terminal obstruction class has a unique squarefree packet word of the declared order, or if the collision fiber is named and repaired.")
    print("n4_transfer_to_current_frontier:")
    print("  two_bit_matching_basis -> finite four-point/unital-block local carrier")
    print("  fixed_C4_filler -> powers-of-two/carry-separated deterministic scaffold")
    print("  SchemeA_S_collision -> quotient debt like HYP-3138 odd leakage or HYP-3139 boundary leakage")
    print("  HYP-3141 edge packets -> attach tail/tip labels before quotienting an edge")
    print("  HYP-3142 U4 exit -> terminal bounded-core sidecar that must not be reached by a lower-order word")
    print()

    print("5. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    print("pairwise_observable=retained class information, injectivity, lower-order exclusion, sidecar repair, LRC transfer, proof readiness, novelty")
    print("binary_gauge=edge A->B iff weighted retained-payload score(A)>score(B); tie Hamiltonian path is lexicographic")
    fp = carrier_tournament()
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print(f"edge_flips_against_name_tie_path={fp['edge_flips_against_name_tie_path']}")
    print("selected_path=" + " -> ".join(fp["selected_path"]))  # type: ignore[index]
    print()

    print("6. NEXT TESTS")
    print("A. Enumerate n=5/n=6 partial assignments and search for minimal class-bases whose fibers separate A000568 classes with lower-order exclusion.")
    print("B. Add a packet-order column to HYP-3141 edge witnesses: which class can first be represented at order 0,1,2,... after fixed filler?")
    print("C. Run the same exact-order audit on HYP-3142 bounded-core U4 sidecars: the terminal k=8 obstruction should not be reachable by a lower-order quotient word.")
    print("D. Treat the q=3 unital four-point block as a C4+matching packet: boundary cycle is filler, matching diagonals are basis bits, and AP/Goddyn-Wong labels are sidecar colors.")


if __name__ == "__main__":
    main()
