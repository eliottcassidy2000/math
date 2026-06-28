#!/usr/bin/env python3
"""n=3 edge-flip kernel, Worpitzky refinement, and ordered function channels.

The prompt asks for the tournament on three vertices from the perspective of
its edges.  This script makes the two-class / three-edge graph exact, then
records how Worpitzky descent data and the functions a*b, a+b, a^b, b^a act
as quotient guardrails for LRC packet work.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations, product


VERTICES = (0, 1, 2)
EDGES = ((0, 1), (0, 2), (1, 2))
# A fixed cyclic reference: 0->1, 1->2, 2->0.  Bits are coin flips against it.
CYCLIC_REFERENCE = {
    (0, 1): 0,  # 0 means smaller -> larger for this edge.
    (0, 2): 1,  # 1 means larger -> smaller, so 2->0.
    (1, 2): 0,
}


def bits_to_orientation(bits: tuple[int, int, int]) -> tuple[int, int, int]:
    """Convert agree/disagree bits relative to the cyclic reference into edge bits.

    Edge bit convention: 0 means i->j for i<j, 1 means j->i.
    Coin bit 0 means keep the cyclic reference; 1 means flip it.
    """
    out = []
    for bit, edge in zip(bits, EDGES):
        out.append(CYCLIC_REFERENCE[edge] ^ bit)
    return tuple(out)


def scores(edge_bits: tuple[int, int, int]) -> tuple[int, int, int]:
    out = [0, 0, 0]
    for bit, (i, j) in zip(edge_bits, EDGES):
        if bit == 0:
            out[i] += 1
        else:
            out[j] += 1
    return tuple(out)


def klass(edge_bits: tuple[int, int, int]) -> str:
    ss = tuple(sorted(scores(edge_bits)))
    if ss == (0, 1, 2):
        return "T"  # transitive
    if ss == (1, 1, 1):
        return "C"  # cyclic
    raise ValueError(ss)


def flip(bits: tuple[int, int, int], idx: int) -> tuple[int, int, int]:
    out = list(bits)
    out[idx] = 1 - out[idx]
    return tuple(out)


def coin_word(bits: tuple[int, int, int]) -> str:
    return "".join("H" if bit == 0 else "T" for bit in bits)


def source_to_sink_order(edge_bits: tuple[int, int, int]) -> tuple[int, int, int] | None:
    if klass(edge_bits) != "T":
        return None
    out = scores(edge_bits)
    # In a transitive 3-tournament the source has outdegree 2, middle 1, sink 0.
    return tuple(vertex for vertex, _ in sorted(enumerate(out), key=lambda item: -item[1]))


def descents(order: tuple[int, int, int]) -> int:
    return sum(1 for a, b in zip(order, order[1:]) if a > b)


def worpitzky_n3_identity_values(limit: int = 8) -> list[tuple[int, int, int]]:
    """Return rows (x, lhs, rhs) for x^3 = C(x+2,3)+4C(x+1,3)+C(x,3)."""

    def choose3(n: int) -> int:
        if n < 3:
            return 0
        return n * (n - 1) * (n - 2) // 6

    rows = []
    for x in range(limit + 1):
        lhs = x**3
        rhs = choose3(x + 2) + 4 * choose3(x + 1) + choose3(x)
        rows.append((x, lhs, rhs))
    return rows


def edge_function_payload(values: tuple[int, int, int] = (1, 2, 3)) -> dict[str, dict[str, int | tuple[int, int]]]:
    payload = {}
    for i, j in combinations(range(3), 2):
        a, b = values[i], values[j]
        payload[f"{a},{b}"] = {
            "a+b": a + b,
            "a*b": a * b,
            "a^b,b^a": (a**b, b**a),
            "ordered_gap": abs(a**b - b**a),
        }
    return payload


def transition_kernel() -> dict[str, object]:
    states = tuple(product((0, 1), repeat=3))
    directed_counts: Counter[tuple[str, str]] = Counter()
    labelled_counts: Counter[tuple[str, str, str]] = Counter()
    per_state = {}
    for bits in states:
        edge_bits = bits_to_orientation(bits)
        cls = klass(edge_bits)
        transitions = []
        for idx, edge in enumerate(EDGES):
            nxt = flip(bits, idx)
            nxt_cls = klass(bits_to_orientation(nxt))
            directed_counts[(cls, nxt_cls)] += 1
            labelled_counts[(cls, f"{edge[0]}{edge[1]}", nxt_cls)] += 1
            transitions.append((f"{edge[0]}{edge[1]}", coin_word(nxt), nxt_cls))
        order = source_to_sink_order(edge_bits)
        per_state[coin_word(bits)] = {
            "coin_weight": sum(bits),
            "class": cls,
            "scores": scores(edge_bits),
            "source_to_sink_order": order,
            "descents": None if order is None else descents(order),
            "transitions": transitions,
        }
    return {
        "states": states,
        "directed_counts": directed_counts,
        "labelled_counts": labelled_counts,
        "per_state": per_state,
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    orientation_payload: int
    quotient_guard: int
    lrc_transfer: int
    worpitzky_signal: int
    function_signal: int
    proof_readiness: int

    def score(self) -> int:
        return (
            3 * self.orientation_payload
            + 3 * self.quotient_guard
            + 3 * self.lrc_transfer
            + 2 * self.worpitzky_signal
            + 2 * self.function_signal
            + 2 * self.proof_readiness
        )


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    def beats(a: Carrier, b: Carrier) -> bool:
        if a.score() != b.score():
            return a.score() > b.score()
        return a.name < b.name

    directed_3cycles = 0
    for a, b, c in combinations(carriers, 3):
        ab, bc, ca = beats(a, b), beats(b, c), beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            directed_3cycles += 1

    # Hamiltonian paths by brute force, since the carrier set is small.
    hpaths = 0
    for order in permutations(carriers):
        if all(beats(order[i], order[i + 1]) for i in range(len(order) - 1)):
            hpaths += 1

    ordered = sorted(carriers, key=lambda item: (-item.score(), item.name))
    return {
        "score_hist": dict(sorted(Counter(c.score() for c in carriers).items())),
        "directed_3cycles": directed_3cycles,
        "hamiltonian_path_count": hpaths,
        "selected_path": [carrier.name for carrier in ordered],
    }


def proof_carrier_tournament() -> dict[str, object]:
    carriers = [
        Carrier("edge_flip_markov_kernel_C_T", 5, 5, 5, 4, 3, 5),
        Carrier("minority_edge_gate_inside_T", 5, 5, 5, 3, 4, 5),
        Carrier("worpitzky_descent_refinement_of_T_fiber", 4, 5, 4, 5, 3, 4),
        Carrier("ordered_exponential_channel_a_to_b", 5, 4, 4, 2, 5, 4),
        Carrier("packet_function_object_not_scalar", 4, 5, 5, 3, 4, 4),
        Carrier("coin_straight_mixed_stationary_measure", 3, 4, 4, 3, 3, 5),
        Carrier("symmetric_sum_product_shadow", 1, 2, 2, 1, 3, 5),
        Carrier("raw_score_class_two_nodes", 2, 2, 2, 1, 1, 5),
    ]
    return tournament_fingerprint(carriers)


def main() -> None:
    tk = transition_kernel()
    print("HYP-3147 n=3 edge-flip Worpitzky/function kernel")
    print("source=codex-2026-06-28-S277")
    print()

    print("1. TWO CLASS / THREE EDGE GRAPH")
    print("coin_reference=cyclic orientation 0->1->2->0; H keeps reference, T flips reference")
    print("class_rule=straight coin words HHH/TTT are cyclic C; 2:1 mixes are transitive T")
    class_counts = Counter(item["class"] for item in tk["per_state"].values())  # type: ignore[union-attr]
    print(f"class_counts={dict(sorted(class_counts.items()))}")
    print("per_state:")
    for word, item in sorted(tk["per_state"].items()):  # type: ignore[union-attr]
        print(
            f"  {word}: class={item['class']} weight={item['coin_weight']} "
            f"scores={item['scores']} order={item['source_to_sink_order']} "
            f"descents={item['descents']}"
        )
        print("    flips=" + ", ".join(f"{edge}->{nxt}/{cls}" for edge, nxt, cls in item["transitions"]))
    print()

    print("2. COLLAPSED FLIP COUNTS AND MARKOV KERNEL")
    directed = tk["directed_counts"]  # type: ignore[assignment]
    print(f"directed_flip_counts={dict(sorted(directed.items()))}")
    print("per_vertex_outdegrees: C->T=3; T->C=1 and T->T=2 on average/exactly")
    print("class_transition_matrix_rows_C_T=[[0,1],[1/3,2/3]]")
    print("stationary_distribution=(C=1/4,T=3/4), matching straight/mixed coin split 2/8 and 6/8")
    print("eigenvalues=(1,-1/3); the negative eigenvalue is the tiny three-edge oscillation left after class collapse")
    print()

    print("3. EDGE-LABEL COUNTS")
    labelled = tk["labelled_counts"]  # type: ignore[assignment]
    for key, count in sorted(labelled.items()):
        print(f"  {key}: {count}")
    print("readout=from a cyclic state every edge is a C->T gate; from a transitive state exactly one minority edge is a T->C gate")
    print()

    print("4. WORPITZKY REFINEMENT")
    descent_counts = Counter(
        item["descents"] for item in tk["per_state"].values() if item["descents"] is not None  # type: ignore[union-attr]
    )
    print(f"transitive_fiber_descent_counts={dict(sorted(descent_counts.items()))}")
    print("Eulerian_A(3,k)=(1,4,1), so Worpitzky splits the six transitive states into 1+4+1")
    print("identity_check=x^3=C(x+2,3)+4*C(x+1,3)+C(x,3)")
    for x, lhs, rhs in worpitzky_n3_identity_values():
        print(f"  x={x}: lhs={lhs} rhs={rhs} ok={lhs == rhs}")
    print("readout=Worpitzky does not replace the C/T kernel; it refines the T fiber by ordered path descents")
    print()

    print("5. FUNCTION QUARTET ON EDGE VALUES")
    payload = edge_function_payload()
    for pair, data in payload.items():
        print(f"  pair={pair}: sum={data['a+b']} product={data['a*b']} exponents={data['a^b,b^a']} ordered_gap={data['ordered_gap']}")
    print("commutative_channels=a+b,a*b; ordered_channels=a^b,b^a")
    print("readout=sum/product see an undirected edge shadow; exponentials remember which endpoint is base and which is exponent")
    print()

    print("6. LRC PACKET TRANSFER")
    print("n3_edge_kernel_packet = class(C/T) + minority_edge_gate + Worpitzky_descent_word + ordered_function_channel")
    print("guardrail=symmetric pair functions alone cannot certify an edge-flip quotient, because they are blind to arc orientation")
    print("HYP3141_update=add edge_flip_class_kernel, minority_edge_gate, and ordered_function_payload to edge witnesses")
    print("HYP3142_update=read the antisymmetric nonmax block as order-sensitive payload, not as disposable noise")
    print("HYP3143_update=packet_order should be computed after separating straight-vs-mixed and then Worpitzky-refining the transitive fiber")
    print()

    print("7. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    print("pairwise_observable=orientation payload, quotient guard strength, LRC transfer, Worpitzky signal, function-order signal, proof readiness")
    print("binary_gauge=edge A->B iff weighted proof-carrier score(A)>score(B); tie Hamiltonian path is lexicographic")
    fp = proof_carrier_tournament()
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("selected_path=" + " -> ".join(fp["selected_path"]))  # type: ignore[index]
    print()

    print("8. NEXT TESTS")
    print("A. Add these fields to the n=4 packet-basis scout: minority/majority flip role, Worpitzky descent word on every transitive subtriangle, and ordered-function payload.")
    print("B. For HYP-3141 edge witnesses, compute the two-state edge kernel on every 3-vertex shadow around a directed proof edge.")
    print("C. For HYP-3142, ask whether the k=8 antisymmetric shell is bounded by summing n=3 ordered-function oscillations with Worpitzky descent weights.")
    print("D. For LRC Rprime rows, treat C/T transitions as a small signed Markov kernel and test whether the -1/3 eigenmode aligns with the HYP-3129 signed SPEC low modes.")


if __name__ == "__main__":
    main()
