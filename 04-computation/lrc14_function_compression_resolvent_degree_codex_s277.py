#!/usr/bin/env python3
"""Function-compression and resolvent-degree guardrail for the LRC14 packet.

This scout sits above HYP-3143..HYP-3149.  It recomputes the small finite
signals from the prompt, then records the proof-facing abstraction:

  * a quotient is legal for a function only when the function is constant on
    quotient fibers or a sidecar reconstructs the lost coordinate;
  * the n=3 edge kernel and n=4 canary/filler tables are preflight tests for
    that rule;
  * the current k=8 bounded-core route remains below degree 5 only because
    these sidecars keep the active dual in a solvable degree <= 4 category.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations, product
from math import prod
from typing import Dict, FrozenSet, Iterable, List, Sequence, Tuple


def choose3(n: int) -> int:
    if n < 3:
        return 0
    return n * (n - 1) * (n - 2) // 6


def pair_function_audit(samples: Sequence[Tuple[int, int]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for a, b in samples:
        rows.append(
            {
                "pair": (a, b),
                "sum_swap_equal": a + b == b + a,
                "product_swap_equal": a * b == b * a,
                "a_pow_b": a**b,
                "b_pow_a": b**a,
                "power_swap_equal": a**b == b**a,
                "ordered_gap": abs(a**b - b**a),
            }
        )
    return rows


K3_EDGES: Tuple[Tuple[int, int], ...] = ((0, 1), (0, 2), (1, 2))
# Edge bit convention: 0 means i->j for i<j, 1 means j->i.
K3_CYCLIC_REFERENCE = {
    (0, 1): 0,
    (0, 2): 1,
    (1, 2): 0,
}


def k3_orientation(bits: Tuple[int, int, int]) -> Tuple[int, int, int]:
    return tuple(K3_CYCLIC_REFERENCE[e] ^ bit for e, bit in zip(K3_EDGES, bits))


def k3_scores(edge_bits: Tuple[int, int, int]) -> Tuple[int, int, int]:
    out = [0, 0, 0]
    for bit, (i, j) in zip(edge_bits, K3_EDGES):
        if bit == 0:
            out[i] += 1
        else:
            out[j] += 1
    return tuple(out)


def k3_class(bits: Tuple[int, int, int]) -> str:
    seq = tuple(sorted(k3_scores(k3_orientation(bits))))
    if seq == (1, 1, 1):
        return "C"
    if seq == (0, 1, 2):
        return "T"
    raise AssertionError(seq)


def flip_bit(bits: Tuple[int, ...], index: int) -> Tuple[int, ...]:
    out = list(bits)
    out[index] = 1 - out[index]
    return tuple(out)


def k3_kernel() -> Dict[str, object]:
    states = tuple(product((0, 1), repeat=3))
    raw: Counter[Tuple[str, str]] = Counter()
    role_counts: Counter[str] = Counter()
    for bits in states:
        source = k3_class(bits)
        weight = sum(bits)
        for i in range(3):
            target = k3_class(flip_bit(bits, i))
            raw[(source, target)] += 1
            if source == "C" and target == "T":
                role_counts["cycle_edge_breaks_to_T"] += 1
            elif source == "T" and target == "C":
                role_counts["minority_edge_closes_cycle"] += 1
            elif source == "T" and target == "T":
                if weight in (1, 2):
                    role_counts["majority_edge_self_flip"] += 1
    sizes = Counter(k3_class(bits) for bits in states)
    quotient = {
        (src, dst): Fraction(raw[(src, dst)], sizes[src]) for src in ("T", "C") for dst in ("T", "C")
    }
    return {"sizes": sizes, "raw": raw, "quotient": quotient, "roles": role_counts}


def coin_kernel() -> Dict[Tuple[str, str], Fraction]:
    states = tuple(product((0, 1), repeat=3))

    def cls(bits: Tuple[int, int, int]) -> str:
        return "same" if bits.count(bits[0]) == 3 else "mix"

    raw: Counter[Tuple[str, str]] = Counter()
    sizes = Counter(cls(bits) for bits in states)
    for bits in states:
        for i in range(3):
            raw[(cls(bits), cls(flip_bit(bits, i)))] += 1
    return {
        (src, dst): Fraction(raw[(src, dst)], sizes[src])
        for src in ("mix", "same")
        for dst in ("mix", "same")
    }


Vertex = int
Arc = Tuple[Vertex, Vertex]
Subset = FrozenSet[str]

N4_VERTICES: Tuple[int, ...] = (0, 1, 2, 3)
N4_PAIRS: Tuple[Arc, ...] = tuple(combinations(N4_VERTICES, 2))
N4_FREE: Dict[str, Arc] = {"a": (0, 2), "b": (1, 3), "c": (0, 3)}
N4_SCORE_CLASS: Dict[Tuple[int, ...], str] = {
    (0, 1, 2, 3): "T",
    (0, 2, 2, 2): "+",
    (1, 1, 1, 3): "-",
    (1, 1, 2, 2): "S",
}


def n4_winners(flips: Iterable[str]) -> Dict[Arc, int]:
    flipped = {N4_FREE[name] for name in flips}
    winners = {}
    for i, j in N4_PAIRS:
        winners[(i, j)] = j if (i, j) in flipped else i
    return winners


def n4_class(flips: Iterable[str]) -> str:
    score = [0, 0, 0, 0]
    for winner in n4_winners(flips).values():
        score[winner] += 1
    return N4_SCORE_CLASS[tuple(sorted(score))]


def n4_fibers() -> Dict[str, List[str]]:
    fibers: Dict[str, List[str]] = defaultdict(list)
    for bits in product((0, 1), repeat=3):
        word = "".join(name for name, bit in zip(("a", "b", "c"), bits) if bit) or "E"
        flips = [] if word == "E" else list(word)
        fibers[n4_class(flips)].append(word)
    return dict(sorted(fibers.items()))


def n4_or_compression(word: str) -> Tuple[int, int]:
    flips = set() if word == "E" else set(word)
    x = int("a" in flips or "c" in flips)
    y = int("b" in flips or "c" in flips)
    return (x, y)


def n4_compression_rows() -> List[Tuple[str, str, Tuple[int, int], str]]:
    rows = []
    for cls, words in n4_fibers().items():
        for word in words:
            x, y = n4_or_compression(word)
            xy_cls = { (0, 0): "T", (1, 0): "+", (0, 1): "-", (1, 1): "S" }[(x, y)]
            rows.append((word, cls, (x, y), xy_cls))
    return sorted(rows, key=lambda row: (row[1], row[0]))


def worpitzky_rows(limit: int = 8) -> List[Tuple[int, int, int]]:
    return [(x, x**3, choose3(x + 2) + 4 * choose3(x + 1) + choose3(x)) for x in range(limit + 1)]


@dataclass(frozen=True)
class DegreeLedgerRow:
    name: str
    degree: int
    preserved: str
    sidecar: str
    wall_status: str


def degree_ledger() -> List[DegreeLedgerRow]:
    return [
        DegreeLedgerRow(
            "symmetric_pair_shadow",
            1,
            "a+b and a*b survive unordered pair fibers",
            "none, if target predicate is symmetric",
            "below Abel-Ruffini wall",
        ),
        DegreeLedgerRow(
            "ordered_pair_channel",
            2,
            "a^b,b^a orientation channel",
            "tail/tip or base/exponent sidecar",
            "solvable by ordered sidecar, not by scalar quotient",
        ),
        DegreeLedgerRow(
            "K3_edge_flip_kernel",
            2,
            "C/T transition law with eigenmode -1/3",
            "minority-edge and Worpitzky descent word",
            "quadratic-sized preflight",
        ),
        DegreeLedgerRow(
            "N4_filler_canary_square",
            2,
            "x,y exact four-class transversal after fixing c",
            "c-canary status plus S-bulk fiber",
            "Klein-square source, no degree growth",
        ),
        DegreeLedgerRow(
            "k8_bounded_core_resolvent",
            4,
            "HYP-3132/HYP-3142 quartic/biquadratic U4 exit",
            "even fold plus odd-coordinate resurrection",
            "deepest known live core; still below degree 5",
        ),
        DegreeLedgerRow(
            "illegal_raw_scalarization",
            5,
            "none: generic quintic-style branch loss",
            "route to named debt before theorem use",
            "Abel-Ruffini wall / not accepted as proof currency",
        ),
    ]


@dataclass(frozen=True)
class Carrier:
    name: str
    order_payload: int
    curve_payload: int
    quotient_legality: int
    lrc_transfer: int
    degree_control: int

    def score(self) -> int:
        return (
            3 * self.order_payload
            + 3 * self.curve_payload
            + 4 * self.quotient_legality
            + 4 * self.lrc_transfer
            + 3 * self.degree_control
        )


def tournament_fingerprint(carriers: Sequence[Carrier]) -> Dict[str, object]:
    def beats(a: Carrier, b: Carrier) -> bool:
        if a.score() != b.score():
            return a.score() > b.score()
        return a.name < b.name

    outdeg = {carrier.name: 0 for carrier in carriers}
    edge_flips = 0
    for a, b in combinations(carriers, 2):
        if beats(a, b):
            outdeg[a.name] += 1
            if a.name > b.name:
                edge_flips += 1
        else:
            outdeg[b.name] += 1
            if b.name > a.name:
                edge_flips += 1

    directed_3cycles = 0
    for a, b, c in combinations(carriers, 3):
        ab, bc, ca = beats(a, b), beats(b, c), beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            directed_3cycles += 1

    hpaths = 0
    for order in permutations(carriers):
        if all(beats(order[i], order[i + 1]) for i in range(len(order) - 1)):
            hpaths += 1

    selected = sorted(carriers, key=lambda item: (-item.score(), item.name))
    return {
        "score_hist": dict(sorted(Counter(carrier.score() for carrier in carriers).items())),
        "outdegrees": dict(sorted(outdeg.items())),
        "directed_3cycles": directed_3cycles,
        "edge_flips_vs_lexical": edge_flips,
        "hamiltonian_path_count": hpaths,
        "selected_path": [carrier.name for carrier in selected],
    }


def proof_carrier_tournament() -> Dict[str, object]:
    carriers = [
        Carrier("function_compression_legality_certificate", 5, 5, 5, 5, 5),
        Carrier("bounded_core_degree4_guardrail", 5, 4, 5, 5, 5),
        Carrier("ordered_pair_tip_tail_sidecar", 5, 3, 5, 5, 4),
        Carrier("fiber_PGF_full_curve_payload", 4, 5, 5, 5, 4),
        Carrier("n4_canary_filler_exact_transversal", 4, 4, 5, 4, 4),
        Carrier("n3_edge_flip_minoritary_gate", 5, 4, 4, 4, 3),
        Carrier("worpitzky_basis_curve", 4, 5, 3, 4, 3),
        Carrier("A000568_n_le_7_tameness_window", 3, 3, 4, 4, 4),
        Carrier("symmetric_sum_product_shadow", 1, 1, 3, 2, 3),
        Carrier("raw_scalar_class_count", 1, 1, 1, 1, 1),
    ]
    return tournament_fingerprint(carriers)


def main() -> None:
    print("HYP-3150 / codex-2026-06-28-S277")
    print("Function-compression resolvent-degree guardrail")
    print("status=exact finite scout plus synthesis; not an LRC14 proof")
    print()

    print("1. PAIR FUNCTIONS AS QUOTIENT TESTS")
    samples = ((2, 3), (2, 4), (4, 2), (3, 5), (5, 3), (4, 4))
    rows = pair_function_audit(samples)
    for row in rows:
        print(
            f"pair={row['pair']} sum_safe={row['sum_swap_equal']} "
            f"product_safe={row['product_swap_equal']} "
            f"a^b={row['a_pow_b']} b^a={row['b_pow_a']} "
            f"power_swap_equal={row['power_swap_equal']} "
            f"ordered_gap={row['ordered_gap']}"
        )
    print("readout: sum/product are quotient-safe on unordered pairs; powers are ordered channels except accidental equalities.")
    print()

    print("2. K3 EDGE-FLIP AND COIN KERNEL")
    kernel = k3_kernel()
    print(f"class_sizes={dict(sorted(kernel['sizes'].items()))}")
    print("raw_labelled_flip_counts:")
    for src in ("T", "C"):
        print(f"  from {src}: to T={kernel['raw'][(src, 'T')]} to C={kernel['raw'][(src, 'C')]}")
    print("quotient_multiplicity_matrix_rows_T_C:")
    for src in ("T", "C"):
        print(f"  from {src}: to T={kernel['quotient'][(src, 'T')]} to C={kernel['quotient'][(src, 'C')]}")
    print("normalized_markov_rows_T_C: from T=(2/3,1/3), from C=(1,0); stationary T=3/4,C=1/4; eigenvalues=1,-1/3")
    print(f"edge_role_counts={dict(sorted(kernel['roles'].items()))}")
    print("coin_kernel_rows_mix_same:")
    ck = coin_kernel()
    for src in ("mix", "same"):
        print(f"  from {src}: to mix={ck[(src, 'mix')]} to same={ck[(src, 'same')]}")
    print()

    print("3. WORPITZKY BASIS WARNING")
    print("identity: x^3 = C(x+2,3)+4*C(x+1,3)+C(x,3)")
    for x, lhs, rhs in worpitzky_rows():
        print(f"  x={x}: lhs={lhs} rhs={rhs} ok={lhs == rhs}")
    print("basis_coefficients=(1,4,1); this is a curve/basis payload, not permission to keep only a scalar.")
    print()

    print("4. N4 CANARY/FILLER COMPRESSION")
    print(f"fixed_path_fibers={n4_fibers()}")
    print("compression=x=(a OR c), y=(b OR c)")
    for word, cls, xy, xy_cls in n4_compression_rows():
        print(f"  word={word:3s} class={cls} -> xy={xy} xy_class={xy_cls}")
    print("readout: the OR map is class-preserving but nonlinear; c is filler when fixed and a canary when erased.")
    print()

    print("5. DEGREE LEDGER")
    for row in degree_ledger():
        print(
            f"  {row.name}: degree={row.degree}; preserved={row.preserved}; "
            f"sidecar={row.sidecar}; status={row.wall_status}"
        )
    print("meta_readout: the current hard core is useful because every accepted carrier stays at degree <= 4; degree 5 is routed to named debt.")
    print()

    print("6. PROOF PACKET FIELDS")
    fields = [
        "function_payload_type",
        "unordered_pair_survival",
        "ordered_pair_sidecar",
        "three_edge_flip_kernel",
        "worpitzky_basis_curve",
        "state_level_pgf_split",
        "compression_map_word",
        "canary_filler_status",
        "resolvent_degree_ceiling",
        "abel_ruffini_wall_status",
        "finite_tameness_window",
        "quotient_legality_status",
        "terminal_exit_or_named_debt",
    ]
    for field in fields:
        print(f"  - {field}")
    print()

    print("7. TOURNAMENT ANALYSIS")
    print("vertices=proof carriers / function payloads, not runners or raw arcs")
    print("pairwise_observable=retained order+curve payload, quotient legality, LRC transfer, and degree-control score")
    print("binary_gauge=A beats B iff weighted retained payload is larger; lexical tie path")
    fp = proof_carrier_tournament()
    print(f"score_hist={fp['score_hist']}")
    print(f"outdegrees={fp['outdegrees']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"edge_flips_vs_lexical={fp['edge_flips_vs_lexical']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("selected_path=")
    for name in fp["selected_path"]:
        print(f"  {name}")


if __name__ == "__main__":
    main()
