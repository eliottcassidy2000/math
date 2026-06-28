#!/usr/bin/env python3
"""HYP-3149 scout: Erdos-870 canary/filler lens for n=4 tournaments.

The prompt gives two models for the four unlabeled tournaments on four
vertices.  This script verifies the tables exactly and records the useful
transfer principle:

  * the fixed-Hamiltonian-path tiling model is a three-bit flip cube;
  * fixing one extra arc gives an order-two source with a deterministic
    filler/canary arc;
  * the second model is the c=0 slice of the first, and flipping c collapses
    the whole slice into the S class.

This refines upstream HYP-3146, HYP-3147, and HYP-3148 by keeping the
fixed-path scaffold and isolating the c=0/c=1 canary-slice test inside it.

That is the same interface shape as the Erdos #870 formalization:
an order-two source, finite fillers, finite-shift/canary control, and a
deletion/nonminimality property that must survive quotienting.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations, product
from typing import Dict, FrozenSet, Iterable, List, Sequence, Tuple


Vertex = int
Arc = Tuple[Vertex, Vertex]
ClassName = str
Subset = FrozenSet[str]

VERTICES: Tuple[Vertex, ...] = (0, 1, 2, 3)
PAIRS: Tuple[Arc, ...] = tuple(combinations(VERTICES, 2))
PATH: Tuple[Arc, ...] = ((0, 1), (1, 2), (2, 3))

# The three non-path arcs in the fixed Hamiltonian path model.  With the
# transitive orientation as the unflipped state, these names match the prompt:
# a -> +, b -> -, c -> S.
FREE3: Dict[str, Arc] = {
    "a": (0, 2),
    "b": (1, 3),
    "c": (0, 3),
}

# The two-source model fixes c=(0,3) as an extra filler/canary arc and leaves
# x=a, y=b free.
FREE2: Dict[str, str] = {"x": "a", "y": "b"}
FIXED_CANARY: str = "c"

SCORE_CLASS: Dict[Tuple[int, ...], ClassName] = {
    (0, 1, 2, 3): "T",
    (0, 2, 2, 2): "+",
    (1, 1, 1, 3): "-",
    (1, 1, 2, 2): "S",
}


def xor_subset(left: Subset, right: Subset) -> Subset:
    return frozenset(left.symmetric_difference(right))


def tournament_from_flips(flips: Iterable[str]) -> Dict[Arc, Vertex]:
    """Return winner for each unordered pair, starting from lower -> higher."""
    flip_arcs = {FREE3[name] for name in flips}
    winners: Dict[Arc, Vertex] = {}
    for i, j in PAIRS:
        winners[(i, j)] = j if (i, j) in flip_arcs else i
    return winners


def scores(winners: Dict[Arc, Vertex]) -> Tuple[int, ...]:
    out = {v: 0 for v in VERTICES}
    for winner in winners.values():
        out[winner] += 1
    return tuple(out[v] for v in VERTICES)


def score_sequence(winners: Dict[Arc, Vertex]) -> Tuple[int, ...]:
    return tuple(sorted(scores(winners)))


def classify(winners: Dict[Arc, Vertex]) -> ClassName:
    seq = score_sequence(winners)
    try:
        return SCORE_CLASS[seq]
    except KeyError as exc:
        raise AssertionError(f"unexpected n=4 score sequence {seq}") from exc


def subset_class(flips: Subset) -> ClassName:
    return classify(tournament_from_flips(flips))


def has_directed_edge(winners: Dict[Arc, Vertex], u: Vertex, v: Vertex) -> bool:
    i, j = sorted((u, v))
    return winners[(i, j)] == u


def contains_path(winners: Dict[Arc, Vertex], path: Sequence[Vertex]) -> bool:
    return all(has_directed_edge(winners, path[i], path[i + 1]) for i in range(len(path) - 1))


def hamiltonian_path_count(winners: Dict[Arc, Vertex]) -> int:
    return sum(1 for path in permutations(VERTICES) if contains_path(winners, path))


def all_labeled_tournaments() -> Iterable[Dict[Arc, Vertex]]:
    names = list(PAIRS)
    for mask in product((0, 1), repeat=len(names)):
        winners: Dict[Arc, Vertex] = {}
        for bit, (i, j) in zip(mask, names):
            winners[(i, j)] = j if bit else i
        yield winners


def table(labels: Sequence[str], label_to_subset: Dict[str, Subset]) -> List[List[str]]:
    rows: List[List[str]] = []
    for r in labels:
        row = [r]
        for c in labels:
            row.append(subset_class(xor_subset(label_to_subset[r], label_to_subset[c])))
        rows.append(row)
    return rows


def format_table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)]
    out = [" ".join(str(x).rjust(w) for x, w in zip(header, widths))]
    for row in rows:
        out.append(" ".join(str(x).rjust(w) for x, w in zip(row, widths)))
    return "\n".join(out)


def partial_scores(fixed_arc_names: Iterable[str]) -> Tuple[int, ...]:
    fixed = set(PATH)
    fixed.update(FREE3[name] for name in fixed_arc_names)
    out = {v: 0 for v in VERTICES}
    for i, j in fixed:
        out[i] += 1
    return tuple(sorted(out.values()))


def class_fibers() -> Dict[ClassName, List[Subset]]:
    fibers: Dict[ClassName, List[Subset]] = defaultdict(list)
    for bits in product((0, 1), repeat=3):
        flips = frozenset(name for name, bit in zip(("a", "b", "c"), bits) if bit)
        fibers[subset_class(flips)].append(flips)
    return dict(fibers)


def slice_readout(c_value: int) -> List[Tuple[str, ClassName]]:
    rows = []
    for x_bit, y_bit in product((0, 1), repeat=2):
        flips = set()
        if x_bit:
            flips.add("a")
        if y_bit:
            flips.add("b")
        if c_value:
            flips.add("c")
        word = f"x={x_bit},y={y_bit},c={c_value}"
        rows.append((word, subset_class(frozenset(flips))))
    return rows


def labeled_class_audit() -> List[Tuple[ClassName, int, int, int]]:
    counts: Counter[ClassName] = Counter()
    h_counts: Dict[ClassName, Counter[int]] = defaultdict(Counter)
    fixed_path_counts: Counter[ClassName] = Counter()
    for winners in all_labeled_tournaments():
        cls = classify(winners)
        counts[cls] += 1
        h_counts[cls][hamiltonian_path_count(winners)] += 1
        if contains_path(winners, (0, 1, 2, 3)):
            fixed_path_counts[cls] += 1
    audit = []
    for cls in ("T", "+", "-", "S"):
        unique_h = list(h_counts[cls])
        assert len(unique_h) == 1
        audit.append((cls, counts[cls], unique_h[0], fixed_path_counts[cls]))
    return audit


@dataclass(frozen=True)
class Carrier:
    name: str
    proof_payload: int
    raw_symmetry: int
    features: Tuple[str, ...]


CARRIERS: Tuple[Carrier, ...] = (
    Carrier(
        "fixed_c_xy_transversal",
        10,
        4,
        (
            "exact_four_class_transversal",
            "order2_xy_source",
            "canary_c_fixed",
            "partial_score_0112",
            "no_S_bulk_collision",
            "deletion_gate_visible",
        ),
    ),
    Carrier(
        "erdos870_order2_plus_filler_interface",
        9,
        2,
        (
            "order2_source",
            "finite_filler",
            "finite_shift_uniformity",
            "canary_exactness",
            "deletion_nonminimality",
            "formal_audit_boundary",
        ),
    ),
    Carrier(
        "tiling_hamiltonian_path_cube",
        7,
        8,
        (
            "fixed_hamiltonian_path",
            "three_arc_flip_cube",
            "S_bulk_exposed",
            "path_conditioned_distribution",
        ),
    ),
    Carrier(
        "edge_tip_tail_information_packet",
        7,
        3,
        (
            "tip_tail_children",
            "observer_orbit",
            "coordinate_resurrection",
            "terminal_exit",
        ),
    ),
    Carrier(
        "S_bulk_collision_fiber",
        6,
        5,
        (
            "five_preimages",
            "lost_coordinate",
            "needs_restoration",
            "scissors_fiber",
        ),
    ),
    Carrier(
        "fiber_pgf_or_distribution_sidecar",
        6,
        4,
        (
            "coefficient_layer",
            "conditional_moment",
            "quotient_legality",
        ),
    ),
    Carrier(
        "raw_score_sequence_scalar",
        2,
        6,
        ("score_sequence", "class_name_only"),
    ),
    Carrier(
        "raw_einheit_group_table_numerology",
        1,
        9,
        ("table_aesthetic", "untyped_group_shadow"),
    ),
)


def edge_key(a: str, b: str) -> Tuple[str, str]:
    ordered = {c.name: i for i, c in enumerate(CARRIERS)}
    return (a, b) if ordered[a] < ordered[b] else (b, a)


def tournament_edges(gauge: str) -> Dict[Tuple[str, str], str]:
    if gauge == "proof":
        key = {c.name: c.proof_payload for c in CARRIERS}
    elif gauge == "raw":
        key = {c.name: c.raw_symmetry for c in CARRIERS}
    else:
        raise ValueError(gauge)
    tie_path = {c.name: i for i, c in enumerate(CARRIERS)}
    edges: Dict[Tuple[str, str], str] = {}
    for a, b in combinations([c.name for c in CARRIERS], 2):
        ka, kb = key[a], key[b]
        if ka > kb:
            winner = a
        elif kb > ka:
            winner = b
        else:
            winner = a if tie_path[a] < tie_path[b] else b
        edges[(a, b)] = winner
    return edges


def score_hist(edges: Dict[Tuple[str, str], str]) -> Dict[int, int]:
    out = Counter({c.name: 0 for c in CARRIERS})
    for (_a, _b), winner in edges.items():
        out[winner] += 1
    return dict(sorted(Counter(out.values()).items()))


def directed_3cycles(edges: Dict[Tuple[str, str], str]) -> int:
    names = [c.name for c in CARRIERS]

    def beats(a: str, b: str) -> bool:
        return edges[edge_key(a, b)] == a

    total = 0
    for a, b, c in combinations(names, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            total += 1
    return total


def hamiltonian_path_count_on_carriers(edges: Dict[Tuple[str, str], str]) -> int:
    names = [c.name for c in CARRIERS]

    def beats(a: str, b: str) -> bool:
        return edges[edge_key(a, b)] == a

    return sum(
        1
        for path in permutations(names)
        if all(beats(path[i], path[i + 1]) for i in range(len(path) - 1))
    )


def selected_path(edges: Dict[Tuple[str, str], str]) -> List[str]:
    names = [c.name for c in CARRIERS]

    def outdegree(name: str) -> int:
        return sum(1 for (_a, _b), winner in edges.items() if winner == name)

    return sorted(names, key=lambda name: (-outdegree(name), name))


def edge_flips_between_gauges() -> int:
    proof = tournament_edges("proof")
    raw = tournament_edges("raw")
    return sum(1 for key in proof if proof[key] != raw[key])


def main() -> None:
    label3 = {
        "E": frozenset(),
        "a": frozenset({"a"}),
        "b": frozenset({"b"}),
        "c": frozenset({"c"}),
    }
    label2 = {
        "E": frozenset(),
        "x": frozenset({FREE2["x"]}),
        "y": frozenset({FREE2["y"]}),
    }

    table3 = table(("*", "E", "a", "b", "c")[1:], label3)
    table2 = table(("*", "E", "x", "y")[1:], label2)

    fibers = class_fibers()
    audit = labeled_class_audit()
    proof_edges = tournament_edges("proof")

    print("HYP-3149 / codex-2026-06-28-erdos870-tournament4-canary-filler")
    print("Tournament-4 two-table canary/filler scout")
    print("namespace=HYP-3149/T1214/LTI-275/LTT-173")
    print()
    print("External source inspected:")
    print("- https://github.com/davidturturean/erdos-870")
    print("- transferable interface: order-two source + deterministic filler + canary/finite-shift deletion property")
    print()
    print("Class dictionary by sorted score sequence:")
    for seq, cls in sorted(SCORE_CLASS.items(), key=lambda item: item[1]):
        print(f"  {cls:>1}: {seq}")
    print()
    print("Fixed Hamiltonian path model: path 0->1->2->3, free arcs a=(0,2), b=(1,3), c=(0,3).")
    print("Partial score sequence with only path fixed:", partial_scores(()))
    print("Prompt table, interpreted as chi3(row xor col) on the visible basis vectors:")
    print(format_table(("*", "E", "a", "b", "c"), table3))
    print()
    print("Four-fixed-arc model: fix c=(0,3) as a deterministic filler/canary, leave x=a and y=b.")
    print("Partial score sequence after fixing path plus c:", partial_scores((FIXED_CANARY,)))
    print("Prompt table, interpreted as the full chi2(row xor col) table on the c=0 slice:")
    print(format_table(("*", "E", "x", "y"), table2))
    print()
    print("Fiber audit of the 3-bit tiling cube:")
    for cls in ("T", "+", "-", "S"):
        words = ["".join(sorted(s)) or "E" for s in fibers[cls]]
        print(f"  {cls}: size={len(words)} words={words}")
    print("  max_fiber_size =", max(len(v) for v in fibers.values()))
    print("  interpretation = the tiling model is not a class group; it is a path-conditioned quotient with an S-bulk fiber.")
    print()
    print("Slice audit:")
    for word, cls in slice_readout(0):
        print(f"  c=0 slice {word} -> {cls}")
    for word, cls in slice_readout(1):
        print(f"  c=1 slice {word} -> {cls}")
    print("  readout = fixing c unflipped gives an exact four-class transversal; flipping c collapses the whole slice to S.")
    print()
    print("All labeled n=4 tournament audit:")
    print("  class labeled_count H_per_labeled_tournament fixed_path_conditioned_count")
    for cls, count, h_count, fixed_count in audit:
        print(f"  {cls:>1} {count:>13} {h_count:>24} {fixed_count:>28}")
    print("  fixed_path counts match the 3-bit cube fibers: T=1, +=1, -=1, S=5.")
    print()
    print("Erdos-870 transfer to the LRC/tournament frontier:")
    print("  order_two_source        = free coordinates x,y after c is fixed")
    print("  deterministic_filler    = fixed canary arc c plus the Hamiltonian path boundary")
    print("  canary_exactness        = c=1 makes every x/y completion land in S, so c cannot be forgotten as raw symmetry")
    print("  finite_shift_uniformity = the same x/y table must survive every legal edge/tip-tail packet, not only this labeling")
    print("  deletion_property       = if one core coordinate is removed, the certificate must name the restoration sidecar or debt")
    print()
    print("Candidate invariant:")
    print("  tournament4_canary_filler_certificate =")
    print("    fixed_path_word + c_canary_status + xy_completion_table")
    print("    + S_bulk_fiber_words + deletion/restoration_sidecar")
    print("    + edge_tip_tail_exit_or_named_debt")
    print()
    print("Tournament Analysis over proof carriers:")
    print("  vertices_are = proof carriers / quotient operators, not runners or raw arcs")
    print("  pairwise_observable = retained proof payload under the order-two-source-plus-filler interface")
    print("  switch_gauge = proof_payload beats raw table symmetry; tie path is the listed carrier order")
    print("  proof_score_hist =", score_hist(proof_edges))
    print("  directed_3cycles =", directed_3cycles(proof_edges))
    print("  hamiltonian_path_count =", hamiltonian_path_count_on_carriers(proof_edges))
    print("  edge_flips_vs_raw_symmetry_gauge =", edge_flips_between_gauges())
    print("  selected_path =", " -> ".join(selected_path(proof_edges)))
    print()
    print("Assumption challenged:")
    print("  The table entries should not be promoted to a group law on isomorphism classes.")
    print("  The useful object is a quotient map with a visible missing-coordinate fiber.")
    print("  In LRC14 terms, c is a canary/filler coordinate: legal if fixed and audited, illegal if erased.")


if __name__ == "__main__":
    main()
