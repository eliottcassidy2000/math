#!/usr/bin/env python3
"""
codex-S107: realizability check for the q=3 unital pair-slot lead.

HYP-2892 observes that the Hermitian unital q=3 has 28 points, matching the
28 pair slots C(8,2) of a k=8 AP layer.  This script tests a sharper question:
can the unital blocks be interpreted as a natural S_8-invariant block system
on the edges of K_8?

The answer is negative for the most canonical interpretation.  No union of
S_8-orbits of 4-edge subsets of K_8 forms a 2-(28,4,1) design on K_8 edges.
Thus the unital remains a useful tight frame, but any LRC application needs a
non-canonical labelling, a weighting, or a separate category-1 realizability
map.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product


def k_edges(n: int) -> list[tuple[int, int]]:
    return list(combinations(range(n), 2))


def component_signatures(block: tuple[int, ...], edges: list[tuple[int, int]]) -> tuple[tuple[int, int, tuple[int, ...]], ...]:
    adj: dict[int, set[int]] = defaultdict(set)
    active = set()
    for idx in block:
        a, b = edges[idx]
        active.add(a)
        active.add(b)
        adj[a].add(b)
        adj[b].add(a)

    seen = set()
    comps = []
    for root in sorted(active):
        if root in seen:
            continue
        stack = [root]
        comp = set()
        seen.add(root)
        while stack:
            v = stack.pop()
            comp.add(v)
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        degs = tuple(sorted((len(adj[v]) for v in comp), reverse=True))
        ecount = sum(len(adj[v]) for v in comp) // 2
        comps.append((len(comp), ecount, degs))
    return tuple(sorted(comps))


def pair_relation(edge_a: tuple[int, int], edge_b: tuple[int, int]) -> str:
    return "adjacent" if set(edge_a) & set(edge_b) else "disjoint"


def orbit_inventory(n: int = 8) -> tuple[list[tuple[int, int]], list[dict[str, object]]]:
    edges = k_edges(n)
    by_sig: dict[tuple[tuple[int, int, tuple[int, ...]], ...], list[tuple[int, ...]]] = defaultdict(list)
    for block in combinations(range(len(edges)), 4):
        by_sig[component_signatures(block, edges)].append(block)

    rows = []
    total_adjacent_pairs = sum(
        1 for a, b in combinations(range(len(edges)), 2)
        if pair_relation(edges[a], edges[b]) == "adjacent"
    )
    total_disjoint_pairs = sum(
        1 for a, b in combinations(range(len(edges)), 2)
        if pair_relation(edges[a], edges[b]) == "disjoint"
    )

    for sig, blocks in sorted(by_sig.items(), key=lambda kv: (len(kv[1]), kv[0])):
        sample = blocks[0]
        rels = Counter(pair_relation(edges[a], edges[b]) for a, b in combinations(sample, 2))
        count = len(blocks)
        rows.append(
            {
                "signature": sig,
                "orbit_size": count,
                "sample_edges": tuple(edges[i] for i in sample),
                "within_adjacent_pairs": rels["adjacent"],
                "within_disjoint_pairs": rels["disjoint"],
                "point_rep_if_taken": Fraction(count * 4, len(edges)),
                "adjacent_lambda_if_taken": Fraction(count * rels["adjacent"], total_adjacent_pairs),
                "disjoint_lambda_if_taken": Fraction(count * rels["disjoint"], total_disjoint_pairs),
            }
        )
    return edges, rows


def invariant_design_solutions(rows: list[dict[str, object]]) -> list[tuple[int, ...]]:
    solutions = []
    for mask in product((0, 1), repeat=len(rows)):
        block_count = sum(rows[i]["orbit_size"] for i, bit in enumerate(mask) if bit)
        if block_count != 63:
            continue
        adj_lam = sum(rows[i]["adjacent_lambda_if_taken"] for i, bit in enumerate(mask) if bit)
        dis_lam = sum(rows[i]["disjoint_lambda_if_taken"] for i, bit in enumerate(mask) if bit)
        if adj_lam == 1 and dis_lam == 1:
            solutions.append(tuple(i for i, bit in enumerate(mask) if bit))
    return solutions


def ap8_packet_stats() -> dict[str, object]:
    edges = k_edges(8)
    by_sum: dict[int, list[tuple[int, int]]] = defaultdict(list)
    by_diff: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for a, b in edges:
        by_sum[a + b].append((a, b))
        by_diff[b - a].append((a, b))
    return {
        "sum_class_size_hist": dict(sorted(Counter(len(v) for v in by_sum.values()).items())),
        "diff_class_size_hist": dict(sorted(Counter(len(v) for v in by_diff.values()).items())),
        "equal_sum_pair_count": sum(len(v) * (len(v) - 1) // 2 for v in by_sum.values()),
        "equal_diff_pair_count": sum(len(v) * (len(v) - 1) // 2 for v in by_diff.values()),
        "sum_classes": {k: v for k, v in sorted(by_sum.items())},
        "diff_classes": {k: v for k, v in sorted(by_diff.items())},
    }


def carrier_tournament() -> list[tuple[str, tuple[int, int, int, int]]]:
    """Proof-carrier ranking: LRC faithfulness before symmetry elegance."""

    rows = [
        ("category1_exact_tiling_AP_plus_GW", (4, 4, 3, 4)),
        ("AP8_sum_difference_packets", (3, 3, 2, 4)),
        ("unital_q3_weighted_pair_frame", (2, 1, 4, 2)),
        ("scalar_additive_energy", (1, 0, 1, 3)),
        ("S8_invariant_unital_block_system", (0, 0, 4, 4)),
    ]
    return sorted(rows, key=lambda x: x[1], reverse=True)


def main() -> None:
    edges, rows = orbit_inventory(8)
    print("S107: q=3 unital pair-slot realizability check")
    print(f"K8 pair slots={len(edges)}; unordered pair-slot pairs={len(edges) * (len(edges) - 1) // 2}")
    print(f"4-edge S8 orbit types={len(rows)}")
    print()
    for i, row in enumerate(rows):
        print(
            f"orbit {i:02d}: size={row['orbit_size']:4d} "
            f"sig={row['signature']} sample={row['sample_edges']} "
            f"point_rep={row['point_rep_if_taken']} "
            f"lambda_adj={row['adjacent_lambda_if_taken']} "
            f"lambda_dis={row['disjoint_lambda_if_taken']}"
        )

    solutions = invariant_design_solutions(rows)
    print()
    print(f"S8-invariant unions realizing 2-(28,4,1): {len(solutions)}")
    if solutions:
        print(f"solutions={solutions}")
    else:
        print("readout: no natural S8-invariant 4-edge block system is the q=3 unital.")

    stats = ap8_packet_stats()
    print()
    print("AP8 natural packet statistics:")
    print(f"  sum class size hist={stats['sum_class_size_hist']}")
    print(f"  diff class size hist={stats['diff_class_size_hist']}")
    print(f"  equal-sum pair count={stats['equal_sum_pair_count']}")
    print(f"  equal-diff pair count={stats['equal_diff_pair_count']}")
    print("  readout: AP pair packets are highly nonuniform, unlike the unital lambda=1 frame.")

    print()
    print("Tournament-analysis carrier ranking:")
    print("  score=(LRC_faithfulness, tight_locus_specificity, incidence_strength, canonicality)")
    for name, score in carrier_tournament():
        print(f"  {name:36s} score={score}")
    print("Hamiltonian path:")
    print("  " + " > ".join(name for name, _ in carrier_tournament()))

    print()
    print("Conclusion:")
    print("  The q=3 unital is a promising weighted tight frame on 28 coordinates,")
    print("  but it is not a canonical AP8 pair-slot block design.  The missing")
    print("  proof obligation is a category-1 labelling/weighting map tied to the")
    print("  AP/Goddyn-Wong tight locus, not an S8-symmetric design identification.")


if __name__ == "__main__":
    main()
