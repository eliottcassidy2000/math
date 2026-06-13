#!/usr/bin/env python3
"""
wiggly_recursion_s20ep.py — How the wiggly fiber structure encodes the Mode B recursion
kind-pasteur-2026-03-23-S20ep

KEY QUESTION: The fiber map phi_B: overlap -> G_n/Z_2 doesn't factor through G_{n-2}.
But HOW MUCH information from G_{n-2} is preserved?

Specifically:
1. If two overlaps are in the same G_{n-2} class, how often are they in the same G_n class?
2. If two overlaps are connected by a wiggly line (arc flip) that's a G_{n-2} edge,
   is it always a G_n edge too? Or can a G_{n-2} edge become a G_n self-loop?
3. What's the "forgetting function": G_n class -> distribution over G_{n-2} classes?

Also: How does the BOUNDARY partition into types?
Two boundaries B1, B2 are "equivalent" if they produce the same fiber map
(up to relabeling the overlap). This equivalence relates to the stabilizer structure.
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WIGGLY RECURSION: Mode B through the fiber lens")
print("  kind-pasteur-2026-03-23-S20ep")
print("=" * 80)

for n in range(4, 7):
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    perms = list(permutations(range(n)))

    overlap_indices = []
    boundary_indices = []
    for k, (i,j) in enumerate(ALL_ARCS):
        if i >= 1 and j <= n-2:
            overlap_indices.append(k)
        else:
            boundary_indices.append(k)

    n_overlap = len(overlap_indices)
    n_boundary = len(boundary_indices)

    sub_arcs = [(i,j) for i in range(1, n-1) for j in range(i+1, n-1)]
    sub_perms = list(permutations(range(1, n-1)))

    def canon(bits):
        best = None
        for p in perms:
            nb = 0
            for k, (i,j) in enumerate(ALL_ARCS):
                pi, pj = p[i], p[j]
                if pi < pj:
                    nk = ALL_ARCS.index((pi, pj))
                    if bits & (1 << k): nb |= (1 << nk)
                else:
                    nk = ALL_ARCS.index((pj, pi))
                    if not (bits & (1 << k)): nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best

    def sub_canon(sub_bits):
        best = None
        for p in sub_perms:
            nb = 0
            for k, (i,j) in enumerate(sub_arcs):
                pi, pj = p[i-1], p[j-1]
                if pi < pj:
                    nk = sub_arcs.index((pi, pj))
                    if sub_bits & (1 << k): nb |= (1 << nk)
                else:
                    nk = sub_arcs.index((pj, pi))
                    if not (sub_bits & (1 << k)): nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best

    canon_map = {}
    for bits in range(1 << m):
        canon_map[bits] = canon(bits)

    # Merged classes
    comp_map = {}
    for cn in set(canon_map.values()):
        comp_bits = ((1 << m) - 1) ^ cn
        comp_map[cn] = canon(comp_bits)

    merged_map = {}
    for cn in set(canon_map.values()):
        merged_map[cn] = min(cn, comp_map[cn])

    merged_classes = set(merged_map.values())
    V_merged = len(merged_classes)

    # Sub-tournament canonical forms
    overlap_to_sub = {}
    for idx, k in enumerate(overlap_indices):
        arc = ALL_ARCS[k]
        if arc in sub_arcs:
            overlap_to_sub[idx] = sub_arcs.index(arc)

    def get_sub_bits(ob):
        sb = 0
        for idx in range(n_overlap):
            if ob & (1 << idx):
                sb |= (1 << overlap_to_sub.get(idx, idx))
        return sb

    sub_canon_map = {}
    for ob in range(1 << n_overlap):
        sub_canon_map[ob] = sub_canon(get_sub_bits(ob))

    sub_classes = sorted(set(sub_canon_map.values()))
    V_sub = len(sub_classes)

    # Sub-metagraph edges (G_{n-2})
    sub_edges = set()
    for ob in range(1 << n_overlap):
        sc1 = sub_canon_map[ob]
        for idx in range(n_overlap):
            ob2 = ob ^ (1 << idx)
            sc2 = sub_canon_map[ob2]
            if sc1 != sc2:
                sub_edges.add((min(sc1,sc2), max(sc1,sc2)))
    E_sub = len(sub_edges)

    print(f"\n{'#'*70}")
    print(f"  n = {n}")
    print(f"  V_merged(n) = {V_merged}, V(n-2) = {V_sub}, E(G_(n-2)) = {E_sub}")
    print(f"{'#'*70}")

    # ================================================================
    # Q1: If two overlaps are in same G_{n-2} class,
    #     how often same G_n class?
    # ================================================================
    same_sub_same_full = 0
    same_sub_diff_full = 0

    for bb in range(1 << n_boundary):
        base_bits = 0
        for idx, k in enumerate(boundary_indices):
            if bb & (1 << idx):
                base_bits |= (1 << k)

        # Group overlaps by sub-class
        sub_to_overlaps = defaultdict(list)
        for ob in range(1 << n_overlap):
            sc = sub_canon_map[ob]
            full_bits = base_bits
            for idx, k in enumerate(overlap_indices):
                if ob & (1 << idx):
                    full_bits |= (1 << k)
            mcn = merged_map[canon_map[full_bits]]
            sub_to_overlaps[sc].append(mcn)

        for sc, full_classes in sub_to_overlaps.items():
            for i in range(len(full_classes)):
                for j in range(i+1, len(full_classes)):
                    if full_classes[i] == full_classes[j]:
                        same_sub_same_full += 1
                    else:
                        same_sub_diff_full += 1

    total_same_sub = same_sub_same_full + same_sub_diff_full
    if total_same_sub > 0:
        print(f"\n  Q1: Same overlap class => same full class?")
        print(f"    Same sub AND same full: {same_sub_same_full} ({same_sub_same_full/total_same_sub*100:.1f}%)")
        print(f"    Same sub BUT diff full: {same_sub_diff_full} ({same_sub_diff_full/total_same_sub*100:.1f}%)")

    # ================================================================
    # Q2: G_{n-2} edge -> always G_n edge?
    #     Can a sub-tournament arc flip that changes the overlap class
    #     preserve the full class?
    # ================================================================
    sub_edge_becomes_gn_edge = 0
    sub_edge_becomes_gn_sl = 0
    sub_sl_becomes_gn_edge = 0
    sub_sl_becomes_gn_sl = 0

    for bb in range(1 << n_boundary):
        base_bits = 0
        for idx, k in enumerate(boundary_indices):
            if bb & (1 << idx):
                base_bits |= (1 << k)

        for ob in range(1 << n_overlap):
            sc1 = sub_canon_map[ob]
            full1 = base_bits
            for idx, k in enumerate(overlap_indices):
                if ob & (1 << idx):
                    full1 |= (1 << k)
            mcn1 = merged_map[canon_map[full1]]

            for idx in range(n_overlap):
                ob2 = ob ^ (1 << idx)
                sc2 = sub_canon_map[ob2]
                full2 = base_bits
                for idx2, k2 in enumerate(overlap_indices):
                    if ob2 & (1 << idx2):
                        full2 |= (1 << k2)
                mcn2 = merged_map[canon_map[full2]]

                sub_is_edge = (sc1 != sc2)
                full_is_edge = (mcn1 != mcn2)

                if sub_is_edge and full_is_edge: sub_edge_becomes_gn_edge += 1
                elif sub_is_edge and not full_is_edge: sub_edge_becomes_gn_sl += 1
                elif not sub_is_edge and full_is_edge: sub_sl_becomes_gn_edge += 1
                else: sub_sl_becomes_gn_sl += 1

    total_transitions = sub_edge_becomes_gn_edge + sub_edge_becomes_gn_sl + sub_sl_becomes_gn_edge + sub_sl_becomes_gn_sl

    print(f"\n  Q2: G_(n-2) edge -> G_n edge? (2x2 table)")
    print(f"    {'':>20} {'G_n edge':>12} {'G_n self-loop':>14}")
    print(f"    {'Sub edge:':>20} {sub_edge_becomes_gn_edge:12d} {sub_edge_becomes_gn_sl:14d}")
    print(f"    {'Sub self-loop:':>20} {sub_sl_becomes_gn_edge:12d} {sub_sl_becomes_gn_sl:14d}")

    sub_edge_total = sub_edge_becomes_gn_edge + sub_edge_becomes_gn_sl
    sub_sl_total = sub_sl_becomes_gn_edge + sub_sl_becomes_gn_sl
    if sub_edge_total > 0:
        print(f"\n    Sub edge -> G_n edge rate: {sub_edge_becomes_gn_edge/sub_edge_total*100:.1f}%")
    if sub_sl_total > 0:
        print(f"    Sub SL -> G_n edge rate:   {sub_sl_becomes_gn_edge/sub_sl_total*100:.1f}%")
    print(f"    Overall G_n edge rate:     {(sub_edge_becomes_gn_edge+sub_sl_becomes_gn_edge)/total_transitions*100:.1f}%")

    # ================================================================
    # Q3: Forgetting function G_n -> G_{n-2}
    #     For each G_n class: what distribution over G_{n-2} classes?
    # ================================================================
    gn_to_sub_dist = defaultdict(Counter)
    for bb in range(1 << n_boundary):
        base_bits = 0
        for idx, k in enumerate(boundary_indices):
            if bb & (1 << idx):
                base_bits |= (1 << k)
        for ob in range(1 << n_overlap):
            sc = sub_canon_map[ob]
            full_bits = base_bits
            for idx, k in enumerate(overlap_indices):
                if ob & (1 << idx):
                    full_bits |= (1 << k)
            mcn = merged_map[canon_map[full_bits]]
            gn_to_sub_dist[mcn][sc] += 1

    print(f"\n  Q3: Forgetting map G_n -> G_(n-2) (distribution of overlap classes per full class):")
    print(f"    {'G_n class':>12} {'Total tilings':>14} {'#sub classes':>12} {'Dominant sub':>14} {'Dominance':>10}")

    for mcn in sorted(gn_to_sub_dist.keys(), key=lambda c: sum(gn_to_sub_dist[c].values())):
        total = sum(gn_to_sub_dist[mcn].values())
        n_sub = len(gn_to_sub_dist[mcn])
        dominant = gn_to_sub_dist[mcn].most_common(1)[0]
        dom_frac = dominant[1] / total
        print(f"    {mcn:12d} {total:14d} {n_sub:12d} {dominant[0]:14d} {dom_frac:10.3f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\n" + "=" * 80)
print("SYNTHESIS")
print("=" * 80)
print("""
THE MODE B RECURSION THROUGH WIGGLY LINES:

1. Same overlap class does NOT guarantee same full class (non-factoring).
   The "extra information" in the boundary-overlap coupling grows with n.

2. The 2x2 table (sub edge vs G_n edge) reveals the MIXING rate:
   A sub-tournament edge (overlap class change) can become a G_n self-loop
   (if the boundary absorbs the change). Conversely, a sub self-loop
   (overlap class preserved) can become a G_n edge (if the boundary
   amplifies an intra-class relabeling into a class change).

3. The forgetting function G_n -> G_{n-2} shows how full classes
   distribute across overlap classes. Classes with more tilings
   tend to use more overlap classes (higher entropy).

4. The fiber bundle is TWISTED: the transition functions depend
   on the boundary. This twist is the obstruction to reducing
   G_n to a product G_{n-2} x boundary.
""")

print("DONE.")
print("=" * 80)
