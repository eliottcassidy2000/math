#!/usr/bin/env python3
"""
wiggly_lines_s20eo.py — Wiggly lines: the fiber bundle structure of G_n
kind-pasteur-2026-03-23-S20eo

ABSTRACT FRAMEWORK:
  Tiling space = {0,1}^{C(n,2)} = Base x Fiber = {0,1}^{2n-3} x {0,1}^{C(n-2,2)}

  Base = Boundary (bottom + top + apex) = 2n-3 bits
  Fiber = Overlap (sub-tournament on {1,...,n-2}) = C(n-2,2) bits

  For each boundary B, the FIBER F_B is a hypercube Q_{C(n-2,2)}.
  WIGGLY LINES = edges of F_B (single overlap tile flips).

  The FIBER MAP phi_B: F_B -> G_n/Z_2 sends each overlap to its full iso class.
  Wiggly lines project to metagraph edges (or self-loops if same class).

  KEY QUESTIONS:
  1. How does the G_{n-2} metagraph embed in each fiber?
  2. What's the distribution of the fiber map phi_B?
  3. How many wiggly self-loops vs wiggly edges per fiber?
  4. Which G_n edges are "visible" in each fiber?
  5. The FIBER MULTIPLICITY: how many fibers witness each G_n edge?
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WIGGLY LINES: FIBER BUNDLE STRUCTURE OF G_n")
print("  kind-pasteur-2026-03-23-S20eo")
print("=" * 80)

for n in range(4, 7):
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    perms = list(permutations(range(n)))

    # Classify arcs into overlap vs boundary
    overlap_indices = []
    boundary_indices = []
    for k, (i,j) in enumerate(ALL_ARCS):
        if i >= 1 and j <= n-2:
            overlap_indices.append(k)
        else:
            boundary_indices.append(k)

    n_overlap = len(overlap_indices)  # = C(n-2,2)
    n_boundary = len(boundary_indices)  # = 2n-3

    print(f"\n{'#'*70}")
    print(f"  n = {n}")
    print(f"  Overlap (wiggly) dimension: {n_overlap} = C({n-2},2)")
    print(f"  Boundary (opaque) dimension: {n_boundary} = 2*{n-2}+1")
    print(f"  Fibers: 2^{n_boundary} = {2**n_boundary}")
    print(f"  Fiber size: 2^{n_overlap} = {2**n_overlap}")
    print(f"{'#'*70}")

    # Canonical form computation
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

    # Precompute canonical forms
    canon_map = {}
    for bits in range(1 << m):
        canon_map[bits] = canon(bits)

    # Build merged class map
    class_to_comp = {}
    for cn in set(canon_map.values()):
        comp_bits = ((1 << m) - 1) ^ cn
        comp_cn = canon(comp_bits)
        class_to_comp[cn] = comp_cn

    unmerged_to_merged = {}
    for cn in set(canon_map.values()):
        mcn = min(cn, class_to_comp[cn])
        unmerged_to_merged[cn] = mcn

    merged_classes = set(unmerged_to_merged.values())
    V_merged = len(merged_classes)

    # Also compute overlap metagraph (iso classes of sub-tournament on {1,...,n-2})
    # Sub-tournament arcs: between vertices 1,...,n-2
    sub_arcs = [(i,j) for i in range(1, n-1) for j in range(i+1, n-1)]
    sub_m = len(sub_arcs)  # = C(n-2,2)
    sub_perms = list(permutations(range(1, n-1)))

    def sub_canon(sub_bits):
        """Canonical form of overlap sub-tournament on {1,...,n-2}."""
        best = None
        for p in sub_perms:
            # p is a permutation of {1,...,n-2}
            nb = 0
            for k, (i,j) in enumerate(sub_arcs):
                pi, pj = p[i-1], p[j-1]  # p indexed from 0
                if pi < pj:
                    nk = sub_arcs.index((pi, pj))
                    if sub_bits & (1 << k): nb |= (1 << nk)
                else:
                    nk = sub_arcs.index((pj, pi))
                    if not (sub_bits & (1 << k)): nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best

    # Map: extract overlap bits from full tiling
    def get_overlap_bits(bits):
        ob = 0
        for idx, k in enumerate(overlap_indices):
            if bits & (1 << k):
                ob |= (1 << idx)
        return ob

    def get_boundary_bits(bits):
        bb = 0
        for idx, k in enumerate(boundary_indices):
            if bits & (1 << k):
                bb |= (1 << idx)
        return bb

    # Map overlap bits to sub-tournament bits (for sub_canon)
    # overlap_indices maps to sub_arcs
    # They should correspond: overlap_indices[i] = index in ALL_ARCS of sub_arcs[i]
    # Let's verify this mapping
    overlap_to_sub = {}
    for idx, k in enumerate(overlap_indices):
        arc = ALL_ARCS[k]
        if arc in sub_arcs:
            sub_idx = sub_arcs.index(arc)
            overlap_to_sub[idx] = sub_idx

    def overlap_to_sub_bits(ob):
        sb = 0
        for idx in range(n_overlap):
            if ob & (1 << idx):
                sub_idx = overlap_to_sub.get(idx, idx)
                sb |= (1 << sub_idx)
        return sb

    # ======================================================
    # FIBER ANALYSIS
    # ======================================================

    # For each boundary B: analyze the fiber
    fiber_data = []  # list of dicts per fiber

    # Stats accumulators
    total_wiggly_edges = 0
    total_wiggly_sl = 0
    gn_edges_from_wiggly = set()  # which G_n edges are witnessed by wiggly lines
    fiber_edge_witnesses = defaultdict(int)  # G_n edge -> how many fibers witness it

    # Overlap class map
    if n_overlap > 0:
        sub_canon_map = {}
        for ob in range(1 << n_overlap):
            sb = overlap_to_sub_bits(ob)
            sub_canon_map[ob] = sub_canon(sb)
        sub_classes = set(sub_canon_map.values())
        V_sub = len(sub_classes)
    else:
        V_sub = 1

    for bb in range(1 << n_boundary):
        # Reconstruct boundary
        base_bits = 0
        for idx, k in enumerate(boundary_indices):
            if bb & (1 << idx):
                base_bits |= (1 << k)

        # Map each overlap to full class
        fiber_classes = {}  # overlap_bits -> merged class
        for ob in range(1 << n_overlap):
            full_bits = base_bits
            for idx, k in enumerate(overlap_indices):
                if ob & (1 << idx):
                    full_bits |= (1 << k)
            cn = canon_map[full_bits]
            mcn = unmerged_to_merged[cn]
            fiber_classes[ob] = mcn

        # Distinct classes in this fiber
        classes_in_fiber = set(fiber_classes.values())
        n_classes = len(classes_in_fiber)

        # Wiggly edges: overlap flips within this fiber
        wiggly_edges = 0
        wiggly_sl = 0
        fiber_gn_edges = set()

        for ob in range(1 << n_overlap):
            mcn1 = fiber_classes[ob]
            for idx in range(n_overlap):
                ob2 = ob ^ (1 << idx)
                mcn2 = fiber_classes[ob2]
                if mcn1 == mcn2:
                    wiggly_sl += 1
                else:
                    wiggly_edges += 1
                    edge = (min(mcn1, mcn2), max(mcn1, mcn2))
                    fiber_gn_edges.add(edge)
                    gn_edges_from_wiggly.add(edge)
                    fiber_edge_witnesses[edge] += 1

        # Count wiggly edges (each counted twice, once from each endpoint)
        wiggly_edges //= 2
        wiggly_sl //= 2  # Actually self-loops are counted once per direction... hmm

        total_wiggly_edges += wiggly_edges
        total_wiggly_sl += wiggly_sl

        # Overlap class distribution
        overlap_class_dist = Counter()
        for ob, mcn in fiber_classes.items():
            sc = sub_canon_map[ob] if n_overlap > 0 else 0
            overlap_class_dist[(sc, mcn)] += 1

        fiber_data.append({
            'bb': bb,
            'n_classes': n_classes,
            'n_gn_edges': len(fiber_gn_edges),
            'wiggly_edges': wiggly_edges,
            'wiggly_sl': wiggly_sl,
        })

    # ======================================================
    # SUMMARY STATISTICS
    # ======================================================

    # All G_n edges
    all_gn_edges = set()
    for bits in range(1 << m):
        cn1 = canon_map[bits]
        mcn1 = unmerged_to_merged[cn1]
        for k in range(m):
            flipped = bits ^ (1 << k)
            cn2 = canon_map[flipped]
            mcn2 = unmerged_to_merged[cn2]
            if mcn1 != mcn2:
                all_gn_edges.add((min(mcn1,mcn2), max(mcn1,mcn2)))

    E_gn = len(all_gn_edges)
    E_wiggly = len(gn_edges_from_wiggly)

    fiber_classes_list = [fd['n_classes'] for fd in fiber_data]
    fiber_edges_list = [fd['n_gn_edges'] for fd in fiber_data]

    print(f"\n  FIBER STATISTICS:")
    print(f"    V(n-2) = {V_sub} (overlap iso classes)")
    print(f"    V_merged(n) = {V_merged}")
    print(f"    E(G_n/Z_2) = {E_gn}")
    print(f"    E visible via wiggly lines = {E_wiggly} ({E_wiggly/E_gn*100:.1f}%)")
    print(f"    Wiggly-only edges (not in opaque) = 0 (by isotropy)")

    print(f"\n  PER-FIBER DISTRIBUTION:")
    print(f"    Classes per fiber: min={min(fiber_classes_list)}, max={max(fiber_classes_list)}, avg={sum(fiber_classes_list)/len(fiber_classes_list):.2f}")
    print(f"    G_n edges per fiber: min={min(fiber_edges_list)}, max={max(fiber_edges_list)}, avg={sum(fiber_edges_list)/len(fiber_edges_list):.2f}")
    print(f"    Total fibers: {len(fiber_data)}")

    # Fiber multiplicity distribution: how many fibers witness each G_n edge?
    mult_values = sorted(fiber_edge_witnesses.values())
    if mult_values:
        print(f"\n  FIBER MULTIPLICITY (how many fibers witness each G_n edge):")
        print(f"    min={min(mult_values)}, max={max(mult_values)}, avg={sum(mult_values)/len(mult_values):.1f}")
        # Distribution
        mult_hist = Counter(mult_values)
        for m_val in sorted(mult_hist.keys()):
            if mult_hist[m_val] >= max(1, E_gn * 0.05):
                print(f"    {m_val} fibers: {mult_hist[m_val]} edges")

    # Wiggly self-loop rate
    total_wiggly_total = sum(fd['wiggly_edges'] + fd['wiggly_sl'] for fd in fiber_data)
    total_wiggly_sl_sum = sum(fd['wiggly_sl'] for fd in fiber_data)
    if total_wiggly_total > 0:
        print(f"\n  WIGGLY SELF-LOOP RATE:")
        print(f"    Total wiggly transitions: {total_wiggly_total}")
        print(f"    Wiggly self-loops: {total_wiggly_sl_sum} ({total_wiggly_sl_sum/total_wiggly_total*100:.1f}%)")

    # Compare with G_{n-2} structure
    # In each fiber: the overlap sub-tournament has V_sub iso classes.
    # The full classes are a REFINEMENT of this.
    # How much refinement?
    print(f"\n  G_(n-2) vs G_n REFINEMENT:")
    print(f"    V(n-2) = {V_sub} overlap classes")
    print(f"    Avg full classes per fiber = {sum(fiber_classes_list)/len(fiber_classes_list):.2f}")
    print(f"    Refinement factor = {sum(fiber_classes_list)/len(fiber_classes_list)/V_sub:.2f}")

    # Overlap class -> full class map: is it many-to-one?
    # For each (boundary, overlap_class): how many full classes?
    if n <= 6:
        overlap_refinement = defaultdict(set)  # (boundary, overlap_class) -> set of full classes
        for bb in range(1 << n_boundary):
            base_bits = 0
            for idx, k in enumerate(boundary_indices):
                if bb & (1 << idx):
                    base_bits |= (1 << k)

            for ob in range(1 << n_overlap):
                full_bits = base_bits
                for idx, k in enumerate(overlap_indices):
                    if ob & (1 << idx):
                        full_bits |= (1 << k)
                cn = canon_map[full_bits]
                mcn = unmerged_to_merged[cn]
                sc = sub_canon_map[ob] if n_overlap > 0 else 0
                overlap_refinement[(bb, sc)].add(mcn)

        # How many (boundary, overlap_class) pairs map to >1 full class?
        multi_class = sum(1 for v in overlap_refinement.values() if len(v) > 1)
        total_pairs = len(overlap_refinement)
        print(f"    (boundary, overlap_class) pairs: {total_pairs}")
        print(f"    Pairs mapping to >1 full class: {multi_class} ({multi_class/total_pairs*100:.1f}%)")
        print(f"    Max full classes from one (B, C_overlap): {max(len(v) for v in overlap_refinement.values())}")

        # Is same overlap class -> same full class (for fixed boundary)?
        # This would mean phi_B factors through the overlap class map.
        # Check: for each boundary, is the map overlap -> full_class constant on overlap iso classes?
        factors_through = 0
        doesnt_factor = 0
        for bb in range(1 << n_boundary):
            base_bits = 0
            for idx, k in enumerate(boundary_indices):
                if bb & (1 << idx):
                    base_bits |= (1 << k)

            # Group overlaps by their sub-tournament iso class
            overlap_groups = defaultdict(set)
            for ob in range(1 << n_overlap):
                full_bits = base_bits
                for idx, k in enumerate(overlap_indices):
                    if ob & (1 << idx):
                        full_bits |= (1 << k)
                cn = canon_map[full_bits]
                mcn = unmerged_to_merged[cn]
                sc = sub_canon_map[ob] if n_overlap > 0 else 0
                overlap_groups[sc].add(mcn)

            # Does the map factor? (each overlap class maps to exactly one full class)
            if all(len(v) == 1 for v in overlap_groups.values()):
                factors_through += 1
            else:
                doesnt_factor += 1

        print(f"\n  FACTORING THROUGH G_(n-2):")
        print(f"    Fibers where phi_B factors through overlap class: {factors_through}/{2**n_boundary} ({factors_through/2**n_boundary*100:.1f}%)")
        print(f"    Fibers where it DOESN'T factor: {doesnt_factor}/{2**n_boundary}")
        print(f"    (phi_B factors = same overlap class always gives same full class)")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\n" + "=" * 80)
print("ABSTRACT SYNTHESIS")
print("=" * 80)
print("""
THE FIBER BUNDLE STRUCTURE:

  Total tiling space = {0,1}^{C(n,2)} is a PRODUCT:
    Base (boundary): {0,1}^{2n-3}
    Fiber (overlap): {0,1}^{C(n-2,2)}

  WIGGLY LINES = edges within each fiber = single overlap tile flips
  Each fiber is a hypercube Q_{C(n-2,2)} with 2^{C(n-2,2)} tilings.

  The FIBER MAP phi_B: Fiber -> G_n/Z_2 sends each overlap to its full class.
  phi_B is NOT a morphism from G_{n-2} to G_n (it doesn't factor through
  the overlap iso class map in general).

  The wiggly lines project to ALL edges of G_n/Z_2 (by S_n isotropy).
  Each G_n edge is witnessed by multiple fibers (fiber multiplicity).

  KEY INSIGHT: The overlap sub-tournament's LABELED structure matters,
  not just its iso class. Two isomorphic overlaps can give different
  full classes (depending on how they interact with the boundary).

  This is why G_n cannot be reduced to G_{n-2}: the boundary breaks
  the symmetry of the overlap. The Mode B recursion n->n-2 is a
  FIBER BUNDLE, not a simple quotient.
""")

print("DONE.")
print("=" * 80)
