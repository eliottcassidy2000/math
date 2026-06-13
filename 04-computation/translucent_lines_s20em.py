#!/usr/bin/env python3
"""
translucent_lines_s20em.py — Translucent lines: the G_{n-2} metagraph inside G_n
kind-pasteur-2026-03-23-S20em

Each tiling of delta_{n-2} (tournament on n vertices) decomposes as:
  boundary (2n-3 tiles) + overlap (C(n-2,2) tiles)

The C(n-2,2) overlap tiles form a sub-tournament on vertices {1,...,n-2}.
Flipping an overlap tile = flipping an arc in this sub-tournament.

TRANSLUCENT LINES: flips within the overlap
OPAQUE LINES: flips of boundary tiles (bottom, top, apex)

Questions:
1. How do translucent vs opaque flips distribute among self-loops vs edges?
2. Does the translucent subgraph reproduce G_{n-2} exactly?
3. What is the fiber structure: for fixed boundary, how many iso classes
   does the overlap range over?
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TRANSLUCENT LINES: THE G_{n-2} METAGRAPH INSIDE G_n")
print("  kind-pasteur-2026-03-23-S20em")
print("=" * 80)

for n in range(3, 8):
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    perms = list(permutations(range(n)))

    # Classify arcs into overlap vs boundary
    # Overlap: arcs between vertices {1, ..., n-2} (the inner n-2 vertices)
    # Bottom wiring: arcs from vertex 0 to {1, ..., n-2}
    # Top wiring: arcs from vertex n-1 to {1, ..., n-2}
    # Apex: arc (0, n-1)
    overlap_arcs = []
    bottom_arcs = []
    top_arcs = []
    apex_arcs = []

    for k, (i,j) in enumerate(ALL_ARCS):
        if i >= 1 and j <= n-2:
            overlap_arcs.append(k)
        elif i == 0 and j == n-1:
            apex_arcs.append(k)
        elif i == 0:
            bottom_arcs.append(k)
        elif j == n-1:
            top_arcs.append(k)
        else:
            # Should not happen
            print(f"  UNEXPECTED ARC: ({i},{j})")

    n_overlap = len(overlap_arcs)
    n_bottom = len(bottom_arcs)
    n_top = len(top_arcs)
    n_apex = len(apex_arcs)

    print(f"\n{'='*60}")
    print(f"  n = {n}, m = {m} arcs")
    print(f"  Overlap: {n_overlap} arcs (C({n-2},2) = {comb(n-2,2)})")
    print(f"  Bottom:  {n_bottom} arcs, Top: {n_top} arcs, Apex: {n_apex} arcs")
    print(f"  Boundary total: {n_bottom+n_top+n_apex} = 2*{n-2}+1 = {2*(n-2)+1}")
    print(f"{'='*60}")

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

    # Compute canonical forms for all tournaments
    canon_map = {}
    for bits in range(1 << m):
        canon_map[bits] = canon(bits)

    # For each tournament: classify each arc flip as translucent or opaque
    # and as self-loop or edge
    trans_sl = 0     # translucent self-loops (labeled)
    trans_edge = 0   # translucent edges (labeled)
    opaque_sl = 0    # opaque self-loops (labeled)
    opaque_edge = 0  # opaque edges (labeled)

    # Track unique edges by type
    trans_edge_pairs = set()
    opaque_edge_pairs = set()

    for bits in range(1 << m):
        cn1 = canon_map[bits]
        for k in range(m):
            flipped = bits ^ (1 << k)
            cn2 = canon_map[flipped]
            is_trans = k in overlap_arcs
            is_sl = (cn1 == cn2)

            if is_trans:
                if is_sl: trans_sl += 1
                else:
                    trans_edge += 1
                    pair = (min(cn1,cn2), max(cn1,cn2))
                    trans_edge_pairs.add(pair)
            else:
                if is_sl: opaque_sl += 1
                else:
                    opaque_edge += 1
                    pair = (min(cn1,cn2), max(cn1,cn2))
                    opaque_edge_pairs.add(pair)

    total_sl = trans_sl + opaque_sl
    total_edge = trans_edge + opaque_edge
    all_edge_pairs = trans_edge_pairs | opaque_edge_pairs

    # Known values
    E_known = {3:1, 4:5, 5:30, 6:290, 7:4086}
    V_known = {3:2, 4:4, 5:12, 6:56, 7:456}

    print(f"\n  LABELED FLIP CLASSIFICATION:")
    print(f"    Translucent self-loops:  {trans_sl:8d}  ({trans_sl/(2**m*m)*100:.1f}%)")
    print(f"    Translucent edges:       {trans_edge:8d}  ({trans_edge/(2**m*m)*100:.1f}%)")
    print(f"    Opaque self-loops:       {opaque_sl:8d}  ({opaque_sl/(2**m*m)*100:.1f}%)")
    print(f"    Opaque edges:            {opaque_edge:8d}  ({opaque_edge/(2**m*m)*100:.1f}%)")
    print(f"    Total:                   {2**m * m:8d}")

    print(f"\n  ORBIT-LEVEL:")
    print(f"    Total metagraph edges:       {len(all_edge_pairs):6d}  (E_known={E_known.get(n,'?')})")
    print(f"    Translucent-reachable edges:  {len(trans_edge_pairs):6d}  ({len(trans_edge_pairs)/len(all_edge_pairs)*100:.1f}% of E)")
    print(f"    Opaque-reachable edges:       {len(opaque_edge_pairs):6d}  ({len(opaque_edge_pairs)/len(all_edge_pairs)*100:.1f}% of E)")
    print(f"    Both (overlap):               {len(trans_edge_pairs & opaque_edge_pairs):6d}")
    print(f"    Trans-only edges:             {len(trans_edge_pairs - opaque_edge_pairs):6d}")
    print(f"    Opaque-only edges:            {len(opaque_edge_pairs - trans_edge_pairs):6d}")

    # Self-loop rates by type
    trans_total = trans_sl + trans_edge
    opaque_total = opaque_sl + opaque_edge
    print(f"\n  SELF-LOOP RATES:")
    print(f"    Translucent: {trans_sl}/{trans_total} = {trans_sl/trans_total:.4f}" if trans_total > 0 else "    Translucent: N/A")
    print(f"    Opaque:      {opaque_sl}/{opaque_total} = {opaque_sl/opaque_total:.4f}" if opaque_total > 0 else "    Opaque: N/A")

    # Fiber structure: for fixed boundary, how many iso classes?
    # Extract boundary bits and overlap bits
    boundary_indices = sorted(bottom_arcs + top_arcs + apex_arcs)
    overlap_indices = sorted(overlap_arcs)

    def get_boundary(bits):
        return tuple((bits >> k) & 1 for k in boundary_indices)

    def get_overlap(bits):
        return tuple((bits >> k) & 1 for k in overlap_indices)

    if n <= 6:  # Only compute fiber structure for small n
        boundary_to_classes = defaultdict(set)
        for bits in range(1 << m):
            b = get_boundary(bits)
            cn = canon_map[bits]
            boundary_to_classes[b].add(cn)

        fiber_sizes = [len(classes) for classes in boundary_to_classes.values()]
        V = V_known.get(n, len(set(canon_map.values())))
        V_sub = V_known.get(n-2, '?')

        print(f"\n  FIBER STRUCTURE (fixed boundary -> overlap iso classes):")
        print(f"    Number of distinct boundaries: {len(boundary_to_classes)}")
        print(f"    = 2^(2n-3) = 2^{2*(n-2)+1} = {2**(2*(n-2)+1)}")
        print(f"    Fiber size (classes per boundary): min={min(fiber_sizes)}, max={max(fiber_sizes)}, avg={sum(fiber_sizes)/len(fiber_sizes):.2f}")
        print(f"    V(n) = {V}, V(n-2) = {V_sub}")

        # How many of V(n-2) overlap classes are reachable from a typical boundary?
        # The overlap on {1,...,n-2} has V(n-2) iso classes.
        # But the full tournament class depends on the ENTIRE tiling, not just overlap.
        # Different overlaps with same boundary can give same or different full classes.

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

    if elapsed > 120:
        print("  Skipping larger n")
        break

print("\nDONE.")
print("=" * 80)
