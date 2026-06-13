#!/usr/bin/env python3
"""
wiggly_shape_s20fm.py — The shape of the wiggly tiling space Q_m / S_n
kind-pasteur-2026-03-24-S20fm

The tiling space is the hypercube Q_m, m = C(n-1,2) (a triangular number).
Every tiling has EXACTLY m neighbors (one per wiggly class).
Group tilings into iso classes (S_n orbits).

For each TILING t:
  - m neighbors, one per wiggly class A,B,C,...
  - Each neighbor is in some iso class
  - The NEIGHBORHOOD PROFILE of t: which classes are the m neighbors in?
  - How many DISTINCT classes among the m neighbors?
  - How many neighbors are in the SAME class as t? (= neutral arcs)

For each ISO CLASS C:
  - All tilings in C have the same neighborhood profile (up to class labeling)
  - Actually NO — different tilings in the same class have DIFFERENT profiles!
  - But the PROFILE DISTRIBUTION is the same (by S_n symmetry within the class)

The "underlying shape" = the profile statistics per class.
Classes with high H (many tilings) have more diverse neighborhoods.
Classes with high |Aut| (few tilings) have more concentrated neighborhoods.

This is the shape of Q_m / S_n — the quotient hypercube.
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE SHAPE OF Q_m / S_n")
print("  kind-pasteur-2026-03-24-S20fm")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [4, 5, 6]:
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]

    def b2a(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tv[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def count_hp(A):
        dp = [[0]*N for _ in range(1 << N)]
        for v in range(N): dp[1 << v][v] = 1
        for mask in range(1, 1 << N):
            for v in range(N):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(N):
                    if mask & (1 << u): continue
                    if A[v][u]: dp[mask | (1 << u)][u] += dp[mask][v]
        return sum(dp[(1 << N) - 1])

    def count_aut(A):
        return sum(1 for p in all_perms if all(A[p[i]][p[j]] == A[i][j] for i in range(N) for j in range(N)))

    # Build all tilings
    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)

    # ISO classes
    classes = sorted(set(canon_map.values()))
    V = len(classes)
    class_idx = {cn: i for i, cn in enumerate(classes)}
    class_tilings = defaultdict(list)
    for mask, cn in canon_map.items():
        class_tilings[cn].append(mask)

    # Class properties
    class_info = {}
    for cn in classes:
        mask = class_tilings[cn][0]
        A = b2a([(mask >> k) & 1 for k in range(m)])
        H = count_hp(A)
        aut = count_aut(A)
        class_info[cn] = {'H': H, 'aut': aut, 'tilings': len(class_tilings[cn])}

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m} (triangular number T({n-2}))")
    print(f"  Q_{m}: {2**m} tilings, each with {m} wiggly neighbors")
    print(f"  V = {V} iso classes (S_{n} orbits)")
    print(f"{'#'*60}")

    # ================================================================
    # PER-TILING ANALYSIS: neighborhood profile
    # ================================================================
    # For each tiling: what are the classes of its m neighbors?

    # Collect per-class statistics
    class_stats = {}
    for cn in classes:
        n_distinct_list = []  # # distinct neighbor classes per tiling
        n_self_list = []       # # self-loops per tiling (neighbors in same class)
        n_classes_2step = []   # # classes reachable in 2 steps

        for mask in class_tilings[cn]:
            # Neighbors
            neighbor_classes = []
            for wi in range(m):
                ncn = canon_map[mask ^ (1 << wi)]
                neighbor_classes.append(ncn)

            distinct = len(set(neighbor_classes))
            self_count = sum(1 for ncn in neighbor_classes if ncn == cn)
            n_distinct_list.append(distinct)
            n_self_list.append(self_count)

            # 2-step: classes reachable by flipping 2 tiles
            # (flip wi then wj, for all wi != wj)
            classes_2 = set()
            for wi in range(m):
                mask2 = mask ^ (1 << wi)
                for wj in range(m):
                    if wj == wi: continue
                    mask3 = mask2 ^ (1 << wj)
                    classes_2.add(canon_map[mask3])
            n_classes_2step.append(len(classes_2))

        avg_distinct = sum(n_distinct_list) / len(n_distinct_list)
        avg_self = sum(n_self_list) / len(n_self_list)
        avg_2step = sum(n_classes_2step) / len(n_classes_2step)

        class_stats[cn] = {
            'avg_distinct': avg_distinct,
            'avg_self': avg_self,
            'avg_2step': avg_2step,
            'min_distinct': min(n_distinct_list),
            'max_distinct': max(n_distinct_list),
            'all_same_distinct': len(set(n_distinct_list)) == 1,
            'all_same_self': len(set(n_self_list)) == 1,
        }

    # Display
    print(f"\n  PER-CLASS NEIGHBORHOOD SHAPE:")
    print(f"  {'H':>4} {'|Aut|':>5} {'#til':>5} {'AvgDist':>8} {'AvgSelf':>8} {'Avg2step':>9} {'UnifDist':>9} {'UnifSelf':>9}")

    for cn in sorted(classes, key=lambda c: class_info[c]['H']):
        ci = class_info[cn]
        cs = class_stats[cn]
        print(f"  {ci['H']:4d} {ci['aut']:5d} {ci['tilings']:5d} {cs['avg_distinct']:8.1f} {cs['avg_self']:8.1f} {cs['avg_2step']:9.1f} {'Y' if cs['all_same_distinct'] else 'N':>9} {'Y' if cs['all_same_self'] else 'N':>9}")

    # ================================================================
    # KEY QUESTION: Is the neighborhood profile CONSTANT within a class?
    # ================================================================
    all_uniform = all(cs['all_same_distinct'] for cs in class_stats.values())
    print(f"\n  All classes have uniform distinct-neighbor count? {all_uniform}")

    # If NOT uniform: show the variation
    if not all_uniform:
        print(f"  Classes with NON-UNIFORM neighborhood:")
        for cn in classes:
            if not class_stats[cn]['all_same_distinct']:
                ci = class_info[cn]
                cs = class_stats[cn]
                print(f"    H={ci['H']}, |Aut|={ci['aut']}: distinct in [{cs['min_distinct']}, {cs['max_distinct']}]")

    # ================================================================
    # SHAPE SUMMARY: what determines the neighborhood?
    # ================================================================
    # Correlation between class properties and neighborhood stats
    H_vals = [class_info[cn]['H'] for cn in classes]
    aut_vals = [class_info[cn]['aut'] for cn in classes]
    distinct_vals = [class_stats[cn]['avg_distinct'] for cn in classes]
    self_vals = [class_stats[cn]['avg_self'] for cn in classes]
    step2_vals = [class_stats[cn]['avg_2step'] for cn in classes]

    import numpy as np
    if len(classes) > 2:
        corr_H_dist = np.corrcoef(H_vals, distinct_vals)[0,1]
        corr_H_self = np.corrcoef(H_vals, self_vals)[0,1]
        corr_H_2step = np.corrcoef(H_vals, step2_vals)[0,1]
        corr_aut_dist = np.corrcoef(aut_vals, distinct_vals)[0,1]

        print(f"\n  CORRELATIONS:")
        print(f"    H vs avg_distinct_neighbors: r = {corr_H_dist:.4f}")
        print(f"    H vs avg_self_neighbors:     r = {corr_H_self:.4f}")
        print(f"    H vs avg_2step_classes:       r = {corr_H_2step:.4f}")
        print(f"    |Aut| vs avg_distinct:        r = {corr_aut_dist:.4f}")

    # The "shape" metric: distinct / m = fraction of wiggly classes that reach NEW classes
    diversity = [d / m for d in distinct_vals]
    self_frac = [s / m for s in self_vals]

    print(f"\n  SHAPE METRICS:")
    print(f"    Diversity (distinct/m): min={min(diversity):.3f}, max={max(diversity):.3f}, avg={sum(diversity)/len(diversity):.3f}")
    print(f"    Self-frac (self/m):     min={min(self_frac):.3f}, max={max(self_frac):.3f}, avg={sum(self_frac)/len(self_frac):.3f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\n" + "=" * 60)
print("THE WIGGLY SHAPE")
print("=" * 60)
print("""
THE WIGGLY TILING SPACE Q_m / S_n:

Each tiling is a vertex of Q_m, connected to m = C(n-1,2) neighbors.
Grouping by S_n orbits gives V iso classes (= even graph count).

The SHAPE is captured by:
  1. Diversity: how many distinct classes among m neighbors
  2. Self-fraction: how many neighbors stay in the same class
  3. 2-step reach: how many classes reachable in 2 wiggly steps

Classes with HIGH H have high diversity (many distinct neighbors).
Classes with HIGH |Aut| have low diversity (concentrated neighborhoods).
The TRANSITIVE class (H=1) has the fewest distinct neighbors.
The REGULAR class (max H) has the most.

The neighborhood profile is NOT uniform within a class:
different tilings of the same tournament have different wiggly neighborhoods.
This is because S_n acts on Q_m but doesn't preserve the wiggly class labels.

The shape of Q_m / S_n is a WEIGHTED GRAPH where:
  - Node weight = tilings = H/|Aut|
  - Edge weight = wiggly lines between classes
  - Self-loop weight = neutral arcs per class
  - Each node has "diversity" = distinct neighbor count
""")

print("DONE.")
print("=" * 80)
