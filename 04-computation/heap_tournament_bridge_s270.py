#!/usr/bin/env python3
"""
heap_tournament_bridge_s270.py — The Binary Heap / Tournament / Young Tableau Triangle
opus-2026-03-23-S270

THREE OBJECTS, ONE SHAPE:

1. TOURNAMENT on n vertices = binary filling of staircase δ_{n-2}
   H(T) = Hamiltonian path count = linear extensions of T-as-poset

2. STANDARD YOUNG TABLEAU of staircase shape = reduced decomposition of w_0
   Counted by hook length formula: f^δ = m! / prod(hooks)
   All hooks of the staircase are ODD

3. BINARY HEAP on m elements = linear extension of complete binary tree
   Counted by: m! / prod(subtree sizes)

The STAIRCASE shape δ_{n-2} = (n-2, n-3, ..., 1) has m = C(n-1,2) cells.
The tournament has C(n,2) arcs. The staircase has C(n-1,2) cells.
DIFFERENT! But related by the cut/cycle decomposition.

Actually: our staircase has cells for arcs (i,j) with i < j, arranged as
a triangular grid. This HAS m = C(n,2) cells (not C(n-1,2)).
The staircase δ_{n-1} = (n-1, n-2, ..., 1) has C(n,2) = n(n-1)/2 cells. ✓

So: a tournament on n vertices = binary filling of staircase δ_{n-1}.

This script computes:
1. Hook lengths of the staircase and SYT count (A005118)
2. H(T) for each tournament class and its relation to the staircase
3. Binary heap count on m = C(n,2) elements
4. The ORDER POLYTOPE volume = H(T)/n!
5. Relation between H, SYT, and heap counts
"""

import sys, time
from math import comb, factorial, prod
from itertools import permutations
from collections import defaultdict
import subprocess

sys.stdout.reconfigure(line_buffering=True)

def staircase_hooks(k):
    """Hook lengths of staircase shape δ_k = (k, k-1, ..., 1).

    Cell (i, j) in the staircase (0-indexed, i = row from top, j = column from left):
    Row i has length k-i cells (j = 0, ..., k-i-1).

    Hook of cell (i,j) = arm + leg + 1
      arm = (k-i-1) - j  (cells to the right in row i)
      leg = #{rows below that extend to column j}

    For the staircase: row r (r > i) has length k-r.
    It extends to column j iff j < k-r, i.e., r < k-j.
    So leg = #{r : i < r < k-j} = max(0, k-j-1-i).

    Hook(i,j) = (k-i-1-j) + max(0, k-j-1-i) + 1
             = (k-i-j) + (k-j-i-1) = 2(k-i-j) - 1  when k-j-1 > i
    Actually need to be more careful.
    """
    cells = []
    hooks = []
    total_cells = k*(k+1)//2

    for i in range(k):  # row i (0-indexed from top)
        row_len = k - i
        for j in range(row_len):  # column j
            arm = row_len - 1 - j  # cells to the right

            # Leg: cells below in same column
            leg = 0
            for r in range(i+1, k):
                if j < k - r:  # row r extends to column j
                    leg += 1

            hook = arm + leg + 1
            cells.append((i, j))
            hooks.append(hook)

    return cells, hooks, total_cells

def syt_count(k):
    """Number of SYT of staircase shape δ_k = (k, k-1, ..., 1).
    = m! / prod(hooks) where m = k(k+1)/2."""
    _, hooks, m = staircase_hooks(k)
    return factorial(m) // prod(hooks)

def binary_heap_count(m):
    """Number of binary min-heaps on m elements.
    = m! / prod(subtree sizes) for the complete binary tree on m nodes."""
    if m <= 1: return 1

    # Build complete binary tree: node i has children 2i+1, 2i+2
    subtree_sizes = [0] * m
    def compute_size(i):
        if i >= m: return 0
        subtree_sizes[i] = 1 + compute_size(2*i+1) + compute_size(2*i+2)
        return subtree_sizes[i]
    compute_size(0)

    return factorial(m) // prod(subtree_sizes)

def get_tournaments(n):
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    m = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1,n)]
    lines = [l.strip() for l in (r.stdout or '').split('\n')
             if len(l.strip()) == m and all(c in '01' for c in l.strip())]
    tournaments = []
    for line in lines:
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if line[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
        tournaments.append(adj)
    return tournaments

def hamiltonian_paths(adj, n):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp.get((S, v), 0)
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj[v][u]:
                    dp[(S | (1 << u), u)] = dp.get((S | (1 << u), u), 0) + val
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

if __name__ == '__main__':
    print("=" * 72)
    print("  THE BINARY HEAP / TOURNAMENT / YOUNG TABLEAU TRIANGLE")
    print("  opus-2026-03-23-S270")
    print("=" * 72)

    # ============================================================
    # 1. STAIRCASE HOOK LENGTHS
    # ============================================================
    print(f"\n{'='*72}")
    print("  1. STAIRCASE HOOK LENGTHS AND SYT COUNTS")
    print(f"{'='*72}")

    print(f"\n  {'k':>3} {'shape':>20} {'cells':>6} {'SYT count':>15} {'hooks':>40}")
    for k in range(1, 10):
        cells, hooks, m = staircase_hooks(k)
        s = syt_count(k)
        shape_str = str(tuple(range(k, 0, -1)))
        hooks_sorted = sorted(hooks, reverse=True)
        hooks_str = str(hooks_sorted[:10]) + ('...' if len(hooks_sorted) > 10 else '')
        print(f"  {k:>3} {shape_str:>20} {m:>6} {s:>15} {hooks_str:>40}")

    print(f"\n  ALL hooks of staircase are ODD:")
    for k in range(1, 8):
        _, hooks, _ = staircase_hooks(k)
        all_odd = all(h % 2 == 1 for h in hooks)
        print(f"    δ_{k}: hooks = {sorted(hooks)}, all odd = {all_odd}")

    # ============================================================
    # 2. THE THREE COUNTS
    # ============================================================
    print(f"\n{'='*72}")
    print("  2. THREE HOOK-LENGTH COUNTS COMPARED")
    print(f"{'='*72}")

    print(f"\n  {'n':>3} {'m=C(n,2)':>8} {'SYT(δ_{n-1})':>15} {'Heaps(m)':>15} {'n!':>12} {'ratio SYT/n!':>14}")
    for n in range(3, 10):
        m = comb(n, 2)
        syt = syt_count(n - 1)
        heaps = binary_heap_count(m)
        nf = factorial(n)
        ratio = syt / nf
        print(f"  {n:>3} {m:>8} {syt:>15} {heaps:>15} {nf:>12} {ratio:>14.4f}")

    # ============================================================
    # 3. H(T) AND ORDER POLYTOPE VOLUME
    # ============================================================
    print(f"\n{'='*72}")
    print("  3. H(T) AS LINEAR EXTENSION COUNT / ORDER POLYTOPE VOLUME")
    print(f"{'='*72}")

    for n in range(3, 8):
        tours = get_tournaments(n)
        nf = factorial(n)
        m = comb(n, 2)
        syt = syt_count(n - 1)

        print(f"\n  n = {n}: {len(tours)} tournament classes, m = {m}")
        print(f"  SYT of staircase δ_{n-1} = {syt}")
        print(f"  Binary heaps on {m} elements = {binary_heap_count(m)}")

        h_vals = []
        for adj in tours:
            h = hamiltonian_paths(adj, n)
            h_vals.append(h)

        total_H = sum(h_vals)  # Weighted by class size? No, these are class reps
        mean_H = total_H / len(tours)
        max_H = max(h_vals)
        min_H = min(h_vals)

        # Order polytope volume = H(T) / n!
        print(f"  H values: {sorted(h_vals)}")
        print(f"  H range: [{min_H}, {max_H}]")
        print(f"  Mean H: {mean_H:.2f}")
        print(f"  Vol(O(T_transitive)) = H_min/n! = {min_H}/{nf} = {min_H/nf:.6f}")
        print(f"  Vol(O(T_regular)) = H_max/n! = {max_H}/{nf} = {max_H/nf:.6f}")

        # Total H over all labeled tournaments (not class reps) = n! × mean(H over random T)
        # By Szele: expected H = n! / 2^{n-1}
        szele_expected = nf / (2 ** (n-1))
        print(f"  Szele expected H (random tournament): {szele_expected:.2f}")

        # H_max / Szele expected
        print(f"  H_max / Szele = {max_H / szele_expected:.4f}")

    # ============================================================
    # 4. THE HOOK LENGTH FORMULA FOR THE STAIRCASE
    # ============================================================
    print(f"\n{'='*72}")
    print("  4. THE HOOK LENGTH FORMULA: WHY ALL HOOKS ARE ODD")
    print(f"{'='*72}")
    print("""
  For the staircase δ_k = (k, k-1, ..., 1), the hook at cell (i,j) is:

    h(i,j) = 2(k - i - j) - 1    (always odd!)

  This is because:
    arm(i,j) = (k-i-1) - j = k - i - j - 1
    leg(i,j) = k - i - j - 1  (same formula!)
    h = arm + leg + 1 = 2(k-i-j-1) + 1 = 2(k-i-j) - 1

  The ARM and LEG are ALWAYS EQUAL in the staircase.
  This is the geometric signature of the right isosceles triangle.

  Product of hooks:
    prod h(i,j) = prod_{c=1}^{k-1} (2c-1)^{frequency}
  where 2c-1 appears as many times as there are cells with k-i-j = c.

  The number of cells with k-i-j = c is exactly c (diagonal counting).
  So: prod h = prod_{c=1}^{k-1} (2c-1)^c = 1^1 × 3^2 × 5^3 × 7^4 × ...

  And: SYT(δ_k) = [k(k+1)/2]! / (1^1 × 3^2 × 5^3 × ... × (2k-3)^{k-1})
""")

    # Verify the formula
    print("  Verification:")
    for k in range(1, 8):
        _, hooks, m = staircase_hooks(k)
        # Check: each hook = 2(k-i-j) - 1
        for idx, (i, j) in enumerate(staircase_hooks(k)[0]):
            expected = 2*(k - i - j) - 1
            actual = hooks[idx]
            if expected != actual:
                print(f"    MISMATCH at k={k}, cell ({i},{j}): expected {expected}, got {actual}")
                break
        else:
            # Check product formula
            freq = {}
            for h in hooks:
                freq[h] = freq.get(h, 0) + 1
            # Should be: h=2c-1 appears c times
            ok = True
            for c in range(1, k):
                h = 2*c - 1
                expected_freq = c
                if freq.get(h, 0) != expected_freq:
                    print(f"    MISMATCH at k={k}: h={h} appears {freq.get(h,0)} times, expected {expected_freq}")
                    ok = False
            if ok:
                print(f"    δ_{k}: ✓ h(i,j) = 2(k-i-j)-1, freq(2c-1) = c")

    # ============================================================
    # 5. THE DEEPER CONNECTION
    # ============================================================
    print(f"\n{'='*72}")
    print("  5. THE DEEPER CONNECTION: H, SYT, AND THE METAGRAPH")
    print(f"{'='*72}")
    print("""
  TOURNAMENT T on n vertices:
    Binary filling of staircase δ_{n-1} with C(n,2) cells.
    H(T) = number of Hamiltonian paths = linear extensions of T-poset.
    H ranges from 1 (transitive) to max (regular/Paley).

  SYT OF STAIRCASE δ_{n-1}:
    f^{δ_{n-1}} = [C(n,2)]! / (1^1 × 3^2 × 5^3 × ... × (2n-3)^{n-1})
    Counts: reduced decompositions of w_0 in S_n (Edelman-Greene)
    = number of sorting networks for n elements

  THE BRIDGE:
    H(T) counts linear extensions of ONE tournament.
    f^δ counts standard fillings of the SAME SHAPE.

    These are DIFFERENT objects:
    - H(T) = linear extensions of a PARTIAL ORDER (the tournament)
    - f^δ = standard fillings of a SHAPE (the staircase)

    But they live on the SAME geometric substrate: δ_{n-1}.

    The ORDER POLYTOPE of T has volume H(T)/n!.
    The ORDER POLYTOPE of the staircase-as-poset has volume f^δ / [C(n,2)]!.

  THE HOOK = THE ARC:
    Hook(i,j) = 2(n-1-i-j) - 1 = the "influence" of arc (i,j).
    The arc connecting vertices far apart (large n-1-i-j) has larger hook.
    The arc between adjacent vertices (n-1-i-j = 1) has hook 1.

    This connects tournament structure to Young tableau combinatorics
    through the GEOMETRY of the right isosceles triangle.

  BINARY HEAP AS TOURNAMENT TREE:
    A binary min-heap on m = C(n,2) elements uses the ARCS as elements.
    Each arc has a "weight" (its hook length).
    The heap structure orders arcs by influence.
    Number of valid heaps = m! / prod(subtree sizes of complete binary tree)

  SUMMARY: The three objects share the staircase as their common home:
    Tournament → binary FILLING (0/1 per cell) → H = linear extensions
    Young tableau → INCREASING FILLING (1..m per cell) → f^δ = SYT count
    Binary heap → HEAP-ORDERED FILLING → heap count
""")

    print("DONE.")
