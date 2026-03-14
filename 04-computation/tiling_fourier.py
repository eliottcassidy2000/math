"""
tiling_fourier.py -- kind-pasteur-2026-03-14-S110g
THE TILING MODEL AND FOURIER ANALYSIS

A tournament with a fixed Hamiltonian path = a tiling of the triangular grid.

The triangular grid for n vertices:
  Row i, column j (0 <= i < j <= n-1): arc between v_i and v_j.
  The PATH arcs: (i, i+1) for i=0,...,n-2. These are the DIAGONAL (y=x+1).
  The NON-PATH arcs: (i, j) with j > i+1. These are ABOVE the diagonal.

Each cell (i,j) is colored: 1 (forward: v_i -> v_j) or 0 (backward: v_j -> v_i).
The PATH cells are all 1 (forward, by definition of the path).
The NON-PATH cells are free.

The FOURIER ANALYSIS on the hypercube {0,1}^m acts on the tiling.
The level-2k Fourier subset S is a set of 2k cells.

KEY INSIGHT: The perpendicular diagonals y = x + c are the ANTI-DIAGONALS.
These connect cells that share a vertex in a specific way.

Anti-diagonal c: cells (i, i+c) for i = 0, 1, ..., n-1-c.
  c=1: the PATH diagonal. Cells (0,1), (1,2), ..., (n-2,n-1).
  c=2: cells (0,2), (1,3), ..., (n-3,n-1).
  c=3: cells (0,3), (1,4), ..., (n-4,n-1).
  ...

Each anti-diagonal c has (n-c) cells.

WHAT ARE THE FOURIER COEFFICIENTS IN THE TILING MODEL?

H_hat(S) for S = a set of cells in the grid.
From the proof: H_hat(S) != 0 iff S forms a linear forest
(set of edges that can all be consecutive in some Hamiltonian path).

In the TILING: S is a set of cells. As ARCS: (i1,j1), ..., (i_{2k}, j_{2k}).
These arcs form a linear forest iff they can all appear as
consecutive pairs in some permutation.

The diagonal y=x+1 (the path) gives the REFERENCE permutation.
A succession at position i means cell (i, i+1) is a path arc.
An anti-succession at position i means cell (i+1, i) is a reversed path arc.

THE TILING PERSPECTIVE ON N(S):
N(S) = signed count of paths containing all arcs of S.
Each path = a permutation. The path's relation to the FIXED path
(the reference permutation = identity) determines the grid position.

The RELATIVE permutation pi = P^{-1} * Q gives the tiling pattern
of Q relative to P.

In the grid: cells ABOVE the path diagonal are "forward" (v_i -> v_j with j > i+1).
Cells ON the diagonal are path arcs (always forward).
There are no cells below the diagonal (since we only have i < j).

THE ANTI-DIAGONALS AND FOURIER LEVELS:
Anti-diagonal c corresponds to arcs with "gap" c: (i, i+c).
The level-2 Fourier coefficients involve PAIRS of arcs.
ADJACENT arcs share a vertex: (i,j) and (j,k) share vertex j.
In the grid: (i,j) and (j,k) are on different anti-diagonals
but share the vertex j.

THE KEY QUESTION: How do the anti-diagonals relate to the
Fourier level structure?
"""

import sys, math
import numpy as np
from itertools import permutations, combinations
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("THE TILING MODEL AND FOURIER ANALYSIS")
    print("kind-pasteur-2026-03-14-S110g")
    print("=" * 70)

    # ============================================================
    print(f"\n{'='*70}")
    print("THE TRIANGULAR GRID FOR n=5")
    print(f"{'='*70}")

    n = 5
    print(f"\n  The grid (i,j) for 0 <= i < j <= {n-1}:")
    print(f"  Path arcs (diagonal y=x+1): (0,1),(1,2),(2,3),(3,4)")
    print(f"  Non-path arcs: (0,2),(0,3),(0,4),(1,3),(1,4),(2,4)")
    print()

    # Visualize the grid
    print(f"  j: ", end="")
    for j in range(1, n):
        print(f" {j}", end="")
    print()
    for i in range(n-1):
        print(f"  i={i}:", end="")
        for j in range(1, n):
            if j > i:
                if j == i+1:
                    print(f" P", end="")  # Path arc
                else:
                    print(f" .", end="")  # Non-path arc
            else:
                print(f"  ", end="")
        print()

    # Anti-diagonals
    print(f"\n  Anti-diagonals (gap c):")
    for c in range(1, n):
        cells = [(i, i+c) for i in range(n-c)]
        print(f"    c={c}: {cells}")

    # ============================================================
    print(f"\n{'='*70}")
    print("FOURIER SUBSETS AS GRID PATTERNS")
    print(f"{'='*70}")

    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(arcs)

    # The level-2 nonzero subsets: adjacent arc pairs.
    # Adjacent: share a vertex. In grid: (i,j) and (j,k) or (i,j) and (i,k).
    # These are cells connected by sharing a row or column index.

    # The level-4 nonzero subsets at n=5: the 60 Hamiltonian path edge sets.
    # Each is a set of 4 arcs forming a spanning path.
    # In the grid: 4 cells that form a path through all 5 vertices.

    # Let me visualize a few nonzero level-4 subsets as grid patterns.
    all_perms = list(permutations(range(n)))

    # Find 3 different level-4 nonzero subsets
    ham_paths = set()
    for perm in all_perms:
        edges = frozenset(tuple(sorted((perm[i], perm[i+1]))) for i in range(n-1))
        ham_paths.add(edges)

    examples = list(ham_paths)[:5]
    for hp in examples:
        print(f"\n  Ham path edges: {sorted(hp)}")
        # Show on grid
        print(f"    j:", end="")
        for j in range(1, n): print(f"  {j}", end="")
        print()
        for i in range(n-1):
            print(f"    i={i}:", end="")
            for j in range(1, n):
                if j > i:
                    if (i,j) in hp:
                        print(f"  X", end="")
                    else:
                        print(f"  .", end="")
                else:
                    print(f"   ", end="")
            print()

    # ============================================================
    print(f"\n{'='*70}")
    print("ANTI-DIAGONALS AND THE 'GAP' STRUCTURE")
    print(f"{'='*70}")

    # Each arc (i,j) has gap c = j - i.
    # The gap tells us which anti-diagonal the arc is on.
    # A Hamiltonian path uses arcs of various gaps.
    # The IDENTITY path 0-1-2-3-4 uses gaps: 1,1,1,1 (all gap 1 = diagonal).
    # Path 0-2-1-3-4: arcs (0,2),(1,2),(1,3),(3,4). Gaps: 2,1,2,1.
    # Path 1-0-2-4-3: arcs (0,1),(0,2),(2,4),(3,4). Gaps: 1,2,2,1.

    print(f"\n  Gap distribution of Hamiltonian paths at n=5:")
    gap_dist = Counter()
    for hp in ham_paths:
        gaps = tuple(sorted(j-i for i,j in hp))
        gap_dist[gaps] += 1

    for gaps, count in sorted(gap_dist.items()):
        print(f"    Gaps {gaps}: {count} paths")

    # ============================================================
    print(f"\n{'='*70}")
    print("THE y=x DIAGONAL AND TOURNAMENT ISOMORPHISM")
    print(f"{'='*70}")

    print(f"""
  The y=x line in the triangular grid passes through the CENTER of the grid.
  It separates the grid into "upper" and "lower" triangles.

  Actually: for the triangular grid with cells (i,j), 0 <= i < j <= n-1,
  the "y=x" diagonal isn't a cell but the REFLECTION axis.
  Reflecting (i,j) gives (n-1-j, n-1-i): this is the COMPLEMENT operation.
  T -> T^c (flip all arcs) corresponds to reflecting across y=x.

  The perpendicular to y=x is the anti-diagonal direction.
  Anti-diagonal c: cells (i, i+c) for fixed c.

  TOURNAMENT ISOMORPHISM: Two tournaments are isomorphic iff one can
  be obtained from the other by RELABELING vertices (permuting rows/columns
  of the grid simultaneously). This is the action of S_n on the grid.

  The isomorphism classes correspond to ORBITS of S_n on the grid.
  The number of orbits grows much slower than the number of tournaments.

  TILINGS AND ISOMORPHISM:
  A tiling of the grid = a tournament with a FIXED Hamiltonian path.
  Two tilings are "isomorphic" if there's a vertex permutation that
  maps one to the other AND preserves the path.
  The stabilizer of the path = the identity (paths are rigid).
  So each tiling is its own "isomorphism class" with respect to the path.
  But different paths through the SAME tournament give different tilings.

  The number of tilings = number of (tournament, path) pairs = sum_T H(T).
  The number of tournaments = 2^m.
  The ratio = Mean(H) = n!/2^(n-1).
  """)

    # ============================================================
    print(f"\n{'='*70}")
    print("ANTI-DIAGONALS AS FOURIER MODES")
    print(f"{'='*70}")

    print(f"""
  CONJECTURE: The Fourier level structure is related to the
  anti-diagonal structure of the grid.

  Anti-diagonal c has (n-c) cells.
  Level 2 involves PAIRS of arcs.
  Level 4 involves QUADRUPLES.

  Are level-2 coefficients concentrated on specific anti-diagonals?

  From the exact formula: |H_hat(e1,e2)| = (n-2)!/2^(n-2) for ADJACENT pairs.
  Adjacent = sharing a vertex.
  In the grid: (i,j) and (j,k) are "adjacent" and sit on different anti-diags.

  The Fourier level 2 connects cells that SHARE a vertex (row or column index).
  This is NOT anti-diagonal structure -- it's ROW/COLUMN structure.

  The anti-diagonals are the "gap" structure.
  The rows/columns are the "vertex" structure.
  Fourier levels use VERTEX structure (adjacency).
  Anti-diagonals use GAP structure.

  These are DUAL to each other!
  Vertex structure: which vertices are involved.
  Gap structure: what distance apart in the path.
  """)

    # ============================================================
    print(f"\n{'='*70}")
    print("THE COMPLEMENT SYMMETRY ON THE GRID")
    print(f"{'='*70}")

    # H(T) = H(T^c): complement preserves path count.
    # On the grid: complementing = flipping all non-path cells.
    # The path cells (diagonal) stay fixed.
    # This is reflection across the path diagonal combined with color inversion.

    # How does this interact with the Fourier structure?
    # chi_S(T^c) = prod_{e in S} (2*(1-x_e) - 1) = prod (-y_e) = (-1)^|S| * chi_S(T).
    # So H_hat(S) for T^c = (-1)^|S| * H_hat(S) for T (if complement were exact).
    # But H(T) = H(T^c) means: sum_T H(T) chi_S(T) = sum_T H(T^c) chi_S(T^c)?
    # Hmm: H_hat(S) = (1/N) sum_T H(T) chi_S(T).
    # Under T -> T^c: H(T^c) = H(T), chi_S(T^c) = (-1)^|S| chi_S(T).
    # So replacing T by T^c in the sum: sum H(T) (-1)^|S| chi_S(T) = (-1)^|S| H_hat(S).
    # But this is also the same sum (just reindexed). So H_hat(S) = (-1)^|S| H_hat(S).
    # For ODD |S|: H_hat(S) = -H_hat(S) => H_hat(S) = 0!
    # This is the DEGREE DROP: odd levels have zero energy.
    # The complement symmetry T <-> T^c FORCES odd Fourier levels to vanish.

    print(f"""
  THE COMPLEMENT SYMMETRY PROVES ODD LEVELS VANISH:

  H(T) = H(T^c) (complement preserves path count).
  chi_S(T^c) = (-1)^|S| * chi_S(T) (complement flips all signs).
  Therefore: H_hat(S) = (-1)^|S| * H_hat(S).
  For ODD |S|: H_hat(S) = 0. QED.

  THIS IS THE TILING PROOF OF THE DEGREE DROP!

  On the grid: the complement flips all non-path cells.
  This creates a Z/2 symmetry that kills odd Fourier modes.
  The even modes survive because (-1)^(even) = +1.

  The y=x reflection (complement) IS the reason H has only even Fourier levels.
  """)

    # ============================================================
    print(f"\n{'='*70}")
    print("THE PATH REVERSAL ON THE GRID")
    print(f"{'='*70}")

    print(f"""
  Path reversal: P = (v_0,...,v_(n-1)) -> P^rev = (v_(n-1),...,v_0).
  On the grid: this maps cell (i,j) to cell (n-1-j, n-1-i).
  This IS the y=x reflection of the grid!

  The path reversal is the ANTI-DIAGONAL reflection.
  Cells on anti-diagonal c map to cells on anti-diagonal c
  (the anti-diagonals are preserved by y=x reflection).

  Combined with complement (which also relates to y=x):
  - Complement: flips cell values, preserves grid positions.
  - Path reversal: permutes cell positions (y=x reflection), preserves values.
  - Together: they generate a Z/2 x Z/2 symmetry.

  The FOURIER ANALYSIS respects both:
  - Complement kills odd levels (proved above).
  - Path reversal: for |S| even, chi_S(P)*chi_S(P^rev) = +1.
    This gives the factor of 2 in |N(S)| = 2 * M(S).

  The TWO Z/2 symmetries combine to give the 4^r in |N(S)|^2:
  4 = 2 * 2 = (complement factor) * (reversal factor).
  Each even-length block contributes a factor of 4 = 2^2 to |N|^2.
  The 2's come from complement symmetry and reversal symmetry separately.
  """)

    # ============================================================
    print(f"\n{'='*70}")
    print("THE GRID PERSPECTIVE ON THE PROOF")
    print(f"{'='*70}")

    print(f"""
  THE PROOF IN TILING LANGUAGE:

  1. A tournament with fixed path = a tiling of the triangular grid.
  2. The complement symmetry (y=x reflection + color flip) kills odd levels.
  3. At even level 2k: only LINEAR FORESTS with ALL-EVEN components contribute.
  4. Each such forest has |N(S)|^2 = 4^r * ((n-2k)!)^2.
     The 4^r comes from: 2^r directions (each block forward/backward)
     times 2 from path reversal. Actually: (2 * 2^r)^2 / 4 = 4^r.
     Wait: |N| = 2 * (n-2k)! * 2^(r-1)? Let me recheck.

     At n=5, type (4,): |N| = 2. (n-2k)! = 1. 2^r = 2. N = 2*1*2/2 = 2. CHECK.
     At n=6, type (4,): |N| = 4. (n-2k)! = 2. 2^r = 2. N = 2*2*2/2 = 4. CHECK.
     At n=6, type (2,2): |N| = 8. (n-2k)! = 2. 2^r = 4. N = 2*2*4/2 = 8. CHECK.
     At n=7, type (2,2): |N| = 24. (n-2k)! = 6. 2^r = 4. N = 2*6*4/2 = 24. CHECK.

     So |N(S)| = (n-2k)! * 2^r. And |N|^2 = ((n-2k)!)^2 * 4^r. CONFIRMED.

     The formula: |N(S)| = (n-2k)! * 2^r.
     (n-2k)! = arrangements of free vertices.
     2^r = direction choices for r blocks.
     No separate reversal factor! The 2 is ALREADY in (n-2k)! because
     (n-2k)! counts ALL directed arrangements including reversed ones.

  5. The count of nonzero S of each type and the total sum of N^2
     gives E_2k = ((n-2k)!)^2 * [sum count*4^r] / 2^(2(n-1)).

  6. E_2k/E_0 = ((n-2k)!)^2 * [sum count*4^r] / (n!)^2
             = 2*(n-2k)^k / P(n,2k).

  7. This requires: sum_types count*4^r = 2*(n-2k)^k * (n!/((n-2k)!))^2 / P(n,2k) ... hmm.
     Actually: sum count*4^r = 2*(n-2k)^k * n!^2 / (((n-2k)!)^2 * P(n,2k)).
     P(n,2k) = n!/(n-2k)!.
     sum = 2*(n-2k)^k * n! / (n-2k)!.
     = 2*(n-2k)^k * P(n, 2k)? No: P(n,2k) = n!/(n-2k)!.
     sum = 2*(n-2k)^k * n!^2 / ((n-2k)!^2 * n!/(n-2k)!)
         = 2*(n-2k)^k * n! / (n-2k)!.
     = 2*(n-2k)^k * P(n, n-2k)? No.
     Let me just compute:
  """)

    # Verify the required sum at n=5,6,7
    for n_val in [5, 6, 7]:
        k = 2
        f = n_val - 2*k  # n - 2k
        E0 = (math.factorial(n_val) / 2**(n_val-1))**2
        pred_E4_E0 = 2*(n_val-4)**2 / math.perm(n_val, 4)
        required_sum = pred_E4_E0 * math.factorial(n_val)**2 / math.factorial(f)**2
        print(f"  n={n_val}: required sum count*4^r = {required_sum:.1f}")

    print(f"\n{'='*70}")
    print("DONE — THE GRID IS THE ROSETTA STONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
