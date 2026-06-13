"""
tiling_theorems.py -- kind-pasteur-2026-03-14-S77
Prove and verify structural theorems about the tiling model.

FROM THE COMPLETE ENUMERATION:
1. GS flip -> GS: ALWAYS (flip of GS tiling is GS)
2. Blueself = 0 at odd n, nonzero at even n
3. Blue line weights are always even
4. H-maximizer is blueself at even n
5. All blue lines connect classes of DIFFERENT t3 parity at odd n

PROOFS AND DEEPER ANALYSIS.
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("TILING STRUCTURE THEOREMS")
    print("kind-pasteur-2026-03-14-S77")
    print("=" * 70)

    # ========================================
    # THEOREM 1: GS flip preserves GS
    # ========================================
    print(f"\n{'='*70}")
    print("THEOREM 1: FLIP OF GS TILING IS ALWAYS GS")
    print(f"{'='*70}")
    print("""
  PROOF:
  A tiling t is GS iff t[i] = t[trans(i)] for all i where trans(i) != i.
  The flip of t is t' where t'[i] = 1 - t[i].
  Then t'[i] = 1 - t[i] = 1 - t[trans(i)] = t'[trans(i)].
  So t' is also GS. QED.

  Key: flip commutes with the GS symmetry because:
  - GS requires paired positions to have EQUAL bits
  - Flip complements ALL bits
  - Complementing equal bits gives equal bits
""")

    # ========================================
    # THEOREM 2: Blue line weights are always even
    # ========================================
    print(f"\n{'='*70}")
    print("THEOREM 2: BLUE LINE WEIGHTS ARE ALWAYS EVEN")
    print(f"{'='*70}")
    print("""
  CLAIM: The number of GS tilings in class A whose flip lands in class B
  is always even (for A != B).

  ARGUMENT: For a GS tiling t in class A with flip(t) in class B:
  - The GS-transpose of t, call it trans(t), is ALSO in class A
    (because GS + transpose = class-preserving for SC classes)
  - flip(trans(t)) = trans(flip(t)) (since flip commutes with trans)
  - trans(flip(t)) is in the SAME class as flip(t) = class B
  - So tilings t and trans(t) are paired, both contributing to the
    A->B count. Unless t = trans(t) (GS fixed point).

  At odd n: the GS map has floor((n-1)/2) fixed points on the grid.
  For paired positions, t and trans(t) are always DIFFERENT tilings.
  But a GS tiling CAN be fixed by an additional symmetry...

  Actually, this needs more careful analysis. Let me verify computationally.
""")

    # Verify from the n=5 and n=6 data
    # Re-run enumeration for n=5,6 and check all blue line weights
    for n in [3, 4, 5, 6]:
        # Quick enumeration
        tiles = []
        for b in range(1, n-1):
            for a in range(b+2, n+1):
                tiles.append((a, b))
        m = len(tiles)
        tile_idx = {t: i for i, t in enumerate(tiles)}

        trans_map = []
        for a, b in tiles:
            a2, b2 = n+1-b, n+1-a
            if a2 < b2:
                a2, b2 = b2, a2
            trans_map.append(tile_idx.get((a2, b2), -1))

        def is_gs(mask):
            for i in range(m):
                j = trans_map[i]
                if j != i and ((mask >> i) & 1) != ((mask >> j) & 1):
                    return False
            return True

        def bits_to_adj(mask):
            A = [[0]*n for _ in range(n)]
            for i in range(1, n):
                A[i][i-1] = 1
            for idx, (a, b) in enumerate(tiles):
                a0, b0 = a-1, b-1
                if (mask >> idx) & 1:
                    A[b0][a0] = 1
                else:
                    A[a0][b0] = 1
            return A

        perms = list(permutations(range(n)))
        def canonicalize(A):
            best = None
            for p in perms:
                s = ''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n))
                if best is None or s < best:
                    best = s
            return best

        # Build class map
        mask_to_class = {}
        for mask in range(1 << m):
            A = bits_to_adj(mask)
            mask_to_class[mask] = canonicalize(A)

        # Check all blue line weights
        flip_mask = (1 << m) - 1
        blue_weights = Counter()  # (classA, classB) -> count

        for mask in range(1 << m):
            if not is_gs(mask):
                continue
            classA = mask_to_class[mask]
            classB = mask_to_class[mask ^ flip_mask]
            if classA != classB:
                edge = tuple(sorted([classA, classB]))
                blue_weights[edge] += 1

        all_even = all(w % 2 == 0 for w in blue_weights.values())
        print(f"  n={n}: {len(blue_weights)} blue edges, all weights even? {all_even}")
        if not all_even:
            for e, w in blue_weights.items():
                if w % 2 != 0:
                    print(f"    ODD weight: {w}")

    # ========================================
    # THEOREM 3: Key formulas
    # ========================================
    print(f"\n{'='*70}")
    print("THEOREM 3: COUNTING FORMULAS")
    print(f"{'='*70}")

    # GS degrees of freedom
    print("\n  GS DOF formula:")
    for n in range(3, 12):
        m = (n-1) * (n-2) // 2
        # Fixed points of GS map: positions (r,c) where c = n-r-c, i.e., c = (n-r)/2
        # This requires n-r to be even. For r from 1 to n-2, n-r ranges from 2 to n-1.
        # Fixed iff n-r is even, i.e., r has same parity as n.
        fixed = sum(1 for r in range(1, n-1) if (n - r) % 2 == 0)
        paired = (m - fixed) // 2
        dof = fixed + paired
        print(f"    n={n}: m={m}, fixed={fixed}, paired={paired}, dof={dof}, #GS=2^{dof}={2**dof}")

    # Formula for fixed points:
    # r has same parity as n, r in [1, n-2]
    # If n even: r even, r in {2, 4, ..., n-2} -> (n-2)/2 values
    # If n odd: r odd, r in {1, 3, ..., n-2} -> (n-1)/2 values
    print("\n  Fixed point formula:")
    print("    n even: #fixed = (n-2)/2")
    print("    n odd:  #fixed = (n-1)/2")
    print("    General: #fixed = floor((n-1)/2)")

    for n in range(3, 12):
        predicted = (n-1) // 2
        actual = sum(1 for r in range(1, n-1) if (n - r) % 2 == 0)
        print(f"    n={n}: predicted={predicted}, actual={actual}, match={predicted == actual}")

    # Formula for GS DOF:
    # dof = fixed + paired = fixed + (m - fixed)/2 = (m + fixed)/2
    # = (C(n-1,2) + floor((n-1)/2)) / 2
    print("\n  GS DOF formula: dof = (C(n-1,2) + floor((n-1)/2)) / 2")
    for n in range(3, 12):
        m = (n-1) * (n-2) // 2
        f = (n-1) // 2
        dof = (m + f) // 2
        actual_dof = f + (m - f) // 2
        print(f"    n={n}: formula={(m+f)//2}, actual={actual_dof}, match={dof == actual_dof}")

    # ========================================
    # THEOREM 4: Certain properties of blue/black lines
    # ========================================
    print(f"\n{'='*70}")
    print("THEOREM 4: PROPERTIES WE CAN BE CERTAIN OF")
    print(f"{'='*70}")
    print("""
  PROVED:
  1. GS flip preserves GS (algebraic, Theorem 1)
  2. Blue line skeleton bipartite at odd n (THM-060, t3 parity)
  3. Blueself requires even n (THM-023, score argument)
  4. GS count = 2^{(m + floor((n-1)/2))/2} (combinatorial)
  5. Weight enumerator = (1+z)^f * (1+z^2)^p (S76 product code)
  6. GS DOF = (C(n-1,2) + floor((n-1)/2)) / 2

  VERIFIED (conjecture status):
  7. All blue line weights are even (verified n=3..6)
  8. H-maximizer is always in a blueself class at even n (verified n=4,6)
  9. At odd n: exactly 0 blueself tilings (verified n=3,5,7)
  10. The number of blackself tilings grows rapidly with n

  OPEN:
  11. Formula for #SC classes as function of n
  12. Formula for #blue edges as function of n
  13. Why is the H-maximizer blueself at even n?
  14. Structure of the black line skeleton
""")

    # Number of SC classes (from OEIS?)
    print("\n  SC class counts (self-converse tournaments):")
    sc_counts = {3: 2, 4: 2, 5: 8, 6: 12, 7: 88}
    for n, count in sorted(sc_counts.items()):
        print(f"    n={n}: {count} SC classes")

    # Total tournament isomorphism classes (OEIS A000568)
    total_counts = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
    print("\n  Total tournament classes (OEIS A000568):")
    for n, count in sorted(total_counts.items()):
        print(f"    n={n}: {count}")

    # NSC pair count = (total - SC) / 2
    for n in sorted(total_counts.keys()):
        nsc_pairs = (total_counts[n] - sc_counts.get(n, 0)) // 2
        print(f"    n={n}: NSC pairs = {nsc_pairs}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
