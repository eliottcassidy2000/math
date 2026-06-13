"""
quiver_dynkin_deep.py -- kind-pasteur-2026-03-14-S68

Corrected root enumeration + deep analysis of:
1. Gabriel numbers = |Phi+| and their connection to forbidden H
2. McKay quiver dimension vectors
3. Quiver representation theory meets tournament invariants
4. Information-theoretic channel interpretation
"""

import numpy as np
from collections import Counter
import math

def cartan_matrix(dtype, rank):
    """Build Cartan matrix for simply-laced Dynkin diagrams."""
    C = np.eye(rank, dtype=int) * 2
    if dtype == 'A':
        for i in range(rank - 1):
            C[i][i+1] = -1; C[i+1][i] = -1
    elif dtype == 'D':
        for i in range(rank - 2):
            C[i][i+1] = -1; C[i+1][i] = -1
        if rank >= 4:
            C[rank-3][rank-1] = -1; C[rank-1][rank-3] = -1
    elif dtype == 'E':
        # E_n Bourbaki: chain 0-2-3-4-..-(n-1), with 1 branching from 3
        C[0][2] = -1; C[2][0] = -1
        C[2][3] = -1; C[3][2] = -1
        C[1][3] = -1; C[3][1] = -1  # branch at node 3
        for i in range(3, rank - 1):
            C[i][i+1] = -1; C[i+1][i] = -1
    return C

def positive_roots_correct(dtype, rank):
    """Enumerate positive roots using root string method.

    Algorithm: start with simple roots. For each root alpha and simple root alpha_i,
    if <alpha, alpha_i^vee> < 0, then alpha + alpha_i is also a positive root.
    """
    C = cartan_matrix(dtype, rank)
    roots = set()

    # Simple roots as standard basis vectors
    for i in range(rank):
        e = tuple([1 if j == i else 0 for j in range(rank)])
        roots.add(e)

    changed = True
    while changed:
        changed = False
        new_roots = set()
        for alpha in roots:
            for i in range(rank):
                # Compute <alpha, alpha_i^vee> = sum_j alpha_j * C[j,i]
                inner = sum(alpha[j] * C[j][i] for j in range(rank))
                if inner < 0:
                    # alpha + alpha_i is a root
                    new = list(alpha)
                    new[i] += 1
                    new_t = tuple(new)
                    if new_t not in roots:
                        new_roots.add(new_t)
                        changed = True
        roots.update(new_roots)

    return sorted(roots)

def main():
    print("=" * 70)
    print("QUIVER-DYNKIN DEEP ANALYSIS (CORRECTED)")
    print("=" * 70)

    KEY1, KEY2 = 2, 3

    # =================================================================
    # PART 1: CORRECTED POSITIVE ROOT COUNTS
    # =================================================================
    print("\n--- PART 1: GABRIEL NUMBERS |Phi+| = nh/2 ---\n")

    dynkin_data = {
        ('A', 1): 2, ('A', 2): 3, ('A', 3): 4, ('A', 4): 5,
        ('A', 5): 6, ('A', 6): 7, ('A', 7): 8,
        ('D', 4): 6, ('D', 5): 8, ('D', 6): 10, ('D', 7): 12,
        ('E', 6): 12, ('E', 7): 18, ('E', 8): 30,
    }

    results = {}
    for (dtype, rank), h in sorted(dynkin_data.items()):
        roots = positive_roots_correct(dtype, rank)
        expected = rank * h // 2
        actual = len(roots)
        results[(dtype, rank)] = (actual, h, roots)

        note = ""
        if (dtype, rank) == ('A', 6): note = " = H_forb_2 = Phi3(4)"
        elif (dtype, rank) == ('E', 6): note = " = C(9,2)"
        elif (dtype, rank) == ('E', 7): note = " = 7*9 = FORBIDDEN at n<=7"
        elif (dtype, rank) == ('E', 8): note = " = |BI| = |S_5|"
        elif (dtype, rank) == ('D', 4): note = " = h(E_6)"
        elif (dtype, rank) == ('A', 7): note = " = C(8,2)"

        status = "OK" if actual == expected else f"FAIL (got {actual})"
        print(f"  {dtype}_{rank}: |Phi+| = {actual:4d}  (nh/2 = {expected}, h={h:2d}) {status}{note}")

    # =================================================================
    # PART 2: THE EXCEPTIONAL ROOT COUNT SEQUENCE
    # =================================================================
    print(f"\n--- PART 2: EXCEPTIONAL ROOT COUNTS ---\n")

    e6_count = results[('E', 6)][0]
    e7_count = results[('E', 7)][0]
    e8_count = results[('E', 8)][0]

    print(f"  |Phi+(E_6)| = {e6_count}")
    print(f"  |Phi+(E_7)| = {e7_count}")
    print(f"  |Phi+(E_8)| = {e8_count}")

    print(f"\n  Total roots (positive + negative):")
    print(f"    |Phi(E_6)| = {2*e6_count} = 2*{e6_count}")
    print(f"    |Phi(E_7)| = {2*e7_count} = 2*{e7_count}")
    print(f"    |Phi(E_8)| = {2*e8_count} = 240 = #roots(E_8)")

    print(f"\n  Connections to tournament theory:")
    print(f"    |Phi+(A_6)| = 21 = H_forb_2 (permanently forbidden)")
    print(f"    |Phi+(E_7)| = 63 = forbidden at n<=7, achievable at n=8=rank(E_8)")
    print(f"    |Phi(E_8)|  = 240 = V(icos)*V(dodec) = 12*20")

    print(f"\n  The sequence |Phi+(E_n)| = {{36, 63, 120}}:")
    print(f"    63/36 = {63/36:.6f} = 7/4 = H_forb_1/KEY_1^2")
    print(f"    120/63 = {120/63:.6f} = 40/21 = 40/H_forb_2")
    print(f"    120/36 = {120/36:.6f} = 10/3 = V(Petersen)/KEY_2")

    # =================================================================
    # PART 3: HEIGHT DISTRIBUTION OF E-TYPE ROOTS
    # =================================================================
    print(f"\n--- PART 3: ROOT HEIGHT DISTRIBUTIONS ---\n")

    for (dtype, rank) in [('E', 6), ('E', 7), ('E', 8)]:
        count, h, roots = results[(dtype, rank)]
        heights = Counter(sum(r) for r in roots)
        highest = max(heights.keys())
        print(f"  {dtype}_{rank} (h={h}, |Phi+|={count}):")
        print(f"    Heights span: 1 to {highest} (should be 1 to {h-1})")
        print(f"    Highest root: {max(roots, key=sum)}")

        # Number of roots at each height
        ht_seq = [heights.get(k, 0) for k in range(1, highest+1)]
        print(f"    Height sequence: {ht_seq}")
        print(f"    Palindromic? {ht_seq == ht_seq[::-1]}")

        # Exponents: heights where there's a "first occurrence" of rank many roots
        print(f"    #roots at height 1: {heights[1]} (= rank = {rank})")
        print(f"    #roots at height {highest}: {heights[highest]} (= 1, highest root)")

    # =================================================================
    # PART 4: McKAY DIMENSIONS AND QUIVER INVARIANTS
    # =================================================================
    print(f"\n--- PART 4: McKAY QUIVER DIMENSIONS ---\n")

    mckay = {
        'BT': {'group': 'BT', 'order': 24, 'dynkin': 'E_6',
               'dims': [1, 1, 1, 2, 2, 2, 3],
               'marks': [1, 1, 2, 3, 2, 1, 1]},  # affine E_6 marks
        'BO': {'group': 'BO', 'order': 48, 'dynkin': 'E_7',
               'dims': [1, 1, 2, 2, 3, 3, 2, 4],
               'marks': [1, 2, 3, 4, 3, 2, 1, 2]},  # affine E_7 marks
        'BI': {'group': 'BI', 'order': 120, 'dynkin': 'E_8',
               'dims': [1, 2, 3, 4, 5, 6, 4, 2, 3],
               'marks': [2, 4, 6, 5, 4, 3, 2, 1, 3]},  # affine E_8 marks
    }

    for name, data in mckay.items():
        dims = data['dims']
        print(f"  {name} ({data['dynkin']}, |G|={data['order']}):")
        print(f"    Irrep dims: {dims}")
        print(f"    Sum(d_i) = {sum(dims)} = h({data['dynkin']})")
        print(f"    Sum(d_i^2) = {sum(d**2 for d in dims)} = |{name}|")
        print(f"    #irreps = {len(dims)} = rank+1")
        print(f"    Max dim = {max(dims)}")

        # The dimension generating function
        dim_counter = Counter(dims)
        print(f"    Dimension spectrum: {dict(sorted(dim_counter.items()))}")

    # =================================================================
    # PART 5: INFORMATION THEORY — TOURNAMENT CHANNEL
    # =================================================================
    print(f"\n--- PART 5: TOURNAMENT INFORMATION CHANNEL ---\n")

    print(f"  Define the TOURNAMENT CHANNEL:")
    print(f"    Input: tournament T on n vertices (C(n,2) bits)")
    print(f"    Output: H(T) (odd integer)")
    print(f"    This is a DETERMINISTIC channel (zero noise)")
    print(f"    But it's MANY-TO-ONE, so information is lost.")

    # H-distribution data
    h_data = {
        3: {1: 6, 3: 2},
        4: {1: 24, 3: 16, 5: 24},
        5: {1: 120, 3: 120, 5: 240, 9: 240, 11: 120, 13: 120, 15: 64},
        6: {1:720, 3:960, 5:2160, 9:2960, 11:1440, 13:1440, 15:2208,
            17:1440, 19:1440, 23:2880, 25:1440, 27:480, 29:2880,
            31:1440, 33:2640, 37:3600, 41:720, 43:1440, 45:480},
    }

    print(f"\n  Channel capacity analysis:")
    for n, dist in sorted(h_data.items()):
        total = sum(dist.values())
        edges = n * (n - 1) // 2
        n_outputs = len(dist)

        # Shannon entropy of H distribution
        probs = [c / total for c in dist.values()]
        H_entropy = -sum(p * math.log2(p) for p in probs if p > 0)

        # Maximum possible entropy (uniform over outputs)
        H_max = math.log2(n_outputs)

        # Input entropy (uniform over tournaments)
        H_input = edges  # log2(2^edges) = edges

        # Mutual information = H_entropy (since channel is deterministic)
        I_mutual = H_entropy

        # Compression ratio
        ratio = I_mutual / H_input

        print(f"    n={n}: input={H_input} bits, I(T;H)={I_mutual:.3f} bits, "
              f"|outputs|={n_outputs}, ratio={ratio:.4f}")

    # =================================================================
    # PART 6: FORBIDDEN H AS EXCLUDED CODEWORDS
    # =================================================================
    print(f"\n--- PART 6: MOAT AS EXCLUDED CODEWORDS ---\n")

    print(f"  The permanent moat {{7, 21}} = excluded H values.")
    print(f"  In coding theory terms: H(T) is a CODE with forbidden codewords.")
    print(f"")
    print(f"  Connection to Hamming code [7,4,3]:")
    print(f"    Code length = 7 = H_forb_1")
    print(f"    Dimension = 4 = rank(F_4) = KEY_1^2")
    print(f"    Distance = 3 = KEY_2")
    print(f"    #codewords = 16 = 2^4")
    print(f"")
    print(f"  The H-CODE:")
    print(f"    'Alphabet' = odd positive integers")
    print(f"    'Forbidden words' = {{7, 21}}")
    print(f"    'Code distance' = ? (minimum gap between achievable H values)")

    # Compute minimum gaps
    for n, dist in sorted(h_data.items()):
        vals = sorted(dist.keys())
        gaps = [vals[i+1] - vals[i] for i in range(len(vals)-1)]
        min_gap = min(gaps) if gaps else 0
        max_gap = max(gaps) if gaps else 0
        print(f"    n={n}: H in {vals}")
        print(f"          min gap={min_gap}, max gap={max_gap}")

    # =================================================================
    # PART 7: QUIVER MUTATION AND THE MOAT
    # =================================================================
    print(f"\n--- PART 7: QUIVER MUTATION INTERPRETATION ---\n")

    print(f"  QUIVER MUTATION (Fomin-Zelevinsky, 2002):")
    print(f"    Mutating at vertex k: reverse arrows at k, apply exchange rule")
    print(f"    Cluster algebras classify finite-type quivers -> Dynkin diagrams!")
    print(f"")
    print(f"  TOURNAMENT MUTATION = ARC FLIP:")
    print(f"    Flipping arc (u,v) in tournament T:")
    print(f"    - Reverses the arrow between u and v")
    print(f"    - Changes H(T) by H(T/e) - H(T'/e')")
    print(f"    - The DELTA is always even (H odd -> H' odd)")
    print(f"")
    print(f"  Arc flip as quiver mutation at a vertex pair:")
    print(f"    Tournament mutations connect ALL n-vertex tournaments")
    print(f"    (the flip graph is connected)")
    print(f"    But cluster mutations preserve Dynkin type!")
    print(f"    Tournament 'type' = H value (H is a mutation invariant class label)")

    # =================================================================
    # PART 8: AUSLANDER-REITEN THEORY
    # =================================================================
    print(f"\n--- PART 8: AUSLANDER-REITEN PERSPECTIVE ---\n")

    print(f"  The AR quiver of a Dynkin algebra kQ:")
    print(f"    Vertices = indecomposable modules (= positive roots by Gabriel)")
    print(f"    Arrows = irreducible morphisms")
    print(f"    Shape: preprojective + regular + preinjective components")
    print(f"")
    print(f"  For Dynkin type, the AR quiver has |Phi+| vertices:")
    print(f"    A_6: 21 = H_forb_2 indecomposables")
    print(f"    E_6: 36 indecomposables")
    print(f"    E_7: 63 indecomposables")
    print(f"    E_8: 120 = |BI| indecomposables")
    print(f"")
    print(f"  CONJECTURE (HYP-new): The permanent forbidden H values")
    print(f"  {{7, 21}} correspond to Gabriel numbers of specific quivers:")
    print(f"    H_forb_1 = 7 = |Phi+(A_3)|... wait, |Phi+(A_3)| = 3*4/2 = 6.")
    print(f"    Actually: 7 = |Phi+(B_2)| = |Phi+(C_2)| = |Phi+(G_2)| (non-simply-laced!)")

    # Check: |Phi+(G_2)| = 2*6/2 = 6. Hmm, G_2 has h=6, rank=2, |Phi+|=6. Not 7.
    # Actually G_2 has 6 positive roots (6 short + 0? No: 3 short + 3 long).
    # B_2 has h=4, rank=2, |Phi+|=4. Not 7.
    # What has |Phi+| = 7?
    # A_n: n(n+1)/2 = 7 -> n=3.something. No integer solution.
    # D_n: n(n-1) = 7. No.
    # So 7 is NOT a Gabriel number for any Dynkin type!

    print(f"\n  CORRECTION: 7 is NOT a Gabriel number for any Dynkin diagram!")
    print(f"  Gabriel numbers (|Phi+|) for ranks 1-8:")
    gabriel_nums = set()
    for dtype in ['A', 'D']:
        for rank in range(1, 9):
            if dtype == 'D' and rank < 4:
                continue
            h = dynkin_data.get((dtype, rank))
            if h:
                g = rank * h // 2
                gabriel_nums.add(g)
                print(f"    {dtype}_{rank}: {g}")
    for rank in [6, 7, 8]:
        h = dynkin_data[('E', rank)]
        g = rank * h // 2
        gabriel_nums.add(g)
        print(f"    E_{rank}: {g}")

    print(f"\n  Gabriel numbers: {sorted(gabriel_nums)}")
    print(f"  Is 7 a Gabriel number? {7 in gabriel_nums}")
    print(f"  Is 21 a Gabriel number? {21 in gabriel_nums} (A_6!)")
    print(f"  Is 63 a Gabriel number? {63 in gabriel_nums} (E_7!)")
    print(f"  Is 120 a Gabriel number? {120 in gabriel_nums} (E_8!)")

    print(f"\n  DEEP OBSERVATION:")
    print(f"    H=21 is the Gabriel number of A_6 (path quiver on 7 vertices)")
    print(f"    7 = the NUMBER OF VERTICES of the A_6 quiver!")
    print(f"    H_forb_1 = #vertices of the quiver whose Gabriel number = H_forb_2")
    print(f"    This is a SELF-REFERENTIAL structure.")

    # =================================================================
    # PART 9: CHANNEL CAPACITY BOUNDS
    # =================================================================
    print(f"\n--- PART 9: ASYMPTOTIC CHANNEL RATE ---\n")

    print(f"  As n -> infinity:")
    print(f"    #distinct H values grows (monotonicity theorem)")
    print(f"    Input entropy = C(n,2) = n(n-1)/2 bits")
    print(f"    H-entropy <= log2(max_H) ~ O(n log n) bits")
    print(f"")
    print(f"  Channel rate R = I(T;H) / C(n,2):")
    for n in [3, 4, 5, 6]:
        dist = h_data[n]
        total = sum(dist.values())
        edges = n * (n - 1) // 2
        probs = [c / total for c in dist.values()]
        H_entropy = -sum(p * math.log2(p) for p in probs if p > 0)
        rate = H_entropy / edges
        print(f"    n={n}: R = {rate:.4f} = {H_entropy:.3f}/{edges}")
    print(f"    Rate -> 0 as n -> inf (H captures vanishing fraction of tournament info)")

    # =================================================================
    # PART 10: SYNTHESIS
    # =================================================================
    print(f"\n{'='*70}")
    print("GRAND SYNTHESIS: QUIVER-DYNKIN-INFORMATION-TOURNAMENT")
    print("=" * 70)
    print(f"""
  THE QUIVER-TOURNAMENT DICTIONARY:
  ==================================

  Tournament theory          |  Quiver/Dynkin theory
  ========================== | =============================
  Tournament T on n vertices |  Complete quiver Q_T (wild type)
  Arc flip (mutation)        |  Quiver mutation
  H(T) = I(Omega(T), 2)     |  Gabriel number |Phi+|
  Forbidden H = 7            |  NOT a Gabriel number (orphan)
  Forbidden H = 21           |  |Phi+(A_6)| = Gabriel number!
  H=63 (n<=7 forbidden)      |  |Phi+(E_7)| = 63
  H=120 (achievable)         |  |Phi+(E_8)| = 120 = |BI|
  Permanent moat {{7, 21}}    |  {{orphan, Phi+(A_6)}}

  INFORMATION-THEORETIC VIEW:
  ============================
  Tournament = C(n,2)-bit codeword
  H(T) = lossy compression (retains ~28% of info at n=6)
  Forbidden H values = excluded codewords
  Hamming [7,4,3] = [H_forb_1, rank(F4), KEY_2]

  McKAY BRIDGE:
  =============
  BT -> affine E_6: sum(dims) = 12 = h(E_6), sum(dims^2) = 24 = |BT|
  BO -> affine E_7: sum(dims) = 18 = h(E_7), sum(dims^2) = 48 = |BO|
  BI -> affine E_8: sum(dims) = 30 = h(E_8), sum(dims^2) = 120 = |BI|

  THE SELF-REFERENTIAL STRUCTURE:
  ================================
  H_forb_2 = 21 = |Phi+(A_6)| where A_6 has 7 = H_forb_1 vertices
  The first forbidden value COUNTS the vertices of the quiver whose
  indecomposable count equals the second forbidden value.
""")

if __name__ == "__main__":
    main()
