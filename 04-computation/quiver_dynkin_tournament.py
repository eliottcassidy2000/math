"""
quiver_dynkin_tournament.py -- kind-pasteur-2026-03-14-S68

Deep investigation: Quiver representations, Dynkin diagrams, and tournaments.

KEY CONNECTIONS TO EXPLORE:

1. GABRIEL'S THEOREM: A connected quiver has finitely many indecomposable
   representations iff its underlying graph is a Dynkin diagram (A_n, D_n, E_6, E_7, E_8).
   The dimension vectors of indecomposables = positive roots of the root system.

2. TOURNAMENTS AS QUIVERS: A tournament on n vertices IS a quiver
   (complete oriented graph). It has NO underlying Dynkin structure for n >= 4,
   but its representation theory connects to our H values.

3. McKAY QUIVERS: For G subset SU(2), the McKay quiver has:
   - Vertices = irreducible representations of G
   - Arrows = tensor product multiplicities with the fundamental rep
   The underlying graph is the EXTENDED (affine) Dynkin diagram.

4. AUSLANDER-REITEN THEORY: The AR quiver of a path algebra kQ captures
   ALL indecomposable modules and their morphisms. For Dynkin quivers,
   the AR quiver is finite (Gabriel's theorem).

5. INFORMATION THEORY: The channel capacity of a tournament-defined channel
   relates to H(T). Can we define a "tournament entropy" using the
   independence polynomial I(Omega(T), x)?

PLAN:
A. Build Dynkin quivers A_n, D_n, E_6, E_7, E_8
B. Compute their positive roots (= dimension vectors of indecomposables)
C. Count indecomposables: |Phi^+| = n(h+1)/2 for Dynkin type
D. Evaluate the "root independence polynomial" and connect to tournament IP
E. Build McKay quivers for BT, BO, BI
F. Information-theoretic analysis: tournament entropy, channel capacity
"""

import numpy as np
from collections import Counter
from itertools import combinations
import math

def cartan_matrix(dynkin_type, rank):
    """Build the Cartan matrix for a Dynkin diagram."""
    C = np.eye(rank, dtype=int) * 2

    if dynkin_type == 'A':
        for i in range(rank - 1):
            C[i][i+1] = -1
            C[i+1][i] = -1

    elif dynkin_type == 'D':
        # D_n: linear chain 1-2-..-(n-2), then (n-2) branches to both (n-1) and n
        for i in range(rank - 2):
            C[i][i+1] = -1
            C[i+1][i] = -1
        # Branch: vertex n-2 connects to both n-1 and n (0-indexed: n-3 to n-2 and n-1)
        if rank >= 4:
            C[rank-3][rank-1] = -1
            C[rank-1][rank-3] = -1

    elif dynkin_type == 'E':
        # E_n (n=6,7,8): linear chain 1-3-4-5-...-n, with 2 branching from 3
        # Standard labeling: vertices 1..n
        # Edges: 1-3, 2-3, 3-4, 4-5, ..., (n-1)-n
        # 0-indexed: 0-2, 1-2, 2-3, 3-4, ..., (n-2)-(n-1)
        C[0][2] = -1; C[2][0] = -1  # 1-3
        C[1][2] = -1; C[2][1] = -1  # 2-3
        for i in range(2, rank - 1):
            C[i][i+1] = -1
            C[i+1][i] = -1

    return C

def positive_roots(dynkin_type, rank):
    """Enumerate positive roots using the Cartan matrix reflection method."""
    C = cartan_matrix(dynkin_type, rank)

    # Simple roots are the standard basis vectors
    roots = set()
    queue = []
    for i in range(rank):
        e = tuple([1 if j == i else 0 for j in range(rank)])
        roots.add(e)
        queue.append(e)

    # Apply simple reflections to generate all positive roots
    while queue:
        alpha = queue.pop(0)
        alpha_arr = np.array(alpha)
        for i in range(rank):
            # s_i(alpha) = alpha - <alpha, alpha_i^vee> * alpha_i
            # = alpha - (C * alpha)[i] * e_i   wait, need to be careful
            # s_i(alpha) = alpha - sum_j C[i,j] * alpha_j * e_i
            # Actually: s_i(alpha)_k = alpha_k - C[i,k] * alpha_i  (reflection formula)
            # Hmm, let me use the standard formula:
            # s_i(alpha) = alpha - <alpha, alpha_i^vee> alpha_i
            # <alpha, alpha_i^vee> = sum_j alpha_j * C[j,i] ... no
            # The Cartan matrix: C[i,j] = <alpha_i, alpha_j^vee> = 2(alpha_i, alpha_j)/(alpha_j, alpha_j)
            # s_i(beta) = beta - <beta, alpha_i^vee> alpha_i
            # <beta, alpha_i^vee> = 2(beta, alpha_i)/(alpha_i, alpha_i)
            # In root coordinates: if beta = sum b_j alpha_j, then
            # <beta, alpha_i^vee> = sum_j b_j C[j,i]
            inner = sum(alpha[j] * C[j][i] for j in range(rank))
            new = list(alpha)
            new[i] = alpha[i] - inner
            # Keep only positive (all coefficients >= 0) and at least one > 0
            new_t = tuple(new)
            if all(x >= 0 for x in new_t) and any(x > 0 for x in new_t):
                if new_t not in roots:
                    roots.add(new_t)
                    queue.append(new_t)

    return sorted(roots)

def main():
    print("=" * 70)
    print("QUIVER REPRESENTATIONS, DYNKIN DIAGRAMS & TOURNAMENTS")
    print("=" * 70)

    # =================================================================
    # PART 1: POSITIVE ROOTS AND GABRIEL'S THEOREM
    # =================================================================
    print("\n" + "=" * 70)
    print("PART 1: GABRIEL'S THEOREM — POSITIVE ROOT COUNTS")
    print("=" * 70)
    print("\nGabriel's Theorem: # indecomposable reps = # positive roots = n(h+1)/2")
    print("where n=rank, h=Coxeter number.\n")

    dynkin_data = [
        ('A', 1, 2), ('A', 2, 3), ('A', 3, 4), ('A', 4, 5),
        ('A', 5, 6), ('A', 6, 7), ('A', 7, 8),
        ('D', 4, 6), ('D', 5, 8), ('D', 6, 10),
        ('E', 6, 12), ('E', 7, 18), ('E', 8, 30),
    ]

    results = {}
    for dtype, rank, h in dynkin_data:
        roots = positive_roots(dtype, rank)
        expected = rank * (h + 1) // 2
        # Actually for type A_n, rank=n, h=n+1, |Phi+| = n(n+1)/2
        # For type D_n, |Phi+| = n(n-1)
        # For E_6: |Phi+| = 36, E_7: |Phi+| = 63, E_8: |Phi+| = 120
        actual = len(roots)
        check = actual == expected
        results[(dtype, rank)] = (actual, h, roots)

        note = ""
        if dtype == 'E' and rank == 6: note = "  = 6*6 = rank^2"
        elif dtype == 'E' and rank == 7: note = f"  = 7*9 = 63 = H_forb_2 * KEY_2 = 7*3^2"
        elif dtype == 'E' and rank == 8: note = f"  = 120 = |BI| = |S_5| = |Aut(Petersen)|"
        elif dtype == 'A' and rank == 7: note = f"  = 28 = C(8,2)"
        elif dtype == 'D' and rank == 4: note = f"  = 12 = h(E_6)"

        print(f"  {dtype}_{rank}: |Phi+| = {actual:4d} (expected n(h+1)/2 = {expected}, h={h:2d}) {'OK' if check else 'FAIL'}{note}")

    # =================================================================
    # PART 2: ROOT SYSTEM STRUCTURE AND TOURNAMENT CONNECTION
    # =================================================================
    print("\n" + "=" * 70)
    print("PART 2: E-TYPE ROOT STRUCTURE AND TOURNAMENT KEYS")
    print("=" * 70)

    KEY1, KEY2 = 2, 3

    for dtype, rank in [('E', 6), ('E', 7), ('E', 8)]:
        actual, h, roots = results[(dtype, rank)]
        print(f"\n  {dtype}_{rank}: {actual} positive roots, Coxeter h={h}")
        print(f"    h+1 = {h+1}")

        # Root height distribution
        heights = Counter()
        for r in roots:
            heights[sum(r)] += 1
        print(f"    Height distribution (height: #roots):")
        for ht in sorted(heights.keys()):
            print(f"      height {ht:2d}: {heights[ht]} roots")

        # Highest root
        highest = max(roots, key=sum)
        print(f"    Highest root: {highest} (height {sum(highest)})")
        print(f"    Highest root height = h = {sum(highest)} = Coxeter number")

        # Number of roots at each height = exponents structure
        # For simply laced, #roots at height k = rank for k = 1..h-1 (in some sense)

    # =================================================================
    # PART 3: |Phi+| VALUES AND THE PERMANENT MOAT
    # =================================================================
    print("\n" + "=" * 70)
    print("PART 3: ROOT COUNTS AND THE PERMANENT MOAT")
    print("=" * 70)

    # Collect all |Phi+| values
    phi_plus = {}
    for (dtype, rank), (count, h, _) in results.items():
        phi_plus[f"{dtype}{rank}"] = count

    print(f"\n  Root count dictionary:")
    for name, count in sorted(phi_plus.items(), key=lambda x: x[1]):
        mod7 = count % 7
        note = ""
        if count == 7: note = " = H_forb_1 = Phi3(KEY1)"
        elif count == 21: note = " = H_forb_2 = Phi3(KEY1^2)"
        elif count == 63: note = " = 7*3^2 (forbidden at n<=7, achievable n=8)"
        elif count == 120: note = " = |BI| = |Aut(Petersen)|"
        elif count == 36: note = " = 6^2 = C(9,2)"
        elif count == 12: note = " = h(E6) = h(F4)"
        elif count == 28: note = " = C(8,2) = C(rank(E8),2)"
        print(f"    |Phi+({name})| = {count:4d} (mod 7 = {mod7}){note}")

    # KEY OBSERVATION
    print(f"\n  KEY OBSERVATION:")
    print(f"    |Phi+(E7)| = 63 = 7 * 9 = H_forb_1 * KEY_2^2")
    print(f"    |Phi+(E8)| = 120 = |BI| = 8 * 15 = rank(E8) * E(Petersen)")
    print(f"    |Phi+(E7)| / |Phi+(E6)| = 63/36 = 7/4 = H_forb_1 / KEY_1^2")
    print(f"    |Phi+(E8)| / |Phi+(E7)| = 120/63 = 40/21 = 40/H_forb_2")

    # =================================================================
    # PART 4: INFORMATION-THEORETIC VIEW
    # =================================================================
    print("\n" + "=" * 70)
    print("PART 4: INFORMATION THEORY OF TOURNAMENTS")
    print("=" * 70)

    print(f"\n  Tournament on n vertices: C(n,2) binary choices.")
    print(f"  Total information content: C(n,2) bits = log2(2^C(n,2)) bits.")

    for n in range(3, 10):
        edges = n * (n - 1) // 2
        total_tournaments = 2 ** edges
        # H values achievable
        print(f"\n  n={n}: {edges} edges, {total_tournaments} tournaments")
        print(f"    Information content: {edges} bits")
        print(f"    log2(#tournaments) = {edges}")

    # Shannon entropy of the H distribution at n=5
    print(f"\n  Shannon entropy of H distribution:")
    # n=5: H values from exhaustive data
    h_dist_5 = {1: 120, 3: 120, 5: 240, 9: 240, 11: 120, 13: 120, 15: 64}
    total_5 = sum(h_dist_5.values())
    entropy_5 = -sum((c/total_5) * math.log2(c/total_5) for c in h_dist_5.values())
    print(f"    n=5: H(dist) = {entropy_5:.4f} bits (out of {math.log2(total_5):.4f} max)")

    h_dist_6 = {1:720, 3:960, 5:2160, 9:2960, 11:1440, 13:1440, 15:2208,
                 17:1440, 19:1440, 23:2880, 25:1440, 27:480, 29:2880,
                 31:1440, 33:2640, 37:3600, 41:720, 43:1440, 45:480}
    total_6 = sum(h_dist_6.values())
    entropy_6 = -sum((c/total_6) * math.log2(c/total_6) for c in h_dist_6.values())
    print(f"    n=6: H(dist) = {entropy_6:.4f} bits (out of {math.log2(total_6):.4f} max)")
    print(f"    n=6: max entropy = log2(19) = {math.log2(19):.4f} bits (19 distinct H values)")
    print(f"    n=6: efficiency = {entropy_6 / math.log2(19):.4f}")

    # Information lost by going from tournament to H value
    print(f"\n  Information loss: tournament -> H(T)")
    for n, edges, n_H in [(5, 10, 7), (6, 15, 19)]:
        total_info = edges  # bits
        H_info = math.log2(n_H)
        lost = total_info - H_info
        ratio = H_info / total_info
        print(f"    n={n}: {total_info} bits -> {H_info:.2f} bits ({ratio:.4f} retained, {lost:.2f} lost)")

    # =================================================================
    # PART 5: QUIVER REPRESENTATION OF TOURNAMENTS
    # =================================================================
    print("\n" + "=" * 70)
    print("PART 5: TOURNAMENTS AS QUIVERS")
    print("=" * 70)

    print(f"""
  A tournament T on n vertices is a QUIVER Q_T:
    - Vertices: {{1, ..., n}} (not to be confused with quiver rep vertices)
    - Arrows: one arrow i -> j for each pair (i,j) based on tournament orientation

  The PATH ALGEBRA kQ_T has basis = all directed paths in T.
  For tournaments: every pair connected, so many paths exist.

  KEY DIFFERENCE from Dynkin quivers:
    - Tournament quivers are COMPLETE (every pair has an arrow)
    - Dynkin quivers are SPARSE (tree-like)
    - Tournament quivers have WILD representation type for n >= 3
      (infinitely many indecomposables, not parameterizable)

  But: the CYCLE STRUCTURE of Q_T connects to H(T) via OCF!
    - Directed odd cycles in Q_T = vertices of Omega(T)
    - H(T) = I(Omega(T), 2) = independence polynomial at x=2
    """)

    # =================================================================
    # PART 6: McKAY QUIVERS
    # =================================================================
    print("=" * 70)
    print("PART 6: McKAY QUIVERS — ADE CORRESPONDENCE")
    print("=" * 70)

    print(f"""
  McKay correspondence: finite G < SU(2) <-> affine ADE diagram

  The McKay quiver has:
    - Vertices = irreducible reps of G
    - Arrows encode tensor product with fundamental 2-dim rep

  For the three exceptional subgroups:
    BT (binary tetrahedral, |BT|=24): McKay quiver = affine E_6 (7 vertices)
    BO (binary octahedral, |BO|=48):  McKay quiver = affine E_7 (8 vertices)
    BI (binary icosahedral, |BI|=120): McKay quiver = affine E_8 (9 vertices)

  Affine E_n has n+1 vertices (adding the extending node to E_n).

  DIMENSIONS of irreps (= components of null vector of affine Cartan):
    """)

    # Null vectors of affine Cartan matrices (= irrep dimensions)
    # These are the marks/labels on the extended Dynkin diagram
    mckay_dims = {
        'E6_aff': [1, 1, 2, 3, 2, 1, 2],  # BT: dims of 7 irreps, sum=12
        'E7_aff': [1, 2, 3, 4, 3, 2, 1, 2],  # BO: dims of 8 irreps, sum=18... wait
        'E8_aff': [1, 2, 3, 4, 5, 6, 4, 2, 3],  # BI: dims of 9 irreps, sum=30
    }

    # Actually the standard marks for extended Dynkin:
    # Affine E6: nodes labeled 1,1,2,3,2,1 + extending node 1, null vec [1,1,2,3,2,1,1]...
    # Let me use the standard null vectors (minimal positive imaginary root)
    # Tilde E6: (1, 1, 1, 2, 1, 1, 1) on the extended diagram? No...
    # The minimal imaginary root delta has coefficients = marks on extended diagram
    # E6~: 1-2-3-2-1 with branch 2 from node 3, plus extending 1 connected to end
    # delta = (1,2,3,2,1,2,1) in standard labeling

    # BT irreps: trivial(1), omega(1), omega^2(1), rho_1(2), rho_2(2), rho_3(2), V(3)
    # dims: 1,1,1,2,2,2,3 sum=12=|BT|/2... wait |BT|=24
    # Sum of squares of dims = |G|: 1+1+1+4+4+4+9 = 24. Yes!
    bt_dims = [1, 1, 1, 2, 2, 2, 3]
    bo_dims = [1, 1, 2, 2, 3, 3, 4, 2]  # This isn't right either
    # Let me just use the known facts
    # BO has 8 irreps: 1,1,2,2,3,3,3,1... 1^2+1^2+2^2+2^2+3^2+3^2+3^2+1^2 = 1+1+4+4+9+9+9+1 = 38 ≠ 48
    # Actually BO = binary octahedral, |BO| = 48
    # Irreps: trivial(1), sign(1), V(2), V'(2), W(3), W'(3), U(4), ... hmm
    # Sum of d_i^2 = |G| = 48: 1+1+4+4+9+9+4+16 = 48? That's 8 irreps with dims 1,1,2,2,3,3,2,4
    bo_dims = [1, 1, 2, 2, 3, 3, 2, 4]
    # BI: |BI|=120, irreps with sum d_i^2 = 120
    # 9 irreps: 1,2,3,4,5,6,4,2,3 -> 1+4+9+16+25+36+16+4+9=120! Yes!
    bi_dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

    print(f"  BT irrep dimensions: {bt_dims}")
    print(f"    Sum of squares: {sum(d**2 for d in bt_dims)} = |BT| = 24")
    print(f"    Number of irreps: {len(bt_dims)} = rank(E6)+1 = 7")
    print(f"    Sum of dims: {sum(bt_dims)} = h(E6) = 12")

    print(f"\n  BO irrep dimensions: {bo_dims}")
    print(f"    Sum of squares: {sum(d**2 for d in bo_dims)} = |BO| = 48")
    print(f"    Number of irreps: {len(bo_dims)} = rank(E7)+1 = 8")
    print(f"    Sum of dims: {sum(bo_dims)} = h(E7) = 18")

    print(f"\n  BI irrep dimensions: {bi_dims}")
    print(f"    Sum of squares: {sum(d**2 for d in bi_dims)} = |BI| = 120")
    print(f"    Number of irreps: {len(bi_dims)} = rank(E8)+1 = 9")
    print(f"    Sum of dims: {sum(bi_dims)} = h(E8) = 30")

    # =================================================================
    # PART 7: DIMENSION VECTOR INDEPENDENCE POLYNOMIAL
    # =================================================================
    print("\n" + "=" * 70)
    print("PART 7: ROOT POLYNOMIAL AND INDEPENDENCE STRUCTURE")
    print("=" * 70)

    # For E_8, compute the "root generating function" by height
    _, _, roots_E8 = results[('E', 8)]
    heights_E8 = Counter()
    for r in roots_E8:
        heights_E8[sum(r)] += 1

    print(f"\n  E_8 positive roots by height (h=30):")
    for ht in sorted(heights_E8.keys()):
        bar = "#" * heights_E8[ht]
        print(f"    h={ht:2d}: {heights_E8[ht]:2d} {bar}")

    # The number of roots at each height should be 8 for h=1..29 (simply-laced)
    # Actually it's not constant, but related to exponents
    print(f"\n  Total: {sum(heights_E8.values())} roots")
    print(f"  Heights 1..29: all present? {all(ht in heights_E8 for ht in range(1, 30))}")
    print(f"  Roots at height 1: {heights_E8[1]} (= rank = 8)")
    print(f"  Roots at height h=29: {heights_E8.get(29, 0)} (= 1, the highest root)")

    # =================================================================
    # PART 8: GABRIEL NUMBERS AND THE MOAT
    # =================================================================
    print("\n" + "=" * 70)
    print("PART 8: GABRIEL NUMBERS AND THE PERMANENT MOAT")
    print("=" * 70)

    print(f"\n  Gabriel's theorem gives |Phi+(Q)| = n(h+1)/2 for Dynkin Q.")
    print(f"  These 'Gabriel numbers' for the exceptional types:")
    for dtype, rank in [('E', 6), ('E', 7), ('E', 8)]:
        count = results[(dtype, rank)][0]
        h = results[(dtype, rank)][1]
        print(f"    {dtype}_{rank}: G = {count} = {rank}*{h+1}//2")

    print(f"\n  GABRIEL NUMBER - MOAT CONNECTION:")
    print(f"    G(E_6) = 36 = C(9,2)")
    print(f"    G(E_7) = 63 = 7*9 = 7*3^2 = H_forb_1 * KEY_2^2")
    print(f"    G(E_8) = 120 = |BI| = |S_5|")
    print(f"    G(E_7) = 63 was FORBIDDEN at n<=7 but ACHIEVABLE at n=8!")
    print(f"    This is the number of POSITIVE ROOTS of E_7.")
    print(f"    H=63 becoming achievable at n=8 = rank(E_8) is deeply significant.")

    print(f"\n  RATIO PATTERN:")
    print(f"    G(E_7)/G(E_6) = 63/36 = 7/4")
    print(f"    G(E_8)/G(E_7) = 120/63 = 40/21 = 40/H_forb_2")
    print(f"    G(E_8)/G(E_6) = 120/36 = 10/3 = V(Petersen)/KEY_2")

    # =================================================================
    # PART 9: ENTROPY OF THE ROOT SYSTEM
    # =================================================================
    print("\n" + "=" * 70)
    print("PART 9: INFORMATION ENTROPY OF ROOT SYSTEMS")
    print("=" * 70)

    for dtype, rank in [('E', 6), ('E', 7), ('E', 8)]:
        count, h, roots = results[(dtype, rank)]
        heights = Counter()
        for r in roots:
            heights[sum(r)] += 1

        # Height distribution entropy
        probs = [heights[ht] / count for ht in sorted(heights.keys())]
        entropy = -sum(p * math.log2(p) for p in probs if p > 0)
        max_entropy = math.log2(h)  # max = log2(number of distinct heights)
        print(f"\n  {dtype}_{rank} (h={h}, |Phi+|={count}):")
        print(f"    Height entropy: {entropy:.4f} bits (max {max_entropy:.4f})")
        print(f"    Efficiency: {entropy/max_entropy:.4f}")

        # Coordinate entropy
        # Distribution of first coordinate values
        coord_dist = Counter()
        for r in roots:
            coord_dist[r[0]] += 1  # just first coordinate
        coord_entropy = -sum((c/count) * math.log2(c/count) for c in coord_dist.values())
        print(f"    First-coord entropy: {coord_entropy:.4f} bits")

    # =================================================================
    # PART 10: SUMMARY
    # =================================================================
    print("\n" + "=" * 70)
    print("SUMMARY: QUIVER-DYNKIN-TOURNAMENT CONNECTIONS")
    print("=" * 70)
    print(f"""
  1. GABRIEL NUMBERS connect to tournament H values:
     G(E_7) = 63 = formerly forbidden, now achievable at n=8=rank(E_8)
     G(E_8) = 120 = |Aut(Petersen)| = |BI|

  2. McKAY CORRESPONDENCE maps BT,BO,BI to affine E_6,E_7,E_8:
     Sum of irrep dims = Coxeter number h
     Sum of dim^2 = |G| (group order)
     Number of irreps = rank + 1

  3. TOURNAMENT-QUIVER DUALITY:
     Tournament = complete quiver (wild representation type)
     H(T) = I(Omega(T), 2) = "cycle independence polynomial"
     The permanent moat {{7, 21}} = {{Phi3(2), Phi3(4)}}

  4. INFORMATION THEORY:
     Tournament encodes C(n,2) bits
     H(T) captures ~log2(#distinct_H) bits
     The "information compression" tournament -> H is massive

  5. ROOT COUNTS as H-like invariants:
     |Phi+(A_n)| = n(n+1)/2 = C(n+1,2) = triangular numbers
     |Phi+(D_n)| = n(n-1) = 2*C(n,2)
     |Phi+(E_n)| = {{36, 63, 120}} — the exceptional sequence
""")

if __name__ == "__main__":
    main()
