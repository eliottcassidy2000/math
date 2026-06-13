"""
lex_maximizer_investigation.py -- kind-pasteur-2026-03-14-S82
THE BIG QUESTION: Is max_H always achieved by a lex product?

FROM S81:
- H(T1 lex T2) = H(T1) * H(T2)^|V1| for |V1|=2 (ALL T1, T2)
- 3-cycle lex transitive-2 gives H=45 = max_H(6)!
- This means max_H(6) = 3 * 15^1... wait, that's 45 = 3 * 15. Hmm.
  Actually: T1 = 3-cycle (H=3, n=3), T2 = transitive (H=1, n=2)
  H(T1 lex T2) should be H(T1) * H(T2)^3 = 3 * 1 = 3. But actual is 45!
  So the formula FAILS here. The 45 comes from a different mechanism.

Wait, let me re-examine. The lex product 3-cycle lex T_2 has 6 vertices.
At n=6, max_H = 45. Is this actually a lex product?

Actually, 45 = 3 * 15 = H(3-cycle) * max_H(5). Is there a recursive
relationship? max_H(6) = H(T_3) * max_H(5)?

max_H values: 1, 1, 3, 5, 15, 45, 189, 661, 3357
Ratios: max_H(n+1)/max_H(n) = 1, 3, 5/3, 3, 3, 4.2, 3.5, 5.08

Note: 45/15 = 3 = H(T_3). And 15/5 = 3 = H(T_3). And 189/45 = 4.2...

Is max_H(n) = 3 * max_H(n-2) for even n? max_H(6) = 3*15 = 45 YES!
And 15 = 3*5 = 3*max_H(4). So max_H(5) = 3*5 = 15? YES (5 = max_H(4))
And max_H(4) = 5. Does 5 = 3*max_H(2)? 3*1 = 3 ≠ 5. NO.

What about: max_H(n) related to max_H(n-1)?
45/15 = 3, 15/5 = 3, 5/3 = 5/3, 3/1 = 3.
Not a constant ratio.

Let me investigate the actual STRUCTURE of the maximizers.
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[j][i] = 1
            else:
                A[i][j] = 1
            idx += 1
    return A

def compute_H(A, n):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = sum(dp.get((pm, u), 0) for u in range(n) if (pm & (1 << u)) and A[u][v])
                if t: dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def lex_product(A1, n1, A2, n2):
    n = n1 * n2
    A = np.zeros((n, n), dtype=int)
    for i1 in range(n1):
        for j1 in range(n2):
            for i2 in range(n1):
                for j2 in range(n2):
                    v1 = i1 * n2 + j1
                    v2 = i2 * n2 + j2
                    if v1 == v2: continue
                    if A1[i1][i2] == 1:
                        A[v1][v2] = 1
                    elif i1 == i2 and A2[j1][j2] == 1:
                        A[v1][v2] = 1
    return A, n

def main():
    print("=" * 70)
    print("LEX PRODUCT MAXIMIZER INVESTIGATION")
    print("kind-pasteur-2026-03-14-S82")
    print("=" * 70)

    maxH = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661}

    # ========================================
    # PART 1: What is the structure of max_H tournaments?
    # ========================================
    print(f"\n{'='*70}")
    print("PART 1: STRUCTURE OF H-MAXIMIZERS")
    print(f"{'='*70}")

    for n in [3, 4, 5, 6]:
        m = n*(n-1)//2
        target_H = maxH[n]

        print(f"\n  n={n}, max_H={target_H}:")

        maximizer_count = 0
        maximizer_scores = set()
        maximizer_c3 = set()

        for bits in range(2**m):
            A = bits_to_adj(bits, n)
            H = compute_H(A, n)
            if H == target_H:
                maximizer_count += 1
                scores = tuple(sorted(A.sum(axis=1).astype(int)))
                c3 = int(np.trace(A @ A @ A)) // 3
                maximizer_scores.add(scores)
                maximizer_c3.add(c3)

        print(f"    Count: {maximizer_count}")
        print(f"    Scores: {sorted(maximizer_scores)}")
        print(f"    c3: {sorted(maximizer_c3)}")

        # What fraction of tournaments are maximizers?
        print(f"    Fraction: {maximizer_count}/{2**m} = {maximizer_count/2**m:.4f}")

    # ========================================
    # PART 2: Lex product construction of maximizers
    # ========================================
    print(f"\n{'='*70}")
    print("PART 2: CAN WE BUILD MAXIMIZERS VIA LEX PRODUCTS?")
    print(f"{'='*70}")

    # The 3-cycle on 3 vertices
    cycle3 = np.array([[0,1,0],[0,0,1],[1,0,0]])
    H_cycle3 = compute_H(cycle3.tolist(), 3)
    print(f"  3-cycle: H = {H_cycle3}")

    # Transitive tournaments
    for nn in [2, 3, 4]:
        trans = np.zeros((nn, nn), dtype=int)
        for i in range(nn):
            for j in range(i+1, nn):
                trans[i][j] = 1
        H_trans = compute_H(trans.tolist(), nn)
        print(f"  Transitive T_{nn}: H = {H_trans}")

    # Build lex products and check H
    print(f"\n  Testing lex products that give n=6:")

    # 2 lex 3
    for b1 in range(2):  # only 2 tournaments on 2 vertices
        A1 = bits_to_adj(b1, 2)
        H1 = compute_H(A1.tolist(), 2)
        for b2 in range(2**(3*(3-1)//2)):
            A2 = bits_to_adj(b2, 3)
            H2 = compute_H(A2.tolist(), 3)
            A_prod, n_p = lex_product(A1, 2, A2, 3)
            H_prod = compute_H(A_prod.tolist(), n_p)
            if H_prod == 45:
                scores = tuple(sorted(A_prod.sum(axis=1).astype(int)))
                print(f"    2 lex 3: H1={H1}, H2={H2}, H_prod={H_prod}, scores={scores}")

    # 3 lex 2
    for b1 in range(2**(3*(3-1)//2)):
        A1 = bits_to_adj(b1, 3)
        H1 = compute_H(A1.tolist(), 3)
        for b2 in range(2):
            A2 = bits_to_adj(b2, 2)
            H2 = compute_H(A2.tolist(), 2)
            A_prod, n_p = lex_product(A1, 3, A2, 2)
            H_prod = compute_H(A_prod.tolist(), n_p)
            if H_prod == 45:
                scores = tuple(sorted(A_prod.sum(axis=1).astype(int)))
                print(f"    3 lex 2: H1={H1}, H2={H2}, H_prod={H_prod}, scores={scores}")

    # 6 lex 1 (trivial)
    # 1 lex 6 (trivial)

    # How many of the 480 maximizers at n=6 are lex products?
    print(f"\n  Checking: which n=6 maximizers (H=45) are lex products?")

    # Build ALL lex products that give 6 vertices
    lex_maximizers = set()

    # 2 lex 3
    for b1 in range(2**(2*(2-1)//2)):
        A1 = bits_to_adj(b1, 2)
        for b2 in range(2**(3*(3-1)//2)):
            A2 = bits_to_adj(b2, 3)
            A_prod, _ = lex_product(A1, 2, A2, 3)
            H_prod = compute_H(A_prod.tolist(), 6)
            if H_prod == 45:
                # Store as frozenset of (row tuples) for comparison
                key = tuple(tuple(int(x) for x in row) for row in A_prod)
                lex_maximizers.add(key)

    # 3 lex 2
    for b1 in range(2**(3*(3-1)//2)):
        A1 = bits_to_adj(b1, 3)
        for b2 in range(2**(2*(2-1)//2)):
            A2 = bits_to_adj(b2, 2)
            A_prod, _ = lex_product(A1, 3, A2, 2)
            H_prod = compute_H(A_prod.tolist(), 6)
            if H_prod == 45:
                key = tuple(tuple(int(x) for x in row) for row in A_prod)
                lex_maximizers.add(key)

    print(f"    Lex products achieving H=45: {len(lex_maximizers)}")

    # Count ALL maximizers
    all_maximizers = set()
    for bits in range(2**(6*5//2)):
        A = bits_to_adj(bits, 6)
        if compute_H(A, 6) == 45:
            key = tuple(tuple(int(x) for x in row) for row in A)
            all_maximizers.add(key)

    print(f"    Total tournaments with H=45: {len(all_maximizers)}")
    print(f"    Lex products / total: {len(lex_maximizers)}/{len(all_maximizers)}")
    overlap = len(lex_maximizers & all_maximizers)
    print(f"    Overlap (lex products that ARE maximizers): {overlap}")
    print(f"    Non-lex maximizers: {len(all_maximizers) - overlap}")

    # ========================================
    # PART 3: The lex product formula — when exactly does it work?
    # ========================================
    print(f"\n{'='*70}")
    print("PART 3: WHEN DOES H(T1 lex T2) = H(T1) * H(T2)^|V1|?")
    print("  From S81: works for |V1|=2 (all T2)")
    print("  Fails for |V1|=3 with 3-cycle T1 and |V2|=2")
    print("  Let me check |V1|=3 more carefully")
    print(f"{'='*70}")

    # For |V1|=3, check all pairs with |V2|=2,3
    for n2 in [2, 3]:
        n1 = 3
        n_prod = n1 * n2

        if n_prod > 9:
            print(f"\n  3 lex {n2}: too large (n={n_prod})")
            continue

        print(f"\n  3 lex {n2} (n_prod={n_prod}):")

        successes = 0
        failures = 0
        failure_details = []

        for b1 in range(2**(n1*(n1-1)//2)):
            A1 = bits_to_adj(b1, n1)
            H1 = compute_H(A1.tolist(), n1)

            for b2 in range(2**(n2*(n2-1)//2)):
                A2 = bits_to_adj(b2, n2)
                H2 = compute_H(A2.tolist(), n2)

                A_prod, _ = lex_product(A1, n1, A2, n2)
                H_prod = compute_H(A_prod.tolist(), n_prod)

                predicted = H1 * (H2 ** n1)

                if H_prod == predicted:
                    successes += 1
                else:
                    failures += 1
                    failure_details.append((H1, H2, H_prod, predicted))

        print(f"    Successes: {successes}, Failures: {failures}")

        if failures > 0:
            # Analyze when it fails
            fail_H1 = Counter(f[0] for f in failure_details)
            fail_H2 = Counter(f[1] for f in failure_details)
            print(f"    Failing H1 values: {dict(fail_H1)}")
            print(f"    Failing H2 values: {dict(fail_H2)}")

            # Key: does it fail ONLY when T1 has cycles (H1 > 1)?
            only_cycle = all(f[0] > 1 for f in failure_details)
            print(f"    Fails only when H1 > 1 (T1 has cycles)? {only_cycle}")

            # Does it fail only when T1 is the 3-cycle (H=3)?
            only_3cycle = all(f[0] == 3 for f in failure_details)
            print(f"    Fails only when H1 = 3 (3-cycle)? {only_3cycle}")

            # What's the actual formula in failure cases?
            for h1, h2, hp, pred in failure_details[:5]:
                print(f"      H1={h1}, H2={h2}: actual={hp}, predicted={pred}, ratio={hp/pred:.4f}")

    # ========================================
    # PART 4: The CORRECT lex product formula
    # ========================================
    print(f"\n{'='*70}")
    print("PART 4: FINDING THE CORRECT LEX PRODUCT FORMULA")
    print("  For T1 lex T2 where T1 is the 3-cycle:")
    print("  H = 45 when H2=1, and the formula H1*H2^3 gives 3.")
    print("  What IS 45 in terms of T1 and T2?")
    print(f"{'='*70}")

    # The 3-cycle lex T_2: 6 vertices, arranged as 3 groups of 2.
    # Within each group: T2 determines internal arcs (1 arc per group).
    # Between groups: 3-cycle determines direction.
    # The 45 ham paths traverse all 6 vertices.

    # Counting directly:
    # A ham path must visit all 3 groups. Order of groups: there are 3! = 6 orderings.
    # But only 3 of these are consistent with the 3-cycle (cyclic orderings).
    # Wait, the 3-cycle gives arcs 0->1, 1->2, 2->0. So valid group orderings
    # are: 0,1,2 or 0,2,1 or 1,0,2... actually ALL group orderings are valid
    # in the lex product because the lex product is a COMPLETE tournament.

    # Within each group of 2: T2 determines the internal order.
    # If T2 is transitive (H=1): 1 internal path.
    # If T2 has both orientations... wait, n2=2 means T2 has 1 arc.
    # Transitive: either 0->1 or 1->0. Either way, H=1.

    # So for T1=3-cycle, T2=transitive-2: the lex product has 6 vertices
    # in 3 groups of 2. H=45 = max_H(6).

    # This is a "blow-up" construction: replace each vertex of the 3-cycle
    # with a pair. The resulting tournament is "doubly regular" (score (2,2,2,3,3,3)).

    print(f"  The 3-cycle lex T_2 is a 'blow-up' of the 3-cycle.")
    print(f"  Each vertex becomes a pair. H=45 because this blow-up")
    print(f"  creates the maximally regular tournament at n=6.")
    print(f"")
    print(f"  GENERAL BLOW-UP CONJECTURE:")
    print(f"  The blow-up of T_p (Paley) by factor k gives:")
    print(f"  n = p*k vertices, and achieves high H.")
    print(f"")
    print(f"  At n=6 = 3*2: blow-up of T_3 by factor 2 -> H=45=max")
    print(f"  At n=7: Paley T_7 directly -> H=189=max")
    print(f"  At n=9 = 3*3: blow-up of T_3 by factor 3 -> H=?")

    # Build T_3 lex T_3 (= blow-up of 3-cycle by factor 3)
    A_3cycle = np.array([[0,1,0],[0,0,1],[1,0,0]])

    # All T2 at n=3
    for b2 in range(2**(3*(3-1)//2)):
        A2 = bits_to_adj(b2, 3)
        H2 = compute_H(A2.tolist(), 3)
        A_prod, n_p = lex_product(A_3cycle, 3, A2, 3)
        H_prod = compute_H(A_prod.tolist(), n_p)
        scores = tuple(sorted(A_prod.sum(axis=1).astype(int)))
        print(f"    3-cycle lex T(H={H2}): n={n_p}, H={H_prod}, scores={scores}")

    # What is max_H(9)? It's 3357.
    print(f"\n  max_H(9) = 3357 (from OEIS)")
    print(f"  3-cycle lex T_3(H=3) gives H=? (computed above)")
    print(f"  3-cycle lex transitive-3(H=1) gives H=? (computed above)")

    # ========================================
    # PART 5: Recursive maximizer structure
    # ========================================
    print(f"\n{'='*70}")
    print("PART 5: RECURSIVE MAXIMIZER STRUCTURE")
    print("  max_H: 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095")
    print("  Looking for recursive patterns")
    print(f"{'='*70}")

    mH = [1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095]

    # Check: max_H(n) = max_H(n-1) * (n-1)?
    print(f"\n  max_H(n) / max_H(n-1):")
    for i in range(1, len(mH)):
        ratio = mH[i] / mH[i-1]
        print(f"    a({i+1})/a({i}) = {mH[i]}/{mH[i-1]} = {ratio:.4f}")

    # Check: max_H(n) = (2n-3) * max_H(n-1) for odd n?
    print(f"\n  Check max_H(n) = (2n-3) * max_H(n-1):")
    for i in range(2, len(mH)):
        n = i + 1
        predicted = (2*n - 3) * mH[i-1]
        print(f"    n={n}: (2n-3)*a(n-1) = {predicted}, actual = {mH[i]}, "
              f"match = {predicted == mH[i]}")

    # Check: max_H(n) = n * max_H(n-1) / something?
    print(f"\n  Check max_H(n) * 2^(n-2) / n!:")
    for i in range(len(mH)):
        n = i + 1
        val = mH[i] * 2**(n-2) / math.factorial(n) if n > 0 else 0
        print(f"    n={n}: max_H * 2^(n-2) / n! = {val:.6f}")

    # This should approach e ≈ 2.718 by Szele's theorem
    print(f"\n  Szele: max_H / (n!/2^(n-1)) -> e = {math.e:.6f}")
    for i in range(len(mH)):
        n = i + 1
        if n > 1:
            ratio = mH[i] / (math.factorial(n) / 2**(n-1))
            print(f"    n={n}: ratio = {ratio:.6f}")

    # ========================================
    # PART 6: Does T_3 lex max(n-3) = max(n) for even n?
    # ========================================
    print(f"\n{'='*70}")
    print("PART 6: BLOW-UP PATTERN")
    print("  Does T_3 lex max_T(k) give max_H(3k) for k = 1,2,...?")
    print(f"{'='*70}")

    # k=1: T_3 lex T_1 = T_3 itself. max_H(3) = 3 = H(T_3). ✓
    # k=2: T_3 lex max_T(2). max_T(2) = transitive(H=1).
    #       H = 45 = max_H(6). ✓ (verified)
    # k=3: T_3 lex max_T(3). max_T(3) = T_3 (3-cycle, H=3).
    #       H = ? (computed above). Does it equal max_H(9) = 3357?

    print(f"  k=1: T_3 lex T_1 = T_3. H=3 = max_H(3). CHECK")
    print(f"  k=2: T_3 lex max_T_2 = T_3 lex trans_2. H=45 = max_H(6). CHECK")
    print(f"  k=3: T_3 lex T_3 = ? H=? vs max_H(9) = 3357")

    # We computed T_3 lex T_3 above. Let me find the H value.
    A_3cycle_list = A_3cycle.tolist()
    A_prod_33, n_33 = lex_product(A_3cycle, 3, A_3cycle, 3)
    H_33 = compute_H(A_prod_33.tolist(), 9)
    print(f"    T_3 lex T_3: H = {H_33}, max_H(9) = 3357")
    print(f"    Match? {H_33 == 3357}")

    if H_33 != 3357:
        print(f"    Ratio: {H_33/3357:.4f}")
        print(f"    The blow-up does NOT give the maximizer at n=9!")
        print(f"    The Paley tournament at n=7 has H=189, not a blow-up.")
        print(f"    So the blow-up construction only works for specific n values.")

    # ========================================
    # PART 7: Iterated lex products
    # ========================================
    print(f"\n{'='*70}")
    print("PART 7: ITERATED LEX PRODUCTS")
    print("  T lex T lex T ... — what happens?")
    print(f"{'='*70}")

    # T_2 lex T_2 = 4-vertex tournament. H?
    trans2 = np.array([[0,1],[0,0]])  # 0->1
    A_22, n22 = lex_product(trans2, 2, trans2, 2)
    H_22 = compute_H(A_22.tolist(), 4)
    print(f"  T_2 lex T_2: n=4, H={H_22}, max_H(4)=5")

    # T_2 lex T_2 lex T_2 = 8-vertex
    A_222, n222 = lex_product(A_22, 4, trans2, 2)
    H_222 = compute_H(A_222.tolist(), 8)
    print(f"  T_2 lex T_2 lex T_2: n=8, H={H_222}, max_H(8)=661")

    # (T_2 lex T_2) lex T_2 vs T_2 lex (T_2 lex T_2)
    A_2_22, n_2_22 = lex_product(trans2, 2, A_22, 4)
    H_2_22 = compute_H(A_2_22.tolist(), 8)
    print(f"  T_2 lex (T_2 lex T_2): n=8, H={H_2_22}")
    print(f"  Associativity: {H_222} vs {H_2_22}, equal? {H_222 == H_2_22}")

    # 3-cycle iterated: T_3 lex T_3 lex T_3 = 27 vertices (too large)
    # But: T_3 lex T_2 lex T_2 = 12 vertices (manageable?)
    # Actually n=12 is too large for brute-force H computation.

    # What about the FORMULA for iterated lex?
    # If H(T lex T) = H(T) * H(T)^|T| = H(T)^{|T|+1}
    # Then H(T^{lex k}) = H(T)^{|T|^{k-1} + |T|^{k-2} + ... + 1}
    # = H(T)^{(|T|^k - 1)/(|T| - 1)}

    # For T_2 (|T|=2, H=1): H(T_2^{lex k}) = 1^anything = 1.
    # But actual H(T_2 lex T_2) = H_22 which we computed.
    # If the formula H(T1 lex T2) = H(T1)*H(T2)^|T1| held universally:
    # H(T_2 lex T_2) = 1 * 1^2 = 1. But actual is 1. OK, trivially works.

    # For T_3 (|T|=3, H=3): H(T_3^{lex 2}) = H(T_3)*H(T_3)^3 = 3*27 = 81.
    # But actual H(T_3 lex T_3) = H_33 = ? (computed above)
    print(f"\n  If formula held: H(T_3 lex T_3) = 3 * 3^3 = 81")
    print(f"  Actual: {H_33}")
    print(f"  Formula {'holds' if H_33 == 81 else 'FAILS'}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
