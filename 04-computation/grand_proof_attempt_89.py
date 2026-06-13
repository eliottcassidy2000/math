#!/usr/bin/env python3
"""
grand_proof_attempt_89.py -- opus-2026-03-14-S89

Attempt to PROVE the Grand Fourier Level Theorem:
  E_{2k}/E_0 = 2(n-2k)^k / P(n,2k)

Strategy: Express E_{2k} in terms of combinatorial sums over tournaments,
then use the symmetry structure of Walsh coefficients.

Key insight: E_{2k} = sum_{|S|=2k} c_hat_S^2
where c_hat_S = (1/2^m) sum_T H(T) chi_S(T)

So E_{2k} = (1/2^{2m}) sum_{|S|=2k} [sum_T H(T) chi_S(T)]^2
           = (1/2^{2m}) sum_{|S|=2k} sum_{T,T'} H(T)H(T') chi_S(T) chi_S(T')
           = (1/2^{2m}) sum_{T,T'} H(T)H(T') sum_{|S|=2k} chi_S(T) chi_S(T')

The inner sum: sum_{|S|=2k} chi_S(T) chi_S(T')
= sum_{|S|=2k} (-1)^{S.(T xor T')}
= sum_{|S|=2k} (-1)^{|S cap D|}  where D = T xor T'

This equals the Krawtchouk polynomial K_2k(|D|, m) = sum_j (-1)^j C(|D|,j) C(m-|D|,2k-j)

And sum over all levels:
sum_{|S|=2k} chi_S(T) chi_S(T') = K_{2k}(d, m) where d = dist(T,T') = |T xor T'|

So: E_{2k} = (1/2^{2m}) sum_{T,T'} H(T)H(T') K_{2k}(d(T,T'), m)

This is a distance-based correlation! The Grand Theorem says this equals
E_0 * 2(n-2k)^k / P(n,2k).
"""

from fractions import Fraction
from math import comb, factorial
from itertools import combinations


def krawtchouk(k, d, m):
    """Krawtchouk polynomial K_k(d, m) = sum_j (-1)^j C(d,j) C(m-d,k-j)."""
    total = 0
    for j in range(min(k, d) + 1):
        total += (-1)**j * comb(d, j) * comb(m - d, k - j)
    return total


def falling_factorial(n, k):
    result = 1
    for i in range(k):
        result *= (n - i)
    return result


def main():
    print("="*70)
    print("GRAND THEOREM PROOF ATTEMPT")
    print("opus-2026-03-14-S89")
    print("="*70)

    # Part 1: Verify the Krawtchouk formulation
    print("\nPart 1: Krawtchouk polynomial values at m=C(n,2)")
    print("="*70)

    for n in range(3, 8):
        m = n * (n - 1) // 2
        print(f"\n  n={n}, m={m}:")
        K_max = (n - 1) // 2
        for k2 in range(0, m + 1, 2):  # level = k2
            if k2 // 2 > K_max:
                break
            kk = k2
            # The Krawtchouk sum for level k2 evaluated at d=0
            # K_{k2}(0, m) = C(m, k2) (all j=0 term)
            val = krawtchouk(kk, 0, m)
            expected = comb(m, kk)
            print(f"    K_{kk}(0, m) = {val}, C(m,{kk}) = {expected}, match = {val == expected}")

    # Part 2: The key identity we need to prove
    print("\nPart 2: The identity we need")
    print("="*70)

    # E_{2k} = sum_{|S|=2k} c_hat_S^2
    # E_0 = Mean^2 = (n!/2^{n-1})^2
    # The Grand Theorem says: E_{2k}/E_0 = 2(n-2k)^k / P(n,2k)

    # Equivalently: sum_{|S|=2k} c_hat_S^2 = 2(n-2k)^k * Mean^2 / P(n,2k)

    # What we KNOW about c_hat_S:
    # For |S|=2: c_hat_S != 0 iff the pair S is adjacent in K_n.
    #   If adjacent: |c_hat_S| = Mean/m = (n!/2^{n-1}) / (n(n-1)/2) = n!/(2^{n-2}*n*(n-1))
    # For |S|=2: there are C(n,2) possible pairs, all adjacent -> C(n,2) = m nonzero.
    #   Wait, not all pairs are adjacent. Adjacent means sharing a vertex.
    #   No -- in the tournament Fourier analysis, S is a SUBSET of edges.
    #   Actually all edge pairs {e} for |S|=1 are "level 1" which vanishes.
    #   For |S|=2 (two edges), "adjacent" means the two edges share a vertex.

    # Let me think about this more carefully.
    # The domain is {0,1}^m where m = C(n,2). Each coordinate is an edge.
    # S is a subset of {1,...,m}, representing a subset of edges of K_n.
    # chi_S(T) = product_{e in S} (-1)^{T_e}

    # The Walsh coefficient c_hat_S = E[H(T) chi_S(T)] (uniform over tournaments)
    # = (1/2^m) sum_T H(T) chi_S(T)

    # The Grand Theorem gives the ENERGY at each level.
    # E_{2k} = sum_{|S|=2k} c_hat_S^2

    # For level 2: |S| = 2, meaning S is a PAIR of edges.
    # c_hat_{e,f} = E[H(T) (-1)^{T_e + T_f}]
    # This measures the correlation between H and the orientations of edges e and f.

    # KEY FACT (verified): c_hat_{e,f} != 0 iff e and f are adjacent (share a vertex).
    # Number of adjacent pairs = n * C(n-1, 2) = n(n-1)(n-2)/2
    # Wait: for each vertex, C(n-1, 2) pairs of edges through it.
    # Total = n * C(n-1,2). But each adjacent pair is counted once (they share exactly one vertex
    # unless they share two, which can't happen for distinct edges).
    # Actually: two edges can share 0 or 1 vertices. If they share 1, they're adjacent.
    # For vertex v, edges through v: (n-1) edges. Pairs: C(n-1, 2).
    # Total adjacent pairs = n * C(n-1, 2) / ... no, each adjacent pair is counted once
    # per shared vertex, and they share exactly one vertex.
    # So total = n * C(n-1, 2) = n * (n-1)(n-2)/2.

    # For n=7: 7 * 15 = 105. Matches!

    print("  Adjacent pair counts:")
    for n in range(3, 10):
        m = n * (n - 1) // 2
        adj_pairs = n * comb(n-1, 2)
        total_pairs = comb(m, 2)
        print(f"    n={n}: m={m}, adj pairs = {adj_pairs}, total pairs = {total_pairs}")

    # E_2 = adj_pairs * |c_hat|^2
    # We know |c_hat_{e,f}| = Mean/m for adjacent pairs.
    # So E_2 = n*C(n-1,2) * (Mean/m)^2 = n*(n-1)(n-2)/2 * Mean^2/(m^2)
    # E_2/E_0 = n*(n-1)(n-2)/2 / m^2 = n*(n-1)(n-2)/2 / [n(n-1)/2]^2
    #         = n*(n-1)*(n-2)/2 * 4/[n^2*(n-1)^2]
    #         = 2(n-2)/(n(n-1))
    # = 2(n-2)/P(n,2). This matches E_{2k}/E_0 = 2(n-2k)^k/P(n,2k) with k=1!

    print("\n  Level-2 proof:")
    for n in range(3, 10):
        m = n * (n - 1) // 2
        adj = n * comb(n - 1, 2)
        e2_over_e0 = Fraction(adj, m * m)
        formula = Fraction(2 * (n - 2), n * (n - 1))
        print(f"    n={n}: adj/m^2 = {e2_over_e0}, formula = {formula}, match = {e2_over_e0 == formula}")

    # Part 3: Level-4 analysis
    print("\nPart 3: Level-4 structure")
    print("="*70)

    # For |S|=4 (four edges forming a subgraph of K_n):
    # c_hat_S != 0 iff the subgraph formed by S has a specific structure.
    # From our earlier analysis at n=7:
    #   Type A: S forms a Hamiltonian path on 5 vertices (P_5 subgraph)
    #   Type B: S forms a "double P_3" (two disjoint paths of length 2)

    # Count Type A: number of P_5 subgraphs of K_n
    # A P_5 on vertices {v1,...,v5} is a path v1-v2-v3-v4-v5.
    # Choose 5 vertices: C(n,5). Number of Hamilton paths on K_5: 5!/2 = 60.
    # So #Type_A = C(n,5) * 60.

    # Count Type B: two disjoint P_3's.
    # A P_3 = path of length 2 on 3 vertices. Choose first P_3: C(n,3)*3 (3 paths in K_3).
    # Choose second P_3 from remaining n-3 vertices: C(n-3,3)*3.
    # Divide by 2 for unordered pair: [C(n,3)*3 * C(n-3,3)*3] / 2.
    # Or: #Type_B = 9 * C(n,3) * C(n-3,3) / 2.

    # Actually P_3 on 3 vertices {a,b,c}: the paths are a-b-c, a-c-b, b-a-c = 3 paths.
    # So #P_3 subgraphs of K_n = C(n,3) * 3.
    # For disjoint pair: C(n,3)*3 * C(n-3,3)*3 / 2 = 9*C(n,3)*C(n-3,3)/2.

    # Wait, I need to be more careful. The "edge subset" S has 4 edges.
    # P_5 has 4 edges. A double P_3 has 2+2=4 edges.

    # For n=7:
    # Type A count: C(7,5)*60 = 21*60 = 1260. Matches!
    # Type B count: 9*C(7,3)*C(4,3)/2 = 9*35*4/2 = 630. Matches!
    # Total nonzero at level 4: 1260 + 630 = 1890.

    print("  Level-4 nonzero subset counts:")
    for n in range(5, 10):
        type_a = comb(n, 5) * 60 if n >= 5 else 0
        type_b = 9 * comb(n, 3) * comb(n - 3, 3) // 2 if n >= 6 else 0
        print(f"    n={n}: Type A (P5) = {type_a}, Type B (2*P3) = {type_b}, total = {type_a + type_b}")

    # From our data: |c_hat_A| = Mean * alpha, |c_hat_B| = Mean * 2*alpha
    # where alpha satisfies:
    # E_4 = type_a * (Mean*alpha)^2 + type_b * (Mean*2*alpha)^2
    # E_4/E_0 = type_a * alpha^2 + type_b * 4*alpha^2

    # Grand Theorem: E_4/E_0 = 2(n-4)^2/P(n,4)

    # So: alpha^2 * (type_a + 4*type_b) = 2(n-4)^2/P(n,4)
    # alpha^2 = 2(n-4)^2 / [P(n,4) * (type_a + 4*type_b)]

    print("\n  Level-4 coefficient analysis:")
    for n in range(5, 10):
        type_a = comb(n, 5) * 60 if n >= 5 else 0
        type_b = 9 * comb(n, 3) * comb(n - 3, 3) // 2 if n >= 6 else 0
        pn4 = falling_factorial(n, 4)
        if pn4 > 0:
            e4_ratio = Fraction(2 * (n-4)**2, pn4)
            count_sum = type_a + 4 * type_b
            if count_sum > 0:
                alpha_sq = e4_ratio / count_sum
                # Compare with 1/(12*C(n,4))^2? No, let's see.
                print(f"    n={n}: E4/E0={e4_ratio}, A+4B={count_sum}, alpha^2={alpha_sq} = {float(alpha_sq):.10f}")
                # alpha = |c_hat_A|/Mean
                # For n=7: alpha_sq should be (3/16 / 78.75)^2... no.
                # Actually alpha = |c_hat_A| / Mean.
                # From data: |c_hat_A| = 3/16 at n=7, Mean = 315/4 = 78.75
                # alpha = (3/16)/(315/4) = 3/(16*315/4) = 3*4/(16*315) = 12/5040 = 1/420
                # alpha^2 = 1/176400
                # Check: alpha^2*(1260+4*630) = 1/176400 * 3780 = 3780/176400 = 1/46.666... = 3/140
                # E4/E0 at n=7 = 2*9/840 = 18/840 = 3/140. YES!
            else:
                print(f"    n={n}: E4/E0={e4_ratio}, only Type A (no B)")
                alpha_sq = e4_ratio / type_a if type_a > 0 else 0
                print(f"    alpha^2 = {alpha_sq}")

    # Part 4: Can we prove the coefficient magnitudes from first principles?
    print("\nPart 4: Coefficient magnitude conjecture")
    print("="*70)

    # For level 2: |c_hat_S| = Mean/m for ALL nonzero S. This is UNIFORM.
    # For level 4: two magnitudes, with |c_hat_B| = 2|c_hat_A|.
    # For level 6: one magnitude (all equal).

    # The uniformity at level 2 comes from the edge-transitive symmetry of K_n.
    # S_n acts on K_n by permuting vertices, hence permuting edges.
    # This action is transitive on edges, hence on single edges.
    # For pairs: S_n is transitive on adjacent pairs (since it's transitive on
    # paths P_3, i.e., ordered triples of vertices).

    # Wait, S_n acts transitively on ORDERED triples, but not on unordered paths.
    # Actually, any P_3 (path a-b-c) can be mapped to any other P_3 (a'-b'-c') by
    # the permutation sending a->a', b->b', c->c' (extended to the other vertices).
    # So S_n is transitive on P_3 subgraphs. Each adjacent edge pair corresponds
    # to a P_3 (the shared vertex is the middle). So all adjacent pairs are equivalent
    # under S_n, hence |c_hat| is the same for all of them.

    # For level 4: S_n is transitive on P_5 subgraphs (Type A) and on double-P_3
    # subgraphs (Type B) separately, but not between types. So two magnitudes.

    # KEY QUESTION: Can we determine |c_hat| for each orbit from the orbit size alone?
    # If E_{2k}/E_0 = 2(n-2k)^k/P(n,2k), and we know the orbit decomposition,
    # we can solve for the magnitudes up to the constraint.

    # For level 2: one orbit (adjacent pairs), size n*C(n-1,2).
    # alpha^2 * n*C(n-1,2) = 2(n-2)/(n(n-1))
    # alpha = sqrt(2(n-2) / (n(n-1) * n*C(n-1,2)))
    # = sqrt(2(n-2) / (n^2*(n-1)*(n-1)(n-2)/2))
    # = sqrt(4 / (n^2*(n-1)^2))
    # = 2/(n(n-1)) = 2/P(n,2) = 1/m

    # So alpha = |c_hat_{e,f}|/Mean = 1/m. Therefore |c_hat| = Mean/m.
    # This is PROVABLE from symmetry + the Grand Theorem!

    print("  Level-2: |c_hat|/Mean = 1/m (follows from symmetry + Grand Theorem)")
    for n in range(3, 8):
        m = n * (n - 1) // 2
        alpha = Fraction(1, m)
        e2 = alpha**2 * n * comb(n-1, 2)
        expected = Fraction(2*(n-2), n*(n-1))
        print(f"    n={n}: alpha=1/{m}, E2/E0 = {e2}, expected = {expected}, match = {e2 == expected}")

    # Part 5: Can we prove the Grand Theorem via the correlation formulation?
    print("\nPart 5: Correlation / Krawtchouk approach")
    print("="*70)

    # E_{2k} = (1/2^{2m}) sum_{T,T'} H(T)H(T') K_{2k}(d(T,T'), m)
    # E_{2k}/E_0 = sum_{T,T'} H(T)H(T') K_{2k}(d, m) / [sum_T H(T)]^2
    # = sum_d f(d) K_{2k}(d, m)
    # where f(d) = sum_{T,T': dist=d} H(T)H(T') / [sum_T H(T)]^2

    # This is a spectral statement about the distance distribution of H.
    # If H were a "nice" function (e.g., a polynomial in the edge variables),
    # f(d) would have known properties.

    # The Grand Theorem is equivalent to:
    # sum_d f(d) K_{2k}(d, m) = 2(n-2k)^k / P(n,2k)
    # for all k = 1, ..., floor((n-1)/2)

    # This is a statement about the "H-weighted distance distribution" of tournaments.
    # It's essentially saying that the Fourier profile of H is determined by a simple formula.

    # Part 6: Check if the Grand Theorem relates to the PERMANENT or IMMANANT
    print("\nPart 6: Connection to permanents and immanants")
    print("="*70)

    # H(T) counts Hamiltonian paths. In matrix terms:
    # H(T) = sum_sigma product_i A_{i,sigma(i)} where sigma ranges over cyclic orderings
    # Actually H(T) = sum over permutations sigma: product_{i=1}^{n-1} A_{sigma(i), sigma(i+1)}
    # This is like a trace of a product, or a permanent-like quantity.

    # The key identity for the variance:
    # E[H^2] = (1/2^m) sum_T H(T)^2
    # E[H]^2 = [(1/2^m) sum_T H(T)]^2

    # Var = E[H^2] - E[H]^2
    # Var/E[H]^2 = E[H^2]/E[H]^2 - 1 = sum_{k>=1} E_{2k}/E_0

    # The Grand Theorem says sum E_{2k}/E_0 = sum 2(n-2k)^k/P(n,2k)
    # = sum (n-2k)^k * (n-2k)! * 2/n!  [after simplifying n!/P(n,2k) = (n-2k)!]

    # So Var/Mean^2 = (2/n!) * sum_{k=1}^K (n-2k)^k * (n-2k)!

    # This is our a(n) formula: Var/Mean^2 = 2*a(n)/n!

    # INSIGHT: Each term (n-2k)^k * (n-2k)! counts something!
    # (n-2k)! is the number of permutations of n-2k objects.
    # (n-2k)^k is a power. So (n-2k)^k * (n-2k)! = number of (k-tuples from [n-2k]) times (n-2k)!
    # Or: it's the number of ways to:
    #   1. Choose a permutation pi of {1,...,n-2k}: (n-2k)! ways
    #   2. For each of k "slots", choose an element from {1,...,n-2k}: (n-2k)^k ways

    # But what does this COUNT in terms of tournament structure?

    # Alternative: (n-2k)^k * (n-2k)! = (n-2k)! * (n-2k)^k
    # = Gamma(n-2k+1) * (n-2k)^k
    # For k=1: (n-2)*(n-2)! = contributions from level-2 Walsh coefficients
    # For k=2: (n-4)^2*(n-4)! = contributions from level-4

    # COMBINATORIAL INTERPRETATION:
    # (n-2k)^k * (n-2k)! = number of words w1...wk from alphabet [n-2k] times perm of [n-2k]
    # = P(n-2k+k, k) * (n-2k)! if the words were distinct... no.

    # Actually (n-2k)^k * (n-2k)! = (n-2k)! * (n-2k)^k
    # For the EGF of E_{2k}/E_0 = 2*(n-2k)^k*(n-2k)! / n!, this is a Poisson-like term.

    # The factor 2/n! * (n-2k)^k * (n-2k)! = 2 * (n-2k)^k / P(n,2k)
    # where P(n,2k) = n!/(n-2k)!

    # This is 2 * [(n-2k)/n] * [(n-2k)/(n-1)] * ... * [(n-2k)/(n-2k+1)]
    # Hmm, not quite. Let's expand:
    # (n-2k)^k / P(n,2k) = (n-2k)^k / [n(n-1)...(n-2k+1)]
    # = prod_{i=0}^{k-1} (n-2k) / prod_{i=0}^{2k-1} (n-i)
    # = (n-2k)^k / prod_{i=0}^{2k-1} (n-i)

    print("  (n-2k)^k / P(n,2k) as product of ratios:")
    for n in [7, 10, 15]:
        K = (n - 1) // 2
        for k in range(1, K + 1):
            j = n - 2 * k
            pn2k = falling_factorial(n, 2 * k)
            ratio = Fraction(j**k, pn2k)
            # Express as product
            factors = []
            for i in range(2 * k):
                factors.append(f"{j}/{n-i}")
            # But that's j^{2k}/P(n,2k), not j^k/P(n,2k).
            # Actually: j^k/P(n,2k) = (j/n) * (j/(n-1)) * ... but only k factors of j
            # vs 2k factors in denominator.
            print(f"    n={n}, k={k}: (n-2k)^k/P(n,2k) = {j}^{k}/{pn2k} = {ratio} = {float(ratio):.8f}")

    # Part 7: The telescoping structure revisited
    print("\nPart 7: Partial fraction decomposition of E_{2k}/E_0")
    print("="*70)

    # E2/E0 = 4/n - 2/(n-1) (proved above)
    # What about E4/E0 = 2(n-4)^2/P(n,4)?
    # 2(n-4)^2 / [n(n-1)(n-2)(n-3)]
    # Partial fractions in n: this is a rational function of n (for fixed k).

    # For k=2: 2(n-4)^2 / [n(n-1)(n-2)(n-3)]
    # = 2(n^2-8n+16) / [n(n-1)(n-2)(n-3)]
    # Partial fractions: A/n + B/(n-1) + C/(n-2) + D/(n-3)
    # A = 2*16/[(-1)(-2)(-3)] = 32/(-6) = -16/3
    # B = 2*9/[1*(-1)(-2)] = 18/2 = 9
    # C = 2*4/[2*1*(-1)] = 8/(-2) = -4
    # D = 2*1/[3*2*1] = 2/6 = 1/3

    # So E4/E0 = -16/(3n) + 9/(n-1) - 4/(n-2) + 1/(3(n-3))

    print("  E4/E0 partial fractions: -16/(3n) + 9/(n-1) - 4/(n-2) + 1/(3(n-3))")
    for n in range(5, 10):
        direct = Fraction(2*(n-4)**2, falling_factorial(n, 4))
        pf = Fraction(-16, 3*n) + Fraction(9, n-1) + Fraction(-4, n-2) + Fraction(1, 3*(n-3))
        print(f"    n={n}: direct={direct}, partial fractions={pf}, match={direct == pf}")

    # General partial fraction for E_{2k}/E_0 = 2(n-2k)^k / prod_{i=0}^{2k-1}(n-i)
    # The poles are at n = 0, 1, 2, ..., 2k-1.
    # Residue at n = r: 2(r-2k)^k / prod_{i!=r, i=0}^{2k-1}(r-i)

    print("\n  General residues for level 2k:")
    for k in range(1, 5):
        print(f"\n  k={k} (level {2*k}):")
        for r in range(2 * k):
            num = 2 * (r - 2 * k) ** k
            denom = 1
            for i in range(2 * k):
                if i != r:
                    denom *= (r - i)
            res = Fraction(num, denom)
            print(f"    Residue at n={r}: {res} = {float(res):.6f}")

    # Part 8: Total Var/Mean^2 as sum of partial fractions
    print("\nPart 8: Total V(n) partial fraction structure")
    print("="*70)

    # V(n) = sum_{k=1}^K E_{2k}/E_0
    # Each E_{2k}/E_0 has partial fractions with poles at 0, 1, ..., 2k-1.
    # As K grows with n, new poles appear.
    # For fixed n, V(n) has contributions from all levels up to 2K.

    # The total residue at n=0 from all levels:
    # sum_{k=1}^K 2(-2k)^k / prod_{i=1}^{2k-1}(-i) = sum_k 2(-2k)^k / ((-1)^{2k-1} (2k-1)!)
    # = sum_k 2(-2k)^k / (-(2k-1)!)
    # = sum_k -2(-2k)^k / (2k-1)!
    # = sum_k -2(-1)^k (2k)^k / (2k-1)!
    # = sum_k -2(-1)^k * 2k * (2k)^{k-1} / (2k-1)!

    # This doesn't simplify obviously. But the DOMINANT contribution is:
    # From k=1: residue at n=0 is 2*(-2)/(-1) = 4 (for E2/E0 = 4/n - 2/(n-1))
    # From k=2: residue at n=0 is -16/3
    # Total at n=0: 4 - 16/3 = -4/3 (for K>=2)

    total_res = {}
    for k in range(1, 8):
        for r in range(2 * k):
            num = 2 * (r - 2 * k) ** k
            denom = 1
            for i in range(2 * k):
                if i != r:
                    denom *= (r - i)
            res = Fraction(num, denom)
            total_res[r] = total_res.get(r, Fraction(0)) + res
        # Print cumulative
        print(f"  After including level {2*k}:")
        for r in sorted(total_res.keys()):
            print(f"    pole at n={r}: cumulative residue = {total_res[r]} = {float(total_res[r]):.6f}")
        print()

    print(f"\n{'='*70}")
    print("DONE -- GRAND THEOREM PROOF ATTEMPT")
    print("="*70)


if __name__ == "__main__":
    main()
