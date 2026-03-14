"""
deep_lead_exploration.py -- kind-pasteur-2026-03-14-S81
DEEP exploration of the most promising leads from S80.
Actually PROVE things and follow connections.

LEAD A: F(T,0) = #{all-descending Ham paths}
LEAD B: The 1/3 ratio — prove Var(H)/Mean(H)^2 → 1/3
LEAD C: Tournament composition and H factorization
LEAD D: The H mod 8 structure — deeper than just "odd"
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

def fwd_poly(A, n):
    """Forward polynomial: coefficients[k] = #{paths with fwd=k}."""
    poly = Counter()
    for perm in permutations(range(n)):
        valid = all(A[perm[i]][perm[i+1]] for i in range(n-1))
        if valid:
            fwd = sum(1 for i in range(n-1) if perm[i] < perm[i+1])
            poly[fwd] += 1
    return poly

def main():
    print("=" * 70)
    print("DEEP LEAD EXPLORATION — kind-pasteur-2026-03-14-S81")
    print("=" * 70)

    # ============================================================
    # LEAD A: F(T,0) AND F(T,n-1) — EXTREME EVALUATIONS
    # ============================================================
    print(f"\n{'='*70}")
    print("LEAD A: F(T,0) AND F(T,n-1) — PATH EXTREMES")
    print("  F(T,0) = #{paths with 0 forward steps = all descending}")
    print("  F(T,n-1) = #{paths with n-1 forward steps = all ascending}")
    print("  CLAIM: F(T,0) + F(T,n-1) = 1 or 2 always.")
    print(f"{'='*70}")

    for n in [3, 4, 5, 6]:
        m = n*(n-1)//2
        N = 2**m

        f0_dist = Counter()
        fn1_dist = Counter()
        sum_dist = Counter()

        count = 0
        for bits in range(N):
            count += 1
            if n >= 6 and count > 10000: break

            A = bits_to_adj(bits, n)
            fp = fwd_poly(A, n)
            f0 = fp.get(0, 0)
            fn1 = fp.get(n-1, 0)

            f0_dist[f0] += 1
            fn1_dist[fn1] += 1
            sum_dist[f0 + fn1] += 1

        print(f"\n  n={n}:")
        print(f"    F(T,0) distribution: {dict(sorted(f0_dist.items()))}")
        print(f"    F(T,n-1) distribution: {dict(sorted(fn1_dist.items()))}")
        print(f"    F(T,0)+F(T,n-1) distribution: {dict(sorted(sum_dist.items()))}")

        # THEOREM: F(T,0) = 1 iff T has a Hamiltonian path that is all-descending
        # An all-descending path: p[0] > p[1] > ... > p[n-1] AND T[p[i]][p[i+1]]=1
        # This means the path visits vertices in DECREASING order, and each arc
        # goes from higher to lower. This is a path in the "upper" subtournament.
        # F(T,0) <= 1 because there's at most one decreasing permutation.
        # WAIT: F(T,0) counts paths where ALL steps are descending.
        # A path (p0,p1,...,p_{n-1}) has fwd=0 iff p0>p1>...>p_{n-1}.
        # There's only ONE such permutation: (n-1, n-2, ..., 0).
        # So F(T,0) = 1 iff T[n-1][n-2]=1, T[n-2][n-3]=1, ..., T[1][0]=1.
        # i.e., F(T,0) = 1 iff T contains the path n-1 -> n-2 -> ... -> 0.

        # Similarly: F(T,n-1) = 1 iff T contains 0 -> 1 -> 2 -> ... -> n-1.

        # Both can be 1 simultaneously (if T contains BOTH paths).
        # Neither can exceed 1.

        # THEOREM: F(T,0), F(T,n-1) in {0,1}, and they're indicator functions
        # for specific Hamiltonian paths.

        # Verify:
        assert all(k <= 1 for k in f0_dist.keys()), f"F(T,0) > 1 at n={n}!"
        assert all(k <= 1 for k in fn1_dist.keys()), f"F(T,n-1) > 1 at n={n}!"
        print(f"    VERIFIED: F(T,0), F(T,n-1) in {{0,1}}")
        print(f"    F(T,0)=1 iff T contains (n-1)->...->0")
        print(f"    F(T,n-1)=1 iff T contains 0->1->...->n-1")

    # ============================================================
    # LEAD A2: THE FULL F-POLYNOMIAL AT EXTREME POINTS
    # ============================================================
    print(f"\n{'='*70}")
    print("LEAD A2: PALINDROME PROPERTY OF F(T,x)")
    print("  THM-030 (proved): the transfer matrix is symmetric")
    print("  Does this imply F(T,x) = x^{n-1} F(T, 1/x)?")
    print(f"{'='*70}")

    for n in [4, 5]:
        m = n*(n-1)//2
        palindrome_count = 0
        total = 0

        for bits in range(2**m):
            A = bits_to_adj(bits, n)
            fp = fwd_poly(A, n)
            total += 1

            # Check palindrome: fp[k] = fp[n-1-k] for all k
            is_palin = all(fp.get(k, 0) == fp.get(n-1-k, 0) for k in range(n))
            if is_palin:
                palindrome_count += 1

        print(f"  n={n}: {palindrome_count}/{total} palindromic F-polynomials")
        print(f"    ALL palindromic? {palindrome_count == total}")

    # ============================================================
    # LEAD B: THE 1/3 RATIO — Var(H)/Mean(H)^2
    # ============================================================
    print(f"\n{'='*70}")
    print("LEAD B: THE 1/3 RATIO — DEEP ANALYSIS")
    print("  Var(H)/Mean(H)^2 ≈ 1/3 universally.")
    print("  From Fourier: H = mean + sum of level-2 terms + level-4 terms + ...")
    print("  Var(H) = sum_{|S|>0} H_hat(S)^2 (Parseval)")
    print("  Mean(H)^2 = H_hat(0)^2")
    print("  So Var/Mean^2 = (sum_{|S|>0} H_hat^2) / H_hat(0)^2")
    print("  = (non-constant energy) / (constant energy)")
    print("  ≈ 0.25/0.75 = 1/3 !!")
    print(f"{'='*70}")

    for n in [3, 4, 5]:
        m = n*(n-1)//2
        N = 2**m

        H_values = np.zeros(N)
        for bits in range(N):
            A = bits_to_adj(bits, n)
            H_values[bits] = compute_H(A, n)

        # Fourier
        H_hat = H_values.copy()
        for i in range(m):
            step = 1 << (i + 1)
            half = 1 << i
            for j in range(0, N, step):
                for k in range(half):
                    u, v = H_hat[j+k], H_hat[j+k+half]
                    H_hat[j+k], H_hat[j+k+half] = u+v, u-v
        H_hat /= N

        E0 = H_hat[0]**2
        E_nonconst = sum(H_hat[S]**2 for S in range(1, N))
        E_total = E0 + E_nonconst

        ratio_fourier = E_nonconst / E0
        ratio_stat = np.var(H_values) / np.mean(H_values)**2

        print(f"\n  n={n}:")
        print(f"    E_nonconst / E_0 = {E_nonconst:.4f} / {E0:.4f} = {ratio_fourier:.6f}")
        print(f"    Var(H)/Mean(H)^2 = {ratio_stat:.6f}")
        print(f"    Match? {abs(ratio_fourier - ratio_stat) < 1e-10}")
        print(f"    Close to 1/3? {abs(ratio_stat - 1/3) < 0.02}")

    # THE KEY INSIGHT:
    print(f"\n  *** THEOREM: Var(H)/Mean(H)^2 = E_nonconst/E_0 EXACTLY ***")
    print(f"  This is a Fourier identity (Parseval): Var = sum of non-constant Fourier energies.")
    print(f"  And E_nonconst/E_0 ≈ 1/3 because the level-2 energy is ~25% of total,")
    print(f"  and E_0 is ~75%, giving ratio 25/75 = 1/3.")
    print(f"")
    print(f"  WHY IS E_0 ≈ 75%?")
    print(f"  E_0 = mean(H)^2. E_total = mean(H^2). E_0/E_total = mean(H)^2/mean(H^2).")
    print(f"  This is 1/(1 + CV^2) where CV = std/mean.")
    print(f"  If CV^2 = 1/3, then E_0/E_total = 3/4 = 75%. Circular but correct!")

    # ============================================================
    # LEAD C: TOURNAMENT COMPOSITION AND H
    # ============================================================
    print(f"\n{'='*70}")
    print("LEAD C: TOURNAMENT COMPOSITION")
    print("  Given tournaments T1 on {1,...,a} and T2 on {1,...,b},")
    print("  define T1[v <- T2]: replace vertex v with T2.")
    print("  All arcs from T2 to T1\\v get direction from v's arcs.")
    print("  QUESTION: H(T1[v <- T2]) = f(H(T1), H(T2), v)?")
    print(f"{'='*70}")

    # Lexicographic composition: T1 x T2
    # Vertices = {(i,j) : i in V(T1), j in V(T2)}
    # Arc (i1,j1) -> (i2,j2) iff T1[i1][i2]=1, OR (i1=i2 and T2[j1][j2]=1)
    # H(T1 x T2) = H(T1) * H(T2)^a where a = |V(T1)|? NO...
    # Actually for lexicographic product: H(T1 lex T2) = H(T1) * H(T2)^{|T1|}

    print(f"\n  Testing LEXICOGRAPHIC PRODUCT: T1 lex T2")
    print(f"  Vertices (i,j), arc rule: i1->i2 in T1, OR (i1=i2 and j1->j2 in T2)")

    for na, nb in [(3, 2), (2, 3), (3, 3)]:
        if na + nb > 6: continue  # too slow for n>6

        # Build lex product for all pairs of small tournaments
        results = []
        for bits_a in range(2**(na*(na-1)//2)):
            A1 = bits_to_adj(bits_a, na)
            H1 = compute_H(A1, na)

            for bits_b in range(2**(nb*(nb-1)//2)):
                A2 = bits_to_adj(bits_b, nb)
                H2 = compute_H(A2, nb)

                # Build lexicographic product
                n_prod = na * nb
                if n_prod > 7: continue

                A_prod = np.zeros((n_prod, n_prod), dtype=int)
                for i1 in range(na):
                    for j1 in range(nb):
                        for i2 in range(na):
                            for j2 in range(nb):
                                v1 = i1 * nb + j1
                                v2 = i2 * nb + j2
                                if v1 == v2: continue
                                if A1[i1][i2] == 1:
                                    A_prod[v1][v2] = 1
                                elif i1 == i2 and A2[j1][j2] == 1:
                                    A_prod[v1][v2] = 1

                H_prod = compute_H(A_prod.tolist(), n_prod)
                results.append((H1, H2, H_prod))

        # Check: H(T1 lex T2) = H(T1) * H(T2)^na ?
        formula_match = all(hp == h1 * h2**na for h1, h2, hp in results)
        # Or: H(T1 lex T2) = H(T1) * (H(T2))^{na} * something?
        if results:
            ratios = set()
            for h1, h2, hp in results:
                if h1 > 0 and h2 > 0:
                    ratios.add(hp / (h1 * h2**na))

            print(f"\n  {na} lex {nb}: {len(results)} products, H_prod / (H1 * H2^{na}) = {sorted(ratios)[:5]}")

            # Try simpler: H_prod = H1^nb * H2^na?
            ratios2 = set()
            for h1, h2, hp in results:
                if h1 > 0 and h2 > 0:
                    ratios2.add(hp / (h1 * h2))
            print(f"  H_prod / (H1 * H2) = {sorted(ratios2)[:5]}")

    # ============================================================
    # LEAD D: H MOD 8 STRUCTURE
    # ============================================================
    print(f"\n{'='*70}")
    print("LEAD D: H MOD 8 STRUCTURE")
    print("  H mod 2 = 1 always (Redei)")
    print("  H mod 4 = 1 + 2*alpha_1 mod 4 (from OCF)")
    print("  H mod 8 = ?")
    print(f"{'='*70}")

    for n in [4, 5, 6]:
        m = n*(n-1)//2
        N = 2**m

        mod8_dist = Counter()
        mod4_dist = Counter()
        count = 0

        for bits in range(N):
            count += 1
            if n >= 6 and count > 10000: break

            A = bits_to_adj(bits, n)
            H = compute_H(A, n)
            mod8_dist[H % 8] += 1
            mod4_dist[H % 4] += 1

        print(f"\n  n={n}:")
        print(f"    H mod 4: {dict(sorted(mod4_dist.items()))}")
        print(f"    H mod 8: {dict(sorted(mod8_dist.items()))}")

        # H mod 4: should be {1, 3} since H is odd and H = 1 + 2*alpha_1
        # alpha_1 even => H mod 4 = 1, alpha_1 odd => H mod 4 = 3
        # H mod 8: H = 1 + 2*a1 + 4*a2, so H mod 8 = (1 + 2*a1 + 4*a2) mod 8
        # = 1 + 2*(a1 mod 4) + 4*(a2 mod 2) mod 8

        # At n=5: alpha_2 = 0, so H mod 8 = (1 + 2*alpha_1) mod 8
        if n == 5:
            print(f"    n=5: alpha_2=0, so H mod 8 = (1 + 2*alpha_1) mod 8")
            print(f"    alpha_1 mod 4 -> H mod 8: 0->1, 1->3, 2->5, 3->7")
            print(f"    H mod 8 distribution counts cycles of specific alpha_1 residues")

    # ============================================================
    # LEAD E: THE FULL DELETION TOWER
    # ============================================================
    print(f"\n{'='*70}")
    print("LEAD E: THE DELETION TOWER — H(T) → H(T-v) → H(T-v-w) → ...")
    print("  Delete vertices one by one. How does H evolve?")
    print("  Connection to the OCF: H(T) - H(T-v) = 2 * sum_C mu(C)")
    print(f"{'='*70}")

    n = 5
    m = n*(n-1)//2

    # For a few tournaments, trace the full deletion tower
    for bits in [0, 42, 341, 682, 1023]:
        A = bits_to_adj(bits, n)
        H = compute_H(A, n)
        print(f"\n  bits={bits}, H(T)={H}:")

        # Delete each vertex
        for v in range(n):
            # Sub-tournament on V \ {v}
            sub_verts = [u for u in range(n) if u != v]
            sub_A = [[A[i][j] for j in sub_verts] for i in sub_verts]
            H_sub = compute_H(sub_A, n-1)
            delta = H - H_sub
            print(f"    Delete vertex {v}: H(T-{v}) = {H_sub}, delta = {delta}")

        # Sum of deltas = sum_v (H(T) - H(T-v)) = n*H - sum_v H(T-v)
        sum_sub = sum(
            compute_H([[A[i][j] for j in range(n) if j != v] for i in range(n) if i != v], n-1)
            for v in range(n)
        )
        R = sum_sub / H if H > 0 else 0
        print(f"    R(T) = sum H(T-v) / H(T) = {sum_sub}/{H} = {R:.4f}")
        print(f"    n - R = {n - R:.4f}")

    # ============================================================
    # LEAD F: H AND THE TUTTE POLYNOMIAL
    # ============================================================
    print(f"\n{'='*70}")
    print("LEAD F: TUTTE POLYNOMIAL EVALUATIONS")
    print("  For the complete graph K_n: T(2,0) = n! (# acyclic orientations)")
    print("  For a tournament T: the Tutte polynomial at (1,0) counts what?")
    print("  Awan-Bernardi B-polynomial might be the right framework")
    print(f"{'='*70}")

    # The complete graph K_n has Tutte polynomial
    # T(x,y) for K_3: T = x^2 + 3x + y (from standard tables)
    # Evaluations: T(2,0) = 4+6+0 = 10... no, K_3 has 3 edges
    # Actually for K_3: T(x,y) = x^2 + x + y
    # T(2,0) = 4+2+0 = 6 = 3! ✓ (acyclic orientations)
    # T(1,0) = 1+1+0 = 2 (= #spanning trees? No, that's T(1,1))
    # T(1,1) = 1+1+1 = 3 = #spanning trees of K_3

    print(f"  K_3: T(2,0) = 6 = 3! (acyclic orientations = transitive tournaments)")
    print(f"  K_4: T(2,0) = 14 (?)")

    # Actually T(2,0) for K_n = number of acyclic orientations
    # = sum over topological sorts = n! (for complete graph, each acyclic
    # orientation = unique topological sort = permutation)
    # Wait: #acyclic orientations of K_n = n! (each corresponds to a total order)
    # And #acyclic orientations = T(2,0) for the Tutte polynomial

    # BUT: acyclic orientations of K_n = TRANSITIVE tournaments
    # And #transitive tournaments on labeled vertices = n!
    # So T_{K_n}(2,0) = n! ✓

    # What about T_{K_n}(x,0)?
    # T(x,0) = chromatic polynomial / x (for connected graphs)
    # chi_{K_n}(k) = k(k-1)(k-2)...(k-n+1) = falling factorial
    # T(x,0) relates to... actually T(1+x, 0) = chi(x+1)/... complicated

    # The question: is there a Tutte evaluation that gives H(T) for a
    # specific tournament? Probably not directly, since the Tutte polynomial
    # is defined for the GRAPH, not for a specific orientation.

    # But the Awan-Bernardi B-polynomial IS defined for digraphs.
    # B(D; x,y,z) satisfies deletion-contraction and specializes to
    # various invariants at different points.

    print(f"\n  The Awan-Bernardi B-polynomial for digraphs")
    print(f"  would be the right framework to connect H to Tutte theory.")
    print(f"  B(D) at some (x,y,z) specialization should give H(D).")
    print(f"  LEAD: Read arXiv:1610.01839 and find the specialization.")

    # ============================================================
    # LEAD G: AUTOCORRELATION FUNCTION OF H
    # ============================================================
    print(f"\n{'='*70}")
    print("LEAD G: AUTOCORRELATION OF H ON THE HYPERCUBE")
    print("  C(d) = E[H(T) * H(T')] where T,T' differ in exactly d arcs")
    print("  This measures how correlated H values are at distance d")
    print(f"{'='*70}")

    n = 5
    m = n*(n-1)//2
    N = 2**m

    # Compute all H
    all_H = np.zeros(N, dtype=float)
    for bits in range(N):
        A = bits_to_adj(bits, n)
        all_H[bits] = compute_H(A, n)

    mean_H = np.mean(all_H)

    # Autocorrelation at distance d
    print(f"\n  n={n}: Autocorrelation C(d) = E[H(T)*H(T')] for Hamming distance d")

    for d in range(m+1):
        # Sample: for each T, pick a random T' at distance d
        # Actually, exact: C(d) = (1/N) * sum_T (1/C(m,d)) * sum_{T' at dist d} H(T)*H(T')
        # This is expensive. Use Fourier instead:
        # C(d) = sum_{|S|=0,2,4,...} H_hat(S)^2 * K_d(|S|) where K is Krawtchouk
        # Actually simpler: C(d) / mean^2 = sum_k P(k) * rho(k)^d where rho(k) relates to level k
        pass

    # Use Fourier: the autocorrelation at distance d is
    # C(d) = sum_S H_hat(S)^2 * (1 - 2d/m)^{|S|}... no, that's for symmetric random walk
    # Exact: C(d) = sum_S H_hat(S)^2 * K_d(|S|, m) where K_d is Krawtchouk polynomial

    # Simpler approach: just compute directly for d=0,1,2
    H_hat = all_H.copy()
    for i in range(m):
        step = 1 << (i + 1)
        half = 1 << i
        for j in range(0, N, step):
            for k in range(half):
                u, v = H_hat[j+k], H_hat[j+k+half]
                H_hat[j+k], H_hat[j+k+half] = u+v, u-v
    H_hat /= N

    # Autocorrelation via Fourier: C(d)/C(0) = sum_k E_k * (1-2k/m)^... no
    # The pairwise correlation for bit flip at position i:
    # E[H(T) * H(T^i)] = sum_S H_hat(S)^2 * (-1)^{[i in S]}
    # Averaged over i: E_i[E[H*H^i]] = sum_S H_hat(S)^2 * (1 - 2|S|/m)

    avg_correlation = 0
    for S in range(N):
        level = bin(S).count('1')
        avg_correlation += H_hat[S]**2 * (1 - 2*level/m)

    print(f"  Average nearest-neighbor correlation: {avg_correlation:.4f}")
    print(f"  Mean(H)^2 = {mean_H**2:.4f}")
    print(f"  Ratio (correlation/mean^2) = {avg_correlation/mean_H**2:.4f}")
    print(f"  This tells us how much H changes under a single arc flip")

    # The correlation structure reveals the "smoothness" of H on the hypercube
    # High correlation = smooth landscape (consistent with unimodality at n=5)

    print(f"\n{'='*70}")
    print("DONE — LEADS A-G EXPLORED")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
