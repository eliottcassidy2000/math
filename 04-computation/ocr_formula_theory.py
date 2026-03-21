#!/usr/bin/env python3
"""
ocr_formula_theory.py — kind-pasteur-2026-03-21-S16

THEORETICAL DERIVATION OF THE EXACT OCR.

Known:
  OCR(3) = 1, OCR(4) = 1, OCR(5) = 18/19, OCR(6) = 12/13.

Key identity: H = 1 + 2*c3 + (higher cycle terms).
Since S2 = -2*c3 + const (Rao's formula), R^2(S2, H) = R^2(c3, H).

For RANDOM tournaments:
  E[H] = n!/2^{n-1}  (linearity of expectation)
  E[c3] = C(n,3)/4    (each triple has P=1/4 of being a 3-cycle)
  Var(c3) = ?          (computable from pair correlations)
  Cov(H, c3) = ?       (the key quantity)

THEOREM (to prove): For random tournaments on n vertices:
  R^2(c3, H) = Var(2*c3) / Var(H) = 4*Var(c3) / Var(H)

This is because H = 1 + 2*c3 + epsilon where epsilon is uncorrelated with c3
(since epsilon involves 5+ cycle terms which are "orthogonal" to c3 contributions).

Wait — is epsilon uncorrelated with c3? Let me check.
If so, Var(H) = Var(2*c3) + Var(epsilon), and
  R^2 = Var(2*c3) / [Var(2*c3) + Var(epsilon)] = 1 / [1 + Var(epsilon)/Var(2*c3)]

At n=5: R^2 = 18/19, so Var(epsilon)/Var(2*c3) = 1/18.
At n=6: R^2 = 12/13, so Var(epsilon)/Var(2*c3) = 1/12.

QUESTIONS:
1. Is epsilon = H - 1 - 2*c3 truly uncorrelated with c3?
2. What is Var(c3) for random tournaments?
3. What is Var(H) for random tournaments?
4. Can we get a closed form for the ratio?
"""

import numpy as np
from itertools import combinations, permutations
from math import comb, factorial, gcd
from fractions import Fraction
from collections import defaultdict

def build_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def ham_paths_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    key = (mask | (1 << w), w)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

print("=" * 72)
print("  OCR FORMULA: THEORETICAL DERIVATION")
print("  kind-pasteur-2026-03-21-S16")
print("=" * 72)

for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    print(f"\n{'='*72}")
    print(f"  n = {n}")
    print(f"{'='*72}")

    all_H = []
    all_c3 = []
    all_c3_sq = []

    for bits in range(total):
        A = build_tournament(n, bits)
        scores = [sum(A[i]) for i in range(n)]
        c3 = comb(n,3) - sum(s*(s-1)//2 for s in scores)
        H = ham_paths_dp(A, n)
        all_H.append(H)
        all_c3.append(c3)

    H = np.array(all_H, dtype=np.int64)
    c3 = np.array(all_c3, dtype=np.int64)
    N = total

    # Exact moments using Python integers
    sum_H = int(np.sum(H))
    sum_c3 = int(np.sum(c3))
    sum_H2 = int(np.sum(H * H))
    sum_c32 = int(np.sum(c3 * c3))
    sum_Hc3 = int(np.sum(H * c3))

    # Exact means and variances
    E_H = Fraction(sum_H, N)
    E_c3 = Fraction(sum_c3, N)
    E_H2 = Fraction(sum_H2, N)
    E_c32 = Fraction(sum_c32, N)
    E_Hc3 = Fraction(sum_Hc3, N)

    Var_H = E_H2 - E_H**2
    Var_c3 = E_c32 - E_c3**2
    Cov_Hc3 = E_Hc3 - E_H * E_c3

    print(f"  E[H] = {E_H} = {float(E_H):.6f}")
    print(f"  E[c3] = {E_c3} = {float(E_c3):.6f}")
    print(f"  Var(H) = {Var_H} = {float(Var_H):.6f}")
    print(f"  Var(c3) = {Var_c3} = {float(Var_c3):.6f}")
    print(f"  Cov(H, c3) = {Cov_Hc3} = {float(Cov_Hc3):.6f}")

    # R^2 = Cov^2 / (Var_H * Var_c3)
    R2 = Cov_Hc3**2 / (Var_H * Var_c3)
    print(f"  R^2(c3, H) = {R2}")

    # Check: Cov(H, c3) = 2*Var(c3)?
    # If H = 1 + 2*c3 + eps with eps independent of c3:
    #   Cov(H, c3) = Cov(1 + 2*c3 + eps, c3) = 2*Var(c3) + Cov(eps, c3)
    #   If Cov(eps, c3) = 0: Cov(H, c3) = 2*Var(c3)
    expected_cov = 2 * Var_c3
    actual_cov = Cov_Hc3
    print(f"\n  Test: Cov(H, c3) = 2*Var(c3)?")
    print(f"    2*Var(c3) = {expected_cov}")
    print(f"    Cov(H,c3) = {actual_cov}")
    print(f"    Equal: {expected_cov == actual_cov}")

    if expected_cov == actual_cov:
        # Then R^2 = (2*Var(c3))^2 / (Var(H) * Var(c3)) = 4*Var(c3)/Var(H)
        print(f"\n  CONFIRMED: epsilon is uncorrelated with c3!")
        print(f"  R^2 = 4*Var(c3) / Var(H) = {4*Var_c3 / Var_H}")

        # Var(epsilon) = Var(H) - 4*Var(c3)
        Var_eps = Var_H - 4 * Var_c3
        print(f"  Var(epsilon) = Var(H) - 4*Var(c3) = {Var_eps}")
        print(f"  Var(epsilon)/Var(H) = {Var_eps / Var_H}")

        # What IS epsilon?
        eps = H - (1 + 2 * c3)
        E_eps = Fraction(int(np.sum(eps)), N)
        Var_eps_check = Fraction(int(np.sum(eps**2)), N) - E_eps**2
        Cov_eps_c3 = Fraction(int(np.sum(eps * c3)), N) - E_eps * E_c3
        print(f"\n  epsilon = H - 1 - 2*c3:")
        print(f"    E[eps] = {E_eps}")
        print(f"    Var(eps) = {Var_eps_check}")
        print(f"    Cov(eps, c3) = {Cov_eps_c3} (should be 0)")
    else:
        diff = actual_cov - expected_cov
        print(f"    Difference: {diff}")

    # Theoretical Var(c3) for random tournaments
    # c3 = sum_{triples} X_T where X_T = indicator(triple T is a 3-cycle)
    # E[X_T] = 1/4 (2 cyclic orientations out of 8 = 2/8 = 1/4)
    # For two triples T1, T2:
    #   disjoint (share 0 verts): E[X_T1 * X_T2] = 1/16
    #   share 1 vertex: E[X_T1 * X_T2] = ?
    #   share 2 vertices (= share an edge): E[X_T1 * X_T2] = ?

    n_triples = comb(n, 3)
    E_c3_theory = Fraction(n_triples, 4)
    print(f"\n  Theoretical E[c3] = C(n,3)/4 = {E_c3_theory}")
    print(f"  Match: {E_c3 == E_c3_theory}")

    # Var(c3) = sum_{T1,T2} Cov(X_T1, X_T2)
    # = n_triples * (1/4 - 1/16) + sum_{pairs sharing k>=1 verts} [E[X1*X2] - 1/16]

    # Two triples sharing 0 vertices: Cov = 0 (independent arcs)
    # Two triples sharing 1 vertex: e.g. {a,b,c} and {a,d,e}
    #   Share 0 arcs. All 5 arcs independent. E[X1*X2] = 1/4 * 1/4 = 1/16
    #   Cov = 0

    # Two triples sharing 2 vertices (= 1 common edge): e.g. {a,b,c} and {a,b,d}
    #   Share arc (a,b). Total arcs: a-b (shared), a-c, b-c, a-d, b-d = 5 arcs.
    #   P(both are 3-cycles) = P(a->b->c->a AND a->b->d->a OR a->c->b->a AND ...)
    #   Each 3-cycle has 2 orientations. For fixed a-b direction:
    #     {a,b,c} is cycle iff either a->b->c->a or a->c->b->a
    #     {a,b,d} is cycle iff either a->b->d->a or a->d->b->a
    #   Condition on direction of a-b:
    #   Case a->b: {a,b,c} cycle iff b->c->a (P=1/4) or a->c,c->b (impossible, b->a not in T)
    #     Wait: with a->b:
    #       a->b->c->a means b->c AND c->a: P = 1/4
    #       a->c->b->a means c->b AND a->c: but a->b is set, need c->b AND a->c: P = 1/4
    #     Total: P(cycle|a->b) = P(b->c,c->a) + P(a->c,c->b) - P(both) = 1/4 + 1/4 - 0 = 1/2
    #     Similarly for {a,b,d}: P(cycle|a->b) = 1/2
    #     P(both cycles | a->b) = P(c3_abc|a->b) * P(c3_abd|a->b | arcs c,d independent)
    #       If c and d don't share an arc: = 1/2 * 1/2 = 1/4
    #     Case a<-b: by symmetry, same = 1/4
    #     Overall: P(both cycles) = 1/4
    #     E[X1*X2] = 1/4
    #     Cov = 1/4 - 1/16 = 3/16

    # How many pairs share exactly 2 vertices (1 edge)?
    # Each edge is in C(n-2, 1) triples. C(n,2) edges. Pairs sharing an edge:
    # For each edge, C(n-2,1) choose 2 = C(n-2, 2) pairs? No.
    # For each edge e, the triples containing e form a set of size (n-2).
    # Pairs of triples sharing e: C(n-2, 2) wait no. C(n-2, 1) triples per edge.
    # Wait: triple {a,b,c} contains edges ab, ac, bc. Edge ab is in triples {a,b,x} for x != a,b.
    # So n-2 triples per edge. Pairs of triples sharing edge ab: C(n-2, 2)?
    # No: each pair (T1,T2) sharing ab corresponds to choosing two different x values.
    # That's C(n-2, 1) * C(n-3, 1) / 2 ... no, just C(n-2, 2) = choosing 2 vertices from remaining n-2.
    # Wait: {a,b,c} and {a,b,d} share edge ab (and vertices a,b). c and d are different.
    # Number of such pairs for edge ab: choose c,d from V \ {a,b}: C(n-2, 2).
    # But this counts each pair once (unordered).
    # Total pairs of triples sharing an edge: C(n,2) * C(n-2, 2)...
    # No wait. We need pairs sharing exactly one edge (2 vertices).
    # Actually pairs sharing 2 vertices always share exactly 1 edge.
    # Number of pairs of triples sharing exactly 2 vertices:
    # Choose the 2 shared vertices: C(n,2). For each, choose the 2 distinct extra vertices: C(n-2, 2).
    # But this overcounts! Pair ({a,b,c}, {a,b,d}) is counted once (shared pair = {a,b}).
    # Wait, what if the shared pair is {a,c}? Then triples are {a,c,b} and {a,c,d}, same thing.
    # Each pair of triples sharing 2 verts has exactly 1 shared edge. Counted exactly once.
    # So: num_edge_sharing_pairs = C(n,2) * C(n-2, 2)? No:
    # Actually for shared vertices {a,b}: triples through them are {a,b,c} for c in V\{a,b}.
    # Pairs of such triples: C(n-2, 2) (choose 2 of the remaining vertices for c1, c2).
    # But wait, each pair ({a,b,c1}, {a,b,c2}) is counted once per shared edge.
    # But these triples share ONLY the edge ab (not ac1 or bc1 with ac2 or bc2).
    # So number of pairs sharing exactly 2 vertices = C(n,2)*C(n-2,2)/1.
    # Hmm, actually no overcounting because we're iterating over the shared edge.
    # Total: C(n,2) * C(n-2,1) * (C(n-2,1)-1) / 2  ... Let me just count it.
    # Shared vertices = {a,b}. Third vertex for T1: any of n-2. Third for T2: any of n-3.
    # Unordered: C(n-2, 2). So total = C(n,2) * C(n-2, 2). But this DOES include pairs
    # that might share more than 2 vertices (impossible for triples: max shared = 3 = identical).
    # So: pairs sharing exactly 2 vertices = C(n,2) * C(n-2, 2).

    n_edge_pairs = comb(n, 2) * comb(n-2, 2) if n >= 4 else 0

    # Var(c3) = sum of variances + 2*sum of covariances
    # = n_triples * (1/4)(3/4)  [variance of each indicator]
    # + 2 * n_edge_pairs * 3/16  [covariance for edge-sharing pairs]
    # + 2 * (remaining pairs) * 0  [independent pairs have 0 covariance]

    Var_c3_theory = Fraction(n_triples * 3, 16) + 2 * Fraction(n_edge_pairs * 3, 16)
    print(f"\n  Theoretical Var(c3):")
    print(f"    n_triples = {n_triples}")
    print(f"    edge-sharing pairs = {n_edge_pairs}")
    print(f"    Var(c3) = {n_triples}*3/16 + 2*{n_edge_pairs}*3/16")
    print(f"            = {Var_c3_theory} = {float(Var_c3_theory):.6f}")
    print(f"    Actual:   {Var_c3} = {float(Var_c3):.6f}")
    print(f"    Match: {Var_c3 == Var_c3_theory}")

print("\n" + "=" * 72)
print("  KEY FORMULA: Var(c3) = 3*(C(n,3) + 2*C(n,2)*C(n-2,2)) / 16")
print("=" * 72)
print()

for n_test in range(3, 10):
    nt = comb(n_test, 3)
    ep = comb(n_test, 2) * comb(n_test - 2, 2) if n_test >= 4 else 0
    vc3 = Fraction(3*(nt + 2*ep), 16)
    # Simplify: C(n,3) + 2*C(n,2)*C(n-2,2)
    #  = n(n-1)(n-2)/6 + 2 * n(n-1)/2 * (n-2)(n-3)/2
    #  = n(n-1)(n-2)/6 + n(n-1)(n-2)(n-3)/2
    #  = n(n-1)(n-2)[1/6 + (n-3)/2]
    #  = n(n-1)(n-2)(3n-8)/12  ... let me check
    inner = Fraction(1, 6) + Fraction(n_test - 3, 2)
    combined = n_test * (n_test-1) * (n_test-2) * inner
    vc3_formula = Fraction(3, 16) * combined

    print(f"  n={n_test}: Var(c3) = {vc3}, 4*Var(c3) = {4*vc3}")
    print(f"    Inner formula: {combined}, 3/16 * inner = {vc3_formula}")
    print(f"    Match: {vc3 == vc3_formula}")

# Now: if Cov(H, c3) = 2*Var(c3) and R^2 = 4*Var(c3)/Var(H):
# Then 1 - R^2 = (Var(H) - 4*Var(c3)) / Var(H) = Var(eps)/Var(H)
# We need Var(H).

print("\n" + "=" * 72)
print("  COMPUTING Var(H) FOR RANDOM TOURNAMENTS")
print("=" * 72)
print()

for n_test in [3, 4, 5, 6]:
    m = n_test*(n_test-1)//2
    total = 1 << m

    sum_H = 0
    sum_H2 = 0
    for bits in range(total):
        A = build_tournament(n_test, bits)
        H = ham_paths_dp(A, n_test)
        sum_H += H
        sum_H2 += H*H

    E_H = Fraction(sum_H, total)
    Var_H = Fraction(sum_H2, total) - E_H**2

    # Theoretical 4*Var(c3)
    nt = comb(n_test, 3)
    ep = comb(n_test, 2) * comb(n_test-2, 2) if n_test >= 4 else 0
    Var4c3 = Fraction(3*(nt + 2*ep), 4)

    R2 = Var4c3 / Var_H if Var_H > 0 else Fraction(1)

    print(f"  n={n_test}: E[H] = {E_H}, Var(H) = {Var_H}")
    print(f"    4*Var(c3) = {Var4c3}")
    print(f"    R^2 = 4*Var(c3)/Var(H) = {R2}")
    print()
