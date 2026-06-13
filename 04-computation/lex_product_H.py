"""
lex_product_H.py -- kind-pasteur-2026-03-14-S81
Deep exploration of the lexicographic product formula for H.

FROM S81: H(T1 lex T2) = H(T1) * H(T2)^|T1| for 2 lex 3.
But NOT for 3 lex 2. WHY?

The lex product T1 lex T2 has vertices V1 x V2 with arcs:
  (i1,j1) -> (i2,j2) iff T1[i1][i2]=1 OR (i1=i2 AND T2[j1][j2]=1)

A Hamiltonian path in T1 lex T2 must visit all |V1|*|V2| vertices.
The key structure: within each "copy" of T2 (for fixed i in V1),
the path must traverse all vertices of T2 before moving to the next i.

THEOREM: H(T1 lex T2) = H(T1) * H(T2)^{|V1|} when |V2| >= 2.

PROOF SKETCH: Each Ham path in T1 lex T2 consists of:
1. A Ham path in T1: (i_0, i_1, ..., i_{a-1})
2. For each i_k: a Ham path in T2 connecting the entry point to exit point

If T2 has a UNIQUE path structure (only 1 Ham path), then each copy
of T2 has exactly 1 way to traverse, giving H(T1 lex T2) = H(T1).
For general T2, each of the a copies of T2 contributes H(T2) paths
(if the entry/exit points are free). But they're NOT free — the
exit of copy k connects to the entry of copy k+1 via a T1 arc.

Actually for lex product: (i_k, j_last) -> (i_{k+1}, j_first) requires
T1[i_k][i_{k+1}]=1. The j values DON'T matter for inter-copy transitions.
So within each copy i_k, we need a Ham path in T2. The choice for each copy
is independent!

Therefore: H(T1 lex T2) = H(T1) * H(T2)^{|V1|}.

Wait, but S81 showed this fails for 3 lex 2. Let me verify more carefully.
"""

import numpy as np
from itertools import permutations
from collections import Counter
import sys

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
    """Lexicographic product T1 lex T2."""
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
    print("LEXICOGRAPHIC PRODUCT AND H — DEEP EXPLORATION")
    print("kind-pasteur-2026-03-14-S81")
    print("=" * 70)

    # Test all small products exhaustively
    for n1, n2 in [(2, 2), (2, 3), (3, 2), (2, 4), (3, 3)]:
        n_prod = n1 * n2
        if n_prod > 8:
            print(f"\n  {n1} lex {n2}: n_prod={n_prod}, too large, skipping")
            continue

        print(f"\n  === {n1} lex {n2} (product size {n_prod}) ===")

        m1 = n1*(n1-1)//2
        m2 = n2*(n2-1)//2

        formula_works = True
        formula_fails = []

        for b1 in range(2**m1):
            A1 = bits_to_adj(b1, n1)
            H1 = compute_H(A1.tolist(), n1)

            for b2 in range(2**m2):
                A2 = bits_to_adj(b2, n2)
                H2 = compute_H(A2.tolist(), n2)

                A_prod, n_p = lex_product(A1, n1, A2, n2)
                H_prod = compute_H(A_prod.tolist(), n_p)

                predicted = H1 * (H2 ** n1)

                if H_prod != predicted:
                    formula_works = False
                    formula_fails.append((H1, H2, H_prod, predicted))

        if formula_works:
            print(f"    H(T1 lex T2) = H(T1) * H(T2)^{n1} : VERIFIED for ALL {2**m1 * 2**m2} pairs!")
        else:
            print(f"    FORMULA FAILS: {len(formula_fails)} / {2**m1 * 2**m2} counterexamples")
            for h1, h2, hp, pred in formula_fails[:5]:
                print(f"      H1={h1}, H2={h2}, H_prod={hp}, predicted={pred}, ratio={hp/pred if pred>0 else 'div0'}")

            # Try alternative formulas
            # H_prod = H1 * H2^n1 doesn't work. Try H1^n2 * H2^n1?
            alt_works = all(hp == h1**n2 * h2**n1 for h1, h2, hp, _ in formula_fails)
            if alt_works:
                print(f"    ALTERNATIVE: H(T1 lex T2) = H(T1)^{n2} * H(T2)^{n1} WORKS!")

            # Try H_prod / (H1 * H2)
            ratios = set(hp / (h1 * h2) if h1*h2 > 0 else None
                         for h1, h2, hp, _ in formula_fails if h1*h2 > 0)
            print(f"    H_prod / (H1 * H2) values: {sorted(ratios)[:10]}")

    # ============================================================
    # DIRECT PRODUCT (tensor product)
    # ============================================================
    print(f"\n{'='*70}")
    print("DIRECT (TENSOR) PRODUCT: T1 x T2")
    print("  (i1,j1) -> (i2,j2) iff T1[i1][i2]=1 AND T2[j1][j2]=1")
    print(f"{'='*70}")

    for n1, n2 in [(2, 2), (2, 3), (3, 2)]:
        n_prod = n1 * n2
        print(f"\n  === {n1} x {n2} (product size {n_prod}) ===")

        m1 = n1*(n1-1)//2
        m2 = n2*(n2-1)//2

        for b1 in range(min(2**m1, 4)):
            A1 = bits_to_adj(b1, n1)
            H1 = compute_H(A1.tolist(), n1)

            for b2 in range(min(2**m2, 4)):
                A2 = bits_to_adj(b2, n2)
                H2 = compute_H(A2.tolist(), n2)

                # Tensor product
                A_prod = np.zeros((n_prod, n_prod), dtype=int)
                for i1 in range(n1):
                    for j1 in range(n2):
                        for i2 in range(n1):
                            for j2 in range(n2):
                                v1 = i1 * n2 + j1
                                v2 = i2 * n2 + j2
                                if v1 != v2 and A1[i1][i2] == 1 and A2[j1][j2] == 1:
                                    A_prod[v1][v2] = 1

                # This is NOT a tournament (may have missing arcs)
                # Check completeness
                is_tourn = True
                for i in range(n_prod):
                    for j in range(i+1, n_prod):
                        if A_prod[i][j] + A_prod[j][i] != 1:
                            is_tourn = False
                            break
                    if not is_tourn: break

                if is_tourn:
                    H_prod = compute_H(A_prod.tolist(), n_prod)
                    print(f"    H1={H1}, H2={H2}: TOURNAMENT, H_prod={H_prod}")
                else:
                    print(f"    H1={H1}, H2={H2}: NOT a tournament (tensor product)")

    # ============================================================
    # SUBSTITUTION PRODUCT
    # ============================================================
    print(f"\n{'='*70}")
    print("SUBSTITUTION PRODUCT: T1[v <- T2]")
    print("  Replace vertex v in T1 with a copy of T2")
    print("  Arcs to/from v in T1 become arcs to/from ALL of T2")
    print(f"{'='*70}")

    # T1 on {0,...,n1-1}, replace vertex 0 with T2 on {0,...,n2-1}
    for n1, n2 in [(3, 2), (3, 3)]:
        n_prod = n1 - 1 + n2
        if n_prod > 7:
            print(f"\n  {n1}[0 <- {n2}]: n_prod={n_prod}, too large, skipping")
            continue

        print(f"\n  === {n1}[0 <- {n2}] (product size {n_prod}) ===")
        m1 = n1*(n1-1)//2
        m2 = n2*(n2-1)//2

        results = []
        for b1 in range(2**m1):
            A1 = bits_to_adj(b1, n1)
            H1 = compute_H(A1.tolist(), n1)

            for b2 in range(2**m2):
                A2 = bits_to_adj(b2, n2)
                H2 = compute_H(A2.tolist(), n2)

                # Build substitution: replace vertex 0 with T2
                # New vertices: T2 vertices (0..n2-1) + T1 vertices 1..n1-1
                A_sub = np.zeros((n_prod, n_prod), dtype=int)

                # T2 arcs among 0..n2-1
                for i in range(n2):
                    for j in range(n2):
                        if i != j:
                            A_sub[i][j] = A2[i][j]

                # T1 arcs among n2..n_prod-1 (= T1 vertices 1..n1-1)
                for i in range(1, n1):
                    for j in range(1, n1):
                        if i != j:
                            A_sub[n2 + i - 1][n2 + j - 1] = A1[i][j]

                # Cross arcs: T2 vertex k to T1 vertex i (i >= 1)
                # Direction = same as vertex 0 in T1 to vertex i
                for k in range(n2):
                    for i in range(1, n1):
                        if A1[0][i] == 1:
                            A_sub[k][n2 + i - 1] = 1
                        else:
                            A_sub[n2 + i - 1][k] = 1

                H_sub = compute_H(A_sub.tolist(), n_prod)
                results.append((H1, H2, H_sub))

        # Analyze
        print(f"    {len(results)} substitutions computed")

        # Does H_sub = H1 * H2? Or H1 * f(H2)?
        for h1, h2, hs in results[:10]:
            ratio = hs / (h1 * h2) if h1 * h2 > 0 else 'N/A'
            print(f"    H1={h1}, H2={h2}: H_sub={hs}, H_sub/(H1*H2)={ratio}")

        # Check if the ratio is constant
        ratios = set()
        for h1, h2, hs in results:
            if h1 > 0 and h2 > 0:
                ratios.add(round(hs / (h1 * h2), 6))

        print(f"    Distinct H_sub/(H1*H2) ratios: {sorted(ratios)[:10]}")
        if len(ratios) == 1:
            print(f"    *** FORMULA: H_sub = {list(ratios)[0]} * H1 * H2 ***")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
