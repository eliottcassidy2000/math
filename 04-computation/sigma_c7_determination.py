#!/usr/bin/env python3
"""
Check: at n=7, does (lambda, sigma) on the full vertex set uniquely determine the tournament?

If so, c7 is trivially (lambda,sigma)-determined.
If not, there might be tournaments with same (lambda,sigma) but different c7.

Also: what is the relationship between sigma(i,j) and the adjacency?
sigma(i,j) = co(i,j) + ci(i,j)
           = #{common out-neighbors} + #{common in-neighbors}
           = (n-2) - P_{ij} - P_{ji}  where P = A²

So sigma(i,j) = (n-2) - (A²)_{ij} - (A²)_{ji}

And lambda(i,j) = A_{ij}·(A²)_{ji} + (1-A_{ij})·(A²)_{ij}

From sigma: (A²)_{ij} + (A²)_{ji} = n-2-sigma(i,j)
From lambda: A_{ij}·(A²)_{ji} + (1-A_{ij})·(A²)_{ij} = lambda(i,j)

Let P = (A²)_{ij}, Q = (A²)_{ji}. Then P+Q = n-2-sigma.
If i→j (A_{ij}=1): lambda = Q.
  So Q = lambda, P = n-2-sigma-lambda.
If j→i (A_{ij}=0): lambda = P.
  So P = lambda, Q = n-2-sigma-lambda.

This means: given (lambda, sigma), and knowing the DIRECTION A_{ij},
we can recover (A²)_{ij} and (A²)_{ji} completely!

But we DON'T know A_{ij} from just (lambda, sigma).
The question: does (lambda, sigma) + score sequence determine A_{ij}?

Actually: P_{ij} = d_i - A_{ij} - co(i,j) (from earlier derivation)
And co(i,j) = (sigma + d_i + d_j - n + 1) / 2

So: P = d_i - A_{ij} - (sigma + d_i + d_j - n + 1)/2
     = d_i/2 - A_{ij} - (sigma + d_j - n + 1)/2

Also: if A_{ij}=1: P = n-2-sigma-lambda.
  → d_i - 1 - (sigma + d_i + d_j - n + 1)/2 = n-2-sigma-lambda
  → d_i - 1 - sigma/2 - d_i/2 - d_j/2 + (n-1)/2 = n-2-sigma-lambda
  → d_i/2 - d_j/2 + (n-3)/2 - sigma/2 = n-2-sigma-lambda
  → d_i - d_j + n - 3 - sigma = 2n - 4 - 2sigma - 2lambda + ... hmm

Let me just verify computationally.

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
from collections import defaultdict
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def lambda_sigma(A, n):
    L = np.zeros((n, n), dtype=int)
    S = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
                if A[u][w] and A[v][w]: S[u][v] += 1; S[v][u] += 1
                if A[w][u] and A[w][v]: S[u][v] += 1; S[v][u] += 1
    return L, S

# At n=7: does (lambda, sigma) uniquely determine the tournament?
n = 7
tb = n*(n-1)//2
np.random.seed(42)

print(f"n={n}: does (lambda, sigma) determine the tournament?")
ls_to_adj = {}
non_unique_count = 0
c7_diff_count = 0

for trial in range(200000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L, S = lambda_sigma(A, n)
    key = (tuple(L[i][j] for i in range(n) for j in range(i+1, n)),
           tuple(S[i][j] for i in range(n) for j in range(i+1, n)))

    if key in ls_to_adj:
        A_prev = ls_to_adj[key]
        if not np.array_equal(A, A_prev):
            non_unique_count += 1
            # Just compare tr7
            tr7_1 = int(np.trace(np.linalg.matrix_power(A, 7)))
            tr7_2 = int(np.trace(np.linalg.matrix_power(A_prev, 7)))
            if tr7_1 != tr7_2:
                c7_diff_count += 1
            if non_unique_count <= 3:
                print(f"  Non-unique pair #{non_unique_count}: tr7={tr7_1} vs {tr7_2} (same={tr7_1==tr7_2})")
    else:
        ls_to_adj[key] = A.copy()

print(f"\n  (lambda,sigma) groups: {len(ls_to_adj)}")
print(f"  Non-unique pairs found: {non_unique_count}")
print(f"  c7 differences: {c7_diff_count}")

if non_unique_count == 0:
    print("  (lambda, sigma) UNIQUELY determines the tournament at n=7!")
    print("  → All cycle counts, and ALL invariants, are (lambda,sigma)-determined.")
else:
    print(f"  (lambda, sigma) does NOT uniquely determine tournaments at n=7.")
    if c7_diff_count == 0:
        print("  But c7 (tr7) is still the same for all same-(lambda,sigma) tournaments!")

# Does the DIRECTION A_{ij} follow from (lambda, sigma, scores)?
# P = (A²)_{ij}, Q = (A²)_{ji}
# P + Q = n-2-sigma(i,j)
# If A_{ij}=1: lambda = Q = (A²)_{ji}. P = n-2-sigma-lambda.
# If A_{ij}=0: lambda = P = (A²)_{ij}. Q = n-2-sigma-lambda.
# So: |P - Q| = |n-2-sigma-2*lambda| (since one of P,Q equals lambda and the other equals n-2-sigma-lambda)
# If lambda ≠ (n-2-sigma)/2: P ≠ Q, so we can determine which is larger.
# A_{ij}=1 ↔ Q = lambda ↔ Q < P (if lambda < n-2-sigma-lambda, i.e., 2*lambda < n-2-sigma)
# A_{ij}=1 ↔ Q = lambda ↔ Q > P (if lambda > n-2-sigma-lambda, i.e., 2*lambda > n-2-sigma)
# Wait, we can't determine P vs Q from just their values without knowing which is (A²)_{ij}.

# But we CAN compute: P_{ij} = d_i - A_{ij} - co(i,j)
# And co(i,j) = (sigma + d_i + d_j - n + 1)/2

# If A_{ij}=1: P = d_i - 1 - co = d_i - 1 - (sigma+d_i+d_j-n+1)/2
# If A_{ij}=0: P = d_i - 0 - co = d_i - (sigma+d_i+d_j-n+1)/2

# So given (lambda, sigma) and the scores, A_{ij}=1 iff P = n-2-sigma-lambda
# i.e., d_i - 1 - (sigma+d_i+d_j-n+1)/2 = n-2-sigma-lambda
# d_i - 1 - sigma/2 - d_i/2 - d_j/2 + n/2 - 1/2 = n - 2 - sigma - lambda
# d_i/2 - d_j/2 - 3/2 - sigma/2 + n/2 = n - 2 - sigma - lambda
# d_i/2 - d_j/2 + (n-3)/2 - sigma/2 = n - 2 - sigma - lambda
# d_i - d_j + n - 3 - sigma = 2n - 4 - 2sigma - 2lambda
# d_i - d_j = 2n - 4 - 2sigma - 2lambda - n + 3 + sigma
# d_i - d_j = n - 1 - sigma - 2lambda

# So: A_{ij} = 1 iff d_i - d_j = n - 1 - sigma(i,j) - 2*lambda(i,j)

# Let's verify this!
print(f"\n{'='*60}")
print("ALGEBRAIC FORMULA: A_ij = 1 iff d_i - d_j = n-1-sigma-2*lambda")
print(f"{'='*60}")

verified = 0
failed = 0
for trial in range(10000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L, S = lambda_sigma(A, n)
    scores = [int(sum(A[i])) for i in range(n)]

    for i in range(n):
        for j in range(i+1, n):
            lam = L[i][j]
            sig = S[i][j]
            di, dj = scores[i], scores[j]
            predicted_diff = n - 1 - sig - 2*lam

            actual_aij = A[i][j]
            actual_diff = di - dj if actual_aij == 1 else dj - di

            # If A_{ij}=1: d_i - d_j should = predicted_diff
            # If A_{ij}=0 (i.e., A_{ji}=1): d_j - d_i should = predicted_diff (by symmetry)

            if actual_aij == 1:
                if di - dj != predicted_diff:
                    failed += 1
                    if failed <= 3:
                        print(f"  FAIL: i={i},j={j}, A_ij={actual_aij}, di-dj={di-dj}, predicted={predicted_diff}")
                else:
                    verified += 1
            else:
                # A_{ji}=1, so d_j - d_i should = n-1-sig-2*lam (swap i↔j)
                if dj - di != predicted_diff:
                    failed += 1
                    if failed <= 3:
                        print(f"  FAIL: i={i},j={j}, A_ij={actual_aij}, dj-di={dj-di}, predicted={predicted_diff}")
                else:
                    verified += 1

print(f"  Verified: {verified}, Failed: {failed}")

if failed > 0:
    print("  Formula FAILS — need to fix the derivation.")
    # Maybe the formula is: A_{ij}=1 iff d_i - d_j = something else.
    # Let me check the correct relationship empirically.
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L, S = lambda_sigma(A, n)
    scores = [int(sum(A[i])) for i in range(n)]
    print(f"\n  Example tournament:")
    for i in range(n):
        for j in range(i+1, n):
            lam = L[i][j]; sig = S[i][j]
            di, dj = scores[i], scores[j]
            aij = A[i][j]
            P = int((np.linalg.matrix_power(A, 2))[i][j])
            Q = int((np.linalg.matrix_power(A, 2))[j][i])
            print(f"    ({i},{j}): A={aij}, d=({di},{dj}), λ={lam}, σ={sig}, P={P}, Q={Q}, d_i-d_j={di-dj}, n-1-σ-2λ={n-1-sig-2*lam}")
else:
    print("  FORMULA VERIFIED! A_{ij}=1 iff d_i - d_j = n-1 - sigma(i,j) - 2*lambda(i,j)")
    print("  This means (lambda, sigma, scores) UNIQUELY determine the tournament!")
    print("  And since scores are determined by lambda (sum of lambda values per vertex),")
    print("  (lambda, sigma) might determine scores too!")

    # Check: do lambda values determine scores?
    print(f"\n  Do lambda values determine scores?")
    # d_i = sum_{j≠i} A_{ij} = (n-1)/2 + (something involving B_{ij})
    # In a tournament: sum_{j} lambda(i,j) = sum_{j} #{3-cycles through (i,j)}
    # = 2 * #{3-cycles containing vertex i} (each cycle contributes to 2 pairs involving i)
    # So per-vertex lambda sum = 2 * t3(i) where t3(i) = #{3-cycles containing i}
    # And t3(i) = d_i * (n-1-d_i) - #{non-edges from i's out-neighbors to i's in-neighbors}
    # Actually: t3(i) = d_i * (d_i - 1)/2 - ? No...
    # t3(i) = #{(j,k): i→j→k→i} = #{k in in-neighbors of i that are out-neighbors of some out-neighbor of i}
    # = (A²)_{ii} = 0 (diagonal of A² is 0 for tournaments)... no, that's wrong.
    # Hmm. Let me compute sum of lambda(i,j) over j and compare to d_i.
    for i in range(n):
        lam_sum = sum(L[i][j] for j in range(n) if j != i)
        print(f"    vertex {i}: d_i={scores[i]}, sum_j lambda(i,j)={lam_sum}")

    # Check: does d_i follow from lambda?
    lam_to_scores = defaultdict(set)
    for trial in range(20000):
        bits = np.random.randint(0, 1 << tb)
        A = bits_to_adj(bits, n)
        L, _ = lambda_sigma(A, n)
        key = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
        scores = tuple(int(sum(A[i])) for i in range(n))
        lam_to_scores[key].add(scores)

    ambig = sum(1 for v in lam_to_scores.values() if len(v) > 1)
    print(f"  Lambda determines scores? Ambiguous: {ambig}/{len(lam_to_scores)}")

print(f"\nDone.")
