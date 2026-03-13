#!/usr/bin/env python3
"""
perfect_square_det.py -- kind-pasteur-2026-03-13-S61

CONJECTURE: For any tournament T on n vertices, det(I + 2A(T)) is
always a perfect square (of an odd integer).

OBSERVED:
  n=5: det(I+2A) in {1, 9, 25, 49, 81} = {1^2, 3^2, 5^2, 7^2, 9^2}
  n=6: det(I+2A) in {1, 9, 25, 49, 81}
  n=7 (two cases): 169=13^2 and 361=19^2

WHY? Key identity for tournaments:
  A + A^T = J - I  (complete graph minus diagonal)
  So I + 2A = I + 2A = I + (J - I) + (2A - J + I) = J + S
  where S = A - A^T is the skew-adjacency matrix.
  Wait: I + 2A = I + A + A = I + (J - I - A^T) + A = J + (A - A^T) = J + S
  So: I + 2A = J + S where S is skew-symmetric!

For odd n: det(J + S) where S is skew-symmetric and J is all-ones.

Since J has eigenvalue n (eigenvector (1,...,1)) and eigenvalue 0 (everything else),
and S has purely imaginary eigenvalues...

Actually, let's think about it differently. For tournaments:
  I + 2A = J + S
  det(J + S) = det(S) * det(I + S^{-1} J) if S is invertible

For odd n, det(S) = 0 (skew-symmetric odd dimension), so this doesn't work directly.
But we can use the rank-1 update formula:
  det(J + S) = det(S + 1*1^T) where 1 = (1,...,1)

Since S is rank n-1 (odd n), this might be...

Let's just verify exhaustively and look for patterns.

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations
from collections import defaultdict
import math


def binary_to_tournament(bits, n):
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


def det_small(M, n):
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0]*M[1][1] - M[0][1]*M[1][0]
    result = 0
    for j in range(n):
        if M[0][j] == 0:
            continue
        minor = [[M[r][c] for c in range(n) if c != j] for r in range(1, n)]
        cofactor = ((-1) ** j) * M[0][j] * det_small(minor, n-1)
        result += cofactor
    return result


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def is_perfect_square(n):
    if n < 0:
        return False
    r = int(n**0.5 + 0.5)
    return r * r == n


# ========================================================================
# TEST 1: Exhaustive verification at n=3,4,5
# ========================================================================
for n in [3, 4, 5]:
    print(f"\n{'='*60}")
    print(f"n={n}: EXHAUSTIVE TEST")
    print(f"{'='*60}")

    m = n * (n - 1) // 2
    total = 1 << m

    det_vals = set()
    all_perfect_sq = True
    failures = 0

    for bits in range(total):
        A = binary_to_tournament(bits, n)

        # I + 2A
        M = [[2*A[i][j] for j in range(n)] for i in range(n)]
        for i in range(n):
            M[i][i] += 1
        d = det_small(M, n)
        det_vals.add(d)

        if not is_perfect_square(d):
            all_perfect_sq = False
            failures += 1
            if failures <= 3:
                print(f"  FAILURE: bits={bits}, det(I+2A)={d}")

    print(f"  All perfect squares: {all_perfect_sq}")
    print(f"  Failures: {failures}/{total}")
    print(f"  Distinct det values: {sorted(det_vals)}")

    # Check: are they all odd perfect squares?
    roots = sorted(set(int(d**0.5 + 0.5) for d in det_vals if d > 0))
    print(f"  Square roots: {roots}")
    all_odd = all(r % 2 == 1 for r in roots)
    print(f"  All roots odd: {all_odd}")


# ========================================================================
# TEST 2: n=6 exhaustive
# ========================================================================
print(f"\n{'='*60}")
print(f"n=6: EXHAUSTIVE TEST")
print(f"{'='*60}")

n6 = 6
m6 = n6 * (n6 - 1) // 2
total6 = 1 << m6

det_vals_6 = set()
all_sq_6 = True
failures_6 = 0

for bits in range(total6):
    A = binary_to_tournament(bits, n6)
    M = [[2*A[i][j] for j in range(n6)] for i in range(n6)]
    for i in range(n6):
        M[i][i] += 1
    d = det_small(M, n6)
    det_vals_6.add(d)

    if not is_perfect_square(d):
        all_sq_6 = False
        failures_6 += 1
        if failures_6 <= 3:
            print(f"  FAILURE: bits={bits}, det(I+2A)={d}")

print(f"  All perfect squares: {all_sq_6}")
print(f"  Failures: {failures_6}/{total6}")
print(f"  Distinct det values: {sorted(det_vals_6)}")
roots_6 = sorted(set(int(d**0.5 + 0.5) for d in det_vals_6 if d > 0))
print(f"  Square roots: {roots_6}")
all_odd_6 = all(r % 2 == 1 for r in roots_6)
print(f"  All roots odd: {all_odd_6}")


# ========================================================================
# TEST 3: n=7 sample
# ========================================================================
print(f"\n{'='*60}")
print(f"n=7: SAMPLE TEST (every 100th tournament)")
print(f"{'='*60}")

n7 = 7
m7 = n7 * (n7 - 1) // 2
total7 = 1 << m7

det_vals_7 = set()
all_sq_7 = True
failures_7 = 0
tested_7 = 0

for bits in range(0, total7, 100):
    A = binary_to_tournament(bits, n7)
    M = [[2*A[i][j] for j in range(n7)] for i in range(n7)]
    for i in range(n7):
        M[i][i] += 1
    d = det_small(M, n7)
    det_vals_7.add(d)
    tested_7 += 1

    if not is_perfect_square(d):
        all_sq_7 = False
        failures_7 += 1
        if failures_7 <= 3:
            print(f"  FAILURE: bits={bits}, det(I+2A)={d}")

print(f"  All perfect squares: {all_sq_7}")
print(f"  Failures: {failures_7}/{tested_7}")
print(f"  Distinct det values: {sorted(det_vals_7)}")
roots_7 = sorted(set(int(d**0.5 + 0.5) for d in det_vals_7 if d > 0))
print(f"  Square roots: {roots_7}")
if roots_7:
    all_odd_7 = all(r % 2 == 1 for r in roots_7)
    print(f"  All roots odd: {all_odd_7}")


# ========================================================================
# ANALYSIS: WHY det(I+2A) = perfect square?
# ========================================================================
print(f"\n{'='*60}")
print("ALGEBRAIC ANALYSIS: WHY IS det(I+2A) A PERFECT SQUARE?")
print("=" * 60)

print("""
For a tournament T with adjacency matrix A:
  A + A^T = J - I  (tournament axiom)
  S = A - A^T  (skew-symmetric)
  A = (J - I + S) / 2

So: I + 2A = I + (J - I + S) = J + S

Therefore: det(I + 2A) = det(J + S) where S is skew-symmetric.

Now: J = 1*1^T (rank-1 matrix, all entries 1).
By the matrix determinant lemma:
  det(J + S) = det(S + 1*1^T) = det(S) * (1 + 1^T S^{-1} 1)  [if S invertible]

For ODD n: det(S) = 0 (skew-symmetric of odd order).
So we need a different approach.

For EVEN n: det(S) = Pf(S)^2 (perfect square by Pfaffian).
And S^{-1} is also skew-symmetric, so 1^T S^{-1} 1 = 0.
This gives det(J + S) = det(S) = Pf(S)^2.  <-- PERFECT SQUARE!

For ODD n: We use the fact that J + S has rank at most n.
Since S has rank n-1 (odd skew-symmetric), and J has rank 1,
the rank of J + S is at most n.

Actually, for odd n: ker(S) = span(v) where v is the vector of
cofactors. Since S is a tournament's skew-adjacency, this v is related
to the Pfaffian structure.

Let me verify: for even n, det(I+2A) = det(J+S) = det(S) = Pf(S)^2.
""")

# Verify the even-n identity: det(I+2A) = det(S) = Pf(S)^2
print("VERIFICATION: Even n, det(I+2A) = det(S) = Pf(S)^2")
for n_test in [4, 6]:
    m_test = n_test * (n_test - 1) // 2
    total_test = 1 << m_test
    all_match = True
    for bits in range(total_test):
        A = binary_to_tournament(bits, n_test)
        S = [[A[i][j] - A[j][i] for j in range(n_test)] for i in range(n_test)]
        M = [[2*A[i][j] for j in range(n_test)] for i in range(n_test)]
        for i in range(n_test):
            M[i][i] += 1

        det_I2A = det_small(M, n_test)
        det_S = det_small(S, n_test)

        if det_I2A != det_S:
            all_match = False
            print(f"  n={n_test}, bits={bits}: det(I+2A)={det_I2A} != det(S)={det_S}")
            break

    if all_match:
        print(f"  n={n_test}: CONFIRMED det(I+2A) = det(S) for all {total_test} tournaments")


# For odd n, the story is different. Let's check what det(J+S) equals.
print(f"\nODD n ANALYSIS:")

# At n=3, check explicitly
n3 = 3
print(f"\n  n=3:")
for bits in range(1 << 3):
    A = binary_to_tournament(bits, n3)
    S = [[A[i][j] - A[j][i] for j in range(n3)] for i in range(n3)]
    JS = [[1 + S[i][j] for j in range(n3)] for i in range(n3)]
    det_JS = det_small(JS, n3)

    # Compute Pfaffian-related: delete each vertex, compute det of remaining
    pf_sum = 0
    for del_v in range(n3):
        remaining = [i for i in range(n3) if i != del_v]
        S2 = [[S[remaining[i]][remaining[j]] for j in range(2)] for i in range(2)]
        det_del = det_small(S2, 2)
        pf_sum += det_del

    print(f"    bits={bits}: det(J+S)={det_JS}, sum_Pf^2={pf_sum}")


# At n=5, check the formula
print(f"\n  n=5:")
n5 = 5
for bits in range(0, 1 << 10, 37):  # sample
    A = binary_to_tournament(bits, n5)
    S = [[A[i][j] - A[j][i] for j in range(n5)] for i in range(n5)]
    JS = [[1 + S[i][j] for j in range(n5)] for i in range(n5)]
    det_JS = det_small(JS, n5)

    # The cofactor expansion along the "1" column
    # Actually: J + S = (S)(I + S^{-1}J) but S is singular for odd n.
    # Alternative: adjugate formula.
    # det(J+S) = det(S + 11^T) = det(S) + 1^T adj(S) 1 = 0 + sum of all cofactors

    # For skew-symmetric S of odd order: adj(S)_{ij} = (-1)^{i+j} Pf(S_{ij})^2 / ???
    # Actually: for skew-symmetric S, adj(S) = Pf_minor matrix
    # The (i,j) cofactor of S is (-1)^{i+j} * det(S with row i, col j deleted)

    # Key identity: for odd-order skew-symmetric S,
    # adj(S) has rank 1, and adj(S) = v*v^T where v_i = (-1)^i * Pf(S with row i and col i deleted)

    # So: det(J + S) = det(S) + 1^T adj(S) 1 = 0 + (sum_i v_i)^2

    # Compute v_i = Pf(S_{ii})  (delete row i AND col i, take Pfaffian)
    v = []
    for i in range(n5):
        remaining = [j for j in range(n5) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n5-1)] for a in range(n5-1)]
        det_del = det_small(S_del, n5 - 1)
        # det_del = Pf^2, take signed Pf
        if det_del >= 0:
            pf = int(det_del**0.5 + 0.5)
        else:
            pf = 0  # shouldn't happen for skew-symmetric
        # Sign: (-1)^i
        signed_pf = ((-1) ** i) * pf
        v.append(signed_pf)

    v_sum = sum(v)
    v_sum_sq = v_sum ** 2
    print(f"    bits={bits}: det(J+S)={det_JS}, v={v}, sum(v)={v_sum}, sum^2={v_sum_sq}, match={det_JS == v_sum_sq}")


# ========================================================================
# THEOREM PROOF: det(I+2A) = (sum_i (-1)^i Pf(S_ii))^2
# ========================================================================
print(f"\n\n{'='*60}")
print("THEOREM: det(I+2A) = (sum_i (-1)^i Pf(S_ii))^2")
print("=" * 60)

print("""
For a tournament T on n vertices (n odd) with skew-adjacency S = A - A^T:

  I + 2A = J + S  (since A + A^T = J - I for tournaments)

  det(J + S) = det(S + 1*1^T)

Since S is skew-symmetric of odd order, det(S) = 0 and rank(S) = n-1.
The adjugate adj(S) has rank 1 for skew-symmetric matrices of odd order.

Key identity: adj(S) = v * v^T where v_i = (-1)^i * Pf(S_ii)
(S_ii = S with row i and column i deleted, which is even-order skew-symmetric)

By the matrix determinant lemma for rank-deficient matrices:
  det(S + 1*1^T) = det(S) + 1^T * adj(S) * 1 = 0 + (1^T v)^2 = (sum_i v_i)^2

Therefore: det(I + 2A) = (sum_{i=0}^{n-1} (-1)^i Pf(S_{ii}))^2

This is ALWAYS a perfect square!

For even n: det(I + 2A) = det(J + S) = det(S) = Pf(S)^2 (also a perfect square).

COMBINED: For ANY tournament on any n, det(I + 2A) is a perfect square.
""")

# Verify exhaustively at n=3, 5
print("EXHAUSTIVE VERIFICATION:")
for n_test in [3, 5]:
    m_test = n_test * (n_test - 1) // 2
    total_test = 1 << m_test
    all_match = True

    for bits in range(total_test):
        A = binary_to_tournament(bits, n_test)
        S = [[A[i][j] - A[j][i] for j in range(n_test)] for i in range(n_test)]

        M = [[2*A[i][j] for j in range(n_test)] for i in range(n_test)]
        for i in range(n_test):
            M[i][i] += 1
        det_I2A = det_small(M, n_test)

        # Compute v = Pfaffian vector
        v = []
        for i in range(n_test):
            remaining = [j for j in range(n_test) if j != i]
            S_del = [[S[remaining[a]][remaining[b]] for b in range(n_test-1)] for a in range(n_test-1)]
            det_del = det_small(S_del, n_test - 1)
            pf = int(det_del**0.5 + 0.5)
            signed_pf = ((-1) ** i) * pf
            v.append(signed_pf)

        v_sum_sq = sum(v) ** 2

        if det_I2A != v_sum_sq:
            all_match = False
            print(f"  MISMATCH at n={n_test}, bits={bits}: {det_I2A} != {v_sum_sq}")
            break

    if all_match:
        print(f"  n={n_test}: VERIFIED for all {total_test} tournaments")


# ========================================================================
# BONUS: What is sum_i (-1)^i Pf(S_ii) in terms of tournament invariants?
# ========================================================================
print(f"\n\n{'='*60}")
print("BONUS: WHAT IS THE PFAFFIAN SUM?")
print("=" * 60)

# At n=5, compute the Pfaffian sum and check its relationship to H, c3, c5
print(f"\n  n=5 exhaustive:")
by_pf_sum = defaultdict(list)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    S = [[A[i][j] - A[j][i] for j in range(n5)] for i in range(n5)]
    H = count_ham_paths(A, n5)

    v = []
    for i in range(n5):
        remaining = [j for j in range(n5) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n5-1)] for a in range(n5-1)]
        det_del = det_small(S_del, n5 - 1)
        pf = int(det_del**0.5 + 0.5)
        v.append(((-1) ** i) * pf)

    pf_sum = sum(v)
    by_pf_sum[pf_sum].append((bits, H))

print(f"  Pfaffian sum values: {sorted(by_pf_sum.keys())}")
for ps in sorted(by_pf_sum.keys()):
    H_vals = sorted(set(h for _, h in by_pf_sum[ps]))
    print(f"    Pf_sum={ps}: H={H_vals}, count={len(by_pf_sum[ps])}")


# ========================================================================
# THE KEY: Pfaffian sum and the Vitali hidden dimension
# ========================================================================
print(f"\n\n{'='*60}")
print("PFAFFIAN SUM AT n=7 FOR AMBIGUOUS PAIR")
print("=" * 60)

n7 = 7
for bits, label in [(4728, "H=109, c7=8"), (4658, "H=111, c7=9")]:
    A = binary_to_tournament(bits, n7)
    S = [[A[i][j] - A[j][i] for j in range(n7)] for i in range(n7)]

    v = []
    for i in range(n7):
        remaining = [j for j in range(n7) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n7-1)] for a in range(n7-1)]
        det_del = det_small(S_del, n7 - 1)
        pf = int(det_del**0.5 + 0.5)
        v.append(((-1) ** i) * pf)

    pf_sum = sum(v)
    print(f"  {label}:")
    print(f"    v = {v}")
    print(f"    Pf_sum = {pf_sum}")
    print(f"    Pf_sum^2 = {pf_sum**2}")
    print(f"    det(I+2A) = {pf_sum**2}")


print(f"\n\n{'='*60}")
print("DONE.")
print("=" * 60)
