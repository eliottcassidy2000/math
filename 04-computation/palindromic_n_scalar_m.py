#!/usr/bin/env python3
"""
PALINDROMIC N AND SCALAR M: A clean characterization.

THEOREM (at odd n):
  N(a,b,j) palindromic for all (a,b) ⟹ M = (H/n)*I

PROOF: At odd n, if N(a,b,j) = N(a,b,n-2-j) for all j, then:
  sum_j (-1)^j N(a,b,j) = sum_j (-1)^j N(a,b,n-2-j)
  = sum_k (-1)^{n-2-k} N(a,b,k) = (-1)^{n} sum_k (-1)^k N(a,b,k)
  (since (-1)^{-2} = 1)
  At odd n: (-1)^n = -1, so alt_sum = -alt_sum ⟹ alt_sum = 0.
  Hence M[a,b] = 0 for all a ≠ b.

QUESTION: Does the converse hold? At n=5:
- Both scalar M classes have palindromic N (verified).
- Do ALL palindromic-N tournaments have scalar M?
- Is palindromicity EQUIVALENT to self-complementarity (T ≅ T^op)?

Also: what happens at EVEN n? Palindromic N doesn't force alt_sum = 0.
"""

from itertools import permutations
from collections import defaultdict
import numpy as np

def ham_paths(A):
    n = len(A)
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if A[perm[k]][perm[k+1]] != 1: ok = False; break
        if ok: paths.append(perm)
    return paths

def compute_N(A):
    n = len(A)
    N = [[[0]*(n-1) for _ in range(n)] for _ in range(n)]
    for p in ham_paths(A):
        for j in range(n-1):
            a, b = p[j], p[j+1]
            N[a][b][j] += 1
            N[b][a][j] += 1
    return N

def is_palindromic_N(N, n):
    for a in range(n):
        for b in range(a+1, n):
            for j in range((n-1)//2 + 1):
                if N[a][b][j] != N[a][b][n-2-j]:
                    return False
    return True

def transfer_matrix_via_consec(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for p in ham_paths(A):
        for a in range(n):
            M[a][a] += (-1)**(list(p).index(a))
        for j in range(n-1):
            a, b = p[j], p[j+1]
            M[a][b] += (-1)**j
            M[b][a] += (-1)**j
    return M

def all_tournaments_canonical(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set()
    results = []
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        key = tuple(tuple(row) for row in A)
        min_key = key
        for perm in permutations(range(n)):
            pkey = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
            if pkey < min_key: min_key = pkey
        if min_key in seen: continue
        seen.add(min_key)
        results.append(A)
    return results

def opposite_tournament(A):
    n = len(A)
    return [[1 - A[i][j] if i != j else 0 for j in range(n)] for i in range(n)]

def is_isomorphic(A, B):
    n = len(A)
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[j]] == B[i][j] for i in range(n) for j in range(n)):
            return True
    return False

# =====================================================================
# n=5: Palindromic N vs Scalar M vs Self-complementary
# =====================================================================
print("=" * 70)
print("n=5: PALINDROMIC N vs SCALAR M vs SELF-COMPLEMENTARY")
print("=" * 70)

n = 5
classes = all_tournaments_canonical(n)
print(f"\n{len(classes)} isomorphism classes")

print(f"\n  {'Class':>5} {'H':>3} {'Scores':>15} {'Palindromic':>12} {'Scalar M':>10} {'Self-comp':>10}")
print("  " + "-" * 60)

for idx, A in enumerate(classes):
    H = len(ham_paths(A))
    scores = sorted([sum(A[i]) for i in range(n)], reverse=True)
    N = compute_N(A)
    M = transfer_matrix_via_consec(A)

    palindromic = is_palindromic_N(N, n)
    scalar = np.array_equal(M, (H//n)*np.eye(n, dtype=int)) if H % n == 0 else False
    A_op = opposite_tournament(A)
    self_comp = is_isomorphic(A, A_op)

    print(f"  {idx+1:5d} {H:3d} {str(scores):>15} {str(palindromic):>12} "
          f"{str(scalar):>10} {str(self_comp):>10}")

# =====================================================================
# ANALYSIS: Palindromic ⟹ Scalar? Self-comp ⟹ Palindromic?
# =====================================================================
print()
print("=" * 70)
print("LOGICAL RELATIONSHIPS")
print("=" * 70)

pal_implies_scalar = True
scalar_implies_pal = True
selfcomp_implies_pal = True
pal_implies_selfcomp = True

for A in classes:
    H = len(ham_paths(A))
    N = compute_N(A)
    M = transfer_matrix_via_consec(A)

    palindromic = is_palindromic_N(N, n)
    scalar = np.array_equal(M, (H//n)*np.eye(n, dtype=int)) if H % n == 0 else False
    A_op = opposite_tournament(A)
    self_comp = is_isomorphic(A, A_op)

    if palindromic and not scalar: pal_implies_scalar = False
    if scalar and not palindromic: scalar_implies_pal = False
    if self_comp and not palindromic: selfcomp_implies_pal = False
    if palindromic and not self_comp: pal_implies_selfcomp = False

print(f"  Palindromic N ⟹ Scalar M: {pal_implies_scalar}")
print(f"  Scalar M ⟹ Palindromic N: {scalar_implies_pal}")
print(f"  Self-complementary ⟹ Palindromic N: {selfcomp_implies_pal}")
print(f"  Palindromic N ⟹ Self-complementary: {pal_implies_selfcomp}")

# =====================================================================
# PATH REVERSAL: N_T(a,b,j) = N_{T^op}(a,b,n-2-j)?
# =====================================================================
print()
print("=" * 70)
print("PATH REVERSAL IDENTITY: N_T(a,b,j) = N_{T^op}(a,b,n-2-j)?")
print("=" * 70)

for A in classes[:3]:
    H = len(ham_paths(A))
    A_op = opposite_tournament(A)
    N_T = compute_N(A)
    N_Top = compute_N(A_op)

    reversal_holds = True
    for a in range(n):
        for b in range(a+1, n):
            for j in range(n-1):
                if N_T[a][b][j] != N_Top[a][b][n-2-j]:
                    reversal_holds = False
                    break

    print(f"  H={H}: N_T(a,b,j) = N_Top(a,b,n-2-j)? {reversal_holds}")

# =====================================================================
# M_T and M_{T^op} relationship
# =====================================================================
print()
print("=" * 70)
print("M_T vs M_{T^op}")
print("=" * 70)

for A in classes[:6]:
    H = len(ham_paths(A))
    A_op = opposite_tournament(A)
    M_T = transfer_matrix_via_consec(A)
    M_Top = transfer_matrix_via_consec(A_op)

    # Check various relationships
    # M_Top = M_T? No (different tournament).
    # M_Top = -M_T? Unlikely.
    # M_Top = M_T with sign changes?

    # From the corrected formula:
    # M_T[a,b] = sum_j (-1)^j N_T(a,b,j)
    # M_Top[a,b] = sum_j (-1)^j N_Top(a,b,j)
    #            = sum_j (-1)^j N_T(a,b,n-2-j)  (by reversal)
    #            = sum_k (-1)^{n-2-k} N_T(a,b,k) = (-1)^n sum_k (-1)^k N_T(a,b,k)
    #            = (-1)^n M_T[a,b]   (for off-diagonal)
    #
    # Hmm wait — this uses the reversal identity which maps (a,b) pairs
    # in T to (a,b) pairs in T^op. But the reversal swaps a,b as well!
    # Actually no: reversing path has {a,b} at {j,j+1} becoming
    # {b,a} at {n-2-j, n-1-j}. So N_T(a,b,j) counts paths where
    # {a,b} is at positions {j,j+1}. Reversal gives {a,b} at
    # {n-2-j, n-1-j} in T^op. This IS {a,b} at positions
    # {n-2-j, n-1-j} = {n-2-j, (n-2-j)+1}. So:
    # N_{T^op}(a,b,n-2-j) = N_T(a,b,j). ✓

    expected = (-1)**n * M_T  # Off-diagonal prediction

    diag_match = all(M_Top[a][a] == M_T[a][a] for a in range(n))
    # Actually diagonal: M[a,a] = sum_j (-1)^j P[a,j]
    # M_Top[a,a] = sum_j (-1)^j P_Top[a,j]
    # Reversal: P_T(a,j) = P_Top(a,n-1-j)
    # M_Top[a,a] = sum_j (-1)^j P_T(a,n-1-j) = sum_k (-1)^{n-1-k} P_T(a,k)
    #            = (-1)^{n-1} M_T[a,a]
    # At odd n: M_Top[a,a] = M_T[a,a]. At even n: M_Top[a,a] = -M_T[a,a].

    offdiag_match = True
    for a in range(n):
        for b in range(n):
            if a != b:
                exp = (-1)**n * M_T[a][b]
                if M_Top[a][b] != exp:
                    offdiag_match = False
                    break

    diag_pred = (-1)**(n-1)
    diag_match2 = all(M_Top[a][a] == diag_pred * M_T[a][a] for a in range(n))

    scores = sorted([sum(A[i]) for i in range(n)], reverse=True)
    print(f"  H={H}, scores={scores}:")
    print(f"    M_Top off-diag = (-1)^n * M_T off-diag? {offdiag_match}")
    print(f"    M_Top diag = (-1)^(n-1) * M_T diag? {diag_match2}")
    if offdiag_match and diag_match2:
        # At odd n: (-1)^n = -1 off-diag, (-1)^{n-1} = 1 diag
        # So M_Top = M_T - 2*diag(M_T) ... no
        # M_Top[a,a] = M_T[a,a], M_Top[a,b] = -M_T[a,b] for a≠b
        # M_Top = 2*diag(M_T) - M_T
        pass

print()
print("=" * 70)
print("THEOREM: M_{T^op} AND M_T RELATIONSHIP")
print("=" * 70)
print("""
From path reversal and the corrected consecutive-position formula:

  M_{T^op}[a,a] = (-1)^{n-1} * M_T[a,a]    (diagonal)
  M_{T^op}[a,b] = (-1)^n * M_T[a,b]         (off-diagonal)

At odd n (both (-1)^{n-1} = 1 and (-1)^n = -1):
  M_{T^op} = diag(M_T) - offdiag(M_T)

At even n ((-1)^{n-1} = -1 and (-1)^n = 1):
  M_{T^op} = -diag(M_T) + offdiag(M_T) = M_T - 2*diag(M_T)

COROLLARY (odd n): If T is self-complementary (T ≅ T^op), then
  M = diag(M) - offdiag(M)  (up to relabeling)
  ⟹ offdiag(M) = 0 ⟹ M is diagonal ⟹ M = diag(M_T[a,a])

But M is SYMMETRIC, so if M is diagonal then M[a,a] = ... but this
doesn't immediately give M = (H/n)*I unless diagonal entries are equal.

Actually, self-complementary + the isomorphism gives:
  For each (a,b): M[sigma(a),sigma(b)] = -M[a,b] where sigma is the
  isomorphism T -> T^op. This forces M[a,b] = 0 for all a≠b in the
  "sigma-symmetric" subspace.
""")
