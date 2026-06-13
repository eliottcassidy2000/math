"""
alpha1_parity_debug.py -- kind-pasteur-2026-03-13-S62

Debug the (H-1)/2 mod 2 = alpha_1 mod 2 identity at n=7.
At n=5 it's perfect (1024/1024), but at n=7 it shows 535/1000.

Since H = 1 + 2*alpha_1 + 4*alpha_2 + ...,
(H-1)/2 = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...,
so (H-1)/2 mod 2 = alpha_1 mod 2 MUST hold.

The failure means either H != I(Omega,2) at n=7, or my alpha_1
computation is wrong. Let me check both.
"""

import numpy as np
from itertools import combinations
from collections import Counter

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

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# ============================================================
# TEST: Check if np.linalg.matrix_power preserves integers
# ============================================================
print("=" * 70)
print("TEST: Integer preservation in matrix_power")
print("=" * 70)

# Create a small test matrix
A_test = np.array([[0, 1, 1, 0, 0],
                    [0, 0, 1, 1, 0],
                    [0, 0, 0, 1, 1],
                    [1, 0, 0, 0, 1],
                    [1, 1, 0, 0, 0]], dtype=int)

print(f"  A dtype: {A_test.dtype}")
A5 = np.linalg.matrix_power(A_test, 5)
print(f"  A^5 dtype: {A5.dtype}")
print(f"  A^5 trace: {np.trace(A5)}, type: {type(np.trace(A5))}")
print(f"  A^5 trace / 5: {np.trace(A5) / 5}")
print(f"  A^5 trace // 5: {np.trace(A5) // 5}")

# Check for n=7 with a larger matrix
A7 = np.random.randint(0, 2, (7, 7), dtype=np.int64)
np.fill_diagonal(A7, 0)
A7_pow = np.linalg.matrix_power(A7, 7)
print(f"\n  Random 7x7, A^7 dtype: {A7_pow.dtype}")
print(f"  Max entry in A^7: {np.max(A7_pow)}")

# The issue: np.linalg.matrix_power might return FLOAT for large powers
# Let me use my own int matrix power
def int_matrix_power(A, k):
    """Integer matrix power, guaranteed to stay in int."""
    n = A.shape[0]
    result = np.eye(n, dtype=np.int64)
    base = A.astype(np.int64)
    while k > 0:
        if k % 2 == 1:
            result = result @ base
        base = base @ base
        k //= 2
    return result

A7_int = int_matrix_power(A7.astype(np.int64), 7)
print(f"  Int power A^7 dtype: {A7_int.dtype}")
print(f"  Match with np version: {np.allclose(A7_pow, A7_int)}")
print(f"  Exact match: {np.array_equal(A7_pow.astype(np.int64), A7_int)}")

# ============================================================
# Fix: use int_matrix_power and recheck alpha_1 parity
# ============================================================
print("\n" + "=" * 70)
print("FIXED: alpha_1 parity check with integer matrix power")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

matches = 0
mismatches_detail = []
n_samples = 200

for trial in range(n_samples):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    alpha_1 = 0
    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            sub = np.zeros((size, size), dtype=np.int64)
            verts = list(combo)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            Ak = int_matrix_power(sub, size)
            c = int(np.trace(Ak)) // size
            alpha_1 += c

    expected = alpha_1 % 2
    actual = ((H - 1) // 2) % 2

    if expected == actual:
        matches += 1
    else:
        mismatches_detail.append({
            'bits': bits, 'H': H, 'alpha_1': alpha_1,
            'expected': expected, 'actual': actual
        })

print(f"  Matches: {matches}/{n_samples}")
if mismatches_detail:
    print(f"  First 5 mismatches:")
    for m in mismatches_detail[:5]:
        print(f"    bits={m['bits']}, H={m['H']}, alpha_1={m['alpha_1']}, "
              f"(H-1)/2 mod 2 = {m['actual']}, alpha_1 mod 2 = {m['expected']}")

# Also check the np version for comparison
print("\n  Checking np.linalg.matrix_power vs int_matrix_power:")
np.random.seed(42)

disagree = 0
for trial in range(100):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    alpha_np = 0
    alpha_int = 0

    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            sub = np.zeros((size, size), dtype=np.int64)
            verts = list(combo)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]

            Ak_np = np.linalg.matrix_power(sub, size)
            Ak_int = int_matrix_power(sub, size)

            c_np = int(np.trace(Ak_np)) // size
            c_int = int(np.trace(Ak_int)) // size

            alpha_np += c_np
            alpha_int += c_int

    if alpha_np != alpha_int:
        disagree += 1
        print(f"    DISAGREE at trial {trial}: np={alpha_np}, int={alpha_int}")

print(f"  Disagreements: {disagree}/100")

# ============================================================
# The tr(Sigma^3) gap and 12-divisibility
# ============================================================
print("\n" + "=" * 70)
print("tr(Sigma^3) gap 12-divisibility: WHY?")
print("=" * 70)

print("""
At n=7, all tr(Sigma^3) gaps between ambiguous simplex profiles are
divisible by 12 = 4 * 3 = 2^2 * 3.

Sigma = (n-2)*J_off - sym(A^2)_off (from THM-179 region)

tr(Sigma^k) involves traces of products of A^2 symmetrization.

The factor 12:
  - Factor 4 = 2^2: comes from sigma being an even function of A
    (sigma(u,v) counts common successors + predecessors)
  - Factor 3: comes from 3-cycle structure

Actually, let me think about this more carefully.

tr(Sigma^3) = sum_{i,j,k} Sigma[i,j] * Sigma[j,k] * Sigma[k,i]

Since Sigma is symmetric: this is 2 * (sum over ordered triples) +
diagonal terms. But Sigma has 0 diagonal (Sigma[i,i] = 0 by def).

So tr(Sigma^3) = sum_{distinct i,j,k} Sigma[i,j]*Sigma[j,k]*Sigma[k,i]
Note this counts DIRECTED triangles in the Sigma graph.

Each undirected triangle {i,j,k} contributes:
  2 * Sigma[i,j]*Sigma[j,k]*Sigma[k,i] (two orientations)

Hmm wait, tr(Sigma^3) = sum_i (Sigma^3)[i,i] = sum_i sum_j sum_k Sigma[i,j]*Sigma[j,k]*Sigma[k,i]This automatically includes all orderings. Each unordered triple {i,j,k}
appears in 6 ordered triples (i,j,k), (i,k,j), (j,i,k), (j,k,i), (k,i,j), (k,j,i).
Due to symmetry of Sigma, the product is the same for:
  (i,j,k) and (k,j,i): S[i,j]S[j,k]S[k,i] = S[k,i]S[i,j]S[j,k] (reordered)
So each unordered triple contributes 6 * S[i,j]*S[j,k]*S[k,i] to tr(Sigma^3).
Wait: 2 * the sum over CYCLIC orderings. There are 2 cyclic orderings
per unordered triple. And 3 rotations per cyclic ordering.
So 2 * 3 = 6 terms per unordered triple, each equal.

tr(Sigma^3) = 6 * sum_{i<j<k} S[i,j]*S[j,k]*S[k,i]
Wait, that's not right either. S[i,j]*S[j,k]*S[k,i] = S[j,i]*S[i,k]*S[k,j]
by symmetry. So the two cyclic orderings give the SAME product.
And 3 rotations also give the same product.
So: tr(Sigma^3) = 6 * sum_{i<j<k} S[i,j]*S[j,k]*S[k,i]? No...

Let me just verify: for Sigma symmetric:
tr(Sigma^3) = sum_i,j,k S[i,j]S[j,k]S[k,i]
For fixed unordered {i,j,k} with i<j<k:
  ordered triples with all distinct: (i,j,k),(i,k,j),(j,i,k),(j,k,i),(k,i,j),(k,j,i)
  (i,j,k): S[i,j]S[j,k]S[k,i] = p
  (i,k,j): S[i,k]S[k,j]S[j,i] = p (by symmetry of S)
  (j,i,k): S[j,i]S[i,k]S[k,j] = p
  (j,k,i): S[j,k]S[k,i]S[i,j] = p
  (k,i,j): S[k,i]S[i,j]S[j,k] = p
  (k,j,i): S[k,j]S[j,i]S[i,k] = p
All 6 are equal to p = S[i,j]*S[j,k]*S[k,i].

So tr(Sigma^3) = 6 * sum_{i<j<k} S[i,j]*S[j,k]*S[k,i]

Therefore tr(Sigma^3) is always divisible by 6!

And if the difference of two tr(Sigma^3) values is divisible by 12,
that means the difference of sum_{i<j<k} S[i,j]*S[j,k]*S[k,i] is
divisible by 2.

Since sigma values are integers, the product S[i,j]*S[j,k]*S[k,i]
is an integer. The question is: why is the DIFFERENCE of sums always even?
""")

# Verify tr(Sigma^3) is always divisible by 6
n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

tr3_mod6 = Counter()
for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    sigma = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            s = 0
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][w] and A[v][w]) or (A[w][u] and A[w][v]):
                    s += 1
            sigma[u][v] = sigma[v][u] = s

    tr3 = int(np.trace(sigma @ sigma @ sigma))
    tr3_mod6[tr3 % 6] += 1

print(f"\n  tr(Sigma^3) mod 6: {dict(sorted(tr3_mod6.items()))}")
print(f"  Always divisible by 6: {all(r == 0 for r in tr3_mod6.keys())}")

# Also check mod 12
tr3_mod12 = Counter()
np.random.seed(42)
for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    sigma = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            s = 0
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][w] and A[v][w]) or (A[w][u] and A[w][v]):
                    s += 1
            sigma[u][v] = sigma[v][u] = s

    tr3 = int(np.trace(sigma @ sigma @ sigma))
    tr3_mod12[tr3 % 12] += 1

print(f"  tr(Sigma^3) mod 12: {dict(sorted(tr3_mod12.items()))}")

print("\n\nDone.")
