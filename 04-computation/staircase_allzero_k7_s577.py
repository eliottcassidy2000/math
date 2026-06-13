#!/usr/bin/env python3
"""
staircase_allzero_k7_s577.py   monad-researcher-2026-06-02-S577

INV-190: compute H(all-0 staircase) at k=7 (n=14) using Held-Karp bitmask DP.
Extend the existing markov_staircase_h.py (which stopped at k=6, n=12)
to k=7 (n=14) and also k=8 (n=16) for a recurrence search.

The all-0 interleaved staircase at n=2k:
  - Pairs: (0,1), (2,3), ..., (2k-2, 2k-1)
  - Global ranking: rank[2p]=p (dominants), rank[2p+1]=k+p (recessives)
  - Within-pair: odd (recessive) beats even (dominant): 2p+1 → 2p
  - Between pairs: lower rank beats higher: rank[i]<rank[j] => i→j
"""

import time


def build_staircase_adj(k):
    """Return adjacency matrix A where A[i][j]=1 iff i→j."""
    n = 2 * k
    rank = {}
    for p in range(k):
        rank[2 * p] = p          # dominant vertices
        rank[2 * p + 1] = k + p  # recessive vertices

    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            same_pair = (i // 2 == j // 2)
            if same_pair:
                # Recessive beats dominant: 2p+1 → 2p
                if i % 2 == 1 and j % 2 == 0:
                    A[i][j] = 1
            else:
                # Lower rank beats higher
                if rank[i] < rank[j]:
                    A[i][j] = 1
    return A, n


def held_karp_H(A, n):
    """
    Count Hamiltonian paths in tournament via Held-Karp bitmask DP.
    dp[mask][v] = #paths ending at v covering vertices in mask.
    Time: O(2^n * n^2). Space: O(2^n * n).
    """
    FULL = (1 << n) - 1
    # Initialize
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            if dp[mask][v] == 0:
                continue
            # Extend: add vertex u not in mask where v→u
            rest = FULL ^ mask
            u = 0
            tmp = rest
            while tmp:
                lsb = tmp & (-tmp)
                u = lsb.bit_length() - 1
                if A[v][u]:
                    dp[mask | lsb][u] += dp[mask][v]
                tmp ^= lsb

    return sum(dp[FULL][v] for v in range(n))


def score_seq(A, n):
    return sorted(sum(A[v]) for v in range(n))


def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if ((A[i][j] and A[j][k] and A[k][i]) or
                        (A[i][k] and A[k][j] and A[j][i])):
                    c3 += 1
    return c3


print("=" * 60)
print("All-0 Interleaved Staircase: H via Held-Karp DP")
print("=" * 60)

known_values = {2: 5, 3: 29, 4: 233, 5: 2489, 6: 33773}

results = []
for k in range(2, 9):  # k=2..8 (n=4..16)
    n = 2 * k
    t0 = time.time()
    A, _ = build_staircase_adj(k)
    H = held_karp_H(A, n)
    elapsed = time.time() - t0

    sc = score_seq(A, n)
    c3 = count_3cycles(A, n)

    expected = known_values.get(k, "?")
    ok = (H == expected) if isinstance(expected, int) else "NEW"

    print(f"\nk={k}, n={n}: H={H}  (expected={expected}, match={ok})  t={elapsed:.3f}s")
    print(f"  score_seq: {sc}")
    print(f"  c3={c3}  (formula k(k-1)={k*(k-1)}, match={c3==k*(k-1)})")
    results.append((k, n, H, c3))

# Recurrence search
print("\n" + "=" * 60)
print("RECURRENCE ANALYSIS")
print("=" * 60)

H_vals = [r[2] for r in results]
c3_vals = [r[3] for r in results]

print(f"\nH(k=2..{1+len(H_vals)}): {H_vals}")
print(f"c3(k=2..{1+len(c3_vals)}): {c3_vals}")

# Growth ratios
print("\nH(k)/H(k-1) growth ratios:")
for i in range(1, len(H_vals)):
    k = results[i][0]
    ratio = H_vals[i] / H_vals[i - 1]
    print(f"  k={k}: {H_vals[i]}/{H_vals[i-1]} = {ratio:.6f}")

# Search for linear recurrences of order 2, 3
print("\nLinear recurrence search (order 2: a[k] = p*a[k-1] + q*a[k-2]):")
for i in range(2, len(H_vals)):
    if H_vals[i-2] == 0 or H_vals[i-1] == 0:
        continue
    # System: a[i] = p*a[i-1] + q*a[i-2]
    #         a[i-1] = p*a[i-2] + q*a[i-3]  (if available)
    if i >= 3:
        # Solve 2x2 system
        a, b = H_vals[i-1], H_vals[i-2]
        c, d = H_vals[i-2], H_vals[i-3]
        det = a * d - b * c
        if det != 0:
            p_num = H_vals[i] * d - b * H_vals[i-1]
            q_num = a * H_vals[i-1] - H_vals[i] * c
            # Check if p and q are rational with small denominators
            from math import gcd
            g1 = gcd(abs(p_num), abs(det))
            g2 = gcd(abs(q_num), abs(det))
            p_frac = (p_num // g1, det // g1)
            q_frac = (q_num // g2, det // g2)
            # Verify
            pred = (p_num * H_vals[i-1] + q_num * H_vals[i-2]) // det
            # Only print if integer p, q
            if det % p_num == 0 and det % q_num == 0 and p_frac[1] in (1,-1) and q_frac[1] in (1,-1):
                p_val = p_frac[0] * p_frac[1]
                q_val = q_frac[0] * q_frac[1]
                verify = p_val * H_vals[i-1] + q_val * H_vals[i-2]
                if verify == H_vals[i]:
                    print(f"  k={results[i][0]}: a[k] = {p_val}*a[k-1] + {q_val}*a[k-2] (from last 3 terms)")

print("\nSearch for p,q,r with a[k] = p*a[k-1] + q*a[k-2] + r*a[k-3] over all k:")
# Use first 4 values to solve for p,q,r
if len(H_vals) >= 4:
    # a[3]=p*a[2]+q*a[1]+r*a[0], a[4]=p*a[3]+q*a[2]+r*a[1], a[5]=p*a[4]+q*a[3]+r*a[2]
    # Build linear system
    import sys
    H = H_vals
    # Row 0: p*H[2]+q*H[1]+r*H[0] = H[3]
    # Row 1: p*H[3]+q*H[2]+r*H[1] = H[4]
    # Row 2: p*H[4]+q*H[3]+r*H[2] = H[5]
    if len(H_vals) >= 6:
        # 3x3 system with exact fractions
        from fractions import Fraction
        M = [
            [Fraction(H[2]), Fraction(H[1]), Fraction(H[0])],
            [Fraction(H[3]), Fraction(H[2]), Fraction(H[1])],
            [Fraction(H[4]), Fraction(H[3]), Fraction(H[2])],
        ]
        b = [Fraction(H[3]), Fraction(H[4]), Fraction(H[5])]
        # Gaussian elimination
        for col in range(3):
            pivot = None
            for row in range(col, 3):
                if M[row][col] != 0:
                    pivot = row
                    break
            if pivot is None:
                print("  Singular system")
                break
            M[col], M[pivot] = M[pivot], M[col]
            b[col], b[pivot] = b[pivot], b[col]
            for row in range(3):
                if row != col and M[row][col] != 0:
                    factor = M[row][col] / M[col][col]
                    for c2 in range(3):
                        M[row][c2] -= factor * M[col][c2]
                    b[row] -= factor * b[col]
        sol = [b[row] / M[row][row] for row in range(3)]
        print(f"  p={sol[0]}, q={sol[1]}, r={sol[2]}")
        # Verify on remaining terms
        ok = True
        for i in range(3, len(H_vals)):
            pred = sol[0] * H[i-1] + sol[1] * H[i-2] + sol[2] * H[i-3]
            if pred != H[i]:
                ok = False
                print(f"  FAILS at k={results[i][0]}: pred={pred}, actual={H[i]}")
        if ok:
            print(f"  Order-3 recurrence VERIFIED for all computed terms!")
            print(f"  a[k] = {sol[0]}*a[k-1] + {sol[1]}*a[k-2] + {sol[2]}*a[k-3]")
        # Use recurrence to predict k=8 if we only computed k=7
        if len(H_vals) >= 6:
            k8_pred = int(sol[0] * H[-1] + sol[1] * H[-2] + sol[2] * H[-3])
            print(f"\n  Predicted H(k=8): {k8_pred}")

# Markov check
print("\nMarkov equation check x²+y²+z² = 3xyz for consecutive H:")
for i in range(len(H_vals) - 2):
    x, y, z = H_vals[i], H_vals[i+1], H_vals[i+2]
    lhs = x**2 + y**2 + z**2
    rhs = 3*x*y*z
    k_vals = (results[i][0], results[i+1][0], results[i+2][0])
    print(f"  k={k_vals}: LHS-RHS = {lhs-rhs}")

# Key result
print("\n" + "=" * 60)
print("KEY RESULT")
print("=" * 60)
new_k = [r for r in results if r[0] >= 7]
for r in new_k:
    k, n, H_val, c3 = r
    print(f"H(all-0 staircase, k={k}, n={n}) = {H_val}")
    print(f"  c3 = {c3} (= k(k-1) = {k*(k-1)}, {'✓' if c3==k*(k-1) else '✗'})")
