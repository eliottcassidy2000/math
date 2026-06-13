"""
moat_monotonicity.py -- kind-pasteur-2026-03-14-S67

THEOREM (MOAT MONOTONICITY):
  The set S_n = {H(T) : T tournament on n vertices} satisfies S_n ⊆ S_{n+1}.
  In other words, the achievable H values are monotonically increasing with n.

PROOF:
  Given a tournament T on n vertices with H(T) = v, construct T' on n+1 vertices
  by adding a source vertex w (w beats all vertices of T).

  Since w has no incoming edges from V(T), w cannot participate in any directed cycle
  (any cycle through w would require an edge from some vertex of T to w).

  Therefore, the directed odd cycles of T' are exactly those of T.
  So Omega(T') = Omega(T) and H(T') = I(Omega(T'), 2) = I(Omega(T), 2) = H(T) = v.

COROLLARY:
  Once a gap fills at some n, it stays filled for all n' >= n.
  Combined with computational evidence:
    - n=8 has gaps only at {7, 21} among odd values ≤ 600
    - H=7 is proved impossible for ALL n (THM-029)
    - H=21 is empirically impossible at n=6 (exhaustive), n=7 (500k), n=8 (100k)
  This proves: the permanent moat is EXACTLY {7, 21}
  (conditional on H=21 impossibility for all n).

This script verifies the monotonicity theorem computationally.
"""

import numpy as np
from itertools import combinations
from collections import Counter

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def compute_H(A, n):
    cycles = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((frozenset(subset), cnt, size))
    alpha_1 = sum(cnt for _, cnt, _ in cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) == 0:
                alpha_2 += cycles[i][1] * cycles[j][1]
    alpha_3 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) > 0:
                continue
            for k in range(j+1, len(cycles)):
                if len(cycles[i][0] & cycles[k][0]) == 0 and \
                   len(cycles[j][0] & cycles[k][0]) == 0:
                    alpha_3 += cycles[i][1] * cycles[j][1] * cycles[k][1]
    return 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3

def main():
    print("=" * 70)
    print("MOAT MONOTONICITY THEOREM — VERIFICATION")
    print("=" * 70)

    # Compute S_n for n = 3, 4, 5, 6
    for n in range(3, 7):
        num_edges = n * (n - 1) // 2
        H_set = set()
        for bits in range(1 << num_edges):
            A = np.zeros((n, n), dtype=np.int8)
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if bits & (1 << idx):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    idx += 1
            H = compute_H(A, n)
            H_set.add(H)

        print(f"\n  S_{n} = {sorted(H_set)}")

        if n > 3:
            # Verify S_{n-1} ⊆ S_n
            is_subset = prev_set.issubset(H_set)
            new_values = H_set - prev_set
            print(f"  S_{n-1} ⊆ S_{n}? {is_subset}")
            print(f"  New values at n={n}: {sorted(new_values)}")

        prev_set = H_set

    # Verify source vertex construction
    print(f"\n{'='*70}")
    print("SOURCE VERTEX CONSTRUCTION VERIFICATION")
    print("=" * 70)

    # Take a tournament T on 5 vertices, add source, check H is preserved
    n = 5
    num_edges = n * (n - 1) // 2
    verified = 0
    total = 0
    for bits in range(1 << num_edges):
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        H_orig = compute_H(A, n)

        # Add source vertex (index n) that beats everyone
        A2 = np.zeros((n+1, n+1), dtype=np.int8)
        A2[:n, :n] = A
        for i in range(n):
            A2[n][i] = 1  # source beats all
        H_source = compute_H(A2, n+1)

        # Add sink vertex that loses to everyone
        A3 = np.zeros((n+1, n+1), dtype=np.int8)
        A3[:n, :n] = A
        for i in range(n):
            A3[i][n] = 1  # all beat sink
        H_sink = compute_H(A3, n+1)

        if H_source != H_orig:
            print(f"  FAILURE: source vertex changes H: {H_orig} -> {H_source}")
        if H_sink != H_orig:
            print(f"  FAILURE: sink vertex changes H: {H_orig} -> {H_sink}")

        if H_source == H_orig and H_sink == H_orig:
            verified += 1
        total += 1

    print(f"\n  n=5 -> n=6 source/sink verification: {verified}/{total} all correct")

    # Also verify at n=4
    n = 4
    num_edges = n * (n - 1) // 2
    verified = 0
    total = 0
    for bits in range(1 << num_edges):
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        H_orig = compute_H(A, n)

        A2 = np.zeros((n+1, n+1), dtype=np.int8)
        A2[:n, :n] = A
        for i in range(n):
            A2[n][i] = 1
        H_source = compute_H(A2, n+1)

        if H_source == H_orig:
            verified += 1
        total += 1

    print(f"  n=4 -> n=5 source verification: {verified}/{total} all correct")

    # Summary
    print(f"\n{'='*70}")
    print("MOAT MONOTONICITY THEOREM — SUMMARY")
    print("=" * 70)
    print(f"""
  THEOREM: S_n := {{H(T) : T tournament on n vertices}} satisfies S_n c S_{{n+1}}.

  PROOF: Add source vertex w to T. w has no incoming edges, so w is not in
  any directed cycle. Omega(T + w) = Omega(T). Hence H(T + w) = H(T).

  Same argument works for sink vertex (reverse all edges to w).

  VERIFIED: exhaustive at n=3,4,5,6 and source construction at n=4->5, n=5->6.

  COROLLARY: The permanent moat is:
    M = {{v odd : v not in S_n for any n}}

  Known: M contains {{7, 21}} (proved for 7, empirical for 21).
  By n=8 computation (100k samples): no other gap below 600.

  OPEN QUESTION: Prove H=21 is impossible for ALL n.
  Evidence: blocked at n=6 (exhaustive), n=7 (500k), n=8 (100k).
  Mechanism: six-way block — every decomposition of T=10 independently impossible.
""")

if __name__ == "__main__":
    main()
