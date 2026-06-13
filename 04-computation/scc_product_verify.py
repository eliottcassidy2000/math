#!/usr/bin/env python3
"""
scc_product_verify.py — opus-2026-03-14-S71g

Quick verification: H = ∏ H(SCCᵢ) for non-strongly-connected tournaments.
Uses Landau's criterion for fast SC checking.
"""

from itertools import permutations, combinations
from collections import defaultdict
from math import comb

def make_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def is_sc(A, n):
    scores = sorted(sum(A[i]) for i in range(n))
    partial = 0
    for k in range(1, n):
        partial += scores[k-1]
        if partial <= comb(k, 2):
            return False
    return True

def find_scc_boundaries(A, n):
    scores = [sum(A[i]) for i in range(n)]
    order = sorted(range(n), key=lambda v: scores[v])
    sorted_scores = [scores[v] for v in order]
    boundaries = []
    partial = 0
    for k in range(1, n+1):
        partial += sorted_scores[k-1]
        if partial == comb(k, 2):
            boundaries.append(k)
    sccs = []
    start = 0
    for end in boundaries:
        sccs.append([order[i] for i in range(start, end)])
        start = end
    return sccs

print("=" * 60)
print("SCC PRODUCT FORMULA VERIFICATION")
print("=" * 60)

for n in range(3, 7):
    total_edges = n * (n - 1) // 2
    num_t = 2 ** total_edges
    sc_count = 0
    non_sc = 0
    mismatches = 0
    scc_patterns = defaultdict(int)

    h_sc = set()
    h_nonsc = set()

    for bits in range(num_t):
        A = make_tournament(bits, n)
        H = count_hp(A, n)

        if is_sc(A, n):
            sc_count += 1
            h_sc.add(H)
            scc_patterns[(n,)] += 1
        else:
            non_sc += 1
            sccs = find_scc_boundaries(A, n)
            pattern = tuple(sorted([len(s) for s in sccs], reverse=True))
            scc_patterns[pattern] += 1

            # Verify product formula
            product = 1
            for scc in sccs:
                m = len(scc)
                sub_A = [[0]*m for _ in range(m)]
                for i in range(m):
                    for j in range(m):
                        if i != j:
                            sub_A[i][j] = A[scc[i]][scc[j]]
                product *= count_hp(sub_A, m)

            if H != product:
                mismatches += 1
            h_nonsc.add(H)

    print(f"\n  n={n}: {sc_count} SC ({100*sc_count/num_t:.1f}%), {non_sc} non-SC")
    print(f"  Product formula mismatches: {mismatches}/{non_sc}")
    print(f"  SC H-values: {sorted(h_sc)}")
    print(f"  Non-SC H-values: {sorted(h_nonsc)}")
    only_sc = h_sc - h_nonsc
    if only_sc:
        print(f"  SC-ONLY H-values (prime): {sorted(only_sc)}")

    print(f"  SCC patterns:")
    for p in sorted(scc_patterns.keys()):
        cnt = scc_patterns[p]
        label = " [SC]" if len(p) == 1 else ""
        print(f"    {'+'.join(str(s) for s in p):10s}: {cnt:6d} ({100*cnt/num_t:5.1f}%){label}")

# Achievability semigroup (using only n≤6 exact data)
print(f"\n{'='*60}")
print("H-SPECTRUM FORBIDDEN VALUES")
print(f"{'='*60}")

all_h = set()
for n in range(3, 7):
    for bits in range(2**(n*(n-1)//2)):
        A = make_tournament(bits, n)
        all_h.add(count_hp(A, n))

# Products of pairs
products = set()
for a in all_h:
    for b in all_h:
        if a * b <= 200:
            products.add(a * b)
all_h.update(products)
# One more round
products2 = set()
for a in all_h:
    if a > 200:
        continue
    for b in [1,3,5,9,11,13,15]:
        if a * b <= 200:
            products2.add(a * b)
all_h.update(products2)

never_odd = sorted(h for h in range(1, 101, 2) if h not in all_h)
print(f"\n  Odd H ≤ 100 NOT achievable (exact n≤6 + products): {never_odd}")

forbidden_7 = [7*3**k for k in range(4) if 7*3**k <= 100]
print(f"  7·3^k: {forbidden_7}")
print(f"  21·3^k: {[21*3**k for k in range(3) if 21*3**k <= 100]}")
unexplained = [h for h in never_odd if h not in {7, 21, 63}]
print(f"  Unexplained: {unexplained}")
