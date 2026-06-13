#!/usr/bin/env python3
"""
WHY is S(T) ≡ 48 mod 64 at n=7?

S(T) = sum_P prod B[P_i][P_{i+1}] where B[i][j] = ±1 for i≠j.
Each term is prod of n-1=6 factors of ±1, so each term is ±1.
There are n!=5040 terms, each ±1.

S = (# positive terms) - (# negative terms) = 2*(# positive) - 5040.
So S ≡ -5040 mod 2, i.e., S is even.
S = 2P - 5040 where P = # positive terms.

S ≡ 48 mod 64 means 2P - 5040 ≡ 48 mod 64, i.e. 2P ≡ 5088 mod 64, i.e. P ≡ 2544 mod 32.
5040/2 = 2520. P = 2520 + S/2. S/2 ≡ 24 mod 32.

Let's verify and look at the contribution structure.

Key insight: reversal preserves sign at odd n. So terms come in pairs {P, P^rev}
with same sign. # paired positive = same for both. # unpaired = 0 (no palindromic perms).
So P is always even, meaning S ≡ 0 mod 4.
Actually S = 2*(sum of signs of canonical reps) where canonical = P < P^rev.
There are 5040/2 = 2520 canonical pairs.
S = 2 * (sum of ±1 for each canonical pair).
So S/2 = P' - (2520 - P') = 2P' - 2520 where P' = # positive canonical reps.
S/2 ≡ 0 - 2520 mod 2 = -2520 mod 2 = 0. So S/2 is even. S ≡ 0 mod 4.

Can we get more? Consider the action of cyclic shifts.
P = (P_0, ..., P_6). A cyclic shift gives (P_1, ..., P_6, P_0).
This changes the product by:
  prod B[P_1][P_2] * ... * B[P_6][P_0] vs prod B[P_0][P_1] * ... * B[P_5][P_6]
  New product = (old product) * B[P_6][P_0] / B[P_0][P_1]

Hmm, this doesn't simplify easily.

Let me try another approach: consider the value of S modulo small powers of 2
by examining the multilinear expansion of S.
"""
from itertools import permutations
import random
from collections import Counter

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

# Approach: Write B[i][j] = 2*A[i][j] - 1. Then:
# prod B[P_i][P_{i+1}] = prod (2*A[P_i][P_{i+1}] - 1)
# = sum_{S subset of edges} (-1)^{|complement(S)|} * 2^{|S|} * prod_{e in S} A_e
# = sum_{k=0}^{n-1} (-1)^{n-1-k} * 2^k * (sum over k-subsets of path-edges of prod A_e)

# So S(T) = sum_P prod (2A_e - 1) = sum_{k=0}^{n-1} (-1)^{n-1-k} * 2^k * C_k
# where C_k = sum_P (sum over k-subsets of P's edges of prod A_e)

# At n=7, n-1=6:
# S = sum_{k=0}^{6} (-1)^{6-k} * 2^k * C_k
# = C_0 - 2*C_1 + 4*C_2 - 8*C_3 + 16*C_4 - 32*C_5 + 64*C_6

# C_0 = sum_P 1 = 7! = 5040
# C_6 = sum_P prod A_e = H(T) (Hamiltonian path count, edges all forward)

# Wait, that's not quite right. C_k should be:
# C_k = sum_P (number of k-subsets of {P_i->P_{i+1}} where all A[P_i][P_{i+1}]=1)

# Actually let me be more careful.
# prod (2A_e - 1) = sum over subsets S of [n-1]:
#   for each edge index i in S, take factor 2A_{e_i}, 
#   for each i not in S, take factor -1
# = sum_S (-1)^{n-1-|S|} * 2^{|S|} * prod_{i in S} A[P_i][P_{i+1}]

# So S(T) = sum_P sum_S (-1)^{n-1-|S|} * 2^{|S|} * prod_{i in S} A[P_i][P_{i+1}]
# = sum_{k=0}^{n-1} (-1)^{n-1-k} * 2^k * D_k
# where D_k = sum_P sum_{S: |S|=k} prod_{i in S} A[P_i][P_{i+1}]

# D_k counts the total number of (permutation, k-edge-subset) pairs where all selected edges exist.

# For n=7: S = D_0 - 2*D_1 + 4*D_2 - 8*D_3 + 16*D_4 - 32*D_5 + 64*D_6

# D_0 = n! * C(n-1, 0) = 5040 * 1 = 5040
# Actually D_0 = sum_P 1 = 5040 (prod of empty set = 1)

# Modulo 64: S ≡ D_0 - 2*D_1 + 4*D_2 - 8*D_3 + 16*D_4 - 32*D_5 (mod 64)

# Let's compute D_k for a specific tournament
n = 7
bits = 0  # transitive tournament
A = tiling_to_tournament(bits, n)

from itertools import combinations

for bits_sample in [0, 7, 42, 1234, 31000]:
    A = tiling_to_tournament(bits_sample % (1 << 15), n)
    
    # Compute D_k for each k
    D = [0] * n
    for P in permutations(range(n)):
        edges = [A[P[i]][P[i+1]] for i in range(n-1)]
        # For each subset size k, count subsets where all selected edges are 1
        for k in range(n):
            count = 0
            for S in combinations(range(n-1), k):
                if all(edges[i] for i in S):
                    count += 1
            D[k] += count
    
    # Verify S
    S_direct = sum(
        (-1 if (n-1-k) % 2 else 1) * (2**k) * D[k]
        for k in range(n)
    )
    
    # Also direct computation
    S_check = 0
    for P in permutations(range(n)):
        prod_val = 1
        for i in range(n-1):
            prod_val *= (2*A[P[i]][P[i+1]] - 1)
        S_check += prod_val
    
    assert S_direct == S_check, f"Mismatch: {S_direct} vs {S_check}"
    
    print(f"bits={bits_sample % (1<<15):>5}: S={S_direct:>6}")
    for k in range(n):
        sign = (-1)**(n-1-k)
        print(f"  D_{k} = {D[k]:>8}, contribution = {sign}*2^{k}*D_{k} = {sign * (2**k) * D[k]:>10}, "
              f"D_{k} mod 2 = {D[k] % 2}, D_{k} mod 4 = {D[k] % 4}")
    
    # The key: S mod 64 = (D_0 - 2D_1 + 4D_2 - 8D_3 + 16D_4 - 32D_5) mod 64
    s_mod64 = sum((-1)**(n-1-k) * (2**k) * D[k] for k in range(6)) % 64
    print(f"  S mod 64 = {s_mod64}, S mod 64 from direct = {S_direct % 64}")
    print()

# What is D_0 mod 64?
print(f"D_0 = n! = {5040}, D_0 mod 64 = {5040 % 64}")

