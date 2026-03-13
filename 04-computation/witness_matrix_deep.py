"""
witness_matrix_deep.py -- kind-pasteur-2026-03-13-S61

Deep analysis of the WITNESS MATRIX and its doubly-sum-preserving property
under Vitali atoms.

The witness matrix W is n x C(n,2), binary.
  W[k][(u,v)] = 1 iff k is a witness for pair {u,v} (k is in a 3-cycle with u,v).
  Column sums = lambda(u,v)
  Row sums = delta(k) = #{3-cycles containing vertex k}

Under a Vitali atom (reversing a (1,1,2,2)-scored 4-subset S):
  - 24 entries of W change (12 become 1, 12 become 0)
  - All column sums (lambda) preserved
  - All row sums (delta) preserved
  - This is a "doubly-stochastic" perturbation

QUESTIONS:
1. Which 24 entries change? Is the pattern universal?
2. Does the change decompose into "swaps" (specific 2x2 submatrices)?
3. Is there a connection to transportation polytope moves?
4. What is the "hidden dimension" — the rank of W or its cokernel?
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter

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

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def witness_matrix(A, n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    W = np.zeros((n, len(pairs)), dtype=int)
    for pi, (u, v) in enumerate(pairs):
        for k in range(n):
            if k == u or k == v:
                continue
            if (A[u][v]*A[v][k]*A[k][u] + A[v][u]*A[u][k]*A[k][v]):
                W[k][pi] = 1
    return W, pairs

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

n = 7
total_bits = n * (n-1) // 2
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
pair_idx = {p: i for i, p in enumerate(pairs)}

print("=" * 60)
print("WITNESS MATRIX CHANGE STRUCTURE UNDER VITALI ATOMS")
print("=" * 60)

np.random.seed(42)
change_patterns = []

for trial in range(3000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        W_A, _ = witness_matrix(A, n)
        W_B, _ = witness_matrix(B, n)
        diff = W_B - W_A

        if np.count_nonzero(diff) == 0:
            continue  # trivial atom

        # Record which (k, (u,v)) entries changed
        changes = []
        for k in range(n):
            for pi in range(len(pairs)):
                if diff[k][pi] != 0:
                    u, v = pairs[pi]
                    changes.append((k, u, v, int(diff[k][pi])))

        # Classify k: in S or not?
        S = set(subset)
        changes_in_S = [(k, u, v, d) for k, u, v, d in changes if k in S]
        changes_out_S = [(k, u, v, d) for k, u, v, d in changes if k not in S]

        # Classify pair (u,v): both in S, one in S, none in S?
        pair_types = {'SS': [], 'SX': [], 'XX': []}
        for k, u, v, d in changes:
            u_in = u in S
            v_in = v in S
            if u_in and v_in:
                pair_types['SS'].append((k, u, v, d))
            elif u_in or v_in:
                pair_types['SX'].append((k, u, v, d))
            else:
                pair_types['XX'].append((k, u, v, d))

        change_patterns.append({
            'subset': subset,
            'total_changes': len(changes),
            'k_in_S': len(changes_in_S),
            'k_out_S': len(changes_out_S),
            'pair_SS': len(pair_types['SS']),
            'pair_SX': len(pair_types['SX']),
            'pair_XX': len(pair_types['XX']),
            'changes': changes,
        })

        break

    if len(change_patterns) >= 100:
        break

print(f"Non-trivial Vitali pairs: {len(change_patterns)}")

# Analyze the change decomposition
print(f"\nChange decomposition by k-membership and pair-type:")
print(f"{'k in S':>8} {'k out S':>8} {'pair SS':>8} {'pair SX':>8} {'pair XX':>8}")
for cp in change_patterns[:10]:
    print(f"{cp['k_in_S']:8d} {cp['k_out_S']:8d} {cp['pair_SS']:8d} {cp['pair_SX']:8d} {cp['pair_XX']:8d}")

# Are the patterns always the same?
k_in_S_dist = Counter(cp['k_in_S'] for cp in change_patterns)
k_out_S_dist = Counter(cp['k_out_S'] for cp in change_patterns)
pair_SS_dist = Counter(cp['pair_SS'] for cp in change_patterns)
pair_SX_dist = Counter(cp['pair_SX'] for cp in change_patterns)
pair_XX_dist = Counter(cp['pair_XX'] for cp in change_patterns)

print(f"\nDistributions:")
print(f"  k in S: {dict(k_in_S_dist)}")
print(f"  k out S: {dict(k_out_S_dist)}")
print(f"  pair SS: {dict(pair_SS_dist)}")
print(f"  pair SX: {dict(pair_SX_dist)}")
print(f"  pair XX: {dict(pair_XX_dist)}")

# Show a detailed example
cp = change_patterns[0]
print(f"\n--- Detailed example: S = {cp['subset']} ---")
print(f"Changes (k, u, v, delta):")
for k, u, v, d in cp['changes']:
    k_type = "IN_S" if k in set(cp['subset']) else "OUT_S"
    u_type = "S" if u in set(cp['subset']) else "X"
    v_type = "S" if v in set(cp['subset']) else "X"
    print(f"  k={k}({k_type}), pair=({u},{v})({u_type}{v_type}), delta={d:+d}")

# Is the change decomposable into 2x2 swaps?
# A "swap" would be: (k1, p1) and (k2, p2) change by +1,
#                      (k1, p2) and (k2, p1) change by -1.
# This preserves both row and column sums.
print(f"\n{'='*60}")
print("SWAP DECOMPOSITION")
print(f"{'='*60}")

# Try to decompose the change into swaps
for cp_idx in range(min(5, len(change_patterns))):
    cp = change_patterns[cp_idx]
    changes = cp['changes']

    # Separate +1 and -1 changes
    plus = [(k, pair_idx[(min(u,v), max(u,v))]) for k, u, v, d in changes if d == 1]
    minus = [(k, pair_idx[(min(u,v), max(u,v))]) for k, u, v, d in changes if d == -1]

    print(f"\nExample {cp_idx}: S={cp['subset']}")
    print(f"  +1 entries: {plus}")
    print(f"  -1 entries: {minus}")
    print(f"  +1 count: {len(plus)}, -1 count: {len(minus)}")

    # Check: can we pair them into 2x2 swaps?
    # A swap is: rows k1,k2 and cols p1,p2 such that
    #   (k1,p1)=+1, (k2,p2)=+1, (k1,p2)=-1, (k2,p1)=-1
    # This requires finding a perfect matching between plus and minus
    # entries that forms such 2x2 blocks.

    # Build the diff matrix
    diff_matrix = np.zeros((n, len(pairs)), dtype=int)
    for k, u, v, d in changes:
        diff_matrix[k][pair_idx[(min(u,v), max(u,v))]] = d

    # The diff matrix has row sums = 0, column sums = 0.
    # It's a "doubly balanced" signed binary matrix.
    row_sums = diff_matrix.sum(axis=1)
    col_sums = diff_matrix.sum(axis=0)
    print(f"  Row sums: {list(row_sums)}")
    print(f"  Col sums: {[int(c) for c in col_sums if c != 0]} (nonzero entries)")
    print(f"  All row sums 0? {all(r == 0 for r in row_sums)}")
    print(f"  All col sums 0? {all(c == 0 for c in col_sums)}")

    # Rank of the diff matrix
    rank_diff = np.linalg.matrix_rank(diff_matrix)
    print(f"  Rank of diff matrix: {rank_diff}")

# ============================================================
# WITNESS MATRIX RANK AND HIDDEN DIMENSION
# ============================================================
print(f"\n{'='*60}")
print("WITNESS MATRIX RANK AND NULL SPACE")
print(f"{'='*60}")

np.random.seed(42)
ranks = []
for trial in range(200):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    W, _ = witness_matrix(A, n)
    r = np.linalg.matrix_rank(W)
    ranks.append(r)

rank_dist = Counter(ranks)
print(f"n={n}: Witness matrix rank distribution: {dict(sorted(rank_dist.items()))}")
print(f"  Matrix size: {n} x {n*(n-1)//2} = {n} x {len(pairs)}")
print(f"  Max possible rank: min({n}, {len(pairs)}) = {min(n, len(pairs))}")

# The null space dimension = C(n,2) - rank = # of "hidden dimensions"
# These represent the degrees of freedom that can change without affecting
# the row/column sums (= lambda and delta).
null_dims = [len(pairs) - r for r in ranks]
print(f"  Null space dimension (hidden DOF): {Counter(null_dims)}")

# At n=5: C(5,2)=10, so max rank = 5 (since n=5 rows)
print(f"\nn=5 check:")
n5 = 5
for bits in range(5):
    A = bits_to_adj(bits, n5)
    W, _ = witness_matrix(A, n5)
    r = np.linalg.matrix_rank(W)
    print(f"  bits={bits}: rank={r}, null dim={10-r}")

# ============================================================
# TRANSPORTATION POLYTOPE CONNECTION
# ============================================================
print(f"\n{'='*60}")
print("TRANSPORTATION POLYTOPE CONNECTION")
print(f"{'='*60}")
print("""
The witness matrix W lives in the TRANSPORTATION POLYTOPE:
  TP(row_sums, col_sums) = {M >= 0 : M*1 = col_sums, 1^T*M = row_sums}

For integer binary matrices, this is the BIRKHOFF POLYTOPE (generalized).
The vertices of TP are specific integer matrices.
The Vitali atom moves between vertices of TP while staying on the polytope.

Key insight: The number of vertices of TP(delta, lambda) = number of
tournaments with the same (delta, lambda) structure = fiber size at Level 1.

The DIMENSION of TP is:
  (n-1)(C(n,2)-1) - (n-1 + C(n,2)-1 - 1) = ... complicated.

But the key observation is: Vitali atoms are EDGES of the transportation
polytope, connecting adjacent vertices.

The diff matrix D = W_B - W_A has rank 1, 2, or 3 (from our analysis).
Rank 1 would mean D = u*v^T — a single swap.
Rank 2 would mean D = u1*v1^T + u2*v2^T — composition of 2 swaps.
""")

# Check if the 24-entry changes are always rank 2
print("Rank of diff matrices in Vitali changes:")
rank_diffs = Counter()
for cp in change_patterns:
    changes = cp['changes']
    diff_matrix = np.zeros((n, len(pairs)), dtype=int)
    for k, u, v, d in changes:
        diff_matrix[k][pair_idx[(min(u,v), max(u,v))]] = d
    r = np.linalg.matrix_rank(diff_matrix)
    rank_diffs[r] += 1
print(f"  {dict(sorted(rank_diffs.items()))}")

# ============================================================
# THE 24-ENTRY STRUCTURE: WHICH ENTRIES CHANGE?
# ============================================================
print(f"\n{'='*60}")
print("THE 24-ENTRY CHANGE: STRUCTURAL ANALYSIS")
print(f"{'='*60}")

# At n=7 with |S|=4, the 4 vertices in S interact with 3 outside vertices.
# C(4,2) = 6 pairs within S. C(4,1)*3 = 12 pairs crossing S-boundary. C(3,2)=3 outside.
# Total: 6 + 12 + 3 = 21 pairs.
# Each of n=7 vertices can witness each pair. But k can only witness (u,v) if k != u,v.
# For k in S (4 vertices): each can witness C(n,2) - C(1,1)*... = many pairs
# For k outside S (3 vertices): same

# The 24 entries that change: let me tally by (k in/out S) x (pair in SS/SX/XX)
print("Tally of changes by type across all Vitali pairs:")
tally = defaultdict(int)
for cp in change_patterns:
    S = set(cp['subset'])
    for k, u, v, d in cp['changes']:
        k_type = "kS" if k in S else "kX"
        u_in, v_in = u in S, v in S
        if u_in and v_in:
            p_type = "SS"
        elif u_in or v_in:
            p_type = "SX"
        else:
            p_type = "XX"
        tally[(k_type, p_type, d)] += 1

print(f"  {'k_type':>6} {'p_type':>6} {'delta':>6} {'count':>6} {'per_pair':>8}")
for key in sorted(tally.keys()):
    count = tally[key]
    per_pair = count / len(change_patterns)
    print(f"  {key[0]:>6} {key[1]:>6} {key[2]:>+6d} {count:>6} {per_pair:>8.1f}")

# Observation: the changes have a specific structure:
# k in S, pair in SS: these are changes WITHIN the reversed subset
# k in S, pair in SX: boundary effects
# k outside S, pair in SS: external witnesses of internal pairs
# k outside S, pair in SX: external witnesses of boundary pairs
# k outside S, pair in XX: should not change (neither pair endpoint is in S)

print("\nDone.")
