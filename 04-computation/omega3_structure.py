#!/usr/bin/env python3
"""
omega3_structure.py — opus-2026-03-13-S70

Deep analysis of Omega_3 at n=5,6:
1. What determines Omega_3? (It's uncorrelated with t_3!)
2. What ARE the 4-cycle Walsh modes that Omega_3 picks out?
3. Can we find an explicit formula like Omega_2 = C(n,3) - t_3?

At n=5: Omega_3 ∈ {3, 5}, Walsh degree {0, 4}, 15 nonzero coefficients.
At n=6: Omega_3 ∈ {9, 10, 11, 12, 15}, Walsh degree {0, 4}, 45 nonzero coefficients.

Omega_3 counts regular 3-paths: (v0,v1,v2,v3) with
v0→v1, v1→v2, v2→v3, AND v0→v2, v1→v3.
This is: {v0,v1,v2} transitive with v0→v1→v2, AND {v1,v2,v3} transitive with v1→v2→v3.
The arc (v0,v3) is NOT constrained.

So Omega_3 = #{ordered 4-tuples (v0,v1,v2,v3) with v0→v1→v2→v3, v0→v2, v1→v3}.
"""

import numpy as np
from math import comb
from itertools import combinations, permutations
from collections import Counter, defaultdict

def adj_matrix(n, T_bits):
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (T_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_omega3(A):
    """Count regular 3-paths directly."""
    n = A.shape[0]
    count = 0
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0: continue
            if not A[v0][v1]: continue
            for v2 in range(n):
                if v2 == v0 or v2 == v1: continue
                if not A[v1][v2]: continue
                if not A[v0][v2]: continue  # regularity: v0→v2
                for v3 in range(n):
                    if v3 == v0 or v3 == v1 or v3 == v2: continue
                    if not A[v2][v3]: continue
                    if not A[v1][v3]: continue  # regularity: v1→v3
                    count += 1
    return count

def count_linked_transitive_pairs(A):
    """Count pairs of overlapping transitive triples sharing an edge.
    A "linked pair" is ({v0,v1,v2}, {v1,v2,v3}) where both are transitive
    with v0→v1→v2 and v1→v2→v3.
    """
    return count_omega3(A)

def count_transitive_triples(A):
    """Count ORDERED transitive triples (v0,v1,v2) with v0→v1→v2, v0→v2."""
    n = A.shape[0]
    count = 0
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0 or not A[v0][v1]: continue
            for v2 in range(n):
                if v2 == v0 or v2 == v1: continue
                if A[v1][v2] and A[v0][v2]:
                    count += 1
    return count

def count_3cycles(A):
    n = A.shape[0]
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    count += 1
    return count

# Count 4-cycles
def count_4cycles(A):
    """Count oriented 4-cycles in tournament.
    An oriented 4-cycle on {a,b,c,d} is a cycle a→b→c→d→a.
    Note: in a tournament, a 4-cycle exists iff the 4-vertex subtournament
    has score sequence (1,1,2,2) = not transitive.
    """
    n = A.shape[0]
    count = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                for d in range(n):
                    if d == a or d == b or d == c: continue
                    if A[c][d] and A[d][a]:
                        count += 1
    return count

def count_4vertex_nontrans(A):
    """Count 4-vertex subsets that are NOT transitive."""
    n = A.shape[0]
    count = 0
    for combo in combinations(range(n), 4):
        # Check if transitive: score sequence should be (0,1,2,3) up to permutation
        scores = []
        for v in combo:
            s = sum(A[v][w] for w in combo if w != v)
            scores.append(s)
        scores.sort()
        if scores != [0, 1, 2, 3]:
            count += 1
    return count

# ============================================================
# n = 5 analysis
# ============================================================
print("="*70)
print("OMEGA_3 STRUCTURE: n = 5")
print("="*70)

n = 5
num_pairs = n*(n-1)//2

omega3_vals = {}
t3_vals = {}
tt_vals = {}  # transitive triples (ordered)
c4_vals = {}  # 4-cycles (oriented)
nt4_vals = {}  # non-transitive 4-vertex subsets

for bits in range(2**num_pairs):
    A = adj_matrix(n, bits)
    omega3_vals[bits] = count_omega3(A)
    t3_vals[bits] = count_3cycles(A)
    tt_vals[bits] = count_transitive_triples(A)
    c4_vals[bits] = count_4cycles(A)
    nt4_vals[bits] = count_4vertex_nontrans(A)

print(f"  Omega_3 values: {sorted(set(omega3_vals.values()))}")
print(f"  t_3 values: {sorted(set(t3_vals.values()))}")
print(f"  Ordered TT values: {sorted(set(tt_vals.values()))}")
print(f"  c_4 (oriented 4-cycles) values: {sorted(set(c4_vals.values()))}")
print(f"  Non-trans 4-sets values: {sorted(set(nt4_vals.values()))}")

# Correlations
o3_list = [omega3_vals[b] for b in range(2**num_pairs)]
t3_list = [t3_vals[b] for b in range(2**num_pairs)]
tt_list = [tt_vals[b] for b in range(2**num_pairs)]
c4_list = [c4_vals[b] for b in range(2**num_pairs)]
nt4_list = [nt4_vals[b] for b in range(2**num_pairs)]

print(f"\n  Correlations with Omega_3:")
print(f"    corr(Omega_3, t_3) = {np.corrcoef(o3_list, t3_list)[0,1]:.6f}")
print(f"    corr(Omega_3, TT) = {np.corrcoef(o3_list, tt_list)[0,1]:.6f}")
print(f"    corr(Omega_3, c_4) = {np.corrcoef(o3_list, c4_list)[0,1]:.6f}")
print(f"    corr(Omega_3, NT4) = {np.corrcoef(o3_list, nt4_list)[0,1]:.6f}")

# Check: is Omega_3 = f(c_4)?
c4_to_o3 = defaultdict(set)
for b in range(2**num_pairs):
    c4_to_o3[c4_vals[b]].add(omega3_vals[b])
print(f"\n  c_4 → Omega_3 mapping:")
for c4, o3_set in sorted(c4_to_o3.items()):
    print(f"    c_4={c4}: Omega_3 ∈ {sorted(o3_set)}")

# Check: is Omega_3 = f(NT4)?
nt4_to_o3 = defaultdict(set)
for b in range(2**num_pairs):
    nt4_to_o3[nt4_vals[b]].add(omega3_vals[b])
print(f"\n  NT4 → Omega_3 mapping:")
for nt4, o3_set in sorted(nt4_to_o3.items()):
    print(f"    NT4={nt4}: Omega_3 ∈ {sorted(o3_set)}")

# Is Omega_3 a function of the number of non-transitive 4-subsets?
# At n=5: C(5,4) = 5 four-vertex subsets.
# NT4 + trans4 = 5. trans4 = 5 - NT4.
# And Omega_3 might count something about these subsets.

# Check: Omega_3 = TT - 3*t_3? or similar
print(f"\n  Looking for formula Omega_3 = a*TT + b:")
for b in range(2**num_pairs):
    # TT = 6*(C(5,3) - t_3) = 6*Omega_2
    # At n=5: Omega_2 = 10 - t_3, TT = 6*Omega_2 = 6*(10-t_3) = 60 - 6*t_3
    pass

# Actually: ordered TT = 6*unordered TT = 6*(C(n,3)-t_3) when we count with multiplicity?
# No: for each unordered transitive triple {a,b,c} with a→b→c, the ordered version is
# (a,b,c) — one order only. But we could also have (a,c,b) as non-ordered, etc.
# Actually for a transitive triple with unique total order a→b→c→ (a→c too),
# the only valid ordered TT is (a,b,c). So ordered_TT = unordered_TT = C(n,3)-t_3.
# Wait, my count function counts (v0,v1,v2) with v0→v1→v2 and v0→v2.
# For each transitive triple, there's exactly ONE such ordering.
# So ordered_TT = unordered_TT = C(n,3)-t_3.

for b in [0, 1, 100, 500]:
    A = adj_matrix(n, b)
    t3 = count_3cycles(A)
    tt = count_transitive_triples(A)
    expected = comb(n,3) - t3
    print(f"    T={b}: TT={tt}, C(n,3)-t3={expected}, match: {tt == expected}")

# So TT = C(n,3) - t_3 = Omega_2. Good.
# And corr(Omega_3, TT) should equal corr(Omega_3, Omega_2) = 0.
# corr(Omega_3, TT) = corr(Omega_3, -t_3) = -corr(Omega_3, t_3) = 0. ✓

# The question is: what 4-VERTEX invariant determines Omega_3?
# Each ordered 4-tuple contributing to Omega_3 uses 4 specific vertices.
# On these 4 vertices, we need: v0→v1→v2→v3, v0→v2, v1→v3.
# That's 5 edges specified. The 6th edge (v0,v3) is free.

# For a given 4-vertex subset {a,b,c,d}: how many ordered 4-tuples
# (v0,v1,v2,v3) from this subset satisfy the Omega_3 condition?

def omega3_on_subset(A, vertices):
    """Count regular 3-paths using exactly these 4 vertices."""
    count = 0
    for perm in permutations(vertices):
        v0, v1, v2, v3 = perm
        if A[v0][v1] and A[v1][v2] and A[v2][v3] and A[v0][v2] and A[v1][v3]:
            count += 1
    return count

print(f"\n  Omega_3 per 4-vertex subset:")
# Check a few tournaments
for b in [0, 100, 500, 900]:
    A = adj_matrix(n, b)
    total = 0
    per_sub = []
    for combo in combinations(range(n), 4):
        c = omega3_on_subset(A, list(combo))
        per_sub.append(c)
        total += c
    print(f"    T={b}: Omega_3={omega3_vals[b]}, sum per-subset={total}, per-sub={per_sub}")

# So Omega_3 = sum over all C(n,4) four-vertex subsets of omega3_on_subset.
# And omega3_on_subset depends on the 4-vertex tournament structure.

# What values can omega3_on_subset take?
print(f"\n  omega3_on_subset values by tournament type:")
per_sub_vals = defaultdict(list)
for b in range(2**num_pairs):
    A = adj_matrix(n, b)
    for combo in combinations(range(n), 4):
        c = omega3_on_subset(A, list(combo))
        # Score sequence of the 4-vertex subtournament
        scores = tuple(sorted(sum(A[v][w] for w in combo if w != v) for v in combo))
        per_sub_vals[scores].append(c)

for scores, vals in sorted(per_sub_vals.items()):
    unique_vals = sorted(set(vals))
    print(f"    scores {scores}: omega3_local ∈ {unique_vals}, counts = {Counter(vals)}")

# A 4-vertex tournament has score sequence either:
# (0,1,2,3) = transitive (1 type)
# (1,1,2,2) = non-transitive (1 type up to isomorphism)
# So:
# - Transitive 4-tour: omega3_local should be deterministic
# - Non-transitive: also deterministic

print("\nDONE.")
