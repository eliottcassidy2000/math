#!/usr/bin/env python3
"""
walsh_glmy_bridge.py — opus-2026-03-13-S70

Bridge between Walsh-Fourier framework and GLMY path homology.

Key question: Can we express GLMY Betti numbers in terms of Walsh coefficients?

For a tournament T on n vertices:
  - H(T) = Hamiltonian path count (Walsh-analyzable via THM-077)
  - chi(T) = Euler characteristic of GLMY complex
  - beta_m(T) = GLMY Betti numbers

For circulant T on Z_p:
  - H = p * H_per (per-eigenspace contribution)
  - chi = p * chi_per where chi_per = sum (-1)^m Omega_m
  - The Omega_m are per-eigenspace dimensions

CONNECTION 1: chi and H
  chi = H - (other terms from Worpitzky)
  Actually: chi = sum (-1)^m dim(Omega_m).
  For a tournament, Omega_0 = n, Omega_1 = #edges = C(n,2) = n(n-1)/2, etc.
  chi = n - C(n,2) + Omega_2 - Omega_3 + ...
  And H = Omega_{n-1} (since n-vertex HP = full allowed path).

  So chi = sum_{m=0}^{n-1} (-1)^m Omega_m, and H = Omega_{n-1}.
  The relationship is: both are linear combinations of the Omega vector.

CONNECTION 2: For circulant tournaments, Walsh of H and chi should be related.
  The Walsh decomposition of H has known structure (THM-077).
  Can we Walsh-decompose chi?

CONNECTION 3: Omega_m for general tournaments via path counting.
  Omega_m = #regular m-paths in T.
  This is a combinatorial invariant of T.

Let's investigate numerically.
"""

import numpy as np
import sys
sys.path.insert(0, '04-computation')

# First: compute Omega for small general tournaments and compare with Walsh data
from itertools import permutations

def adj_matrix(n, T):
    """Build adjacency matrix from tournament encoding.
    T encodes pairs (i,j) with i<j: T[k]=1 means i→j, T[k]=0 means j→i.
    """
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if T[idx]:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_allowed_paths(A, length):
    """Count allowed paths of given length in tournament with adjacency matrix A."""
    n = A.shape[0]
    count = 0

    def dfs(path, depth):
        nonlocal count
        if depth == length:
            count += 1
            return
        last = path[-1]
        for v in range(n):
            if v not in path and A[last][v] == 1:
                path.append(v)
                dfs(path, depth + 1)
                path.pop()

    for start in range(n):
        dfs([start], 0)

    return count

def count_regular_paths(A, length):
    """Count REGULAR m-paths (= dim Omega_m) in tournament.
    A regular path is one where all faces are also allowed paths.
    For path (v0,...,vm), face i = (v0,...,v_{i-1},v_{i+1},...,vm).
    Face i is allowed iff v_{i-1} → v_{i+1} (for 0 < i < m).
    Faces 0 and m are always allowed (they're subpaths).

    So a path (v0,...,vm) is regular iff for all 0 < i < m: A[v_{i-1}][v_{i+1}] = 1.
    """
    n = A.shape[0]
    if length <= 1:
        return count_allowed_paths(A, length)

    count = 0

    def dfs(path, depth):
        nonlocal count
        if depth == length:
            # Check regularity: for each interior vertex i (1 <= i <= depth-1),
            # we need A[path[i-1]][path[i+1]] = 1
            regular = True
            for i in range(1, depth):
                if A[path[i-1]][path[i+1]] != 1:
                    regular = False
                    break
            if regular:
                count += 1
            return
        last = path[-1]
        for v in range(n):
            if v not in path and A[last][v] == 1:
                path.append(v)
                dfs(path, depth + 1)
                path.pop()

    for start in range(n):
        dfs([start], 0)

    return count

def omega_vector(A):
    """Compute Omega_m for all m."""
    n = A.shape[0]
    result = []
    for m in range(n):
        result.append(count_regular_paths(A, m))
    return result

def hamilton_count(A):
    """Count Hamiltonian paths."""
    n = A.shape[0]
    count = 0
    for perm in permutations(range(n)):
        valid = all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1))
        if valid:
            count += 1
    return count

def walsh_transform(f_vals, n):
    """Walsh-Hadamard transform of a function on tournaments.
    f_vals: dict mapping tournament encoding (tuple of bits) to value.
    Returns: dict mapping subset S (frozenset) to Walsh coefficient.
    """
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((i,j))
    num_pairs = len(pairs)

    # Walsh transform: hat{f}[S] = (1/2^num_pairs) * sum_T f(T) * (-1)^{|S ∩ T|}
    # where T is identified with its arc set
    results = {}

    # Only compute for small |S|
    from itertools import combinations
    for size in range(min(5, num_pairs+1)):
        for combo in combinations(range(num_pairs), size):
            S = frozenset(combo)
            val = 0
            for T_bits, f_val in f_vals.items():
                sign = 1
                for idx in S:
                    if T_bits[idx]:
                        sign *= -1
                val += sign * f_val
            results[S] = val / (2**num_pairs)

    return results

# Test at n=5
print("="*70)
print("WALSH-GLMY BRIDGE: n=5 ANALYSIS")
print("="*70)

n = 5
num_pairs = n*(n-1)//2

# Collect data for all tournaments
H_vals = {}
chi_vals = {}
omega_vecs = {}
omega_m_vals = {m: {} for m in range(n)}

print(f"\nComputing Omega for all 2^{num_pairs} = {2**num_pairs} tournaments on {n} vertices...")

for bits in range(2**num_pairs):
    T = tuple((bits >> i) & 1 for i in range(num_pairs))
    A = adj_matrix(n, T)

    omega = omega_vector(A)
    H = omega[-1]  # H = Omega_{n-1}
    chi = sum((-1)**m * omega[m] for m in range(n))

    H_vals[T] = H
    chi_vals[T] = chi
    omega_vecs[T] = tuple(omega)
    for m in range(n):
        omega_m_vals[m][T] = omega[m]

print(f"Done. H range: [{min(H_vals.values())}, {max(H_vals.values())}]")
print(f"chi range: [{min(chi_vals.values())}, {max(chi_vals.values())}]")

# Check: chi = H - ... what?
# chi = Omega_0 - Omega_1 + Omega_2 - Omega_3 + Omega_4
# Omega_0 = n = 5, Omega_1 = C(n,2) = 10 for tournaments
# So chi = 5 - 10 + Omega_2 - Omega_3 + H
# chi - H = -5 + Omega_2 - Omega_3

print(f"\nOmega_0 = {set(omega_m_vals[0].values())} (should be {{{n}}})")
print(f"Omega_1 = {set(omega_m_vals[1].values())} (should be {{{n*(n-1)//2}}} for tournaments)")
print(f"Omega_2 range: [{min(omega_m_vals[2].values())}, {max(omega_m_vals[2].values())}]")
print(f"Omega_3 range: [{min(omega_m_vals[3].values())}, {max(omega_m_vals[3].values())}]")
print(f"Omega_4 = H range: [{min(omega_m_vals[4].values())}, {max(omega_m_vals[4].values())}]")

# Walsh transform of chi
print(f"\nWalsh transform of chi (low degrees):")
chi_hat = walsh_transform(chi_vals, n)
for S, val in sorted(chi_hat.items(), key=lambda x: (len(x[0]), x[0])):
    if abs(val) > 1e-10 and len(S) <= 4:
        print(f"  S={set(S)}: hat_chi = {val:.6f}")

# Walsh transform of H
print(f"\nWalsh transform of H (low degrees):")
H_hat = walsh_transform(H_vals, n)
for S, val in sorted(H_hat.items(), key=lambda x: (len(x[0]), x[0])):
    if abs(val) > 1e-10 and len(S) <= 4:
        print(f"  S={set(S)}: hat_H = {val:.6f}")

# Walsh transform of Omega_2
print(f"\nWalsh transform of Omega_2 (low degrees):")
o2_hat = walsh_transform(omega_m_vals[2], n)
for S, val in sorted(o2_hat.items(), key=lambda x: (len(x[0]), x[0])):
    if abs(val) > 1e-10 and len(S) <= 4:
        print(f"  S={set(S)}: hat_Omega2 = {val:.6f}")

# Walsh transform of Omega_3
print(f"\nWalsh transform of Omega_3 (low degrees):")
o3_hat = walsh_transform(omega_m_vals[3], n)
for S, val in sorted(o3_hat.items(), key=lambda x: (len(x[0]), x[0])):
    if abs(val) > 1e-10 and len(S) <= 4:
        print(f"  S={set(S)}: hat_Omega3 = {val:.6f}")

# Key question: is chi a WALSH POLYNOMIAL in the arc variables?
# i.e., can chi(T) be written as a polynomial in the arc signs?
# Since chi = Omega_0 - Omega_1 + Omega_2 - Omega_3 + Omega_4
# and each Omega_m is a function of T, chi is a function of T.
# But is it a LOW-DEGREE Walsh polynomial?

print(f"\nDegree analysis of Walsh transforms:")
for name, hat_vals in [("H", H_hat), ("chi", chi_hat), ("Omega_2", o2_hat), ("Omega_3", o3_hat)]:
    nonzero_degrees = set()
    for S, val in hat_vals.items():
        if abs(val) > 1e-10:
            nonzero_degrees.add(len(S))
    print(f"  {name}: nonzero Walsh degrees = {sorted(nonzero_degrees)}")

# chi correlation with H and cycle counts
print(f"\nCorrelation analysis:")
H_list = [H_vals[T] for T in sorted(H_vals.keys())]
chi_list = [chi_vals[T] for T in sorted(chi_vals.keys())]
corr = np.corrcoef(H_list, chi_list)[0,1]
print(f"  corr(H, chi) = {corr:.6f}")

# Check: is chi = H + const for some tournaments?
diffs = set()
for T in H_vals:
    diffs.add(chi_vals[T] - H_vals[T])
print(f"  chi - H takes values: {sorted(diffs)}")

# Omega_m correlation with H
for m in [2, 3]:
    om_list = [omega_m_vals[m][T] for T in sorted(H_vals.keys())]
    corr = np.corrcoef(H_list, om_list)[0,1]
    print(f"  corr(H, Omega_{m}) = {corr:.6f}")

# Check: does Omega_2 have a simple formula?
# For tournaments, every 3-subset {a,b,c} has edges a→b→c, a→c or some permutation.
# The 3-cycle condition: {a,b,c} is a 3-cycle if it has a cycle.
# A 2-path (a,b,c) is allowed iff a→b and b→c.
# It's regular iff additionally a→c (face 1 = edge (a,c) must be allowed).
# So Omega_2 = #{ordered triples (a,b,c) with a→b, b→c, a→c} = #{transitive triples} × 3!/?
# No: (a,b,c) with a→b, b→c, a→c is a TRANSITIVE triple with the path going a→b→c and a→c.
# The regular condition: face 0 = (b,c) [always allowed], face 1 = (a,c) [need a→c], face 2 = (a,b) [always allowed].
# So Omega_2 = #{ordered (a,b,c): a→b, b→c, AND a→c}.
# This is exactly 2 × (number of transitive triples), since for each transitive triple {a<b<c}
# with unique linear order a→b→c, there's exactly one such ordered triple...
# Actually no, (a,b,c) is an ORDERED triple. For a transitive triple where the order is a→b→c→...,
# there are 1 path a→b→c. Hmm.
# Each transitive triple has a unique source and sink. If source=a, middle=b, sink=c,
# then a→b, a→c, b→c. The 2-paths within this triple are:
# (a,b,c): a→b ✓, b→c ✓, regularity: a→c ✓ → contributes
# (a,c,b): a→c ✓, c→b ✗ → not a path
# etc. So each transitive triple contributes exactly 1 regular 2-path.
# For a 3-cycle triple: say a→b, b→c, c→a. Then:
# (a,b,c): a→b ✓, b→c ✓, but a→c? We have c→a, so a→c ✗ → not regular.
# (b,c,a): b→c ✓, c→a ✓, but b→a? We have a→b, so b→a ✗ → not regular.
# (c,a,b): c→a ✓, a→b ✓, but c→b? We have b→c, so c→b ✗ → not regular.
# So 3-cycles contribute 0 regular 2-paths!

# Therefore: Omega_2 = #transitive triples = C(n,3) - t_3 where t_3 = #3-cycles.
# Let's verify:

def count_3cycles(A):
    n = A.shape[0]
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if {i,j,k} form a 3-cycle
                has_cycle = (
                    (A[i][j] and A[j][k] and A[k][i]) or
                    (A[i][k] and A[k][j] and A[j][i])
                )
                if has_cycle:
                    count += 1
    return count

print(f"\n{'='*70}")
print(f"OMEGA_2 = C(n,3) - t_3 VERIFICATION (n={n})")
print(f"{'='*70}")

mismatches = 0
for bits in range(min(100, 2**num_pairs)):
    T = tuple((bits >> i) & 1 for i in range(num_pairs))
    A = adj_matrix(n, T)
    omega2 = omega_m_vals[2][T]
    t3 = count_3cycles(A)
    from math import comb as comb_
    expected = comb_(n, 3) - t3
    if omega2 != expected:
        mismatches += 1
        print(f"  MISMATCH: T={T}, Omega_2={omega2}, C(5,3)-t3={expected}")

from math import comb
print(f"  Checked {min(100, 2**num_pairs)} tournaments: {mismatches} mismatches")
if mismatches == 0:
    print(f"  CONFIRMED: Omega_2 = C(n,3) - t_3 for all tested tournaments")

# Now check if chi = f(H, t3, ...) has a SIMPLE formula
# chi = 5 - 10 + Omega_2 - Omega_3 + H
# = -5 + (C(5,3)-t3) - Omega_3 + H
# = -5 + 10 - t3 - Omega_3 + H
# = 5 - t3 - Omega_3 + H

# Is Omega_3 = f(t3, t4, ...)?
# Or more fundamentally: what determines Omega_3?

print(f"\n{'='*70}")
print(f"CHI = 5 - t3 - Omega_3 + H VERIFICATION")
print(f"{'='*70}")

sample = list(H_vals.keys())[:100]
mismatches = 0
for T in sample:
    A = adj_matrix(n, T)
    t3 = count_3cycles(A)
    chi = chi_vals[T]
    H = H_vals[T]
    o3 = omega_m_vals[3][T]
    expected_chi = 5 - t3 - o3 + H
    if chi != expected_chi:
        mismatches += 1

print(f"  Verified chi = 5 - t3 - Omega_3 + H for {len(sample)} tournaments: {mismatches} mismatches")

# What is Omega_3 in terms of tournament invariants?
print(f"\n{'='*70}")
print(f"OMEGA_3 ANALYSIS")
print(f"{'='*70}")

# Omega_3 counts regular 3-paths (v0,v1,v2,v3) where:
# v0→v1, v1→v2, v2→v3 (path condition)
# v0→v2 (face 1 regularity)
# v1→v3 (face 2 regularity)
# So it counts ordered 4-tuples (v0,v1,v2,v3) all distinct with
# v0→v1, v0→v2, v1→v2, v1→v3, v2→v3
# i.e., {v0,v1,v2} is transitive with v0→v1→v2, and {v1,v2,v3} is transitive with v1→v2→v3.
# The overlap: v1→v2 is shared.
# What about v0→v3? Not required for regularity!

# Check: does v0→v3 matter?
count_with_v03 = 0
count_without_v03 = 0
T_example = tuple((0 >> i) & 1 for i in range(num_pairs))
A_ex = adj_matrix(n, T_example)

# Actually let me just analyze the distribution of Omega_3
o3_vals = sorted(set(omega_m_vals[3].values()))
print(f"  Omega_3 takes values: {o3_vals}")

# Distribution
from collections import Counter
o3_dist = Counter(omega_m_vals[3].values())
print(f"  Distribution: {dict(sorted(o3_dist.items()))}")

# Correlation with t3
t3_list = []
o3_list = []
for T in sorted(H_vals.keys()):
    A = adj_matrix(n, T)
    t3_list.append(count_3cycles(A))
    o3_list.append(omega_m_vals[3][T])

corr = np.corrcoef(t3_list, o3_list)[0,1]
print(f"  corr(t3, Omega_3) = {corr:.6f}")
print(f"  corr(H, Omega_3) = {np.corrcoef(H_list, o3_list)[0,1]:.6f}")

# Check if Omega_3 = f(t3) exactly
t3_to_o3 = defaultdict(set)
for T in H_vals:
    A = adj_matrix(n, T)
    t3 = count_3cycles(A)
    o3 = omega_m_vals[3][T]
    t3_to_o3[t3].add(o3)

print(f"\n  t3 → Omega_3 (is it a function?):")
for t3, o3_set in sorted(t3_to_o3.items()):
    print(f"    t3={t3}: Omega_3 ∈ {sorted(o3_set)}")

# Check chi distribution
chi_dist = Counter(chi_vals.values())
print(f"\n  chi distribution: {dict(sorted(chi_dist.items()))}")
print(f"  Unique chi values: {sorted(set(chi_vals.values()))}")

print("\nDONE.")
