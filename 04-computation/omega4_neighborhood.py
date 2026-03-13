#!/usr/bin/env python3
"""
omega4_neighborhood.py — opus-2026-03-13-S70

THM-151: Omega_3 = C(n,4) - Σ_v[c3(N_in(v)) + c3(N_out(v))]
Can we find a similar formula for Omega_4?

Omega_4 uses 5-vertex paths. At n=6, Omega_4 decomposes over 5-vertex subsets.
Each 5-vertex subtournament's omega4_local is determined by its isomorphism type.

But for a neighborhood formula, we need to express Omega_4 in terms of
vertex neighborhoods (in-neighborhoods, out-neighborhoods).

Strategy: at n=6, each 5-vertex subset is "6 minus one vertex".
So Omega_4 = Σ_v omega4_local(T - v) where T - v is the 5-tournament
obtained by deleting vertex v.

Can omega4_local(T') be expressed in terms of cycle counts and neighborhood
structure of the DELETED vertex?
"""

import numpy as np
from math import comb
from itertools import combinations, permutations
from collections import defaultdict

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

def count_regular_m_paths(A, m):
    n = A.shape[0]
    count = 0
    def dfs(path_set, last, prev, depth):
        nonlocal count
        if depth == m:
            count += 1
            return
        for v in range(n):
            if v in path_set: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path_set.add(v)
            dfs(path_set, v, last, depth + 1)
            path_set.remove(v)
    for start in range(n):
        dfs({start}, start, -1, 0)
    return count

def count_3cycles(A, vertices):
    count = 0
    vlist = list(vertices)
    for i in range(len(vlist)):
        for j in range(i+1, len(vlist)):
            for k in range(j+1, len(vlist)):
                a, b, c = vlist[i], vlist[j], vlist[k]
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    count += 1
    return count

def count_omega3_local(A, vertices):
    """Count regular 3-paths on exactly these vertices."""
    count = 0
    for perm in permutations(vertices):
        v0, v1, v2, v3 = perm
        if (A[v0][v1] and A[v1][v2] and A[v2][v3] and
            A[v0][v2] and A[v1][v3]):
            count += 1
    return count

n = 6
num_pairs = n*(n-1)//2
total = 2**num_pairs

print("="*70)
print(f"OMEGA_4 NEIGHBORHOOD FORMULA: n = {n}")
print("="*70)

# For each tournament, compute Omega_4 via vertex deletion
# and relate to neighborhood invariants
print("  Computing for all tournaments...")

data = []
for bits in range(total):
    A = adj_matrix(n, bits)

    omega4 = count_regular_m_paths(A, 4)

    # For each vertex v, compute omega4_local(T - v)
    per_vertex = []
    for v in range(n):
        others = [w for w in range(n) if w != v]
        # Count regular 4-paths on others
        count = 0
        for perm in permutations(others):
            v0, v1, v2, v3, v4 = perm
            if (A[v0][v1] and A[v1][v2] and A[v2][v3] and A[v3][v4] and
                A[v0][v2] and A[v1][v3] and A[v2][v4]):
                count += 1
        per_vertex.append(count)

    # Neighborhood invariants
    # For each vertex v: in-neighbors, out-neighbors
    # Compute c3 and omega3 of neighborhoods
    c3_in = []
    c3_out = []
    omega3_neigh_in = []
    omega3_neigh_out = []

    for v in range(n):
        n_in = [w for w in range(n) if w != v and A[w][v]]
        n_out = [w for w in range(n) if w != v and A[v][w]]

        c3_in.append(count_3cycles(A, n_in) if len(n_in) >= 3 else 0)
        c3_out.append(count_3cycles(A, n_out) if len(n_out) >= 3 else 0)

        # Omega_3 on the in-neighborhood (if 4+ vertices)
        if len(n_in) >= 4:
            omega3_neigh_in.append(sum(count_omega3_local(A, list(combo))
                                      for combo in combinations(n_in, 4)))
        else:
            omega3_neigh_in.append(0)

        if len(n_out) >= 4:
            omega3_neigh_out.append(sum(count_omega3_local(A, list(combo))
                                       for combo in combinations(n_out, 4)))
        else:
            omega3_neigh_out.append(0)

    data.append({
        'bits': bits,
        'omega4': omega4,
        'per_vertex': per_vertex,
        'c3_in': c3_in,
        'c3_out': c3_out,
        'omega3_in': omega3_neigh_in,
        'omega3_out': omega3_neigh_out,
    })

    if bits % 5000 == 0 and bits > 0:
        print(f"    {bits}/{total}...")

print("  Done computing.")

# Verify: Omega_4 = sum per_vertex
mismatches = 0
for d in data:
    if d['omega4'] != sum(d['per_vertex']):
        mismatches += 1
print(f"\n  Omega_4 = sum per_vertex: {mismatches} mismatches")

# Check: is Omega_4 = C(n,5) - Σ_v[correction(v)]?
# Where correction(v) = something about v's neighborhoods
print(f"\n  C(n,5) = {comb(n,5)}")

# Try: Omega_4 = C(n,5) - Σ_v [c3_in(v)*something + c3_out(v)*something]
# At n=6: each vertex has in-degree and out-degree summing to 5
# c3_in, c3_out count 3-cycles in neighborhoods of 5 vertices

# Collect summary statistics
total_c3_in = [sum(d['c3_in']) for d in data]
total_c3_out = [sum(d['c3_out']) for d in data]
total_o3_in = [sum(d['omega3_in']) for d in data]
total_o3_out = [sum(d['omega3_out']) for d in data]
omega4_list = [d['omega4'] for d in data]

# Check correlations
print(f"\n  Correlations with Omega_4:")
print(f"    corr(sum_c3_in, Omega_4) = {np.corrcoef(total_c3_in, omega4_list)[0,1]:.6f}")
print(f"    corr(sum_c3_out, Omega_4) = {np.corrcoef(total_c3_out, omega4_list)[0,1]:.6f}")
print(f"    corr(sum_o3_in, Omega_4) = {np.corrcoef(total_o3_in, omega4_list)[0,1]:.6f}")
print(f"    corr(sum_o3_out, Omega_4) = {np.corrcoef(total_o3_out, omega4_list)[0,1]:.6f}")

# Try: Omega_4 = a*Σc3_in + b*Σc3_out + c*Σo3_in + d*Σo3_out + e
from numpy.linalg import lstsq
X = np.column_stack([total_c3_in, total_c3_out, total_o3_in, total_o3_out,
                     np.ones(total)])
y = np.array(omega4_list, dtype=float)
coeffs, residuals, _, _ = lstsq(X, y, rcond=None)
print(f"\n  Linear fit: Omega_4 = {coeffs[0]:.4f}*Σc3_in + {coeffs[1]:.4f}*Σc3_out "
      f"+ {coeffs[2]:.4f}*Σo3_in + {coeffs[3]:.4f}*Σo3_out + {coeffs[4]:.4f}")
predicted = X @ coeffs
max_error = np.max(np.abs(y - predicted))
print(f"  Max error: {max_error:.6f}")
if max_error < 0.01:
    print(f"  *** EXACT FIT! ***")

# Also try simpler: Omega_4 = a*Σ(c3_in+c3_out) + b*Σ(o3_in+o3_out) + c
total_c3 = [a+b for a,b in zip(total_c3_in, total_c3_out)]
total_o3 = [a+b for a,b in zip(total_o3_in, total_o3_out)]

X2 = np.column_stack([total_c3, total_o3, np.ones(total)])
coeffs2, residuals2, _, _ = lstsq(X2, y, rcond=None)
print(f"\n  Simpler fit: Omega_4 = {coeffs2[0]:.4f}*Σc3 + {coeffs2[1]:.4f}*Σo3 + {coeffs2[2]:.4f}")
predicted2 = X2 @ coeffs2
max_error2 = np.max(np.abs(y - predicted2))
print(f"  Max error: {max_error2:.6f}")
if max_error2 < 0.01:
    print(f"  *** EXACT FIT! ***")

# Try with global t3 as well
t3_list = [count_3cycles(adj_matrix(n, d['bits']), range(n)) for d in data]
X3 = np.column_stack([total_c3, total_o3, t3_list, np.ones(total)])
coeffs3, _, _, _ = lstsq(X3, y, rcond=None)
print(f"\n  With t3: Omega_4 = {coeffs3[0]:.4f}*Σc3 + {coeffs3[1]:.4f}*Σo3 + "
      f"{coeffs3[2]:.4f}*t3 + {coeffs3[3]:.4f}")
predicted3 = X3 @ coeffs3
max_error3 = np.max(np.abs(y - predicted3))
print(f"  Max error: {max_error3:.6f}")

print("\nDONE.")
