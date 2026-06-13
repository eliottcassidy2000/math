#!/usr/bin/env python3
"""
Even-r polynomial analysis of the non-VT scalar-M tournament at n=5.

Since M = 3*I for both the VT (Paley/C_5) and non-VT class,
they must have the SAME even-r polynomial:
  M(r) = c_0 + c_2*r^2 + c_4*r^4
with c_0 + c_2/4 + c_4/16 = 3*I.

But c_4 = 24*I universally, so c_0 + c_2/4 = 3*I - (24/16)*I = 3*I - 1.5*I = 1.5*I.

Question: do the non-VT scalar-M tournaments have the same c_0, c_2 as the VT ones?
At n=5, there are 2 VT iso classes (Paley: c_0 = 0.5*I, C_5: c_0 = -2*I).

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
import numpy as np

def count_paths_weighted(A, verts, r_val, start=None, end=None):
    total = 0.0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        w = 1.0
        for i in range(len(p)-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

def transfer_matrix_r(A, r_val):
    n = len(A)
    M = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0.0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_weighted(A, S_verts, r_val, end=a)
                    bb = count_paths_weighted(A, R_verts, r_val, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

# The non-VT scalar-M tournament (from our analysis)
# bits=001100 in tiling encoding, but let's build it directly
# Adjacency: 0->4, 1->0,3, 2->0,1, 3->0,2, 4->1,2,3
A_nonVT = [
    [0, 0, 0, 0, 1],
    [1, 0, 0, 1, 0],
    [1, 1, 0, 0, 0],
    [1, 0, 1, 0, 0],
    [0, 1, 1, 1, 0],
]

# Paley tournament on Z/5: gen = {1,2} (QR mod 5)
A_paley = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(5):
        if i != j and (j-i) % 5 in {1, 2}:
            A_paley[i][j] = 1

# C_5 tournament: gen = {1,4} (just the cycle)
A_C5 = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(5):
        if i != j and (j-i) % 5 in {1, 4}:
            A_C5[i][j] = 1

n = 5
# Sample points for polynomial fitting in u = r^2
u_samples = np.array([0.0, 0.04, 0.16, 0.36, 0.64])
r_samples = np.sqrt(u_samples)

print("=" * 70)
print("EVEN-r POLYNOMIAL FOR NON-VT SCALAR-M vs VT TOURNAMENTS AT n=5")
print("=" * 70)

for name, A in [("Non-VT (scores 1,2,2,2,3)", A_nonVT),
                ("Paley (regular)", A_paley),
                ("C_5 (regular)", A_C5)]:
    print(f"\n  {name}:")

    # Compute full M(r) at sample points
    M_samples = [transfer_matrix_r(A, rv) for rv in r_samples]

    # For scalar-M tournaments, M = f(r)*I, so just track M[0,0]
    diag_vals = [M[0][0] for M in M_samples]
    offdiag_vals = [M[0][1] for M in M_samples]

    # Fit degree-2 polynomial in u = r^2 (since c_4 = 24*I, max degree is 2)
    poly_diag = np.polyfit(u_samples, diag_vals, 2)
    poly_off = np.polyfit(u_samples, offdiag_vals, 2)

    c4_diag = poly_diag[0]
    c2_diag = poly_diag[1]
    c0_diag = poly_diag[2]

    print(f"    M[0,0](r) = {c0_diag:.4f} + {c2_diag:.4f}*r^2 + {c4_diag:.4f}*r^4")
    print(f"    M[0,1](r) = {poly_off[2]:.4f} + {poly_off[1]:.4f}*r^2 + {poly_off[0]:.4f}*r^4")

    # Full matrix c_0, c_2
    # Since M is scalar: all entries are the same
    # But let's check the FULL matrix structure at r=0
    M0 = M_samples[0]
    print(f"    M(r=0) = c_0:")
    for row in M0:
        print(f"      {[f'{x:.2f}' for x in row]}")

    # Check all diagonal entries
    print(f"    Diagonal c_k: c0={c0_diag:.4f}, c2={c2_diag:.4f}, c4={c4_diag:.4f}")
    print(f"    At r=0.5: M[0,0] = {c0_diag + c2_diag*0.25 + c4_diag*0.0625:.4f}")

    # Count 3-cycles
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                t3 += A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i]
    print(f"    3-cycles: {t3}")
    print(f"    tr(c_2) = n*c2_diag = {n*c2_diag:.4f} (expected 12*t3 - 30 = {12*t3 - 30})")

# =====================================================================
print("\n" + "=" * 70)
print("WHAT MAKES THE NON-VT TOURNAMENT POSITION-UNIFORM?")
print("=" * 70)

# The anti-automorphism gives N[v,j] = N[sigma(v), n-1-j]
# The Z/3Z automorphism gives N[1,j] = N[2,j] = N[3,j]
# Together: N[0,j] = N[4, 4-j]
# And: N[0,j] + N[4,j] + 3*N[1,j] = 15
# Let's verify N is indeed constant = 3

A = A_nonVT
ham = [p for p in permutations(range(n)) if all(A[p[i]][p[i+1]] == 1 for i in range(n-1))]
print(f"\n  H = {len(ham)} Hamiltonian paths")

N = np.zeros((n, n), dtype=int)
for p in ham:
    for j, v in enumerate(p):
        N[v][j] += 1

print(f"  N[v,j] matrix:")
for v in range(n):
    print(f"    vertex {v} (deg={sum(A_nonVT[v])}): {list(N[v])}")

# The paths themselves
print(f"\n  All {len(ham)} Hamiltonian paths:")
for p in ham:
    edges = [f"{p[i]}->{p[i+1]}" for i in range(n-1)]
    print(f"    {p}: {', '.join(edges)}")

# Check: position of vertex 0 in each path
pos_0 = [p.index(0) for p in ham]
pos_4 = [p.index(4) for p in ham]
print(f"\n  Vertex 0 positions: {sorted(pos_0)} (histogram: {[pos_0.count(j) for j in range(n)]})")
print(f"  Vertex 4 positions: {sorted(pos_4)} (histogram: {[pos_4.count(j) for j in range(n)]})")

# =====================================================================
# Why is vertex 0 (out-degree 1!) at position 0 exactly 3 times?
# Vertex 0 has only ONE outgoing edge (to vertex 4).
# If vertex 0 is first (position 0), the path must start 0->4->...
# That means there are exactly 3 paths starting with 0->4.
# =====================================================================

paths_starting_0 = [p for p in ham if p[0] == 0]
print(f"\n  Paths starting with vertex 0: {len(paths_starting_0)}")
for p in paths_starting_0:
    print(f"    {p}")

paths_ending_0 = [p for p in ham if p[4] == 0]
print(f"\n  Paths ending with vertex 0: {len(paths_ending_0)}")
for p in paths_ending_0:
    print(f"    {p}")

# Note: vertex 0 has in-degree 3 (from 1,2,3) and out-degree 1 (to 4).
# Starting paths need 0's only out-neighbor (4) to be second.
# Ending paths need one of 0's in-neighbors (1,2,3) to be second-to-last.

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
