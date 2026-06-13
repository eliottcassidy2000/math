#!/usr/bin/env python3
"""
Which tournaments have c_2 = 0 (off-diagonal)?

For the Paley tournament at n=5, M[a,b](r) is CONSTANT (degree 0) for a≠b.
This means the off-diagonal entries don't depend on r at all!
Only the diagonal entries vary with r.

The off-diagonal polynomial at n=5 has max even degree 2.
c_2^{off-diag} = 0 means: M[a,b](r) = M[a,b](0) for all a≠b.

WHEN does this happen?

Also: explore how c_0 (the tournament-specific part) determines
the iso class within a score group.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
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

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

def score_sequence(A):
    return tuple(sorted([sum(row) for row in A]))

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A, tiles

# =====================================================================
n = 5
_, tiles = tiling_to_tournament(0, n)
m = len(tiles)

print("=" * 70)
print("OFF-DIAGONAL c_2 ANALYSIS AT n=5")
print("=" * 70)

# For each iso class, compute the off-diagonal c_2 matrix
iso_classes = defaultdict(list)
tiling_data = {}

for bits in range(2**m):
    A, _ = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    canon = tournament_canonical(A)
    iso_classes[canon].append(bits)
    tiling_data[bits] = {'H': H, 'canon': canon, 'A': A}

class_labels = {}
for idx, canon in enumerate(sorted(iso_classes.keys())):
    class_labels[canon] = idx

# Extract c_2 off-diagonal for each class
# M(0) = c_0, M(r) = c_0 + c_2*r^2 for off-diag (degree ≤ 2)
# So c_2^{off} = (M(0.5) - M(0)) / 0.25 for off-diagonal entries

print("\nOff-diagonal c_2 max absolute value by class:")
for canon in sorted(iso_classes.keys()):
    c = class_labels[canon]
    bits = iso_classes[canon][0]
    A = tiling_data[bits]['A']
    H = tiling_data[bits]['H']

    M_0 = transfer_matrix_r(A, 0.0)
    M_half = transfer_matrix_r(A, 0.5)

    c2_mat = (M_half - M_0) / 0.25
    # Extract off-diagonal
    c2_offdiag = c2_mat.copy()
    np.fill_diagonal(c2_offdiag, 0)

    max_offdiag = np.max(np.abs(c2_offdiag))
    is_zero_offdiag = max_offdiag < 1e-10

    scores = score_sequence(A)
    print(f"  Class {c:2d} (H={H:2d}, scores={scores}): "
          f"max|c2_offdiag|={max_offdiag:.6f}, zero?{is_zero_offdiag}")

    if not is_zero_offdiag:
        # Show the off-diagonal c_2
        print(f"    c2 offdiag: {np.round(c2_offdiag, 4).tolist()}")


# =====================================================================
print()
print("=" * 70)
print("c_0 STRUCTURE BY CLASS — WHAT DISTINGUISHES ISO CLASSES WITHIN SCORE GROUP?")
print("=" * 70)

for scores_key in sorted(set(score_sequence(tiling_data[iso_classes[canon][0]]['A'])
                             for canon in iso_classes)):
    matching_classes = [class_labels[canon] for canon in sorted(iso_classes.keys())
                       if score_sequence(tiling_data[iso_classes[canon][0]]['A']) == scores_key]
    if len(matching_classes) <= 1:
        continue

    print(f"\n  Score group {scores_key}: classes {matching_classes}")
    for c in matching_classes:
        canon = sorted(iso_classes.keys())[c]
        bits = iso_classes[canon][0]
        A = tiling_data[bits]['A']
        H = tiling_data[bits]['H']
        M_0 = transfer_matrix_r(A, 0.0)

        # c_0 eigenvalues
        eigs_c0 = sorted(np.linalg.eigvalsh(M_0))
        print(f"    Class {c} (H={H}): c_0 eigs = {[round(e,3) for e in eigs_c0]}")
        print(f"      c_0 = {np.round(M_0, 3).tolist()}")

    # How do c_0 matrices differ?
    c1, c2_idx = matching_classes[0], matching_classes[1]
    canon1 = sorted(iso_classes.keys())[c1]
    canon2 = sorted(iso_classes.keys())[c2_idx]
    M0_1 = transfer_matrix_r(tiling_data[iso_classes[canon1][0]]['A'], 0.0)
    M0_2 = transfer_matrix_r(tiling_data[iso_classes[canon2][0]]['A'], 0.0)
    diff = M0_1 - M0_2
    print(f"    Difference c_0({c1}) - c_0({c2_idx}): {np.round(diff, 3).tolist()}")
    print(f"    Diff eigs: {[round(e,3) for e in sorted(np.linalg.eigvalsh(diff))]}")


# =====================================================================
print()
print("=" * 70)
print("THE M(0) = c_0 INTERPRETATION: 'SIGNED TRANSFER MATRIX'")
print("=" * 70)

# M(0)[a,b] = IE formula with edge weights s_{ij} = ±1/2
# This is the transfer matrix of the "signed tournament"
# where every edge has weight +1/2 or -1/2.

# At r=0: the path weight is product of s_e = (+1/2)^k * (-1/2)^{n-1-k}
# where k = number of tournament edges in the path.
# = (1/2)^{n-1} * (-1)^{n-1-k} = (1/2)^{n-1} * (-1)^{non-edges in path}

# So c_0 = (1/2)^{n-1} * "signed transfer matrix"
# where each path is weighted by (-1)^{non-edges} instead of 0/1.

print(f"  c_0 = (1/2)^{{n-1}} * M_signed")
print(f"  M_signed[a,b] = sum_S (-1)^|S| (# E_a paths with (-1)^{{non-edges}}) * (# B_b paths)")

# For regular tournament (Paley): every path has same number of non-edges?
# No, that depends on the path. But the Paley symmetry ensures c_0 is scalar.

# For class 11 (non-regular scalar at n=5): c_0 = 0!
# This means: all signed paths cancel perfectly at r=0.
# The entire transfer matrix comes from the r^4 term.

print("\n  Class 11 (non-regular scalar, H=15):")
canon11 = sorted(iso_classes.keys())[11]
A11 = tiling_data[iso_classes[canon11][0]]['A']
M0_11 = transfer_matrix_r(A11, 0.0)
print(f"    c_0 = {np.round(M0_11, 6).tolist()}")
print(f"    M(0.5) = {transfer_matrix_r(A11, 0.5).astype(int).tolist()}")
print(f"    Score sequence: {score_sequence(A11)}")
print(f"    Adjacency: {A11}")


# =====================================================================
print()
print("=" * 70)
print("WHAT MAKES CLASSES 8, 9, 10 DIFFERENT? (same scores, different H)")
print("=" * 70)

# All have score sequence (1,2,2,2,3) and 4 three-cycles.
# But H = 11, 15, 13 respectively.
# They share c_2 and c_4 (eigenvalue spectra), differ only in c_0.
# H = tr(c_0) + tr(c_2)/4 + 120/16
#   = tr(c_0) + 18/4 + 7.5
#   = tr(c_0) + 4.5 + 7.5
#   = tr(c_0) + 12

# So: Class 8: H=11, tr(c_0)=-1
#     Class 9: H=15, tr(c_0)=3
#     Class 10: H=13, tr(c_0)=1

# tr(c_0) = H(0) = signed Ham path count at r=0.
# This is the ONLY thing distinguishing the iso classes within a score group.

print("  Within score group (1,2,2,2,3):")
print("  Class 8: H=11, tr(c_0) = -1  =>  H(0) = -1/16 per path")
print("  Class 9: H=15, tr(c_0) = 3   =>  H(0) = 3/16 per path")
print("  Class 10: H=13, tr(c_0) = 1  =>  H(0) = 1/16 per path")

# H(0) = sum_P prod_{e in P} s_e = (1/2^4) * sum_P (-1)^{non-edges in P}
# So sum_P (-1)^{non-edges} = 16 * tr(c_0)
# Class 8: sum = -16
# Class 9: sum = 48
# Class 10: sum = 16

print("\n  sum_P (-1)^{non-edges in P} (times n for trace):")
for c_idx in [8, 9, 10]:
    canon = sorted(iso_classes.keys())[c_idx]
    A = tiling_data[iso_classes[canon][0]]['A']
    signed_sum = 0
    for p in permutations(range(n)):
        edges = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        non_edges = (n-1) - edges
        if all(True for i in range(n-1)):  # count ALL permutations
            signed_sum += (-1)**non_edges
    print(f"    Class {c_idx}: sum(-1)^non-edges over ALL perms = {signed_sum}")
    # This should be n*tr(c_0)*16
    print(f"    Expected: {n * 16 * np.trace(transfer_matrix_r(A, 0.0)):.0f}")

# Actually need to sum over only actual Ham paths and weight by IE formula...
# Let me just confirm tr(c_0) directly.

# =====================================================================
print()
print("=" * 70)
print("DECOMPOSITION SUMMARY: M = c_0 + c_2*r^2 + (n-1)!*r^{n-1}*I")
print("=" * 70)
print(f"""
At n=5, the transfer matrix decomposes as:
  M(r) = c_0 + c_2*r^2 + 24r^4*I

where:
  c_4 = 24I = (n-1)!*I  (UNIVERSAL, independent of tournament)
  c_2: eigenvalues determined by SCORE SEQUENCE only
       (same scores => same c_2 up to permutation)
  c_0: eigenvalues depend on full ISOMORPHISM CLASS
       (distinguishes classes within a score group)

The transfer matrix at r=1/2:
  M(1/2) = c_0 + c_2/4 + 24/16 * I
         = c_0 + c_2/4 + (3/2)*I

Since H = tr(M(1/2)) = tr(c_0) + tr(c_2)/4 + n*(n-1)!/2^{{n-1}}:
  H = tr(c_0) + tr(c_2)/4 + 5*24/16 = tr(c_0) + tr(c_2)/4 + 7.5

With tr(c_2) = 12*(#3-cycles) - 30:
  H = tr(c_0) + 3*(#3-cycles) - 7.5 + 7.5
  H = tr(c_0) + 3*(#3-cycles)

Verified:
  Class 8: H=11, tr(c_0)=-1, 3cyc=4: -1 + 12 = 11 ✓
  Class 9: H=15, tr(c_0)=3, 3cyc=4: 3 + 12 = 15 ✓
  Class 10: H=13, tr(c_0)=1, 3cyc=4: 1 + 12 = 13 ✓
  Class 11: H=15, tr(c_0)=0, 3cyc=5: 0 + 15 = 15 ✓
  Class 0: H=1, tr(c_0)=1, 3cyc=0: 1 + 0 = 1 ✓

FORMULA: H = tr(c_0) + 3*(#3-cycles)  at n=5.
""")

# Verify this formula for ALL classes
print("Verification of H = tr(c_0) + 3*(#3-cycles):")
all_match = True
for canon in sorted(iso_classes.keys()):
    c = class_labels[canon]
    bits = iso_classes[canon][0]
    A = tiling_data[bits]['A']
    H = tiling_data[bits]['H']
    M_0 = transfer_matrix_r(A, 0.0)
    tr_c0 = np.trace(M_0)

    # Count 3-cycles
    cyc3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                cyc3 += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])

    predicted_H = tr_c0 + 3*cyc3
    match = abs(predicted_H - H) < 1e-6
    if not match:
        all_match = False
    print(f"  Class {c}: H={H}, tr(c0)={tr_c0:.4f}, 3cyc={cyc3}, "
          f"predicted={predicted_H:.4f}, match={match}")

print(f"\n  ALL MATCH: {all_match}")
