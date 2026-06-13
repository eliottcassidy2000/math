#!/usr/bin/env python3
"""
Full polynomial extraction of M(r) = c_0 + c_2*r^2 + c_4*r^4 at n=5.

Previous analysis used only M(0) and M(1/2) to extract c_0 and "c_2",
but M[0,0] has degree 4 in r (with only even powers: r^0 and r^4).
So we need at LEAST 3 sample points per entry to extract c_0, c_2, c_4.

Extracting the CORRECT even-power polynomial for each iso class
and checking:
1. tr(c_2) = 12 * (#3-cycles)?
2. Score sequence determines c_k eigenvalues?
3. c_4 structure?

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

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check all 3! orderings for 3-cycles
                cyc = (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
                count += cyc
    return count

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
print(f"FULL EVEN-r POLYNOMIAL M(r) = c_0 + c_2*r^2 + c_4*r^4 AT n={n}")
print("=" * 70)

# Collect iso classes
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

# For even-power polynomial of degree 4: M(r) = c_0 + c_2*r^2 + c_4*r^4
# We need M at r^2 = 0, t1, t2, i.e. r = 0, sqrt(t1), sqrt(t2)
# Easier: sample at r = 0, 0.2, 0.4, 0.6, 0.8, 1.0 and fit even polynomial
# by substituting u = r^2 and fitting c_0 + c_2*u + c_4*u^2

r_samples = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
u_samples = np.array([r**2 for r in r_samples])

class_polys = {}

for canon in sorted(iso_classes.keys()):
    c = class_labels[canon]
    bits = iso_classes[canon][0]
    A = tiling_data[bits]['A']
    H = tiling_data[bits]['H']

    # Sample M at each r
    M_samples = [transfer_matrix_r(A, rv) for rv in r_samples]

    # Extract c_0, c_2, c_4 for each entry
    c0 = np.zeros((n, n))
    c2 = np.zeros((n, n))
    c4 = np.zeros((n, n))

    for a in range(n):
        for b in range(n):
            vals = np.array([M_samples[i][a][b] for i in range(len(r_samples))])
            # Fit degree-2 polynomial in u = r^2
            coeffs = np.polyfit(u_samples, vals, 2)  # coeffs[0]*u^2 + coeffs[1]*u + coeffs[2]
            c4[a][b] = coeffs[0]
            c2[a][b] = coeffs[1]
            c0[a][b] = coeffs[2]

    # Verify reconstruction
    M_check = transfer_matrix_r(A, 0.35)
    M_recon = c0 + c2 * 0.35**2 + c4 * 0.35**4
    err = np.max(np.abs(M_check - M_recon))

    scores = score_sequence(A)
    cyc3 = count_3cycles(A)

    eigs_c0 = sorted(np.linalg.eigvalsh(c0))
    eigs_c2 = sorted(np.linalg.eigvalsh(c2))
    eigs_c4 = sorted(np.linalg.eigvalsh(c4))

    class_polys[c] = {
        'H': H, 'c0': c0, 'c2': c2, 'c4': c4,
        'scores': scores, 'cyc3': cyc3,
        'eigs_c0': eigs_c0, 'eigs_c2': eigs_c2, 'eigs_c4': eigs_c4,
        'A': A
    }

    print(f"\n  Class {c:2d} (H={H:2d}, scores={scores}, 3-cyc={cyc3}, recon_err={err:.2e}):")
    print(f"    tr(c0)={np.trace(c0):.4f}, tr(c2)={np.trace(c2):.4f}, tr(c4)={np.trace(c4):.4f}")
    print(f"    c0 eigs: {[round(e,3) for e in eigs_c0]}")
    print(f"    c2 eigs: {[round(e,3) for e in eigs_c2]}")
    print(f"    c4 eigs: {[round(e,3) for e in eigs_c4]}")

# =====================================================================
print()
print("=" * 70)
print("tr(c_k) FORMULAS")
print("=" * 70)

print("\n  tr(c_0) vs H(0), tr(c_2) vs 3-cycles, tr(c_4) pattern:")
for c in sorted(class_polys.keys()):
    d = class_polys[c]
    print(f"  Class {c:2d}: H={d['H']:2d}, 3cyc={d['cyc3']}, "
          f"tr(c0)={np.trace(d['c0']):.4f}, "
          f"tr(c2)={np.trace(d['c2']):.4f}, "
          f"tr(c4)={np.trace(d['c4']):.4f}")

# Check tr(c_2) = 12 * (#3-cycles)?
print("\n  tr(c_2) / (#3-cycles):")
for c in sorted(class_polys.keys()):
    d = class_polys[c]
    if d['cyc3'] > 0:
        ratio = np.trace(d['c2']) / d['cyc3']
        print(f"    Class {c}: ratio = {ratio:.4f}")

# Check: tr(c_0) + 0.25*tr(c_2) + 0.0625*tr(c_4) = H at odd n
print("\n  Verification: tr(c_0) + 0.25*tr(c_2) + 0.0625*tr(c_4) = H?")
for c in sorted(class_polys.keys()):
    d = class_polys[c]
    lhs = np.trace(d['c0']) + 0.25*np.trace(d['c2']) + 0.0625*np.trace(d['c4'])
    print(f"    Class {c}: LHS={lhs:.4f}, H={d['H']}, match={abs(lhs - d['H']) < 1e-6}")

# =====================================================================
print()
print("=" * 70)
print("c_k EIGENVALUES BY SCORE SEQUENCE")
print("=" * 70)

# Group by score sequence
score_groups = defaultdict(list)
for c, d in class_polys.items():
    score_groups[d['scores']].append(c)

for scores, classes in sorted(score_groups.items()):
    print(f"\n  Scores {scores}: classes {classes}")
    for c in classes:
        d = class_polys[c]
        print(f"    Class {c} (H={d['H']}): "
              f"c0 eigs={[round(e,3) for e in d['eigs_c0']]}, "
              f"c2 eigs={[round(e,3) for e in d['eigs_c2']]}, "
              f"c4 eigs={[round(e,3) for e in d['eigs_c4']]}")

# =====================================================================
# Check if c_k matrices within a score group are permutation-equivalent
# =====================================================================
print()
print("=" * 70)
print("PERMUTATION EQUIVALENCE WITHIN SCORE GROUPS")
print("=" * 70)

for scores, classes in sorted(score_groups.items()):
    if len(classes) <= 1:
        continue
    print(f"\n  Scores {scores}: classes {classes}")
    for k_name, k_key in [('c0', 'c0'), ('c2', 'c2'), ('c4', 'c4')]:
        c1, c2_idx = classes[0], classes[1]
        mat1 = class_polys[c1][k_key]
        mat2 = class_polys[c2_idx][k_key]

        found = False
        for perm in permutations(range(n)):
            P = np.zeros((n, n))
            for i in range(n):
                P[i][perm[i]] = 1
            rotated = P @ mat2 @ P.T
            if np.allclose(mat1, rotated, atol=1e-6):
                found = True
                break

        print(f"    {k_name}(class {c1}) ~ {k_name}(class {c2_idx}) by perm? {found}")

# =====================================================================
# The exact c_2 formula: c_2 = alpha*I + beta*(J - A - A^T)?
# =====================================================================
print()
print("=" * 70)
print("EXACT c_2 FORMULA SEARCH")
print("=" * 70)

I_mat = np.eye(n)
J_mat = np.ones((n, n))

for c in sorted(class_polys.keys()):
    d = class_polys[c]
    A_np = np.array(d['A'], dtype=float)
    S = A_np - A_np.T  # Skew adjacency
    c2_mat = d['c2']

    # Try: c_2 = a*I + b*J + c*S + d*(A+A^T) + e*S^2
    S2 = S @ S
    AAt = A_np + A_np.T

    basis = [I_mat.flatten(), J_mat.flatten(), S.flatten(), AAt.flatten(), S2.flatten()]
    X = np.column_stack(basis)
    coeffs, residual, _, _ = np.linalg.lstsq(X, c2_mat.flatten(), rcond=None)

    recon = sum(coeffs[i] * np.reshape(basis[i], (n,n)) for i in range(len(basis)))
    err = np.max(np.abs(c2_mat - recon))

    print(f"  Class {c} (H={d['H']}): c2 ≈ {coeffs[0]:.3f}*I + {coeffs[1]:.3f}*J + "
          f"{coeffs[2]:.3f}*S + {coeffs[3]:.3f}*(A+A^T) + {coeffs[4]:.3f}*S^2, err={err:.6f}")

# Try simpler: c_2 = a*I + b*S^2
print("\n  Simpler: c_2 = a*I + b*S^2?")
for c in sorted(class_polys.keys()):
    d = class_polys[c]
    A_np = np.array(d['A'], dtype=float)
    S = A_np - A_np.T
    S2 = S @ S
    c2_mat = d['c2']

    basis = [I_mat.flatten(), S2.flatten()]
    X = np.column_stack(basis)
    coeffs, _, _, _ = np.linalg.lstsq(X, c2_mat.flatten(), rcond=None)

    recon = coeffs[0] * I_mat + coeffs[1] * S2
    err = np.max(np.abs(c2_mat - recon))

    print(f"    Class {c}: c2 ≈ {coeffs[0]:.3f}*I + {coeffs[1]:.3f}*S^2, err={err:.6f}")


print()
print("=" * 70)
print("DONE")
print("=" * 70)
