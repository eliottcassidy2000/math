#!/usr/bin/env python3
"""
Investigating the shared c_2 spectrum among isomorphism classes.

Classes 8 (H=11), 9 (H=15), 10 (H=13) at n=5 all share identical c_2 eigenvalues.
What is c_2 algebraically? Can we understand why some classes share it?

c_2 is extracted from M(r) = c_0 + c_2*r^2 via:
  c_2 = (M(1/2) - M(0)) / (1/4)

At the algebraic level, M(r) is the IE transfer matrix where each edge
weight is t_{ij} = r + s_{ij} with s_{ij} = A[i][j] - 1/2 in {-1/2, +1/2}.

Each M[a,b] = sum over IE-decompositions of products of (r + s_e) terms.
Expanding: (r + s_e1)(r + s_e2)...(r + s_ek) = sum_{j=0}^k C(k,j) r^j * (product of k-j s-values)

The r^2 coefficient picks terms where exactly 2 of the edge-weights contribute r,
and the rest contribute their s-values.

For M[a,b] at n=5: each term in IE involves paths of various lengths.
The r^2 term comes from paths where exactly 2 edges contribute r.

QUESTION: Is c_2 determined by some simple combinatorial invariant
(score sequence? number of 3-cycles? etc.)?

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

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
    count = 0
    for p in permutations(range(n)):
        if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)):
            count += 1
    return count

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
        for j in range(n):
            for k in range(n):
                if i != j and j != k and i != k:
                    if A[i][j] == 1 and A[j][k] == 1 and A[k][i] == 1:
                        count += 1
    return count // 3  # each 3-cycle counted 3 times

# =====================================================================
n = 5
_, tiles = tiling_to_tournament(0, n)
m = len(tiles)

print("=" * 70)
print("c_2 SPECTRUM SHARING AT n=5")
print("=" * 70)

# Collect data by class
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

# Compute c_0, c_2 for representative of each class
class_data = {}
for canon in sorted(iso_classes.keys()):
    c = class_labels[canon]
    bits = iso_classes[canon][0]
    A = tiling_data[bits]['A']
    H = tiling_data[bits]['H']

    M_0 = transfer_matrix_r(A, 0.0)
    M_half = transfer_matrix_r(A, 0.5)

    c0 = M_0
    c2 = (M_half - M_0) / 0.25

    scores = score_sequence(A)
    cyc3 = count_3cycles(A)

    eigs_c2 = sorted(np.linalg.eigvalsh(c2))

    class_data[c] = {
        'H': H, 'c0': c0, 'c2': c2, 'scores': scores,
        'cyc3': cyc3, 'eigs_c2': eigs_c2, 'tr_c2': np.trace(c2),
        'bits': bits, 'A': A
    }

    print(f"  Class {c:2d} (H={H:2d}): scores={scores}, 3-cycles={cyc3}, "
          f"tr(c2)={np.trace(c2):.0f}, c2 eigs={[round(e,2) for e in eigs_c2]}")

# Group classes by c_2 spectrum
print("\nGrouping by c_2 eigenvalue spectrum:")
spectrum_groups = defaultdict(list)
for c, d in class_data.items():
    key = tuple(round(e, 2) for e in d['eigs_c2'])
    spectrum_groups[key].append(c)

for spec, classes in sorted(spectrum_groups.items()):
    h_vals = [class_data[c]['H'] for c in classes]
    print(f"  spectrum {spec}: classes {classes}, H values {h_vals}")

# =====================================================================
# Is c_2 determined by score sequence?
# =====================================================================
print()
print("=" * 70)
print("IS c_2 DETERMINED BY SCORE SEQUENCE?")
print("=" * 70)

score_groups = defaultdict(list)
for c, d in class_data.items():
    score_groups[d['scores']].append(c)

for scores, classes in sorted(score_groups.items()):
    if len(classes) > 1:
        c2_match = np.allclose(class_data[classes[0]]['c2'], class_data[classes[1]]['c2'])
        print(f"  scores={scores}: classes {classes}, c2 match? {c2_match}")
        for c in classes:
            print(f"    class {c}: tr(c2)={class_data[c]['tr_c2']:.0f}, "
                  f"c2 eigs={[round(e,2) for e in class_data[c]['eigs_c2']]}")
    else:
        print(f"  scores={scores}: class {classes[0]} (unique)")

# =====================================================================
# tr(c_2) = 4 * H? Check relationship
# =====================================================================
print()
print("=" * 70)
print("tr(c_2) vs H RELATIONSHIP")
print("=" * 70)

for c, d in sorted(class_data.items()):
    ratio = d['tr_c2'] / d['H'] if d['H'] != 0 else float('inf')
    tr_c0 = np.trace(d['c0'])
    print(f"  Class {c}: H={d['H']}, tr(c0)={tr_c0:.2f}, tr(c2)={d['tr_c2']:.0f}, "
          f"tr(c2)/H={ratio:.4f}, tr(c0)+0.25*tr(c2)={tr_c0 + 0.25*d['tr_c2']:.2f}")
    # At r=1/2: tr(M(1/2)) = H at odd n
    # tr(M(r)) = tr(c0) + tr(c2)*r^2
    # tr(M(1/2)) = tr(c0) + tr(c2)/4 = H
    # So tr(c2) = 4*(H - tr(c0))

# =====================================================================
# The actual c_2 matrices for the shared-spectrum classes
# =====================================================================
print()
print("=" * 70)
print("c_2 MATRICES FOR CLASSES WITH SHARED SPECTRUM")
print("=" * 70)

for c in [8, 9, 10]:
    d = class_data[c]
    print(f"\n  Class {c} (H={d['H']}):")
    print(f"    c_2 = {np.round(d['c2'], 2).tolist()}")

# Check if they are related by permutation
print("\nAre the c_2 matrices related by vertex permutation?")
for c1, c2_idx in [(8, 9), (8, 10), (9, 10)]:
    c2_1 = class_data[c1]['c2']
    c2_2 = class_data[c2_idx]['c2']

    found_perm = False
    for perm in permutations(range(n)):
        P = np.zeros((n, n))
        for i in range(n):
            P[i][perm[i]] = 1
        rotated = P @ c2_2 @ P.T
        if np.allclose(c2_1, rotated, atol=1e-6):
            found_perm = True
            break

    print(f"  c2({c1}) ~ c2({c2_idx}) by permutation? {found_perm}")

# =====================================================================
# What about the characteristic polynomial of c_2?
# =====================================================================
print()
print("=" * 70)
print("CHARACTERISTIC POLYNOMIAL OF c_2")
print("=" * 70)

for c, d in sorted(class_data.items()):
    charpoly = np.round(np.poly(d['c2']), 4)
    print(f"  Class {c}: charpoly(c2) = {charpoly.tolist()}")


# =====================================================================
# Is there a formula c_2 = f(J, D_s, ...) where J is all-ones, D_s is score diagonal?
# =====================================================================
print()
print("=" * 70)
print("DECOMPOSING c_2 IN TERMS OF J, I, AND TOURNAMENT STRUCTURE")
print("=" * 70)

J = np.ones((n, n))
I = np.eye(n)

for c, d in sorted(class_data.items()):
    A = d['A']
    A_np = np.array(A, dtype=float)
    # Score matrix
    S = A_np - A_np.T  # skew-adjacency
    scores_arr = np.sum(A_np, axis=1)
    D_s = np.diag(scores_arr)

    c2 = d['c2']
    tr_c2 = np.trace(c2)

    # Try c_2 = alpha*I + beta*J + gamma*A + delta*A^T
    # Solve using least squares
    basis = [I.flatten(), J.flatten(), A_np.flatten(), A_np.T.flatten()]
    X = np.column_stack(basis)
    coeffs, residual, _, _ = np.linalg.lstsq(X, c2.flatten(), rcond=None)

    recon = coeffs[0]*I + coeffs[1]*J + coeffs[2]*A_np + coeffs[3]*A_np.T
    err = np.max(np.abs(c2 - recon))

    print(f"  Class {c} (H={d['H']}): c2 ≈ {coeffs[0]:.2f}*I + {coeffs[1]:.2f}*J + "
          f"{coeffs[2]:.2f}*A + {coeffs[3]:.2f}*A^T, max_err={err:.4f}")


# =====================================================================
# tr(c_0) pattern
# =====================================================================
print()
print("=" * 70)
print("tr(c_0) PATTERN")
print("=" * 70)

for c, d in sorted(class_data.items()):
    tr_c0 = np.trace(d['c0'])
    cyc3 = d['cyc3']
    print(f"  Class {c}: H={d['H']}, 3-cycles={cyc3}, tr(c0)={tr_c0:.4f}")

# At n=5: tr(c0) = H(r=0). What is H(r=0)?
# H(0) = sum_P prod_e s_e = sum_P prod_e (A[e]-1/2) = sum_P (-1)^{#non-edges in P} / 2^4
# This is the "signed" Ham path count divided by 16.
print("\nH(r=0) directly computed:")
for c, d in sorted(class_data.items()):
    A = d['A']
    total = 0.0
    for p in permutations(range(n)):
        w = 1.0
        for i in range(n-1):
            w *= (A[p[i]][p[i+1]] - 0.5)
        total += w
    tr_c0 = np.trace(d['c0'])
    print(f"  Class {c}: H(0)={total:.6f}, tr(c0)/n={tr_c0/n:.6f}, "
          f"n*H(0)={n*total:.6f}, tr(c0)={tr_c0:.6f}")
    # tr(c0) should equal n * M[0,0](0) if scalar, otherwise it's the trace
    # At odd n: tr(M(r)) = sum_a M[a,a](r). And H(r) = ?
    # Actually H(r) = total weighted ham paths, not tr(M(r)).
    # tr(M(r)) = H(r) at odd n? Let me check.
    # H(r) = sum_P w(P). tr(M(r)) = sum_a M[a,a](r).
    # M[a,a](r) = sum_P (-1)^{pos(a,P)} w(P).
    # So tr(M(r)) = sum_P w(P) sum_a (-1)^{pos(a,P)} = sum_P w(P) * 0 = 0 (at even n)
    # At odd n: sum_a (-1)^{pos(a,P)} = 1 for all P.
    # Wait, that's only at n odd. sum_{a} (-1)^{pos(a)} where pos ranges over 0..n-1.
    # Each a appears at some unique position, so it's sum_{j=0}^{n-1} (-1)^j.
    # At odd n: this sum = 1. At even n: this sum = 0.
    # So tr(M(r)) = H(r) at odd n. YES.

print("\n  => At odd n=5: tr(M(r)) = H(r)")
print("  => tr(c0) = H(0) and tr(c2) = 4*(H(1/2) - H(0))")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
