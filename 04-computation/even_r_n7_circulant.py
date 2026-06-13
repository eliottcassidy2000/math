#!/usr/bin/env python3
"""
Even-r polynomial decomposition at n=7 for circulant tournaments.

At n=7: M(r) = c_0 + c_2*r^2 + c_4*r^4 + c_6*r^6
With c_6 = (n-1)!*I = 720*I (proved at odd n).

Since all circulant tournaments have scalar M = (H/n)*I:
  M(r) = [c_0(r=0)/n + c_2(r=0)/n * r^2 + c_4/n * r^4 + 720/n * r^6] * I

So c_k is scalar for each k (when T is circulant, all c_k are scalar).

Questions:
1. What is tr(c_2) as a function of tournament invariants?
2. Does tr(c_4) have a nice formula?
3. How does H decompose into even-r contributions?

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
import numpy as np

def circulant_tournament(n, gen_set):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in gen_set:
                A[i][j] = 1
    return A

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

def transfer_matrix_entry_r(A, a, b, r_val):
    n = len(A)
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
    return total

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

# =====================================================================
n = 7
print("=" * 70)
print(f"EVEN-r POLYNOMIAL AT n={n}: CIRCULANT TOURNAMENTS")
print("=" * 70)

# All circulant gen sets
half = list(range(1, (n+1)//2))
gen_sets = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.add(frozenset(gs))

# For scalar M: only need M[0,0](r) since M = (M[0,0])*I
# M[0,0](r) = c_0^00 + c_2^00 * r^2 + c_4^00 * r^4 + c_6^00 * r^6
# Need 4 sample points in u = r^2

u_samples = np.array([0.0, 0.04, 0.16, 0.25, 0.36])  # r = 0, 0.2, 0.4, 0.5, 0.6
r_samples = np.sqrt(u_samples)

for gs in sorted(gen_sets):
    A = circulant_tournament(n, gs)
    H = ham_path_count(A)
    t3 = count_3cycles(A)

    # Compute M[0,0](r) at sample points
    vals = [transfer_matrix_entry_r(A, 0, 0, rv) for rv in r_samples]

    # Fit polynomial in u = r^2 of degree 3
    poly_coeffs = np.polyfit(u_samples, vals, 3)
    # poly_coeffs[0]*u^3 + poly_coeffs[1]*u^2 + poly_coeffs[2]*u + poly_coeffs[3]
    c6_diag = poly_coeffs[0]
    c4_diag = poly_coeffs[1]
    c2_diag = poly_coeffs[2]
    c0_diag = poly_coeffs[3]

    # Since M is scalar: tr(c_k) = n * c_k^{00}
    tr_c0 = n * c0_diag
    tr_c2 = n * c2_diag
    tr_c4 = n * c4_diag
    tr_c6 = n * c6_diag

    # Verify reconstruction
    M_check = transfer_matrix_entry_r(A, 0, 0, 0.3)
    M_recon = c0_diag + c2_diag * 0.09 + c4_diag * 0.0081 + c6_diag * 0.000729
    err = abs(M_check - M_recon)

    print(f"\n  gen={sorted(gs)}: H={H}, t3={t3}, err={err:.2e}")
    print(f"    tr(c_0)={tr_c0:.4f}, tr(c_2)={tr_c2:.4f}, tr(c_4)={tr_c4:.4f}, tr(c_6)={tr_c6:.4f}")
    print(f"    c_6 = {c6_diag:.1f}*I (expected 720)")

    # Verify: tr(c_0) + tr(c_2)/4 + tr(c_4)/16 + tr(c_6)/64 = H?
    H_recon = tr_c0 + tr_c2/4 + tr_c4/16 + tr_c6/64
    print(f"    H = tr(c0) + tr(c2)/4 + tr(c4)/16 + tr(c6)/64 = {H_recon:.4f} (expected {H})")


# =====================================================================
print()
print("=" * 70)
print("FORMULAS FOR tr(c_k) AT n=7")
print("=" * 70)

# Collect data
data = []
for gs in sorted(gen_sets):
    A = circulant_tournament(n, gs)
    H = ham_path_count(A)
    t3 = count_3cycles(A)

    vals = [transfer_matrix_entry_r(A, 0, 0, rv) for rv in r_samples]
    poly_coeffs = np.polyfit(u_samples, vals, 3)

    tr_c0 = n * poly_coeffs[3]
    tr_c2 = n * poly_coeffs[2]
    tr_c4 = n * poly_coeffs[1]
    tr_c6 = n * poly_coeffs[0]

    data.append((H, t3, tr_c0, tr_c2, tr_c4, tr_c6))

print("\n  H  | t3 | tr(c0) | tr(c2) | tr(c4) | tr(c6)")
print("  " + "-" * 55)
for H, t3, c0, c2, c4, c6 in sorted(data):
    print(f"  {H:3d} | {t3:2d} | {c0:6.1f} | {c2:7.1f} | {c4:7.1f} | {c6:7.1f}")

# Check tr(c_2) = f(t3)?
print("\n  tr(c_2) vs 3-cycles:")
for H, t3, c0, c2, c4, c6 in sorted(set(data)):
    print(f"    t3={t3}, tr(c2)={c2:.1f}")

# Check tr(c_4) = f(t3)?
print("\n  tr(c_4) vs 3-cycles:")
for H, t3, c0, c2, c4, c6 in sorted(set(data)):
    print(f"    t3={t3}, tr(c4)={c4:.1f}")

# Count 5-cycles
def count_5cycles(A):
    n = len(A)
    count = 0
    for perm in permutations(range(n), 5):
        if (A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and
            A[perm[2]][perm[3]] and A[perm[3]][perm[4]] and A[perm[4]][perm[0]]):
            count += 1
    return count // 5  # each 5-cycle counted 5 times

print("\n  tr(c_4) vs 3-cycles and 5-cycles:")
for gs in sorted(gen_sets):
    A = circulant_tournament(n, gs)
    H = ham_path_count(A)
    t3 = count_3cycles(A)
    t5 = count_5cycles(A)

    vals = [transfer_matrix_entry_r(A, 0, 0, rv) for rv in r_samples]
    poly_coeffs = np.polyfit(u_samples, vals, 3)
    tr_c4 = n * poly_coeffs[1]

    print(f"    gen={sorted(gs)}: t3={t3}, t5={t5}, tr(c4)={tr_c4:.1f}")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
