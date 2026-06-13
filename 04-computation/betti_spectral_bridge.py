#!/usr/bin/env python3
"""
betti_spectral_bridge.py — opus-2026-03-13-S71

Investigate connections between:
1. Spectral data (Q_k, F_p, amplification) from S67e/f
2. GLMY path homology (Betti numbers, Omega profiles) from S70

Key questions:
- Does chi/n = ±1 always hold for circulant tournaments? What determines the sign?
- Is there a formula for total Betti sum in terms of spectral data?
- How do boundary map ranks relate to Q_k values?
- Does the Fibonacci product identity (prod(1+Q_k) = F_p) connect to Betti?
"""

import numpy as np
from itertools import permutations

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def get_regular_paths(A, m):
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def boundary_matrix(paths_m, paths_m_minus_1):
    if not paths_m or not paths_m_minus_1:
        return np.zeros((len(paths_m_minus_1) if paths_m_minus_1 else 0,
                         len(paths_m) if paths_m else 0), dtype=int)
    path_to_idx = {p: i for i, p in enumerate(paths_m_minus_1)}
    m = len(paths_m[0]) - 1
    matrix = np.zeros((len(paths_m_minus_1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(1, m):
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in path_to_idx:
                matrix[path_to_idx[face], j] += sign
    return matrix

def compute_betti(A, max_dim=None):
    n = A.shape[0]
    if max_dim is None:
        max_dim = n - 1
    all_paths = {}
    for m in range(max_dim + 1):
        all_paths[m] = get_regular_paths(A, m)
    ranks = {}
    for m in range(1, max_dim + 1):
        if all_paths[m] and all_paths[m-1]:
            B = boundary_matrix(all_paths[m], all_paths[m-1])
            ranks[m] = np.linalg.matrix_rank(B)
        else:
            ranks[m] = 0
    betti = []
    for m in range(max_dim + 1):
        omega_m = len(all_paths[m])
        rank_dm = ranks.get(m, 0)
        rank_dm_plus_1 = ranks.get(m+1, 0)
        beta_m = omega_m - rank_dm - rank_dm_plus_1
        betti.append(beta_m)
    omegas = [len(all_paths[m]) for m in range(max_dim + 1)]
    return betti, omegas, ranks

def spectral_data(n, S):
    """Compute Q_k = |S_hat(k)|^2 for all k."""
    omega = np.exp(2j * np.pi / n)
    Q = []
    for k in range(n):
        S_hat = sum(omega**(k*s) for s in S)
        Q.append(abs(S_hat)**2)
    return Q

# ============================================================
print("="*70)
print("BETTI-SPECTRAL BRIDGE")
print("="*70)

# Test cases: interval and Paley tournaments
test_cases = []

for n in [5, 7, 9]:
    m = (n-1)//2
    S = set(range(1, m+1))
    test_cases.append((f"Interval n={n}", n, S))

for p in [7, 11]:
    if p % 4 != 3: continue
    QR = {a % p for a in range(1, p) if pow(a, (p-1)//2, p) == 1}
    test_cases.append((f"Paley p={p}", p, QR))

# Also test C_7^{1,3,5}
test_cases.append(("C7_135", 7, {1, 3, 5}))

for name, n, S in test_cases:
    print(f"\n{'='*50}")
    print(f"  {name}  (n={n}, S={sorted(S)})")
    print(f"{'='*50}")

    A = circulant_tournament(n, S)
    Q = spectral_data(n, S)
    betti, omegas, ranks = compute_betti(A)

    print(f"  Q_k: {[f'{q:.4f}' for q in Q]}")
    print(f"  prod(1+Q_k, k≠0): {np.prod([1+Q[k] for k in range(1,n)]):.4f}")

    chi = sum((-1)**m * betti[m] for m in range(len(betti)))
    total_beta = sum(betti)

    print(f"  Omega:  {omegas}")
    print(f"  Betti:  {betti}")
    print(f"  chi = {chi}, chi/n = {chi/n}")
    print(f"  total_beta = {total_beta}, total_beta/n = {total_beta/n}")

    # Spectral invariants
    mean_Q = np.mean(Q[1:])
    var_Q = np.var(Q[1:])
    prod_Q = np.prod([Q[k] for k in range(1, n)])
    sum_Q = sum(Q[k] for k in range(1, n))

    print(f"  mean(Q_k, k≠0) = {mean_Q:.4f}")
    print(f"  var(Q_k, k≠0) = {var_Q:.4f}")
    print(f"  sum(Q_k, k≠0) = {sum_Q:.4f}")
    print(f"  prod(Q_k, k≠0) = {prod_Q:.4f}")

    # Key ratios
    if n > 3:
        print(f"  total_Omega = {sum(omegas)}")
        print(f"  total_Omega/n = {sum(omegas)/n:.4f}")

    # Boundary map ranks
    print(f"  ranks(∂_m): {[ranks.get(m,0) for m in range(1, n)]}")

    # Per-eigenspace (divide by n)
    per_eig_betti = [b//n for b in betti]
    per_eig_omega = [o//n for o in omegas]
    print(f"  per_eig_betti: {per_eig_betti}")
    print(f"  per_eig_omega: {per_eig_omega}")

    # Poincare polynomial
    print(f"  Poincare poly P(t) = {' + '.join(f'{b}t^{m}' for m, b in enumerate(betti) if b > 0)}")

# ============================================================
print(f"\n{'='*70}")
print("CHI ANALYSIS: SIGN PATTERN")
print("="*70)

# What determines sign(chi)?
print("\n  For circulant tournaments on Z_n:")
print("  chi > 0 means more even-dim homology; chi < 0 means more odd-dim")
print()

# For interval tournaments:
for n in [3, 5, 7, 9]:
    m = (n-1)//2
    S = set(range(1, m+1))
    A = circulant_tournament(n, S)
    betti, omegas, _ = compute_betti(A)
    chi = sum((-1)**m * betti[m] for m in range(len(betti)))
    print(f"  Interval n={n}: chi={chi:+d}, chi/n={chi//n:+d}, (n-1)/2={m}")

print()
for p in [3, 7]:
    QR = {a % p for a in range(1, p) if pow(a, (p-1)//2, p) == 1}
    A = circulant_tournament(p, QR)
    betti, omegas, _ = compute_betti(A)
    chi = sum((-1)**m * betti[m] for m in range(len(betti)))
    print(f"  Paley p={p}: chi={chi:+d}, chi/p={chi//p:+d}")

# ============================================================
print(f"\n{'='*70}")
print("FIBONACCI CONNECTION")
print("="*70)

# For Interval tournaments, prod(1+Q_k) = F_p (Fibonacci)
# Is total_Omega related to F_p?

def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a

for n in [5, 7, 9]:
    m = (n-1)//2
    S = set(range(1, m+1))
    A = circulant_tournament(n, S)
    Q = spectral_data(n, S)
    betti, omegas, _ = compute_betti(A)

    F_n = fib(n)
    prod_1_Q = np.prod([1+Q[k] for k in range(1, n)])
    total_omega = sum(omegas)
    total_beta = sum(betti)

    print(f"\n  Interval n={n}:")
    print(f"    F_{n} = {F_n}")
    print(f"    prod(1+Q_k) = {prod_1_Q:.1f}")
    print(f"    total_Omega = {total_omega}")
    print(f"    total_Omega / n = {total_omega/n:.4f}")
    print(f"    total_beta = {total_beta}")
    print(f"    total_beta / n = {total_beta/n:.4f}")
    print(f"    F_n / total_Omega_per_eig = {F_n / (total_omega/n):.6f}")

    # Alternating sum of Omega (= chi via ranks)
    alt_omega = sum((-1)**m * omegas[m] for m in range(len(omegas)))
    print(f"    alternating sum of Omega = {alt_omega}")
    print(f"    ... / n = {alt_omega / n:.4f}")

# ============================================================
print(f"\n{'='*70}")
print("BETTI CONCENTRATION: WHERE IS HOMOLOGY CONCENTRATED?")
print("="*70)

# For Paley p=7: all Betti are in dims 0 and 3-6
# For Interval n=7: Betti in dims 0,1,3,4,5
# Is there a pattern related to Q_k?

for name, n, S in test_cases:
    A = circulant_tournament(n, S)
    betti, omegas, _ = compute_betti(A)
    nonzero_dims = [m for m, b in enumerate(betti) if b > 0]
    zero_dims = [m for m, b in enumerate(betti) if b == 0]
    print(f"  {name}: nonzero β at dims {nonzero_dims}, zero at {zero_dims}")

# ============================================================
print(f"\n{'='*70}")
print("HOMOLOGICAL DIMENSION vs SPECTRAL FLATNESS")
print("="*70)

# Is the top nonzero Betti dimension related to spectral properties?
for name, n, S in test_cases:
    A = circulant_tournament(n, S)
    Q = spectral_data(n, S)
    betti, _, _ = compute_betti(A)
    top_dim = max(m for m, b in enumerate(betti) if b > 0)
    flatness = 1 - np.std(Q[1:]) / np.mean(Q[1:])  # 1=flat, 0=not
    print(f"  {name}: top_dim={top_dim}/{n-1}, flatness={flatness:.4f}, "
          f"m/(n-1)={n-1}")

print("\nDONE.")
