#!/usr/bin/env python3
"""
counterexample_anatomy.py — opus-2026-03-06-S1

Comprehensive spectral and structural anatomy of the THM-025 counterexample:
the n=9 tournament T where I(Omega(T), x) = 1 + 94x + 10x^2 + x^3 is NOT real-rooted.

Computes:
1. Eigenvalues of A (adjacency matrix)
2. Eigenvalues of S = A - A^T (skew-symmetric)
3. tr(A^k) for k=1..9
4. Irving-Omar generating function: sum_{k odd} tr(A^k)/k
5. Signed adjacency B = 2A - J eigenvalues (J = all-ones minus identity)
6. AA^T eigenvalues (correlation with H via kind-pasteur)
7. Comparison with 5 random real-rooted n=9 tournaments
8. Cayley transform C = (I+A)(I-A)^{-1} eigenvalues
9. Pfaffian of S: det(S) = Pf(S)^2
"""

import numpy as np
from numpy.linalg import eig, eigvals, det, inv
from itertools import combinations
import random

np.set_printoptions(precision=6, suppress=True, linewidth=120)

# ============================================================
# The counterexample tournament (THM-025)
# ============================================================
T_adj = np.array([
    [0, 1, 0, 1, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 1, 1, 1, 1, 0],
    [0, 0, 1, 0, 0, 1, 0, 1, 0],
    [1, 1, 0, 1, 0, 0, 0, 1, 0],
    [1, 1, 0, 0, 1, 0, 1, 1, 1],
    [0, 1, 0, 1, 1, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 1, 1, 0],
], dtype=float)

n = 9

def score_sequence(A):
    return sorted(A.sum(axis=1).astype(int))

def find_directed_odd_cycles(A, n):
    """Find all directed odd cycles (length 3,5,7,9) up to rotation normalization."""
    cycles = set()
    N = A.shape[0]

    def dfs(path, start):
        cur = path[-1]
        length = len(path)
        if length > 1 and A[cur, start] == 1 and length % 2 == 1:
            # Normalize by rotation: start from minimum vertex
            min_idx = path.index(min(path))
            normalized = tuple(path[min_idx:] + path[:min_idx])
            cycles.add(normalized)
        if length < N:
            for nxt in range(N):
                if nxt != start and nxt not in path and A[cur, nxt] == 1:
                    dfs(path + [nxt], start)

    for v in range(N):
        dfs([v], v)

    return cycles

def build_omega_graph(cycles):
    """Build the conflict graph Omega: vertices=odd directed cycles, edges=shared vertices."""
    cycle_list = list(cycles)
    vertex_sets = [set(c) for c in cycle_list]
    m = len(cycle_list)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if vertex_sets[i] & vertex_sets[j]:
                adj[i][j] = True
                adj[j][i] = True
    return cycle_list, adj

def independence_polynomial(adj, m):
    """Compute independence polynomial by inclusion over all independent sets."""
    # For moderate m, enumerate all independent sets via backtracking
    coeffs = [0] * (m + 1)

    def backtrack(idx, current_set, excluded):
        coeffs[len(current_set)] += 1
        for i in range(idx, m):
            if i not in excluded:
                new_excluded = excluded | {j for j in range(m) if adj[i][j]}
                new_excluded.add(i)
                backtrack(i + 1, current_set | {i}, new_excluded)

    backtrack(0, set(), set())
    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def is_real_rooted(coeffs):
    """Check if polynomial with given coefficients is real-rooted."""
    if len(coeffs) <= 2:
        return True
    # numpy poly uses highest degree first
    poly = list(reversed(coeffs))
    roots = np.roots(poly)
    return all(abs(r.imag) < 1e-8 for r in roots)

def compute_H(A):
    """Compute H(T) = number of Hamiltonian paths in tournament with adjacency matrix A."""
    N = A.shape[0]
    # DP over bitmask
    dp = {}
    for v in range(N):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, N + 1):
        for mask in range(1 << N):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(N):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                count = 0
                for u in range(N):
                    if (prev_mask & (1 << u)) and A[u, v] == 1:
                        count += dp.get((prev_mask, u), 0)
                if count > 0:
                    dp[(mask, v)] = count
    full_mask = (1 << N) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(N))

def random_tournament(n):
    """Generate a random n-tournament."""
    A = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A

# ============================================================
print("=" * 70)
print("COUNTEREXAMPLE ANATOMY: THM-025 (n=9, non-real-rooted I(Omega,x))")
print("=" * 70)

scores = score_sequence(T_adj)
print(f"\nScore sequence: {scores}")
H_val = compute_H(T_adj)
print(f"H(T) = {H_val}")

# ============================================================
# 1. Eigenvalues of A
# ============================================================
print("\n" + "=" * 70)
print("1. EIGENVALUES OF A (adjacency matrix)")
print("=" * 70)
evals_A = eigvals(T_adj)
evals_A_sorted = sorted(evals_A, key=lambda x: -abs(x))
print("Eigenvalues (sorted by |lambda|):")
for i, e in enumerate(evals_A_sorted):
    if abs(e.imag) < 1e-10:
        print(f"  lambda_{i+1} = {e.real:.6f}")
    else:
        print(f"  lambda_{i+1} = {e.real:.6f} + {e.imag:.6f}i  (|lambda| = {abs(e):.6f})")
print(f"\nSpectral radius rho(A) = {max(abs(e) for e in evals_A):.6f}")

# ============================================================
# 2. Eigenvalues of S = A - A^T
# ============================================================
print("\n" + "=" * 70)
print("2. EIGENVALUES OF S = A - A^T (skew-symmetric)")
print("=" * 70)
S = T_adj - T_adj.T
evals_S = eigvals(S)
evals_S_sorted = sorted(evals_S, key=lambda x: -abs(x))
print("Eigenvalues (should be pure imaginary or zero):")
for i, e in enumerate(evals_S_sorted):
    print(f"  mu_{i+1} = {e.imag:.6f}i  (real part: {e.real:.2e})")

# ============================================================
# 3. tr(A^k) for k=1..9
# ============================================================
print("\n" + "=" * 70)
print("3. TRACE OF A^k FOR k=1,...,9")
print("=" * 70)
Ak = np.eye(n)
traces = []
for k in range(1, 10):
    Ak = Ak @ T_adj
    tr = np.trace(Ak)
    traces.append(tr)
    print(f"  tr(A^{k}) = {tr:.0f}")

print("\nInterpretation: tr(A^k) = number of closed directed walks of length k")
print(f"  tr(A^3) = {traces[2]:.0f}  (3*[directed 3-cycles])")

# ============================================================
# 4. Irving-Omar generating function
# ============================================================
print("\n" + "=" * 70)
print("4. IRVING-OMAR: sum_{k odd} tr(A^k)/k for k=1,3,5,7,9")
print("=" * 70)
io_sum = 0
for k in [1, 3, 5, 7, 9]:
    val = traces[k-1] / k
    io_sum += val
    print(f"  tr(A^{k})/{k} = {traces[k-1]:.0f}/{k} = {val:.4f}")
print(f"\n  Total sum = {io_sum:.4f}")
print(f"  This relates to sum of arctanh(eigenvalues) of A")

# Also compute via eigenvalues
print("\n  Verification via eigenvalues: sum_i sum_{k odd} lambda_i^k / k")
eig_sum = 0
for e in evals_A:
    for k in [1, 3, 5, 7, 9]:
        eig_sum += e**k / k
print(f"  Eigenvalue-based sum = {eig_sum.real:.4f} + {eig_sum.imag:.2e}i")

# ============================================================
# 5. Signed adjacency B = 2A - J
# ============================================================
print("\n" + "=" * 70)
print("5. SIGNED ADJACENCY B = 2A - J (J = all-ones minus identity)")
print("=" * 70)
J = np.ones((n, n)) - np.eye(n)
B = 2 * T_adj - J
evals_B = eigvals(B)
evals_B_sorted = sorted(evals_B, key=lambda x: -abs(x))
print("Eigenvalues of B:")
for i, e in enumerate(evals_B_sorted):
    if abs(e.imag) < 1e-10:
        print(f"  beta_{i+1} = {e.real:.6f}")
    else:
        print(f"  beta_{i+1} = {e.real:.6f} + {e.imag:.6f}i  (|beta| = {abs(e):.6f})")
print(f"\nNote: B is skew-symmetric (B = S = A - A^T when J = ones-I)")
print(f"  Actually B = 2A - (J+I) + I = 2A - ones + 2I? Let's verify:")
print(f"  B is skew-symmetric: {np.allclose(B, -B.T)}")
print(f"  B = S: {np.allclose(B, S)}")

# ============================================================
# 6. AA^T eigenvalues
# ============================================================
print("\n" + "=" * 70)
print("6. EIGENVALUES OF AA^T")
print("=" * 70)
AAT = T_adj @ T_adj.T
evals_AAT = np.sort(np.real(eigvals(AAT)))[::-1]
print("Eigenvalues of AA^T (sorted descending):")
for i, e in enumerate(evals_AAT):
    print(f"  sigma_{i+1}^2 = {e:.6f}")
print(f"\n  lambda_1(AA^T) = {evals_AAT[0]:.6f}")
print(f"  Singular values of A: {np.sqrt(np.maximum(evals_AAT, 0))}")
print(f"\n  kind-pasteur found corr(H, lambda_1(AA^T)) = -0.97 for H-maximizers")
print(f"  H(T) = {H_val}, lambda_1(AA^T) = {evals_AAT[0]:.4f}")

# ============================================================
# 7. Comparison with real-rooted tournaments
# ============================================================
print("\n" + "=" * 70)
print("7. COMPARISON WITH RANDOM REAL-ROOTED n=9 TOURNAMENTS")
print("=" * 70)

random.seed(42)
np.random.seed(42)

real_rooted_data = []
non_real_rooted_data = []
attempts = 0
total_computable = 0

print("Sampling random n=9 tournaments (classifying real-rooted vs not)...")
while (len(real_rooted_data) < 5 or len(non_real_rooted_data) < 5) and attempts < 2000:
    attempts += 1
    A_rand = random_tournament(n)

    # Find odd cycles
    cycles = find_directed_odd_cycles(A_rand, n)
    if len(cycles) == 0:
        continue

    cycle_list, omega_adj = build_omega_graph(cycles)
    m = len(cycle_list)

    if m > 50:  # Skip if too many cycles
        continue

    total_computable += 1
    ip_coeffs = independence_polynomial(omega_adj, m)
    rr = is_real_rooted(ip_coeffs)

    # Compute spectral data
    H_rand = compute_H(A_rand)
    evals_rand = eigvals(A_rand)
    S_rand = A_rand - A_rand.T
    evals_S_rand = eigvals(S_rand)
    AAT_rand = A_rand @ A_rand.T
    evals_AAT_rand = np.sort(np.real(eigvals(AAT_rand)))[::-1]

    traces_rand = []
    Ak_rand = np.eye(n)
    for k in range(1, 10):
        Ak_rand = Ak_rand @ A_rand
        traces_rand.append(np.trace(Ak_rand))

    entry = {
        'scores': score_sequence(A_rand),
        'H': H_rand,
        'num_cycles': m,
        'ip_coeffs': ip_coeffs,
        'real_rooted': rr,
        'rho_A': max(abs(e) for e in evals_rand),
        'evals_A': sorted(evals_rand, key=lambda x: -abs(x)),
        'max_imag_S': max(abs(e.imag) for e in evals_S_rand),
        'lambda1_AAT': evals_AAT_rand[0],
        'evals_AAT': evals_AAT_rand,
        'traces': traces_rand,
    }

    if rr and len(real_rooted_data) < 5:
        real_rooted_data.append(entry)
    elif not rr and len(non_real_rooted_data) < 5:
        non_real_rooted_data.append(entry)

comparison_data = real_rooted_data  # for backward compat below

print(f"\nSampled {total_computable} computable tournaments in {attempts} attempts")
print(f"Found {len(real_rooted_data)} real-rooted, {len(non_real_rooted_data)} non-real-rooted\n")

# Print counterexample data for comparison
cycles_ce = find_directed_odd_cycles(T_adj, n)
cycle_list_ce, omega_adj_ce = build_omega_graph(cycles_ce)
m_ce = len(cycle_list_ce)
ip_coeffs_ce = independence_polynomial(omega_adj_ce, m_ce)

print(f"COUNTEREXAMPLE:")
print(f"  Scores: {scores}")
print(f"  H(T) = {H_val}")
print(f"  |Omega| = {m_ce} odd cycles")
print(f"  I(Omega,x) = {' + '.join(f'{c}x^{i}' if i > 0 else str(c) for i, c in enumerate(ip_coeffs_ce) if c != 0)}")
print(f"  Real-rooted: {is_real_rooted(ip_coeffs_ce)}")
print(f"  rho(A) = {max(abs(e) for e in evals_A):.6f}")
print(f"  max |Im(mu)| (S) = {max(abs(e.imag) for e in evals_S_sorted):.6f}")
print(f"  lambda_1(AA^T) = {evals_AAT[0]:.6f}")
print(f"  tr(A^3) = {traces[2]:.0f}")

def print_tournament_summary(label, d):
    print(f"\n{label}:")
    print(f"  Scores: {d['scores']}")
    print(f"  H(T) = {d['H']}")
    print(f"  |Omega| = {d['num_cycles']} odd cycles")
    ip_str = ' + '.join(f'{c}x^{i}' if i > 0 else str(c) for i, c in enumerate(d['ip_coeffs']) if c != 0)
    print(f"  I(Omega,x) = {ip_str}")
    print(f"  rho(A) = {d['rho_A']:.6f}")
    print(f"  max |Im(mu)| (S) = {d['max_imag_S']:.6f}")
    print(f"  lambda_1(AA^T) = {d['lambda1_AAT']:.6f}")
    print(f"  tr(A^3) = {d['traces'][2]:.0f}")

for idx, d in enumerate(real_rooted_data):
    print_tournament_summary(f"REAL-ROOTED #{idx+1}", d)

for idx, d in enumerate(non_real_rooted_data):
    print_tournament_summary(f"NON-REAL-ROOTED #{idx+1}", d)

# Summary statistics
def print_summary_block(label, data_list):
    rhos = [d['rho_A'] for d in data_list]
    imags = [d['max_imag_S'] for d in data_list]
    aats = [d['lambda1_AAT'] for d in data_list]
    hs = [d['H'] for d in data_list]
    tr3s = [d['traces'][2] for d in data_list]
    ncycs = [d['num_cycles'] for d in data_list]
    return {
        'rho_A': (np.mean(rhos), min(rhos), max(rhos)),
        'max_imag_S': (np.mean(imags), min(imags), max(imags)),
        'lambda1_AAT': (np.mean(aats), min(aats), max(aats)),
        'H': (np.mean(hs), min(hs), max(hs)),
        'tr3': (np.mean(tr3s), min(tr3s), max(tr3s)),
        'ncyc': (np.mean(ncycs), min(ncycs), max(ncycs)),
    }

ce_rho = max(abs(e) for e in evals_A)
ce_imag = max(abs(e.imag) for e in evals_S_sorted)

print(f"\n--- SUMMARY COMPARISON ---")
header = f"{'Metric':<22} {'Counterex.':>12}"
if real_rooted_data:
    header += f" {'RR avg':>10} {'RR range':>22}"
if non_real_rooted_data:
    header += f" {'NRR avg':>10} {'NRR range':>22}"
print(header)
print("-" * len(header))

metrics = [
    ('rho(A)', ce_rho, 'rho_A'),
    ('max|Im(mu_S)|', ce_imag, 'max_imag_S'),
    ('lambda_1(AAT)', evals_AAT[0], 'lambda1_AAT'),
    ('H(T)', H_val, 'H'),
    ('tr(A^3)', traces[2], 'tr3'),
    ('|Omega|', m_ce, 'ncyc'),
]

rr_stats = print_summary_block("RR", real_rooted_data) if real_rooted_data else None
nrr_stats = print_summary_block("NRR", non_real_rooted_data) if non_real_rooted_data else None

for name, ce_val, key in metrics:
    row = f"{name:<22} {ce_val:>12.4f}" if isinstance(ce_val, float) else f"{name:<22} {ce_val:>12}"
    if rr_stats:
        avg, lo, hi = rr_stats[key]
        row += f" {avg:>10.2f} {f'[{lo:.2f}, {hi:.2f}]':>22}"
    if nrr_stats:
        avg, lo, hi = nrr_stats[key]
        row += f" {avg:>10.2f} {f'[{lo:.2f}, {hi:.2f}]':>22}"
    print(row)

# ============================================================
# 8. Cayley transform C = (I+A)(I-A)^{-1}
# ============================================================
print("\n" + "=" * 70)
print("8. CAYLEY TRANSFORM C = (I + A)(I - A)^{-1}")
print("=" * 70)
I_mat = np.eye(n)
I_minus_A = I_mat - T_adj
det_IminusA = det(I_minus_A)
print(f"det(I - A) = {det_IminusA:.6f}")

if abs(det_IminusA) > 1e-10:
    C = (I_mat + T_adj) @ inv(I_minus_A)
    evals_C = eigvals(C)
    evals_C_sorted = sorted(evals_C, key=lambda x: -abs(x))
    print("I - A is invertible. Cayley transform eigenvalues:")
    for i, e in enumerate(evals_C_sorted):
        if abs(e.imag) < 1e-10:
            print(f"  gamma_{i+1} = {e.real:.6f}")
        else:
            print(f"  gamma_{i+1} = {e.real:.6f} + {e.imag:.6f}i  (|gamma| = {abs(e):.6f})")
    print(f"\n  Spectral radius rho(C) = {max(abs(e) for e in evals_C):.6f}")
else:
    print("I - A is SINGULAR. Cayley transform does not exist.")

# ============================================================
# 9. Pfaffian of S
# ============================================================
print("\n" + "=" * 70)
print("9. PFAFFIAN OF S = A - A^T")
print("=" * 70)
det_S = det(S)
print(f"det(S) = {det_S:.6f}")
print(f"Since n=9 is odd, S is 9x9 skew-symmetric => det(S) = 0 (always for odd-dimensional skew-symmetric)")
print(f"Verification: det(S) ~ 0? {abs(det_S) < 1e-6}")

# For the even-dimensional principal submatrices, compute Pfaffians
print(f"\nPfaffian is only defined for even-dimensional skew-symmetric matrices.")
print(f"Computing det(S[I,I]) for all 8x8 principal subminors (removing vertex v):")

for v in range(n):
    idx = [i for i in range(n) if i != v]
    S_sub = S[np.ix_(idx, idx)]
    det_sub = det(S_sub)
    pf_sq = det_sub
    pf_sign = "+" if det_sub >= 0 else "-"
    pf_val = np.sqrt(abs(det_sub))
    print(f"  Remove v={v}: det(S_8x8) = {det_sub:>12.2f}, |Pf| = {pf_val:>8.2f}")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)

# Roots of the independence polynomial
print(f"\nRoots of I(Omega, x) = 1 + {ip_coeffs_ce[1]}x + {ip_coeffs_ce[2]}x^2 + {ip_coeffs_ce[3]}x^3:")
poly_reversed = list(reversed(ip_coeffs_ce))
roots = np.roots(poly_reversed)
for i, r in enumerate(roots):
    if abs(r.imag) < 1e-10:
        print(f"  root_{i+1} = {r.real:.6f} (real)")
    else:
        print(f"  root_{i+1} = {r.real:.6f} + {r.imag:.6f}i (COMPLEX)")
print(f"\nDiscriminant check: the polynomial has complex roots => NOT real-rooted.")
