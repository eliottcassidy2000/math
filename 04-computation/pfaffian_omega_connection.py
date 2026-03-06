#!/usr/bin/env python3
"""
PFAFFIAN-OMEGA CONNECTION EXPLORER
===================================

Deep investigation of the relationship between:
1. Pfaffian subminors Pf(S_v) of the skew-adjacency S = A - A^T
2. The Omega conflict graph structure
3. H(T) = I(Omega(T), 2)

Key observations from counterexample anatomy:
- |Pf(S_v)| = [5, 21, 3, 9, 5, 15, 5, 3, 9] for v=0,...,8
- Vertex 1 (out-deg=1, "sink") has largest |Pf| = 21
- H(T) = 237

Classical results:
- For a tournament, det(S_v) = (# spanning arborescences rooted at v)^2 [Matrix-Tree thm]
  Actually: det(L_v) = # arborescences, where L = diag(d^-) - A^T (Laplacian)
- Pf(S_v) is related but different (skew-symmetric, not Laplacian)

This script explores:
A. Whether H(T) can be expressed via Pfaffians
B. The Irving-Omar determinantal formula connection
C. Relationship between arborescences and Hamiltonian paths
D. Whether Omega structure is visible in the spectral/Pfaffian data

Author: opus-2026-03-06-S19
"""

import numpy as np
from numpy.linalg import det, eigvals
from itertools import combinations, permutations
from collections import Counter, defaultdict
import sys, os

# ============================================================
# Tournament data
# ============================================================
T_CE = np.array([
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

def compute_H(A):
    """DP Hamiltonian path count."""
    N = A.shape[0]
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
                count = sum(dp.get((prev_mask, u), 0) for u in range(N)
                           if (prev_mask & (1 << u)) and A[u, v] == 1)
                if count > 0:
                    dp[(mask, v)] = count
    full = (1 << N) - 1
    return sum(dp.get((full, v), 0) for v in range(N))

def random_tournament(n):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            if np.random.random() < 0.5:
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A

def find_directed_odd_cycles(A, N):
    cycles = set()
    def dfs(path, start):
        cur = path[-1]
        L = len(path)
        if L > 1 and A[cur, start] == 1 and L % 2 == 1:
            min_idx = path.index(min(path))
            normalized = tuple(path[min_idx:] + path[:min_idx])
            cycles.add(normalized)
        if L < N:
            for nxt in range(N):
                if nxt != start and nxt not in path and A[cur, nxt] == 1:
                    dfs(path + [nxt], start)
    for v in range(N):
        dfs([v], v)
    return cycles

def independence_poly(cycles, N):
    cycle_list = list(cycles)
    vsets = [set(c) for c in cycle_list]
    m = len(cycle_list)
    adj = {}
    for i in range(m):
        adj[i] = frozenset(j for j in range(m) if j != i and vsets[i] & vsets[j])
    # DC
    memo = {}
    def solve(verts):
        if verts in memo:
            return memo[verts]
        if not verts:
            return [1]
        v = max(verts, key=lambda u: len(adj[u] & verts))
        p1 = solve(verts - {v})
        p2 = solve(verts - (adj[v] & verts) - {v})
        maxlen = max(len(p1), len(p2) + 1)
        result = [0] * maxlen
        for i in range(len(p1)):
            result[i] += p1[i]
        for i in range(len(p2)):
            result[i + 1] += p2[i]
        memo[verts] = result
        return result
    return solve(frozenset(range(m)))

# ============================================================
print("=" * 70)
print("PFAFFIAN-OMEGA CONNECTION EXPLORER")
print("=" * 70)

S = T_CE - T_CE.T
H_val = compute_H(T_CE)
print(f"H(T) = {H_val}")

# ============================================================
# A. Pfaffian subminors and their combinatorial meaning
# ============================================================
print("\n" + "=" * 70)
print("A. PFAFFIAN SUBMINORS OF S = A - A^T")
print("=" * 70)

# For each vertex v, compute det(S_v) and |Pf(S_v)|
pf_vals = {}
for v in range(n):
    idx = [i for i in range(n) if i != v]
    S_sub = S[np.ix_(idx, idx)]
    d = det(S_sub)
    pf = np.sqrt(abs(d))
    pf_vals[v] = int(round(pf))
    out_deg = int(sum(T_CE[v]))
    in_deg = n - 1 - out_deg
    print(f"  v={v}: out-deg={out_deg}, in-deg={in_deg}, det(S_v)={d:>8.0f}, |Pf(S_v)|={pf_vals[v]}")

print(f"\n  Sum of |Pf|^2 = {sum(pf**2 for pf in pf_vals.values())}")
print(f"  Sum of |Pf| = {sum(pf_vals.values())}")
print(f"  H(T) = {H_val}")
print(f"  Ratio H(T) / sum|Pf| = {H_val / sum(pf_vals.values()):.6f}")

# ============================================================
# B. Arborescence counts (Matrix-Tree theorem for tournaments)
# ============================================================
print("\n" + "=" * 70)
print("B. ARBORESCENCES VIA LAPLACIAN")
print("=" * 70)

# In-degree Laplacian: L = diag(d^in) - A^T
# det(L_v) = # arborescences rooted at v
d_in = T_CE.T.sum(axis=1)
L = np.diag(d_in) - T_CE.T

print("In-degree Laplacian eigenvalues:")
eigs_L = sorted(eigvals(L).real, reverse=True)
for i, e in enumerate(eigs_L):
    print(f"  lambda_{i+1} = {e:.4f}")

arb_counts = {}
for v in range(n):
    idx = [i for i in range(n) if i != v]
    L_sub = L[np.ix_(idx, idx)]
    arb = det(L_sub)
    arb_counts[v] = int(round(arb.real))
    print(f"  Arborescences rooted at v={v}: {arb_counts[v]}")

total_arb = sum(arb_counts.values())
print(f"\n  Total arborescences = {total_arb}")
print(f"  H(T) = {H_val}")
print(f"  Ratio H/Arb = {H_val / total_arb:.6f}")

# ============================================================
# C. Irving-Omar: W(z) = det(I + zA^T) / det(I - zA)
# ============================================================
print("\n" + "=" * 70)
print("C. IRVING-OMAR DETERMINANTAL FORM")
print("=" * 70)

# W(z) encodes the odd cycle structure: log W(z) = 2 * sum_{k odd} tr(A^k) z^k / k
# Evaluate W at small z values and compare to I(Omega, z)

# First compute I(Omega, x) for the counterexample
cycles_ce = find_directed_odd_cycles(T_CE, n)
ip_ce = independence_poly(cycles_ce, n)
print(f"  |Omega| = {len(cycles_ce)} directed odd cycles")
print(f"  I(Omega, x) = {ip_ce}")
print(f"  I(Omega, 2) = {sum(c * 2**k for k, c in enumerate(ip_ce))}")

I = np.eye(n)
for z_val in [0.1, 0.5, 1.0, 2.0]:
    num = det(I + z_val * T_CE.T)
    den = det(I - z_val * T_CE)
    if abs(den) > 1e-10:
        W = num / den
        I_omega_z = sum(c * z_val**k for k, c in enumerate(ip_ce))
        print(f"  z={z_val}: W(z)={W:.4f}, I(Omega,z)={I_omega_z:.4f}, ratio={W/I_omega_z if abs(I_omega_z) > 1e-10 else 'inf':.4f}")
    else:
        print(f"  z={z_val}: det(I-zA) = 0 (singular)")

# ============================================================
# D. H(T-v) and Pfaffian connection
# ============================================================
print("\n" + "=" * 70)
print("D. H(T-v) VS |Pf(S_v)| — VERTEX DELETION HAMILTONIANS")
print("=" * 70)

for v in range(n):
    idx = [i for i in range(n) if i != v]
    A_sub = T_CE[np.ix_(idx, idx)]
    H_sub = compute_H(A_sub)
    pf = pf_vals[v]
    print(f"  v={v}: H(T-v) = {H_sub:>4d}, |Pf(S_v)| = {pf:>3d}, "
          f"ratio = {H_sub/pf if pf > 0 else 'inf':>8.4f}, "
          f"H(T-v)/pf^2 = {H_sub/pf**2 if pf > 0 else 'inf':>8.4f}")

# Claim A check: H(T) = H(T-v) + 2 * sum_C mu(C)
print(f"\n  Sum of H(T-v) = {sum(compute_H(T_CE[np.ix_([i for i in range(n) if i != v], [i for i in range(n) if i != v])]) for v in range(n))}")
print(f"  H(T) * (n-1) if proportional would be {H_val * (n-1)}")

# ============================================================
# E. Pattern across multiple tournaments
# ============================================================
print("\n" + "=" * 70)
print("E. PFAFFIAN-H CORRELATION ACROSS 50 RANDOM TOURNAMENTS")
print("=" * 70)

np.random.seed(42)
pf_h_data = []

for trial in range(50):
    A = random_tournament(n)
    H_r = compute_H(A)
    S_r = A - A.T

    # Get Pfaffian values
    pfs = []
    for v in range(n):
        idx = [i for i in range(n) if i != v]
        S_sub = S_r[np.ix_(idx, idx)]
        d = det(S_sub)
        pfs.append(int(round(np.sqrt(abs(d)))))

    sum_pf = sum(pfs)
    sum_pf2 = sum(p**2 for p in pfs)
    max_pf = max(pfs)

    pf_h_data.append((H_r, sum_pf, sum_pf2, max_pf, pfs))

# Correlations
hs = np.array([d[0] for d in pf_h_data], dtype=float)
sum_pfs = np.array([d[1] for d in pf_h_data], dtype=float)
sum_pf2s = np.array([d[2] for d in pf_h_data], dtype=float)
max_pfs = np.array([d[3] for d in pf_h_data], dtype=float)

def corr(x, y):
    return np.corrcoef(x, y)[0, 1]

print(f"  Corr(H, sum|Pf|)   = {corr(hs, sum_pfs):.4f}")
print(f"  Corr(H, sum|Pf|^2) = {corr(hs, sum_pf2s):.4f}")
print(f"  Corr(H, max|Pf|)   = {corr(hs, max_pfs):.4f}")

# Linear regression: H ~ a * sum_pf + b
from numpy.polynomial import polynomial as P
coeffs_fit = np.polyfit(sum_pfs, hs, 1)
print(f"\n  Best fit: H ~ {coeffs_fit[0]:.4f} * sum|Pf| + {coeffs_fit[1]:.4f}")
residuals = hs - (coeffs_fit[0] * sum_pfs + coeffs_fit[1])
print(f"  Mean residual: {np.mean(residuals):.2f}, std: {np.std(residuals):.2f}")

# Counterexample
ce_sum_pf = sum(pf_vals.values())
ce_pred = coeffs_fit[0] * ce_sum_pf + coeffs_fit[1]
print(f"\n  Counterexample: predicted H = {ce_pred:.0f}, actual = {H_val}")

# ============================================================
# F. DEEPER: det(I + xA^T) as a polynomial in x
# ============================================================
print("\n" + "=" * 70)
print("F. det(I + xA^T) AND det(I - xA) AS POLYNOMIALS")
print("=" * 70)

# det(I + xA^T) = product (1 + x * lambda_i(A^T)) = product (1 + x * conj(lambda_i(A)))
# For a tournament, these are the eigenvalues of A

evals_A = eigvals(T_CE)
print("Eigenvalues of A:")
for i, e in enumerate(sorted(evals_A, key=lambda x: -abs(x))):
    print(f"  {e.real:+.6f} {e.imag:+.6f}i")

# Compute det(I + xA^T) as polynomial by expanding
# For a 9x9 matrix, det(I + xB) = sum_{k=0}^{9} e_k(B) * x^k
# where e_k(B) = k-th elementary symmetric function of eigenvalues of B
evals_AT = eigvals(T_CE.T)

from numpy.polynomial import polynomial as NPoly
# Build det(I + xA^T) = product_i (1 + x * lambda_i)
poly_num = np.array([1.0])
for e in evals_AT:
    poly_num = np.convolve(poly_num, [1.0, e])

# Build det(I - xA) = product_i (1 - x * lambda_i)
poly_den = np.array([1.0])
for e in evals_A:
    poly_den = np.convolve(poly_den, [1.0, -e])

print(f"\nCoefficients of det(I + xA^T) [degree 0 to 9]:")
for k in range(len(poly_num)):
    c = poly_num[k]
    print(f"  x^{k}: {c.real:>12.2f}{' + ' + str(round(c.imag, 6)) + 'i' if abs(c.imag) > 0.01 else ''}")

print(f"\nCoefficients of det(I - xA) [degree 0 to 9]:")
for k in range(len(poly_den)):
    c = poly_den[k]
    print(f"  x^{k}: {c.real:>12.2f}{' + ' + str(round(c.imag, 6)) + 'i' if abs(c.imag) > 0.01 else ''}")

# These should have real coefficients (eigenvalues come in conjugate pairs + real ones)
print(f"\nVerification — all coefficients real:")
print(f"  det(I+xA^T) max|imag coeff|: {max(abs(c.imag) for c in poly_num):.2e}")
print(f"  det(I-xA) max|imag coeff|: {max(abs(c.imag) for c in poly_den):.2e}")

# Evaluate W(z) = det(I+zA^T)/det(I-zA) as ratio of polynomials
# At z=2: these relate to...?
print(f"\n  det(I + 2A^T) = {np.polyval(poly_num[::-1].real, 2):.0f}")
print(f"  det(I - 2A) = {np.polyval(poly_den[::-1].real, 2):.0f}")
num_at_2 = sum(poly_num[k].real * 2**k for k in range(len(poly_num)))
den_at_2 = sum(poly_den[k].real * 2**k for k in range(len(poly_den)))
print(f"  W(2) = {num_at_2:.0f} / {den_at_2:.0f} = {num_at_2/den_at_2:.4f}")

# ============================================================
# G. CREATIVE: Can we extract I(Omega, x) from det(I + xA^T)?
# ============================================================
print("\n" + "=" * 70)
print("G. EXTRACTING OMEGA FROM DETERMINANTAL DATA")
print("=" * 70)

# Key identity: det(I + xA) = sum_{k=0}^n sigma_k(A) * x^k
# where sigma_k = k-th elementary symmetric function of eigenvalues
# = sum of all k x k principal minors of A

# For tournament A, the k-th elementary symmetric function of eigenvalues
# relates to k-vertex sub-tournaments

# Compute principal minors of A^T (= all k-subsets)
print("Principal minors of A^T by size:")
for k in range(1, n+1):
    total = 0
    count = 0
    for subset in combinations(range(n), k):
        idx = list(subset)
        minor = det(T_CE.T[np.ix_(idx, idx)])
        total += minor
        count += 1
    print(f"  k={k}: sum of {count} minors = {total:.2f}  (e_{k}(A^T) from poly: {poly_num[k].real:.2f})")

# Compare with Omega independence polynomial
print(f"\nI(Omega, x) = {ip_ce}")
print(f"det(I + xA^T) coefficients = {[round(c.real, 2) for c in poly_num]}")
print(f"\nThese are NOT the same polynomial — det has degree {n}, I(Omega) has degree {len(ip_ce)-1}")
print(f"But both evaluate to related quantities at x=2:")
print(f"  I(Omega, 2) = {H_val}")
print(f"  W(2) = det(I+2A^T)/det(I-2A) = {num_at_2/den_at_2:.4f}")

# ============================================================
# H. CREATIVE: Permanent vs Determinant approach
# ============================================================
print("\n" + "=" * 70)
print("H. PERMANENT OF (I + xA) — does it relate to H(T)?")
print("=" * 70)

# perm(I + xA) = sum over permutations product_i (I + xA)[i, sigma(i)]
# At x=0: perm(I) = 1
# At x=1: perm(I + A) = sum_sigma product_i (delta_{i,sigma(i)} + A[i,sigma(i)])
# Expanding: sum over subsets S of [n], sum over derangements of complement,
#   product of A entries for the non-fixed points

# H(T) = perm(A) restricted to derangements that form a single cycle (Hamiltonian)
# Actually H(T) is the sum over permutations that are a SINGLE n-cycle... no.
# H(T) = number of orderings (v1,...,vn) with v_i -> v_{i+1} for all i.
# This equals permanent of A? No, perm(A) counts perfect matchings in bipartite(V_left, V_right).

# Actually: H(T) = sum over (v1,...,vn) prod A[v_i, v_{i+1}]
# = sum over (sigma in S_n) prod A[sigma(i), sigma(i+1 mod n)]... no, that's Hamiltonian cycles.
# H(T) is Hamiltonian PATHS: = sum_{sigma} prod_{i=1}^{n-1} A[sigma(i), sigma(i+1)]

# Let's compute perm(A) for comparison
def permanent(M):
    """Compute permanent via Ryser's formula."""
    n = M.shape[0]
    result = 0
    for subset_mask in range(1, 1 << n):
        cols = [j for j in range(n) if subset_mask & (1 << j)]
        k = len(cols)
        prod = 1
        for i in range(n):
            s = sum(M[i, j] for j in cols)
            prod *= s
        result += ((-1) ** (n - k)) * prod
    return result

perm_A = permanent(T_CE)
perm_I_plus_A = permanent(np.eye(n) + T_CE)
print(f"  perm(A) = {perm_A}")
print(f"  perm(I + A) = {perm_I_plus_A}")
print(f"  H(T) = {H_val}")
print(f"  perm(A) relates to cycle covers, not Hamiltonian paths directly")

# ============================================================
# I. THE KEY QUESTION: What is special about n=9?
# ============================================================
print("\n" + "=" * 70)
print("I. STRUCTURAL TRANSITION AT n=9")
print("=" * 70)

# At n <= 8, Omega(T) is always claw-free => I(Omega, x) always real-rooted
# At n = 9, claw K_{1,3} can appear in Omega
# THE QUESTION: What property of the tournament creates a claw in Omega?

# A claw in Omega = 4 directed odd cycles: one center cycle C0 that shares
# vertices with each of C1, C2, C3, while C1, C2, C3 are pairwise vertex-disjoint.

# For n=9 counterexample, let's find all claws in Omega
print("\nSearching for claws (K_{1,3}) in Omega of counterexample...")
cycle_list = list(cycles_ce)
vsets = [set(c) for c in cycle_list]
m = len(cycle_list)

# Build adjacency
adj_omega = [[False]*m for _ in range(m)]
for i in range(m):
    for j in range(i+1, m):
        if vsets[i] & vsets[j]:
            adj_omega[i][j] = adj_omega[j][i] = True

# Find independent triples
ind_triples = []
for i in range(m):
    for j in range(i+1, m):
        if adj_omega[i][j]:
            continue
        for k in range(j+1, m):
            if not adj_omega[i][k] and not adj_omega[j][k]:
                ind_triples.append((i, j, k))

print(f"  Independent triples in Omega: {len(ind_triples)}")

# For each independent triple, check if there's a cycle adjacent to all three => claw
claws = []
for trip in ind_triples:
    i, j, k = trip
    for c in range(m):
        if c in trip:
            continue
        if adj_omega[c][i] and adj_omega[c][j] and adj_omega[c][k]:
            claws.append((c, i, j, k))

print(f"  Claws (K_{{1,3}}) found: {len(claws)}")
if claws:
    # Show first few
    for claw_idx, (center, a, b, c) in enumerate(claws[:5]):
        print(f"\n  Claw #{claw_idx+1}:")
        print(f"    Center: cycle {center} = {cycle_list[center]} (len={len(cycle_list[center])})")
        print(f"    Leaf 1: cycle {a} = {cycle_list[a]} (len={len(cycle_list[a])})")
        print(f"    Leaf 2: cycle {b} = {cycle_list[b]} (len={len(cycle_list[b])})")
        print(f"    Leaf 3: cycle {c} = {cycle_list[c]} (len={len(cycle_list[c])})")
        print(f"    Center shares with leaf 1: {vsets[center] & vsets[a]}")
        print(f"    Center shares with leaf 2: {vsets[center] & vsets[b]}")
        print(f"    Center shares with leaf 3: {vsets[center] & vsets[c]}")

# ============================================================
# J. CREATIVE: Is there a weight function making I(Omega, x) always real-rooted?
# ============================================================
print("\n" + "=" * 70)
print("J. VERTEX WEIGHTING TO RESTORE REAL-ROOTEDNESS")
print("=" * 70)

# The idea: assign weight w(C) to each cycle C in Omega based on its structure.
# The weighted independence polynomial I_w(G, x) = sum_S prod_{C in S} w(C) * x^|S|
# At x=2 with appropriate weights, I_w(Omega, 2) = H(T).
# Can we find weights such that I_w is ALWAYS real-rooted?

# For the counterexample, I = 1 + 94x + 10x^2 + x^3
# The coefficients are: a0=1, a1=94, a2=10, a3=1
# Newton fails at k=2: 10^2 = 100 < 94*1*3/2 = 141

# To make Newton hold, we need a2^2 >= a1 * a3 * 3/2
# So a2 >= sqrt(a1 * a3 * 3/2)
# If a1=94, a3=1: a2 >= sqrt(141) ~ 11.87, so a2 >= 12

# What if we increase a2? We need more independent pairs.
# Currently: 10 independent pairs. All involve cycles (0,4,6) and (2,5,8).
# These are the ONLY 2 cycles avoiding vertex 3 (the hub).

# What if we SPLIT vertex 3 into a "gadget" that creates more non-hub cycles?
# This is a structural operation on the tournament, not a weight change.

# Instead, try: weight each cycle by 1/|C| or 1/2^|C| or H(complement)
print("  Weighting schemes on the counterexample:")

cycle_lens = [len(c) for c in cycle_list]
len_counts = Counter(cycle_lens)
print(f"  Cycle length distribution: {dict(len_counts)}")

# Weight by 1 (unweighted) — baseline
baseline_at_2 = sum(c * 2**k for k, c in enumerate(ip_ce))
print(f"  Unweighted: I(Omega, 2) = {baseline_at_2}")

# Weight by complement Hamiltonian count
print("\n  Computing H(T[complement(C)]) for each cycle...")
comp_H = []
for i, cyc in enumerate(cycle_list):
    comp_verts = [v for v in range(n) if v not in vsets[i]]
    if comp_verts:
        A_comp = T_CE[np.ix_(comp_verts, comp_verts)]
        h_comp = compute_H(A_comp)
    else:
        h_comp = 1
    comp_H.append(h_comp)

comp_H_dist = Counter(comp_H)
print(f"  H(complement) distribution: {dict(comp_H_dist)}")

# Sum of complement H's
print(f"  Sum of H(complement(C)): {sum(comp_H)}")
print(f"  Sum of 2*H(complement(C)) + 1 for each independent set size:")

# Weighted beta_1 = sum w(C) * 2 where w(C) = H(complement(C))
# This gives sum_C H(comp(C)) * x, which at x=2 gives...
# But the OCF already DOES this via the recursion.

# Actually, OCF says H(T) = sum_{indep S in Omega} 2^|S| * PRODUCT_{C in S} ...?
# No. OCF says H(T) = I(Omega, 2) = sum_k alpha_k * 2^k.
# Each alpha_k counts independent k-sets equally.
# The mu-weighted version is the RECURSIVE formulation, which when unrolled gives plain I(Omega, 2).

print("\n  The OCF recursion assigns effective weight mu(C) = H(T[V\\V(C)]) to each cycle.")
print("  When unrolled, this yields the UNWEIGHTED I(Omega, 2).")
print("  So the question is: can we define a DIFFERENT graph G(T) such that")
print("  I(G(T), 2) = H(T) and I(G(T), x) is always real-rooted?")

# ============================================================
# K. CREATIVE: Omega restricted to 3-cycles only
# ============================================================
print("\n" + "=" * 70)
print("K. WHAT HAPPENS WITH OMEGA_3 (3-cycles only)?")
print("=" * 70)

cycles_3 = [c for c in cycle_list if len(c) == 3]
print(f"  3-cycles: {len(cycles_3)}")
vsets_3 = [set(c) for c in cycles_3]
m3 = len(cycles_3)

adj3 = [[False]*m3 for _ in range(m3)]
for i in range(m3):
    for j in range(i+1, m3):
        if vsets_3[i] & vsets_3[j]:
            adj3[i][j] = adj3[j][i] = True

# Independence poly of Omega_3
ip3 = independence_poly(set(cycles_3), n)
print(f"  I(Omega_3, x) = {ip3}")
print(f"  I(Omega_3, 2) = {sum(c * 2**k for k, c in enumerate(ip3))}")
print(f"  H(T) = {H_val}")
print(f"  Difference: {H_val - sum(c * 2**k for k, c in enumerate(ip3))}")

# Check real-rootedness of I(Omega_3)
if len(ip3) > 2:
    roots_3 = np.roots(list(reversed(ip3)))
    all_real_3 = all(abs(r.imag) < 1e-8 for r in roots_3)
    print(f"  Real-rooted: {all_real_3}")
    for r in roots_3:
        im_part = f"+{r.imag:.4f}i" if abs(r.imag) > 1e-8 else ""
        print(f"    root: {r.real:.4f}{im_part}")

# ============================================================
# L. MOST CREATIVE: OCF via multilinear algebra
# ============================================================
print("\n" + "=" * 70)
print("L. MULTILINEAR PERSPECTIVE: H(T) AS HYPERDETERMINANT?")
print("=" * 70)

# H(T) counts Hamiltonian PATHS, which are sequences (v1,...,v_n) with
# v_sigma(1) -> v_sigma(2) -> ... -> v_sigma(n)
# = sum_{sigma in S_n} prod_{i=1}^{n-1} A[sigma(i), sigma(i+1)]

# This is the CHAIN permanent: C-perm(A) = sum_{sigma} prod A[sigma(i), sigma(i+1)]
# Different from standard permanent (bipartite matching) or hafnian.

# Can we relate this to a determinant of a LARGER matrix?
# For an n-vertex tournament, define M_{ij} = A[i][j] for i != j, M_{ii} = 0.
# Then the transfer matrix T_v approach uses M(i,j) entries.

# Actually: the classical result is that the number of Hamiltonian paths
# in a tournament of order n is ALWAYS odd (Rédei's theorem).
# More specifically, H(T) ≡ 1 (mod 2).

# The OCF says H(T) = I(Omega(T), 2). So I(Omega(T), 2) is always odd.
# Is there an analogous statement about I(Omega(T), x) modulo 2?

# Actually, I(Omega, 2) = 1 + a_1*2 + a_2*4 + a_3*8
# For this to be odd, we need 1 + 0 (mod 2) = 1, which is automatic
# since a_0 = 1 and all other terms are even.
# So Rédei's theorem is trivially satisfied by OCF!

print("  H(T) = I(Omega(T), 2) = 1 + 2*a_1 + 4*a_2 + 8*a_3 + ...")
print("  This is always 1 (mod 2), confirming Rédei's theorem trivially!")
print(f"  For counterexample: 1 + 2*94 + 4*10 + 8*1 = {1 + 2*94 + 4*10 + 8*1}")
print(f"  = {H_val} ≡ {H_val % 2} (mod 2)")

# More: H(T) mod 4?
# H(T) = 1 + 2*a_1 + 4*(a_2 + 2*a_3 + ...)
# H(T) mod 4 depends on a_1 mod 2
# If a_1 is even: H ≡ 1 (mod 4)
# If a_1 is odd: H ≡ 3 (mod 4)
a1 = ip_ce[1] if len(ip_ce) > 1 else 0
print(f"\n  a_1 = {a1} ({'even' if a1 % 2 == 0 else 'odd'})")
print(f"  H(T) mod 4 = {H_val % 4} (predicted: {'1' if a1 % 2 == 0 else '3'})")
print(f"  a_1 = number of directed odd cycles in T")

# Rédei strengthening: H(T) ≡ 1 (mod 2) for all T.
# OCF strengthening: is a_1 (# directed odd cycles) always EVEN?
# That would give H(T) ≡ 1 (mod 4).
print(f"\n  Is the number of directed odd cycles always even?")
# For the counterexample: a_1 = 94, which is even.
# Check a few random tournaments:
for trial in range(10):
    A = random_tournament(n)
    cycs = find_directed_odd_cycles(A, n)
    print(f"    Trial {trial}: {len(cycs)} directed odd cycles ({'even' if len(cycs) % 2 == 0 else 'ODD'})")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
