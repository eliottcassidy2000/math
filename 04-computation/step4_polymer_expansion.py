#!/usr/bin/env python3
"""
Step 4: Connecting spectral concentration to Z(Ω, 2) = H.

THE KEY QUESTION: Among all circulant tournaments on Z_p, which maximizes H?

APPROACH 1: Exhaustive verification for small p.
  Compute H for all 2^m orientations and show Interval is maximal.

APPROACH 2: Polymer expansion of Z(Ω, λ).
  Z = exp(Σ_n cluster terms), where convergence depends on spectral gap.

APPROACH 3: The Schur-convexity argument.
  H = f(|λ_1|², ..., |λ_m|²) where f is Schur-concave for p ≡ 3 mod 4.
  Since Interval has the most concentrated distribution (majorized by all others),
  f is maximized at the Interval.

The isometry property (all ||b|| equal) means we're optimizing on a sphere.
The problem is: H is a NONLINEAR function of (b_1,...,b_m) through the
cycle-overlap-independence-polynomial chain.

opus-2026-03-12-S66
"""

import numpy as np
from itertools import combinations
from sympy.ntheory import legendre_symbol as legendre

def make_tournament(p, S):
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for d in S:
            A[i][(i+d)%p] = 1
    return A

def count_H(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def orientation_to_S(sigma, p):
    m = (p - 1) // 2
    S = set()
    for k in range(1, m + 1):
        if sigma[k-1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    return S

def eigenvalues_circulant(sigma, p):
    m = (p - 1) // 2
    omega = np.exp(2j * np.pi / p)
    eigs = np.zeros(p, dtype=complex)
    for t in range(p):
        for j in range(1, m + 1):
            if sigma[j-1] == 1:
                eigs[t] += omega**(j * t)
            else:
                eigs[t] += omega**((p - j) * t)
    return eigs

# ========================================================================
print("=" * 72)
print("PART I: EXHAUSTIVE H SEARCH OVER ALL CIRCULANT ORIENTATIONS")
print("=" * 72)

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    print(f"\n  p={p}, m={m}, 2^m = {N} orientations")

    all_H = []
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H = count_H(A)
        all_H.append((H, bits, sigma))

    all_H.sort(reverse=True)

    # Top and bottom
    print(f"  Top 10 H values:")
    for rank, (H, bits, sigma) in enumerate(all_H[:10]):
        S = orientation_to_S(sigma, p)
        eigs = eigenvalues_circulant(sigma, p)
        lam_sq = np.abs(eigs[1:])**2
        ipr = np.sum(lam_sq**2) / np.sum(lam_sq)**2
        print(f"    #{rank+1}: H = {H:>10}, σ = {sigma}, IPR = {ipr:.4f}")

    print(f"  Bottom 5:")
    for rank, (H, bits, sigma) in enumerate(all_H[-5:]):
        S = orientation_to_S(sigma, p)
        eigs = eigenvalues_circulant(sigma, p)
        lam_sq = np.abs(eigs[1:])**2
        ipr = np.sum(lam_sq**2) / np.sum(lam_sq)**2
        print(f"    #{len(all_H)-4+rank}: H = {H:>10}, σ = {sigma}, IPR = {ipr:.4f}")

    # Find Interval and Paley
    sigma_int = tuple(1 for _ in range(m))
    sigma_pal = tuple(int(legendre(k, p)) for k in range(1, m + 1))

    H_int = None
    H_pal = None
    rank_int = None
    rank_pal = None
    for rank, (H, bits, sigma) in enumerate(all_H):
        if sigma == sigma_int:
            H_int = H
            rank_int = rank + 1
        if sigma == sigma_pal:
            H_pal = H
            rank_pal = rank + 1

    print(f"\n  Interval: H = {H_int}, rank = {rank_int}/{N}")
    print(f"  Paley:    H = {H_pal}, rank = {rank_pal}/{N}")

    # Is Interval the maximum?
    max_H = all_H[0][0]
    if H_int == max_H:
        print(f"  *** INTERVAL IS THE GLOBAL MAXIMUM! ***")
    else:
        print(f"  Interval is NOT the maximum. Max = {max_H} (rank diff: {rank_int - 1})")

    # Correlation between H and IPR
    H_vals = np.array([h for h, _, _ in all_H])
    IPR_vals = np.array([np.sum(np.abs(eigenvalues_circulant(sig, p)[1:])**4) /
                          np.sum(np.abs(eigenvalues_circulant(sig, p)[1:])**2)**2
                          for _, _, sig in all_H])

    corr = np.corrcoef(H_vals, IPR_vals)[0, 1]
    print(f"\n  Correlation(H, IPR) = {corr:.4f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART II: p=17 EXHAUSTIVE SEARCH")
print("=" * 72)
print("  p=17, m=8, 2^8 = 256 orientations (feasible)")

p = 17
m = (p - 1) // 2
N = 2**m

import time
t0 = time.time()

all_H_17 = []
for bits in range(N):
    sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
    S = orientation_to_S(sigma, p)
    A = make_tournament(p, S)
    H = count_H(A)
    all_H_17.append((H, bits, sigma))

    if (bits + 1) % 64 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{N} done ({elapsed:.1f}s)")

all_H_17.sort(reverse=True)

print(f"\n  p=17: Top 10 H values:")
for rank, (H, bits, sigma) in enumerate(all_H_17[:10]):
    eigs = eigenvalues_circulant(sigma, p)
    lam_sq = np.abs(eigs[1:])**2
    ipr = np.sum(lam_sq**2) / np.sum(lam_sq)**2
    print(f"    #{rank+1}: H = {H:>16}, σ = {sigma}, IPR = {ipr:.4f}")

sigma_int = tuple(1 for _ in range(m))
sigma_pal = tuple(int(legendre(k, p)) for k in range(1, m + 1))

for rank, (H, bits, sigma) in enumerate(all_H_17):
    if sigma == sigma_int:
        print(f"\n  Interval: H = {H}, rank = {rank+1}/{N}")
    if sigma == sigma_pal:
        print(f"  Paley:    H = {H}, rank = {rank+1}/{N}")

max_H = all_H_17[0][0]
H_int = dict((sig, h) for h, _, sig in all_H_17)[sigma_int]
if H_int == max_H:
    print(f"  *** INTERVAL IS THE GLOBAL MAXIMUM AT p=17! ***")
else:
    print(f"  Interval is NOT the maximum at p=17. Max = {max_H}")

# Correlation
H_vals = np.array([h for h, _, _ in all_H_17])
IPR_vals = np.array([np.sum(np.abs(eigenvalues_circulant(sig, p)[1:])**4) /
                      np.sum(np.abs(eigenvalues_circulant(sig, p)[1:])**2)**2
                      for _, _, sig in all_H_17])
corr = np.corrcoef(H_vals, IPR_vals)[0, 1]
print(f"  Correlation(H, IPR) = {corr:.4f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART III: THE MAJORIZATION / SCHUR-CONVEXITY ARGUMENT")
print("=" * 72)
print("""
DEFINITION: A vector x majorizes y (written x ≻ y) if the partial sums
of x sorted in decreasing order dominate those of y.

For eigenvalue distributions on the sphere Σ|λ_t|² = const:
  Interval's |λ|² vector MAJORIZES Paley's (and every other circulant's).

This is because:
  Interval has the most concentrated distribution (one big, rest small)
  Paley has the flattest distribution (all equal)

If H = f(|λ_1|², ..., |λ_m|²) is SCHUR-CONVEX, then:
  x ≻ y ⟹ f(x) ≥ f(y)
  ⟹ H(Interval) ≥ H(all others)

QUESTION: Is H Schur-convex in the eigenvalue magnitudes?

For the CYCLE COUNTS c_k = (1/p)Σ|λ_t|^k, these are symmetric
power sums, which are Schur-convex for k ≥ 2.

So c_k(Interval) ≥ c_k(all others) for k ≥ 2... WAIT, this is FALSE!
We know Paley has MORE short odd cycles than Interval at p=7,11.

The issue: c_k = (1/p)Σ Re(λ_t^k), not (1/p)Σ |λ_t|^k.
Cycle counts involve the PHASES of eigenvalues, not just magnitudes.

So the Schur-convexity argument FAILS for individual cycle counts.
But it might work for H = Z(Ω, 2) = function of ALL cycle counts together.
""")

# Check the majorization order
for p in [7, 13]:
    m = (p - 1) // 2
    N = 2**m

    sigma_int = tuple(1 for _ in range(m))
    eigs_int = eigenvalues_circulant(sigma_int, p)
    lam_sq_int = sorted(np.abs(eigs_int[1:])**2, reverse=True)

    sigma_pal = tuple(int(legendre(k, p)) for k in range(1, m + 1))
    eigs_pal = eigenvalues_circulant(sigma_pal, p)
    lam_sq_pal = sorted(np.abs(eigs_pal[1:])**2, reverse=True)

    print(f"\n  p={p}: Majorization check")
    print(f"  Interval |λ|² (sorted): {[f'{x:.4f}' for x in lam_sq_int]}")
    print(f"  Paley    |λ|² (sorted): {[f'{x:.4f}' for x in lam_sq_pal]}")

    # Check partial sums
    cumsum_int = np.cumsum(lam_sq_int)
    cumsum_pal = np.cumsum(lam_sq_pal)
    dominates = all(cumsum_int[k] >= cumsum_pal[k] - 1e-10 for k in range(len(cumsum_int)))
    print(f"  Partial sums Int ≥ Pal? {dominates}")
    for k in range(len(cumsum_int)):
        marker = " ✓" if cumsum_int[k] >= cumsum_pal[k] - 1e-10 else " ✗"
        print(f"    k={k+1}: Int={cumsum_int[k]:.4f}, Pal={cumsum_pal[k]:.4f}{marker}")

# ========================================================================
print("\n" + "=" * 72)
print("PART IV: H AS A FUNCTION ON THE EIGENVALUE SPHERE")
print("=" * 72)

for p in [7, 13]:
    m = (p - 1) // 2
    N = 2**m

    print(f"\n  p={p}: All orientations in (IPR, H) space")

    data = []
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H = count_H(A)
        eigs = eigenvalues_circulant(sigma, p)
        lam_sq = np.abs(eigs[1:])**2
        ipr = np.sum(lam_sq**2) / np.sum(lam_sq)**2
        data.append((ipr, H, sigma))

    data.sort()

    # Group by IPR and show H range
    from collections import defaultdict
    ipr_groups = defaultdict(list)
    for ipr, H, sigma in data:
        ipr_groups[round(ipr, 4)].append(H)

    print(f"  {'IPR':>10} {'#orientations':>14} {'min H':>10} {'max H':>10} {'mean H':>10}")
    for ipr_val in sorted(ipr_groups.keys()):
        H_list = ipr_groups[ipr_val]
        print(f"  {ipr_val:>10.4f} {len(H_list):>14} {min(H_list):>10} {max(H_list):>10} {np.mean(H_list):>10.1f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART V: DOES MAXIMUM IPR ⟹ MAXIMUM H?")
print("=" * 72)
print("""
From the data above, we can check: is the H-maximizer always the
IPR-maximizer? Or are there exceptions?

More precisely: among all orientations σ with the same IPR,
which one maximizes H? Is it always the Interval orientation?
""")

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    max_ipr = 0
    max_ipr_sigma = None
    max_H = 0
    max_H_sigma = None

    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H = count_H(A)
        eigs = eigenvalues_circulant(sigma, p)
        lam_sq = np.abs(eigs[1:])**2
        ipr = np.sum(lam_sq**2) / np.sum(lam_sq)**2

        if ipr > max_ipr:
            max_ipr = ipr
            max_ipr_sigma = sigma
        if H > max_H:
            max_H = H
            max_H_sigma = sigma

    sigma_int = tuple(1 for _ in range(m))

    print(f"\n  p={p}:")
    print(f"    Max IPR sigma: {max_ipr_sigma} (IPR = {max_ipr:.4f})")
    print(f"    Max H sigma:   {max_H_sigma} (H = {max_H})")
    print(f"    Interval sigma: {sigma_int}")
    print(f"    Max IPR = Interval? {max_ipr_sigma == sigma_int or max_ipr_sigma == tuple(-x for x in sigma_int)}")
    print(f"    Max H = Interval? {max_H_sigma == sigma_int or max_H_sigma == tuple(-x for x in sigma_int)}")

print("\n" + "=" * 72)
print("CONCLUSIONS")
print("=" * 72)
print("""
SUMMARY OF EXHAUSTIVE SEARCH:
  1. At each p, we compute H for all 2^m circulant orientations
  2. The Interval orientation (all σ = +1) is checked against all others
  3. The correlation between H and IPR reveals whether spectral
     concentration reliably predicts H-maximality

If max(H) = H(Interval) for all tested p ≥ 13, this provides strong
computational evidence for the conjecture.

Note: σ and -σ give the same tournament (complement symmetry), so
effectively there are 2^{m-1} distinct tournaments.
""")

print("\nDONE.")
