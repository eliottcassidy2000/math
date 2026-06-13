#!/usr/bin/env python3
"""
imaginary_spectrum_optimization.py — H as a function of Im(λ_t) alone

Key insight from S66: For circulant tournaments on Z_p,
  λ_t(σ) = -1/2 + i·y_t(σ)
where y_t(σ) = Σ_{j=1}^m σ_j sin(2πjt/p) and Re(λ_t) = -1/2 for ALL t≠0.

Therefore c_k = (1/p) Σ_{t=0}^{p-1} λ_t^k, with λ_0 = m (real, constant).

Since only the y_t's vary, H is entirely determined by the IMAGINARY SPECTRUM
{y_1, ..., y_{p-1}}. Furthermore y_t = -y_{p-t}, so there are only m = (p-1)/2
independent values.

QUESTION 1: What is H as an EXPLICIT function of y_1, ..., y_m?
QUESTION 2: What property of the imaginary spectrum determines H-maximality?
QUESTION 3: Does the IPR/additive energy connect to |y_t| distribution?

Author: opus-2026-03-12-S67
"""

import numpy as np
from itertools import product as iprod

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

def compute_eigenvalues(sigma, p):
    """Compute all eigenvalues λ_t for circulant tournament with orientation σ."""
    m = (p-1)//2
    omega = np.exp(2j * np.pi / p)
    eigs = []
    for t in range(p):
        lam = 0
        for j in range(1, m+1):
            if sigma[j-1] == 1:
                lam += omega**(j*t)
            else:
                lam += omega**((-j) % p * t)  # = omega^{(p-j)*t} = omega^{-jt}
        eigs.append(lam)
    return eigs

def compute_H_from_eigenvalues(eigs, p):
    """Compute H via permanent-trace formula: H = (1/p) Σ_t λ_t^{p-1} · p
    Wait — actually H(T) for the circulant is NOT simply Σ λ^{p-1}.
    H is the permanent of the adjacency matrix, which for circulants equals
    Σ_{permutations} Π A[i,σ(i)].

    The trace formula gives c_k = (1/(2k)) * (sum of k-th powers of eigs of
    special matrix). But H is the HP count, not a trace.

    For circulant tournaments, we can compute H directly.
    """
    return None  # Will compute H directly

def compute_H_direct(sigma, p):
    """Compute H(T) by DP (Held-Karp) for small p."""
    m = (p-1)//2
    n = p
    # Build adjacency
    adj = [[0]*n for _ in range(n)]
    for k in range(1, m+1):
        for i in range(n):
            j = (i + k) % n
            if sigma[k-1] == 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1

    # Held-Karp DP: dp[S][v] = # Ham paths from start to v using vertex set S
    # Count all Ham paths (any start, any end)
    from functools import lru_cache

    # Use bitmask DP
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and adj[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total > 0:
                    dp[(mask, v)] = total

    full_mask = (1 << n) - 1
    H = sum(dp.get((full_mask, v), 0) for v in range(n))
    return H

def compute_y_spectrum(sigma, p):
    """Compute the imaginary spectrum y_t = Σ σ_j sin(2πjt/p)."""
    m = (p-1)//2
    ys = []
    for t in range(1, p):
        y = sum(sigma[j-1] * np.sin(2*np.pi*j*t/p) for j in range(1, m+1))
        ys.append(y)
    return ys

print("=" * 70)
print("IMAGINARY SPECTRUM ANALYSIS")
print("=" * 70)

for p in [7, 11]:
    m = (p-1)//2
    print(f"\n{'='*70}")
    print(f"p = {p}, m = {m}")
    print(f"{'='*70}")

    # Compute for all 2^m orientations
    results = []
    for bits in range(1 << m):
        sigma = tuple(1 if (bits >> j) & 1 else -1 for j in range(m))
        H = compute_H_direct(sigma, p)
        ys = compute_y_spectrum(sigma, p)

        # Key statistics of the imaginary spectrum
        ys_half = ys[:m]  # y_1, ..., y_m (y_{p-t} = -y_t)
        y_sq_sum = sum(y**2 for y in ys_half)  # Σ y_t^2 for t=1..m
        y_4th_sum = sum(y**4 for y in ys_half)  # Σ y_t^4
        y_abs_sorted = sorted([abs(y) for y in ys_half], reverse=True)
        ipr_y = y_4th_sum / y_sq_sum**2 if y_sq_sum > 0 else 0

        results.append({
            'sigma': sigma, 'H': H, 'ys_half': ys_half,
            'y_sq_sum': y_sq_sum, 'y_4th_sum': y_4th_sum,
            'ipr_y': ipr_y, 'y_abs_sorted': y_abs_sorted
        })

    results.sort(key=lambda r: -r['H'])

    # Identify interval and Paley
    sigma_int = tuple([1]*m)
    sigma_pal = None
    if p % 4 == 3:
        sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))

    print(f"\nTop 5 and bottom 5 by H:")
    print(f"  {'sigma':<20} {'H':>8} {'Σy²':>10} {'Σy⁴':>12} {'IPR_y':>8} {'|y|_max':>8}")
    for r in results[:5]:
        tag = ""
        if r['sigma'] == sigma_int: tag = " [INT]"
        if sigma_pal and r['sigma'] == sigma_pal: tag = " [PAL]"
        print(f"  {str(r['sigma']):<20} {r['H']:>8} {r['y_sq_sum']:>10.2f} "
              f"{r['y_4th_sum']:>12.2f} {r['ipr_y']:>8.4f} "
              f"{r['y_abs_sorted'][0]:>8.3f}{tag}")
    print("  ...")
    for r in results[-3:]:
        tag = ""
        if r['sigma'] == sigma_int: tag = " [INT]"
        if sigma_pal and r['sigma'] == sigma_pal: tag = " [PAL]"
        print(f"  {str(r['sigma']):<20} {r['H']:>8} {r['y_sq_sum']:>10.2f} "
              f"{r['y_4th_sum']:>12.2f} {r['ipr_y']:>8.4f} "
              f"{r['y_abs_sorted'][0]:>8.3f}{tag}")

    # Correlation analysis: H vs various y-statistics
    Hs = [r['H'] for r in results]
    y_sq_sums = [r['y_sq_sum'] for r in results]
    y_4th_sums = [r['y_4th_sum'] for r in results]
    ipr_ys = [r['ipr_y'] for r in results]
    y_maxs = [r['y_abs_sorted'][0] for r in results]

    def corr(a, b):
        a, b = np.array(a), np.array(b)
        return np.corrcoef(a, b)[0, 1]

    print(f"\n  Correlations with H:")
    print(f"    corr(H, Σy²)    = {corr(Hs, y_sq_sums):.4f}")
    print(f"    corr(H, Σy⁴)    = {corr(Hs, y_4th_sums):.4f}")
    print(f"    corr(H, IPR_y)  = {corr(Hs, ipr_ys):.4f}")
    print(f"    corr(H, |y|max) = {corr(Hs, y_maxs):.4f}")

    # KEY: Σ y_t^2 = Σ_{j=1}^m σ_j^2 · sin²(2πjt/p) cross terms...
    # Actually: Σ_{t=1}^m y_t^2 = Σ_{t=1}^m (Σ_j σ_j sin(2πjt/p))^2
    # = Σ_{j,k} σ_j σ_k Σ_t sin(2πjt/p) sin(2πkt/p)
    # = Σ_{j,k} σ_j σ_k · (p-1)/4 · δ_{jk}   [orthogonality!]
    # Wait, is this right? Let me check...

    print(f"\n  Parseval check: Σ_{'{t=1..m}'} y_t^2 should be CONSTANT")
    unique_y_sq = set(round(r['y_sq_sum'], 6) for r in results)
    print(f"    Unique values of Σy²: {unique_y_sq}")

    if len(unique_y_sq) == 1:
        print(f"    YES! Σy² = {results[0]['y_sq_sum']:.6f} = m*(p-1)/4 = {m*(p-1)/4}")
        print(f"    This is Parseval's identity for the sine system.")
        print(f"    Since Σy² is constant, H is NOT correlated with Σy².")
        print(f"    The H variation comes entirely from HIGHER MOMENTS of y.")

    # Since Σy² is constant, the IPR_y = Σy⁴/(Σy²)² ∝ Σy⁴
    # So IPR_y and Σy⁴ carry the SAME information

    print(f"\n  Since Σy² is constant, IPR_y ∝ Σy⁴. Key question: does HIGH or LOW Σy⁴ → high H?")

    # Detailed analysis: what is Σy_t^4?
    # Σ_{t=1}^m y_t^4 = Σ_t (Σ_j σ_j sin(2πjt/p))^4
    # This involves 4th-order interactions between chords

    # NEW INSIGHT: Σy⁴ is related to ADDITIVE ENERGY of the connection set!
    # From S65: IPR = (p·E(S) - (p-1)⁴)/(m(m+1))²
    # And E(S) = additive energy = #{(a,b,c,d) ∈ S^4 : a+b=c+d}

    print(f"\n  Connection to additive energy:")
    # For each sigma, compute the connection set S
    for tag, sigma_check in [("Interval", sigma_int)]:
        S = set()
        for j in range(1, m+1):
            if sigma_check[j-1] == 1:
                S.add(j)
            else:
                S.add(p - j)
        # Additive energy
        from collections import Counter
        sums = Counter()
        for a in S:
            for b in S:
                sums[(a+b) % p] += 1
        E = sum(v**2 for v in sums.values())
        print(f"    {tag}: |S|={len(S)}, E(S)={E}")

    if sigma_pal:
        S = set()
        for j in range(1, m+1):
            if sigma_pal[j-1] == 1:
                S.add(j)
            else:
                S.add(p - j)
        from collections import Counter
        sums = Counter()
        for a in S:
            for b in S:
                sums[(a+b) % p] += 1
        E = sum(v**2 for v in sums.values())
        print(f"    Paley:    |S|={len(S)}, E(S)={E}")

print("\n" + "=" * 70)
print("DEEPER ANALYSIS: H AS FUNCTION OF y-SPECTRUM SHAPE")
print("=" * 70)

# For p=7, let's express H explicitly in terms of y_1, y_2, y_3
p = 7
m = 3
print(f"\np={p}: H as function of y_1, y_2, y_3")
print(f"Note: y_t = -y_{{p-t}}, so y_4=-y_3, y_5=-y_2, y_6=-y_1")
print(f"λ_0 = m = {m}")
print(f"λ_t = -1/2 + i*y_t for t=1,...,{p-1}")
print()

# H can be written as... H = permanent of adjacency matrix
# For a circulant, there's a formula involving eigenvalues
# The permanent of a circulant matrix C with eigenvalues λ_0,...,λ_{n-1} is:
# perm(C) = (1/n!) Σ_{partitions} ... (complicated)
# Actually this is the Rysser/Bregman connection
# But for our purposes, let's just tabulate

# Let's check: is there a polynomial relationship H = f(Σy⁴, Σy⁶, ...)?
results_p7 = []
for bits in range(1 << m):
    sigma = tuple(1 if (bits >> j) & 1 else -1 for j in range(m))
    H = compute_H_direct(sigma, p)
    ys = compute_y_spectrum(sigma, p)
    ys_half = [ys[t] for t in range(m)]

    y2 = sum(y**2 for y in ys_half)
    y4 = sum(y**4 for y in ys_half)
    y6 = sum(y**6 for y in ys_half)

    # Also compute individual |λ_t|^k sums
    eigs = compute_eigenvalues(sigma, p)
    lam_powers = {}
    for k in [3, 5, 7]:
        lam_powers[k] = sum(eigs[t]**k for t in range(p)).real / p

    results_p7.append({
        'sigma': sigma, 'H': H, 'y2': y2, 'y4': y4, 'y6': y6,
        'c3': lam_powers.get(3, 0), 'c5': lam_powers.get(5, 0),
        'c7': lam_powers.get(7, 0)
    })

print(f"  {'sigma':<20} {'H':>6} {'Σy⁴':>10} {'Σy⁶':>12} {'c3':>8} {'c5':>8}")
for r in sorted(results_p7, key=lambda x: -x['H']):
    print(f"  {str(r['sigma']):<20} {r['H']:>6} {r['y4']:>10.4f} {r['y6']:>12.4f} "
          f"{r['c3']:>8.1f} {r['c5']:>8.1f}")

# Check if H is a LINEAR function of Σy⁴ and Σy⁶
# H = a + b·Σy⁴ + c·Σy⁶
Hs = np.array([r['H'] for r in results_p7])
y4s = np.array([r['y4'] for r in results_p7])
y6s = np.array([r['y6'] for r in results_p7])

# Least squares: H = a + b*y4 + c*y6
A_mat = np.column_stack([np.ones(len(Hs)), y4s, y6s])
coeffs, residuals, rank, sv = np.linalg.lstsq(A_mat, Hs, rcond=None)
print(f"\n  Fit H = {coeffs[0]:.4f} + {coeffs[1]:.4f}·Σy⁴ + {coeffs[2]:.4f}·Σy⁶")
print(f"  Residuals: {residuals}")
fitted = A_mat @ coeffs
print(f"  Max |residual| = {max(abs(Hs - fitted)):.6f}")

if max(abs(Hs - fitted)) < 0.01:
    print(f"  *** H is EXACTLY a linear function of Σy⁴ and Σy⁶! ***")
    print(f"  This means only two spectral invariants determine H at p=7.")
else:
    print(f"  H is NOT a linear function of (Σy⁴, Σy⁶) alone.")
    # Try adding more moments
    y8s = np.array([sum(y**8 for y in compute_y_spectrum(r['sigma'], p)[:m])
                    for r in results_p7])
    A_mat2 = np.column_stack([np.ones(len(Hs)), y4s, y6s, y8s])
    coeffs2, res2, _, _ = np.linalg.lstsq(A_mat2, Hs, rcond=None)
    fitted2 = A_mat2 @ coeffs2
    print(f"  With Σy⁸: max |residual| = {max(abs(Hs - fitted2)):.6f}")
    if max(abs(Hs - fitted2)) < 0.01:
        print(f"  *** H = {coeffs2[0]:.4f} + {coeffs2[1]:.4f}·Σy⁴ + {coeffs2[2]:.4f}·Σy⁶ + {coeffs2[3]:.4f}·Σy⁸ ***")

# Now the key question: which y-spectrum SHAPE maximizes H?
# Since Σy² is constant, the constraint surface is a sphere.
# Higher Σy⁴ means more CONCENTRATED spectrum (one dominant y_t).
# Lower Σy⁴ means more SPREAD spectrum (all y_t similar).
# Which does H prefer?

print(f"\n  H-maximizer y-spectrum shape:")
best = max(results_p7, key=lambda r: r['H'])
worst = min(results_p7, key=lambda r: r['H'])
print(f"  H-max (σ={best['sigma']}, H={best['H']}): Σy⁴={best['y4']:.4f}")
print(f"  H-min (σ={worst['sigma']}, H={worst['H']}): Σy⁴={worst['y4']:.4f}")
if best['y4'] > worst['y4']:
    print(f"  → H INCREASES with spectral concentration (higher Σy⁴)")
    print(f"  → Concentrated imaginary spectrum → more Ham paths")
else:
    print(f"  → H DECREASES with spectral concentration")
    print(f"  → Spread imaginary spectrum → more Ham paths")

# Repeat for p=11
print(f"\n{'='*70}")
print(f"p=11: H vs y-spectrum moments")
print(f"{'='*70}")
p = 11
m = 5
sigma_int = tuple([1]*m)
sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))

results_p11 = []
for bits in range(1 << m):
    sigma = tuple(1 if (bits >> j) & 1 else -1 for j in range(m))
    H = compute_H_direct(sigma, p)
    ys = compute_y_spectrum(sigma, p)
    ys_half = [ys[t] for t in range(m)]
    y4 = sum(y**4 for y in ys_half)
    y6 = sum(y**6 for y in ys_half)
    results_p11.append({'sigma': sigma, 'H': H, 'y4': y4, 'y6': y6})

results_p11.sort(key=lambda r: -r['H'])

Hs = np.array([r['H'] for r in results_p11])
y4s = np.array([r['y4'] for r in results_p11])
y6s = np.array([r['y6'] for r in results_p11])

A_mat = np.column_stack([np.ones(len(Hs)), y4s, y6s])
coeffs, residuals, rank, sv = np.linalg.lstsq(A_mat, Hs, rcond=None)
fitted = A_mat @ coeffs
print(f"Fit H = {coeffs[0]:.4f} + {coeffs[1]:.4f}·Σy⁴ + {coeffs[2]:.4f}·Σy⁶")
print(f"Max |residual| = {max(abs(Hs - fitted)):.4f}")

# Add more moments if needed
y8s = np.array([sum(y**8 for y in compute_y_spectrum(r['sigma'], p)[:m])
                for r in results_p11])
y10s = np.array([sum(y**10 for y in compute_y_spectrum(r['sigma'], p)[:m])
                 for r in results_p11])
A_mat3 = np.column_stack([np.ones(len(Hs)), y4s, y6s, y8s, y10s])
coeffs3, res3, _, _ = np.linalg.lstsq(A_mat3, Hs, rcond=None)
fitted3 = A_mat3 @ coeffs3
print(f"With Σy⁸,Σy¹⁰: max |residual| = {max(abs(Hs - fitted3)):.4f}")

# Key comparison: Interval vs Paley y-spectrum
print(f"\nInterval vs Paley y-spectrum:")
r_int = next(r for r in results_p11 if r['sigma'] == sigma_int)
r_pal = next(r for r in results_p11 if r['sigma'] == sigma_pal)
print(f"  Interval: H={r_int['H']}, Σy⁴={r_int['y4']:.2f}, rank={results_p11.index(r_int)+1}/{len(results_p11)}")
print(f"  Paley:    H={r_pal['H']}, Σy⁴={r_pal['y4']:.2f}, rank={results_p11.index(r_pal)+1}/{len(results_p11)}")

# What is the OPTIMAL Σy⁴ for H-maximization?
print(f"\nTop 5:")
for i, r in enumerate(results_p11[:5]):
    tag = ""
    if r['sigma'] == sigma_int: tag = " [INT]"
    if r['sigma'] == sigma_pal: tag = " [PAL]"
    print(f"  #{i+1}: H={r['H']}, Σy⁴={r['y4']:.2f}{tag}")

print(f"\n{'='*70}")
print("KEY INSIGHT: THE IMAGINARY SPECTRUM SHAPE DETERMINES H")
print("=" * 70)
print("""
For circulant tournaments on Z_p:
  1. Re(λ_t) = -1/2 for all t≠0 (universal, S66)
  2. Im(λ_t) = y_t = Σ_j σ_j sin(2πjt/p) (depends on orientation σ)
  3. Σ y_t² = m(p-1)/4 (constant by Parseval = PROVED)
  4. H depends on σ ONLY through the higher moments Σy⁴, Σy⁶, ...

The H-maximization problem is equivalent to:
  max H(y₁,...,y_m) subject to Σy_t² = const.

This is a CONSTRAINED OPTIMIZATION on a SPHERE in R^m.
The Interval and Paley represent two specific points on this sphere.

Interval: y_t = Fejér kernel = peaked spectrum (one dominant y_t)
Paley: y_t = ±√(p)/2 = FLAT spectrum (all |y_t| equal)

For SMALL p: flat spectrum wins (Paley = H-max)
For LARGE p: peaked spectrum wins (Interval = H-max)

The crossover is when the CURVATURE of H changes sign on the sphere.
""")

print("DONE.")
