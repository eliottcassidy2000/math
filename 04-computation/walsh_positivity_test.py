#!/usr/bin/env python3
"""
Walsh positivity test — the key to proving Interval maximizes H.

CENTRAL QUESTION: For p ≥ 13, are ALL Walsh coefficients ĥ[S] ≥ 0?

If yes, then H(σ) = Σ_S ĥ[S] χ_S(σ) ≤ Σ_S ĥ[S] = H(all-ones) = H(Interval)
since χ_S(σ) ≤ 1 with equality at all-ones. ONE-LINE PROOF.

This script computes H once per orientation (the bottleneck), then does
all Walsh analysis from the cached H values.

opus-2026-03-12-S67
"""

import numpy as np
import time
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

def compute_walsh_spectrum(H_vec, m):
    """Compute full Walsh-Hadamard transform of H_vec."""
    N = len(H_vec)
    h_hat = np.zeros(N)
    for S_bits in range(N):
        val = 0.0
        for bits in range(N):
            chi = 1
            # χ_S(σ) = Π_{k∈S} σ_k
            # σ_k = +1 if bit k is set, -1 otherwise
            # So χ_S(σ) = (-1)^(number of k∈S where bit k is 0)
            # = (-1)^(|S| - popcount(S & bits))
            common = bin(S_bits & bits).count('1')
            deg = bin(S_bits).count('1')
            chi = (-1)**(deg - common)
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N
    return h_hat

print("=" * 72)
print("WALSH POSITIVITY TEST — COMPUTING H VALUES")
print("=" * 72)

# Cache H values for all primes
cached = {}

for p in [7, 11, 13, 17]:
    m = (p - 1) // 2
    N = 2**m

    print(f"\n  Computing H for p={p} ({N} orientations)...", flush=True)
    t0 = time.time()

    H_vec = np.zeros(N, dtype=np.int64)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)
        if (bits + 1) % 64 == 0:
            elapsed = time.time() - t0
            print(f"    {bits+1}/{N} ({elapsed:.1f}s)", flush=True)

    elapsed = time.time() - t0
    print(f"    Done in {elapsed:.1f}s")

    cached[p] = (m, N, H_vec)

# ========================================================================
print("\n" + "=" * 72)
print("WALSH TRANSFORM AND POSITIVITY CHECK")
print("=" * 72)

for p in [7, 11, 13, 17]:
    m, N, H_vec = cached[p]

    print(f"\n  p={p}, m={m}:")

    # Compute Walsh spectrum
    h_hat = compute_walsh_spectrum(H_vec.astype(float), m)

    # Positivity check
    all_positive = True
    by_degree = {}

    for S_bits in range(N):
        deg = bin(S_bits).count('1')
        if deg not in by_degree:
            by_degree[deg] = {'pos': 0, 'neg': 0, 'zero': 0,
                               'pos_sum': 0.0, 'neg_sum': 0.0,
                               'coeffs': []}
        if abs(h_hat[S_bits]) < 0.001:
            by_degree[deg]['zero'] += 1
        elif h_hat[S_bits] > 0:
            by_degree[deg]['pos'] += 1
            by_degree[deg]['pos_sum'] += h_hat[S_bits]
        else:
            by_degree[deg]['neg'] += 1
            by_degree[deg]['neg_sum'] += h_hat[S_bits]
            if deg > 0:
                all_positive = False

        if abs(h_hat[S_bits]) > 0.001 and deg > 0:
            S_set = sorted(k+1 for k in range(m) if (S_bits >> k) & 1)
            by_degree[deg]['coeffs'].append((S_set, h_hat[S_bits]))

    # Print summary by degree
    for deg in sorted(by_degree.keys()):
        d = by_degree[deg]
        n_total = d['pos'] + d['neg'] + d['zero']
        if d['pos'] + d['neg'] == 0:
            continue
        sign_str = "ALL+" if d['neg'] == 0 else "ALL-" if d['pos'] == 0 else "MIXED"
        print(f"    Degree {deg}: {d['pos']}+, {d['neg']}-, {d['zero']}=0 [{sign_str}]"
              f"  sum(+)={d['pos_sum']:.2f}, sum(-)={d['neg_sum']:.2f}")

        # Print individual coefficients for small cases
        if len(d['coeffs']) <= 20:
            for S_set, coeff in sorted(d['coeffs'], key=lambda x: -abs(x[1])):
                sign = "+" if coeff > 0 else "-"
                print(f"      [{sign}] ĥ{S_set} = {coeff:.4f}")

    if all_positive:
        print(f"    ★★★ ALL NONZERO WALSH COEFFICIENTS ARE POSITIVE! ★★★")
        print(f"    ⟹ H(all-ones) = H(Interval) is the UNIQUE MAXIMUM!")
    else:
        # Check: is f₂ maximized at all-ones? Is f₄?
        bits_all_ones = N - 1
        for deg in [2, 4, 6]:
            if deg not in by_degree or not by_degree[deg]['coeffs']:
                continue
            # f_deg at all-ones = sum of all ĥ[S] with |S|=deg (since χ_S(1) = 1)
            f_at_ones = sum(c for _, c in by_degree[deg]['coeffs'])
            # f_deg maximum over all σ
            f_max = -float('inf')
            f_max_sigma = None
            for bits in range(N):
                f_val = 0.0
                for S_set, coeff in by_degree[deg]['coeffs']:
                    chi = 1
                    for k in S_set:
                        if not ((bits >> (k-1)) & 1):
                            chi *= -1
                    f_val += coeff * chi
                if f_val > f_max:
                    f_max = f_val
                    sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
                    f_max_sigma = sigma

            is_max = abs(f_at_ones - f_max) < 0.01
            print(f"    f_{deg}(all-ones) = {f_at_ones:.4f}, max = {f_max:.4f} "
                  f"{'← ALL-ONES IS MAX' if is_max else '← NOT max, at ' + str(f_max_sigma)}")

    # The big picture: H(all-ones) vs max H
    H_int = H_vec[N-1]
    H_max = max(H_vec)
    H_max_idx = np.argmax(H_vec)
    H_max_sigma = tuple(1 if (H_max_idx >> k) & 1 else -1 for k in range(m))
    print(f"\n    H(Interval) = {H_int}")
    print(f"    H(max) = {H_max} at σ = {H_max_sigma}")
    print(f"    Interval is max? {H_int == H_max}")

# ========================================================================
print("\n" + "=" * 72)
print("THE TRANSITION PATTERN")
print("=" * 72)

print("""
Summary of Walsh coefficient signs:

  p=7 (m=3): Degree 2 has MIXED signs → Paley optimizes degree 2 → Paley wins
  p=11 (m=5): Both degree 2 and 4 have coefficients → degree 2 dominates
  p=13 (m=6): Degree 4 dominates (81.6% of variance)
  p=17 (m=8): ???

KEY OBSERVATION from previous analysis:
  - Degree-1 Walsh coefficients are ALL ZERO (complement symmetry)
  - Only EVEN degrees appear
  - The dominant degree INCREASES with p

PROOF STRATEGY:
  If ALL degree-4 coefficients are positive for p ≥ 13, and degree-4
  dominates, then all-ones maximizes the dominant term.
  The degree-2 "loss" (where Paley or other σ might be better)
  is overcome by the degree-4 "gain".

  This is analogous to the SECOND-ORDER phase transition in physics:
  the quartic term (degree 4) overcomes the quadratic term (degree 2)
  above the critical temperature.
""")

# Final comparison: degree-2 vs degree-4 energy fraction
for p in [7, 11, 13, 17]:
    m, N, H_vec = cached[p]
    h_hat = compute_walsh_spectrum(H_vec.astype(float), m)

    energy = {}
    for S_bits in range(1, N):
        deg = bin(S_bits).count('1')
        energy[deg] = energy.get(deg, 0.0) + h_hat[S_bits]**2

    total_e = sum(energy.values())
    e2 = energy.get(2, 0)
    e4 = energy.get(4, 0)
    print(f"  p={p}: E₂/E_total = {e2/total_e:.4f}, E₄/E_total = {e4/total_e:.4f}, "
          f"E₄/E₂ = {e4/max(e2,1e-10):.4f}")

print("\nDONE.")
