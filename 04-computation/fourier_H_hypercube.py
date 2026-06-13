"""
fourier_H_hypercube.py -- kind-pasteur-2026-03-14-S73
Fourier analysis of H(T) on the Boolean hypercube {0,1}^m.

H(T) is a multilinear polynomial of degree d = 2*floor((n-1)/2) in m = C(n,2)
arc variables. Its Fourier expansion on the hypercube reveals the
interaction structure between arcs.

FOURIER ON {-1,+1}^m (Walsh-Hadamard):
Let y_e = 2*x_e - 1 in {-1,+1}. Then:
  H(T) = sum_{S subset [m]} H_hat(S) * prod_{e in S} y_e

where H_hat(S) = (1/2^m) * sum_T H(T) * prod_{e in S} y_e(T)

KEY PROPERTIES:
- H_hat(empty) = mean(H) = n!/2^{n-1}
- sum_S H_hat(S)^2 = (1/2^m) * sum_T H(T)^2 (Parseval)
- Degree d polynomial => H_hat(S) = 0 for |S| > d
- The Degree Drop means: H_hat(S) = 0 for |S| > 2*floor((n-1)/2)

QUESTIONS:
1. What is the Fourier spectrum of H? Which S have large H_hat(S)?
2. How much "energy" (sum of squares) is at each level?
3. Does the Fourier structure explain H=7 impossibility?
4. Connection to influence: Inf_e(H) = sum_{S ni e} H_hat(S)^2
5. Total influence = sum of all influences = noise sensitivity
6. Is H "noise stable" or "noise sensitive"?
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def compute_all_H(n):
    """Compute H(T) for all tournaments on n vertices."""
    m = n * (n - 1) // 2
    H_values = np.zeros(2**m, dtype=float)
    for bits in range(2**m):
        A = bits_to_adj(bits, n)
        H_values[bits] = compute_H_dp(A, n)
    return H_values

def fourier_transform(f, m):
    """
    Compute all 2^m Fourier coefficients of f: {0,1}^m -> R.
    Uses the Walsh-Hadamard transform.
    f_hat(S) = (1/2^m) * sum_x f(x) * chi_S(x)
    where chi_S(x) = (-1)^{sum_{i in S} x_i} = prod_{i in S} (1-2*x_i)
    """
    N = 2**m
    # Convert to {-1,+1} basis
    f_hat = np.zeros(N)

    for S in range(N):
        total = 0.0
        for x in range(N):
            # chi_S(x) = (-1)^{popcount(S & x)}
            parity = bin(S & x).count('1') % 2
            chi = 1 - 2 * parity
            total += f[x] * chi
        f_hat[S] = total / N

    return f_hat

def fourier_transform_fast(f, m):
    """Fast Walsh-Hadamard transform."""
    N = 2**m
    f_hat = f.copy()

    for i in range(m):
        step = 1 << (i + 1)
        half = 1 << i
        for j in range(0, N, step):
            for k in range(half):
                u = f_hat[j + k]
                v = f_hat[j + k + half]
                f_hat[j + k] = u + v
                f_hat[j + k + half] = u - v

    f_hat /= N
    return f_hat

def main():
    print("=" * 70)
    print("FOURIER ANALYSIS OF H ON THE BOOLEAN HYPERCUBE")
    print("kind-pasteur-2026-03-14-S73")
    print("=" * 70)

    for n in [3, 4, 5]:
        m = n * (n - 1) // 2
        N = 2**m
        print(f"\n{'='*70}")
        print(f"n = {n}, m = {m} arcs, N = {N} tournaments")
        print(f"{'='*70}")

        # Compute all H values
        H_values = compute_all_H(n)

        # Fourier transform
        H_hat = fourier_transform_fast(H_values, m)

        # Analyze by level (= |S| = popcount of index)
        level_energy = defaultdict(float)
        level_count_nonzero = defaultdict(int)
        level_max = defaultdict(float)

        for S in range(N):
            level = bin(S).count('1')
            energy = H_hat[S]**2
            level_energy[level] += energy
            if abs(H_hat[S]) > 1e-10:
                level_count_nonzero[level] += 1
            level_max[level] = max(level_max[level], abs(H_hat[S]))

        total_energy = sum(level_energy.values())

        print(f"\n  FOURIER SPECTRUM BY LEVEL:")
        print(f"  {'Level':>6} {'#nonzero':>10} {'Energy':>14} {'%Energy':>10} {'Max |coeff|':>14}")
        for level in sorted(level_energy.keys()):
            pct = 100 * level_energy[level] / total_energy if total_energy > 0 else 0
            print(f"  {level:6d} {level_count_nonzero[level]:10d} "
                  f"{level_energy[level]:14.4f} {pct:9.2f}% "
                  f"{level_max[level]:14.6f}")

        # Parseval check
        parseval_lhs = sum(H_hat[S]**2 for S in range(N))
        parseval_rhs = np.mean(H_values**2)
        print(f"\n  Parseval: sum H_hat^2 = {parseval_lhs:.4f}, mean H^2 = {parseval_rhs:.4f}, match = {abs(parseval_lhs - parseval_rhs) < 0.01}")

        # Key Fourier coefficients
        print(f"\n  H_hat(empty set) = mean(H) = {H_hat[0]:.6f}")
        print(f"  Expected: n!/2^(n-1) = {math.factorial(n)/2**(n-1):.6f}")

        # Level-1 coefficients (individual arc influences)
        print(f"\n  LEVEL-1 FOURIER COEFFICIENTS (single arc):")
        arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
        for idx in range(m):
            S = 1 << idx
            print(f"    arc {arcs[idx]}: H_hat = {H_hat[S]:+.6f}")

        # Influence of each arc
        print(f"\n  ARC INFLUENCES (Inf_e = sum_{{S ni e}} H_hat(S)^2):")
        for idx in range(m):
            inf = 0
            for S in range(N):
                if S & (1 << idx):
                    inf += H_hat[S]**2
            print(f"    arc {arcs[idx]}: influence = {inf:.6f}")

        # Total influence
        total_inf = 0
        for idx in range(m):
            for S in range(N):
                if S & (1 << idx):
                    total_inf += H_hat[S]**2
        print(f"\n  Total influence = {total_inf:.6f}")
        print(f"  Average influence per arc = {total_inf/m:.6f}")

        # Noise stability: rho-correlated noise
        # Stab_rho(f) = sum_S rho^|S| * H_hat(S)^2
        for rho in [0.5, 0.9, 0.99]:
            stab = sum(rho**bin(S).count('1') * H_hat[S]**2 for S in range(N))
            print(f"  Noise stability at rho={rho}: {stab:.6f} (normalized: {stab/H_hat[0]**2:.6f})")

        # Top 10 largest Fourier coefficients
        sorted_coeffs = sorted(range(N), key=lambda S: abs(H_hat[S]), reverse=True)
        print(f"\n  TOP 15 FOURIER COEFFICIENTS:")
        for rank, S in enumerate(sorted_coeffs[:15]):
            level = bin(S).count('1')
            arc_set = [arcs[i] for i in range(m) if S & (1 << i)]
            print(f"    rank {rank+1}: S={arc_set}, |S|={level}, H_hat={H_hat[S]:+.6f}")

    # ========================================
    # SPECIAL ANALYSIS: Why H=7 is impossible
    # ========================================
    print(f"\n{'='*70}")
    print("FOURIER PERSPECTIVE ON H=7 IMPOSSIBILITY")
    print("  H(T) = sum_S H_hat(S) * chi_S(T)")
    print("  For H = 7, we need sum_S H_hat(S) * chi_S(T) = 7")
    print("  The Fourier coefficients constrain which H values are achievable")
    print(f"{'='*70}")

    n = 5
    m = n * (n - 1) // 2
    H_values = compute_all_H(n)
    H_hat = fourier_transform_fast(H_values, m)

    # H range and spectrum
    H_set = sorted(set(int(h) for h in H_values))
    print(f"\n  n={n}: achievable H values = {H_set}")
    print(f"  Missing odd values in [1,{max(H_set)}]: {[h for h in range(1, max(H_set)+1, 2) if h not in H_set]}")

    # The constraint: H must be an odd integer that's achievable
    # From Fourier: H = H_hat(0) + sum_{|S|>=1} H_hat(S) * chi_S
    # H_hat(0) = 7.5 (mean). So deviations from mean must be half-integers.

    print(f"\n  H_hat(0) = {H_hat[0]:.4f} (= n!/2^(n-1))")
    print(f"  For H=7: deviation from mean = 7 - {H_hat[0]:.1f} = {7 - H_hat[0]:.1f}")
    print(f"  For H=9: deviation from mean = 9 - {H_hat[0]:.1f} = {9 - H_hat[0]:.1f}")

    # Can the Fourier coefficients produce deviation = -0.5?
    # This requires sum_{|S|>=1} H_hat(S) * chi_S = -0.5
    # where chi_S in {-1,+1} for each nonzero S

    # What is the range of sum H_hat(S) * chi_S over all sign choices?
    nonzero_coeffs = [(S, H_hat[S]) for S in range(1, 2**m) if abs(H_hat[S]) > 1e-10]
    max_deviation = sum(abs(c) for _, c in nonzero_coeffs)
    print(f"\n  Max possible |deviation| = {max_deviation:.4f}")
    print(f"  This gives H range [{H_hat[0] - max_deviation:.1f}, {H_hat[0] + max_deviation:.1f}]")

    # The deviation -0.5 for H=7 is small. Why can't it be achieved?
    # Answer: the chi_S values are NOT independent! They're determined by
    # the tournament T, and a valid tournament constrains which sign patterns
    # (chi_S(T)) can occur.

    # Compute the actual achieved deviations
    deviations = sorted(set(int(h) - H_hat[0] for h in H_values))
    print(f"  Achieved deviations from mean: {[f'{d:+.1f}' for d in deviations]}")

    # ========================================
    # SPECTRAL GAP / CONCENTRATION
    # ========================================
    print(f"\n{'='*70}")
    print("SPECTRAL CONCENTRATION: WHERE IS THE ENERGY?")
    print(f"{'='*70}")

    for n in [3, 4, 5]:
        m = n * (n - 1) // 2
        H_values = compute_all_H(n)
        H_hat = fourier_transform_fast(H_values, m)

        # Energy at each level
        total = sum(H_hat[S]**2 for S in range(2**m))
        level_pct = defaultdict(float)
        for S in range(2**m):
            level_pct[bin(S).count('1')] += H_hat[S]**2 / total * 100

        print(f"\n  n={n}: Energy distribution by Fourier level:")
        for k in sorted(level_pct.keys()):
            bar = '#' * int(level_pct[k] / 2)
            print(f"    level {k}: {level_pct[k]:6.2f}% {bar}")

    # ========================================
    # DEGREE DROP IN FOURIER LANGUAGE
    # ========================================
    print(f"\n{'='*70}")
    print("DEGREE DROP IN FOURIER LANGUAGE")
    print("  Degree Drop Theorem: H_hat(S) = 0 for |S| > 2*floor((n-1)/2)")
    print("  This means the TOP Fourier levels are EMPTY")
    print(f"{'='*70}")

    for n in [3, 4, 5]:
        m = n * (n - 1) // 2
        d = 2 * ((n-1) // 2)  # degree = 2*floor((n-1)/2)
        H_values = compute_all_H(n)
        H_hat = fourier_transform_fast(H_values, m)

        max_nonzero_level = max(bin(S).count('1') for S in range(2**m) if abs(H_hat[S]) > 1e-10)

        print(f"\n  n={n}: predicted max level = {d}, actual max level = {max_nonzero_level}")
        print(f"    m = {m}, n-1 = {n-1}, degree drop = {(n-1) - d}")

        # How many levels are "wasted" (all zero)?
        empty_levels = [k for k in range(d+1, m+1) if all(
            abs(H_hat[S]) < 1e-10 for S in range(2**m) if bin(S).count('1') == k)]
        print(f"    Empty levels above degree: {empty_levels}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
