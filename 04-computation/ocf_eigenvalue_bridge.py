"""
ocf_eigenvalue_bridge.py — Bridge between eigenvalue power sums and I(Ω, 2)

The OCF gives H(T) = I(Ω(T), 2) where Ω(T) is the odd-cycle complex.
For circulant tournaments, eigenvalues encode cycle counts:
  C_ℓ = Tr(A^ℓ)/ℓ  for ℓ = 3, 5 (exact in tournaments)

The bridge question: Can we express I(Ω, 2) directly in terms of eigenvalues?

Strategy: For circulant tournaments on Z_p, the eigenvalue multiset
{|λ_k|} (equivalently {y_k²} where λ_k = -1/2 + iy_k) parametrizes
the tournament. We know:
  - C_3 depends on Σy_k² (universal → C_3 = C_3(p))
  - C_5 depends on Σy_k⁴ (Paley minimizes → max C_5)
  - C_7, C_9, ... depend on higher power sums

If H = I(Ω, 2) can be written as a polynomial in the Newton power sums
σ_j = Σ y_k^{2j}, then the spectral flatness argument extends.

This script tests whether H is a polynomial in σ_1, σ_2, σ_3, ...
for circulant tournaments on Z_p.

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def circulant_adj(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1
    return A


def circulant_eigenvalues(n, S):
    omega = np.exp(2j * np.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def circulant_tournaments(n):
    k = (n - 1) // 2
    pairs = [(d, n - d) for d in range(1, k + 1)]
    for bits in range(1 << k):
        S = set()
        for i in range(k):
            if bits & (1 << i):
                S.add(pairs[i][0])
            else:
                S.add(pairs[i][1])
        yield frozenset(S)


def is_paley(n, S):
    qr = set(pow(x, 2, n) for x in range(1, n))
    return set(S) == qr or set(S) == set(range(1, n)) - qr


def main():
    print("OCF-EIGENVALUE BRIDGE: IS H A POLYNOMIAL IN POWER SUMS?")
    print("=" * 75)

    for p in [7, 11, 13]:
        print(f"\n{'='*75}")
        print(f"p = {p} (p mod 4 = {p % 4})")
        print(f"{'='*75}")

        tournaments = list(circulant_tournaments(p))
        m = (p - 1) // 2  # number of free y² values

        data = []
        for S in tournaments:
            evals = circulant_eigenvalues(p, S)
            y_vals = [evals[k].imag for k in range(1, p)]
            y2_vals = [y**2 for y in y_vals]

            # Newton power sums of y²: σ_j = Σ (y_k²)^j
            sigma = {}
            for j in range(1, m + 1):
                sigma[j] = sum(y2**j for y2 in y2_vals)

            A = circulant_adj(p, S)
            H = hamiltonian_paths_dp(A, p)

            paley = is_paley(p, S)
            data.append({'S': S, 'H': H, 'sigma': sigma, 'paley': paley,
                         'y2': sorted(set(round(y2, 8) for y2 in y2_vals))})

        # Check: how many distinct (σ_1, σ_2, ...) tuples are there?
        sigma_tuples = set()
        for d in data:
            t = tuple(round(d['sigma'][j], 6) for j in range(1, m + 1))
            sigma_tuples.add(t)

        H_values = set(d['H'] for d in data)

        print(f"  {len(tournaments)} circulants, {len(H_values)} distinct H values, {len(sigma_tuples)} distinct sigma tuples")

        # Is H determined by σ_1, σ_2?
        sigma12_to_H = {}
        for d in data:
            key = (round(d['sigma'][1], 6), round(d['sigma'][2], 6))
            if key in sigma12_to_H:
                if sigma12_to_H[key] != d['H']:
                    print(f"  σ_1, σ_2 do NOT determine H")
                    break
            else:
                sigma12_to_H[key] = d['H']
        else:
            print(f"  σ_1, σ_2 DETERMINE H uniquely!")

        # Print the data
        print(f"\n  {'H':>10} {'σ_1':>10} {'σ_2':>10} {'σ_3':>10} {'Paley':>6}")
        seen = set()
        for d in sorted(data, key=lambda x: -x['H']):
            key = (d['H'], round(d['sigma'][2], 4))
            if key in seen:
                continue
            seen.add(key)
            s = d['sigma']
            label = "PALEY" if d['paley'] else ""
            print(f"  {d['H']:>10} {s[1]:>10.4f} {s[2]:>10.4f} {s.get(3, 0):>10.4f} {label:>6}")

        # Try to find a polynomial H = a₀ + a₁σ₁ + a₂σ₂ + a₃σ₁² + a₄σ₃ + ...
        # Since σ₁ is universal, the dependence is really on σ₂, σ₃, ...
        #
        # Set up linear system: H = c_0 + c_1*σ_2 + c_2*σ_3 + c_3*σ_2² + ...
        # We have as many equations as distinct (H, sigma) pairs.

        unique_data = []
        seen_sigma = set()
        for d in data:
            key = tuple(round(d['sigma'][j], 6) for j in range(1, m + 1))
            if key not in seen_sigma:
                seen_sigma.add(key)
                unique_data.append(d)

        n_pts = len(unique_data)
        print(f"\n  Fitting polynomial H(σ_2, σ_3, ...) with {n_pts} data points")

        if n_pts <= 10:
            # Try linear in σ_2 first
            sigmas = np.array([d['sigma'][2] for d in unique_data])
            Hs = np.array([d['H'] for d in unique_data])

            if n_pts >= 2:
                A_mat = np.column_stack([np.ones(n_pts), sigmas])
                coeffs, residuals, _, _ = np.linalg.lstsq(A_mat, Hs, rcond=None)
                H_pred = A_mat @ coeffs
                max_err = max(abs(H_pred - Hs))
                print(f"  Linear fit H = {coeffs[0]:.4f} + {coeffs[1]:.4f}*σ_2: max error = {max_err:.4f}")

            # Try quadratic in σ_2
            if n_pts >= 3:
                A_mat = np.column_stack([np.ones(n_pts), sigmas, sigmas**2])
                coeffs, residuals, _, _ = np.linalg.lstsq(A_mat, Hs, rcond=None)
                H_pred = A_mat @ coeffs
                max_err = max(abs(H_pred - Hs))
                print(f"  Quadratic fit H = {coeffs[0]:.4f} + {coeffs[1]:.4f}*σ_2 + {coeffs[2]:.6f}*σ_2²: max error = {max_err:.4f}")

            # Try linear in σ_2 and σ_3
            if n_pts >= 3 and 3 in unique_data[0]['sigma']:
                sigma3s = np.array([d['sigma'][3] for d in unique_data])
                A_mat = np.column_stack([np.ones(n_pts), sigmas, sigma3s])
                coeffs, residuals, _, _ = np.linalg.lstsq(A_mat, Hs, rcond=None)
                H_pred = A_mat @ coeffs
                max_err = max(abs(H_pred - Hs))
                print(f"  Linear fit H = {coeffs[0]:.4f} + {coeffs[1]:.4f}*σ_2 + {coeffs[2]:.6f}*σ_3: max error = {max_err:.4f}")

    # DEEP ANALYSIS: The H-σ relationship at p=7
    print(f"\n{'='*75}")
    print("DEEP ANALYSIS: H vs σ_2 at p=7")
    print(f"{'='*75}")

    p = 7
    tournaments = list(circulant_tournaments(p))
    for S in tournaments:
        evals = circulant_eigenvalues(p, S)
        y_vals = [evals[k].imag for k in range(1, p)]
        y2_vals = [y**2 for y in y_vals]
        sigma2 = sum(y2**2 for y2 in y2_vals)

        A = circulant_adj(p, S)
        H = hamiltonian_paths_dp(A, p)
        paley = is_paley(p, S)
        label = " PALEY" if paley else ""
        print(f"  S={sorted(S)}, σ_2 = Σy⁴ = {sigma2:.4f}, H = {H}{label}")

    # At p=7, there are only 2 distinct σ_2 values (18.375 and 46.375)
    # and 2 distinct H values (189 and 175).
    # H = 189 when σ_2 = 18.375, H = 175 when σ_2 = 46.375
    # Linear: H = 189 - (14/28)*(σ_2 - 18.375) = 189 - 0.5*(σ_2 - 18.375) = 198.1875 - 0.5*σ_2
    # Check: 198.1875 - 0.5*18.375 = 198.1875 - 9.1875 = 189 ✓
    # Check: 198.1875 - 0.5*46.375 = 198.1875 - 23.1875 = 175 ✓
    print(f"\n  H = 198.1875 - 0.5 * Σy⁴")
    print(f"  = 198.1875 - 0.5 * σ_2")
    print(f"  Equivalently: H = c₀(p) - (1/2) * σ_2")
    print(f"  where c₀(7) = 198.1875 = 1587/8")

    # This is EXACT! H is a LINEAR function of σ_2 at p=7!
    # Does this hold at p=11?
    print(f"\n{'='*75}")
    print("DEEP ANALYSIS: H vs σ_2 at p=11")
    print(f"{'='*75}")

    p = 11
    data_11 = []
    tournaments = list(circulant_tournaments(p))
    for S in tournaments:
        evals = circulant_eigenvalues(p, S)
        y_vals = [evals[k].imag for k in range(1, p)]
        sigma2 = sum(y**4 for y in y_vals)
        sigma3 = sum(y**6 for y in y_vals)

        A = circulant_adj(p, S)
        H = hamiltonian_paths_dp(A, p)
        paley = is_paley(p, S)
        data_11.append({'sigma2': sigma2, 'sigma3': sigma3, 'H': H, 'paley': paley, 'S': S})

    # Unique values
    seen = set()
    unique = []
    for d in sorted(data_11, key=lambda x: x['sigma2']):
        key = round(d['sigma2'], 4)
        if key not in seen:
            seen.add(key)
            unique.append(d)
            label = " PALEY" if d['paley'] else ""
            print(f"  σ_2={d['sigma2']:10.4f}, σ_3={d['sigma3']:12.4f}, H={d['H']:>8}{label}")

    # Check linearity in σ_2
    if len(unique) >= 2:
        sigmas = np.array([d['sigma2'] for d in unique])
        Hs = np.array([d['H'] for d in unique])

        # Linear fit
        A_mat = np.column_stack([np.ones(len(unique)), sigmas])
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, Hs, rcond=None)
        H_pred = A_mat @ coeffs
        errors = abs(H_pred - Hs)
        print(f"\n  Linear fit: H = {coeffs[0]:.4f} + {coeffs[1]:.6f} * σ_2")
        print(f"  Max error: {max(errors):.4f}")
        print(f"  Is linear: {max(errors) < 0.01}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
