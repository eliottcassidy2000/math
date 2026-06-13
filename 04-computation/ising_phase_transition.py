"""
ising_phase_transition.py — The tournament Ising model and its phase transition

DEEP CONNECTION: The H-maximization problem on circulant tournaments is
EXACTLY a ground state problem for a generalized Ising model.

Setup:
  - m = (p-1)/2 "spins" σ_i ∈ {±1}
  - σ_i = +1 means chord type i is "forward" (i ∈ S)
  - Energy function E(σ) = -H(σ) (we minimize E = maximize H)
  - Walsh expansion: E = -H₀ - Σ_{|S|=2} Ĥ(S)χ_S - Σ_{|S|=4} Ĥ(S)χ_S - ...

Phase transition:
  - At p=7 (m=3): ONLY degree-2 interactions → Paley is ground state
  - At p=11 (m=5): degree-2 + degree-4, still Paley wins
  - At p≥19 (m≥9): higher-degree terms dominate → Interval is ground state

KEY INSIGHT: The "temperature" in this analogy is 1/m (or 1/log(p)).
As m increases, the energy landscape becomes more complex (exponentially
more local minima), and the globally smooth (Paley) ground state is
replaced by a locally structured (Interval) one.

This is analogous to:
  - Spin glass transition in random field Ising models
  - Phase transition in MAX-CUT on random graphs
  - The replica symmetry breaking in Sherrington-Kirkpatrick model

NEW TECHNIQUE: Can we use the Ising model machinery
(transfer matrix, correlation functions, susceptibility)
to PROVE the crossover and possibly close Alon's n^{3/2} gap?

Author: opus-2026-03-12-S60b
"""
import sys
import time
import math
import numpy as np
from collections import defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def circulant_adj(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1
    return A


def paley_set(p):
    return frozenset(pow(x, 2, p) for x in range(1, p))


def interval_set(p):
    return frozenset(range(1, (p-1)//2 + 1))


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


def full_walsh_expansion(p):
    """Compute the complete Walsh-Fourier expansion of H on {±1}^m."""
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    N = 1 << m

    # Compute H for all orientations
    all_H = {}
    for bits in range(N):
        S = set()
        sigma = []
        for i in range(m):
            if bits & (1 << i):
                S.add(pairs[i][0])
                sigma.append(1)
            else:
                S.add(pairs[i][1])
                sigma.append(-1)
        A = circulant_adj(p, S)
        H = hamiltonian_paths_dp(A, p)
        all_H[tuple(sigma)] = H

    # Walsh transform: Ĥ(S) = (1/N) Σ_σ H(σ) χ_S(σ)
    walsh = {}
    for size in range(m + 1):
        for subset in combinations(range(m), size):
            coeff = 0.0
            for sigma, H in all_H.items():
                prod = 1
                for idx in subset:
                    prod *= sigma[idx]
                coeff += prod * H
            coeff /= N
            if abs(coeff) > 1e-6:
                walsh[subset] = coeff

    return walsh, all_H


def ising_analysis(p, walsh, all_H):
    """Analyze the Ising model structure."""
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]

    # Energy by degree
    energy_by_deg = defaultdict(float)
    coeffs_by_deg = defaultdict(list)
    for subset, coeff in walsh.items():
        deg = len(subset)
        energy_by_deg[deg] += coeff**2
        coeffs_by_deg[deg].append((subset, coeff))

    # Interaction matrix J (degree 2)
    J = np.zeros((m, m))
    for subset, coeff in walsh.items():
        if len(subset) == 2:
            i, j = subset
            J[i][j] = J[j][i] = coeff

    # Eigendecomposition of J
    eigenvalues, eigenvectors = np.linalg.eigh(J)

    # Paley and Interval orientations
    qr = paley_set(p)
    sigma_P = tuple(1 if pairs[i][0] in qr else -1 for i in range(m))
    sigma_I = tuple(1 if pairs[i][0] <= m else -1 for i in range(m))

    # Quadratic form values
    Q_P = sum(J[i][j] * sigma_P[i] * sigma_P[j] for i in range(m) for j in range(m))
    Q_I = sum(J[i][j] * sigma_I[i] * sigma_I[j] for i in range(m) for j in range(m))

    # Degree-4 contribution
    D4_P = sum(coeff * np.prod([sigma_P[i] for i in subset])
               for subset, coeff in walsh.items() if len(subset) == 4)
    D4_I = sum(coeff * np.prod([sigma_I[i] for i in subset])
               for subset, coeff in walsh.items() if len(subset) == 4)

    # Total H decomposition
    H0 = walsh.get((), 0)
    H_P = all_H[sigma_P]
    H_I = all_H[sigma_I]

    return {
        'energy_by_deg': dict(energy_by_deg),
        'J': J,
        'eigenvalues': eigenvalues,
        'Q_P': Q_P, 'Q_I': Q_I,
        'D4_P': D4_P, 'D4_I': D4_I,
        'H0': H0, 'H_P': H_P, 'H_I': H_I,
        'sigma_P': sigma_P, 'sigma_I': sigma_I,
        'coeffs_by_deg': dict(coeffs_by_deg),
    }


def main():
    print("ISING PHASE TRANSITION IN TOURNAMENT H-MAXIMIZATION")
    print("=" * 75)

    results = {}
    for p in [7, 11, 13]:
        m = (p - 1) // 2
        print(f"\n{'='*75}")
        print(f"p = {p}, m = {m} spins")
        print(f"{'='*75}")

        t0 = time.time()
        walsh, all_H = full_walsh_expansion(p)
        elapsed = time.time() - t0
        print(f"  Walsh expansion computed in {elapsed:.1f}s ({1 << m} configurations)")

        analysis = ising_analysis(p, walsh, all_H)
        results[p] = analysis

        # Energy spectrum
        total_var = sum(v for k, v in analysis['energy_by_deg'].items() if k > 0)
        print(f"\n  Energy spectrum (variance decomposition):")
        for deg in sorted(analysis['energy_by_deg'].keys()):
            E = analysis['energy_by_deg'][deg]
            if deg > 0 and total_var > 0:
                frac = E / total_var * 100
                print(f"    Degree {deg}: {E:>14.2f} ({frac:>5.1f}% of variance)")
            elif deg == 0:
                print(f"    Degree {deg}: {E:>14.2f} (mean²)")

        # Phase transition diagnostic
        print(f"\n  Ising diagnostics:")
        print(f"    J eigenvalues: {analysis['eigenvalues']}")
        print(f"    Q(Paley) = {analysis['Q_P']:.4f}")
        print(f"    Q(Interval) = {analysis['Q_I']:.4f}")
        print(f"    Q(P) > Q(I)? {analysis['Q_P'] > analysis['Q_I']}")

        print(f"    D4(Paley) = {analysis['D4_P']:.4f}")
        print(f"    D4(Interval) = {analysis['D4_I']:.4f}")

        print(f"\n    H decomposition:")
        print(f"    Paley: H₀={analysis['H0']:.2f}, +Q₂={analysis['Q_P']:.2f}, "
              f"+D₄={analysis['D4_P']:.2f} → H={analysis['H_P']}")
        print(f"    Interval: H₀={analysis['H0']:.2f}, +Q₂={analysis['Q_I']:.2f}, "
              f"+D₄={analysis['D4_I']:.2f} → H={analysis['H_I']}")

        winner = "PALEY" if analysis['H_P'] > analysis['H_I'] else "INTERVAL"
        margin = analysis['H_P'] - analysis['H_I']
        print(f"    Winner: {winner} by {abs(margin)}")

    # Phase transition summary
    print(f"\n{'='*75}")
    print("PHASE TRANSITION SUMMARY")
    print("=" * 75)

    print(f"\n  {'p':>4} {'m':>4} {'deg-2 %':>10} {'Q(P)':>12} {'Q(I)':>12} {'D4(P)':>10} {'D4(I)':>10} {'Winner':>8}")
    for p in sorted(results.keys()):
        r = results[p]
        total_var = sum(v for k, v in r['energy_by_deg'].items() if k > 0)
        deg2_pct = r['energy_by_deg'].get(2, 0) / total_var * 100 if total_var > 0 else 0
        winner = "P" if r['H_P'] > r['H_I'] else "I"
        print(f"  {p:>4} {(p-1)//2:>4} {deg2_pct:>10.1f}% {r['Q_P']:>12.2f} {r['Q_I']:>12.2f} "
              f"{r['D4_P']:>10.2f} {r['D4_I']:>10.2f} {winner:>8}")

    # The key prediction
    print(f"\n  PREDICTION: The phase transition occurs when degree-2 energy")
    print(f"  fraction drops below ~80% and degree-4+ terms begin favoring Interval.")
    print(f"  At p=7: 100% degree-2 → Paley wins (pure Ising)")
    print(f"  At p=11: 84% degree-2 → Paley still wins")
    print(f"  At p=13: ??? → check (p ≡ 1 mod 4, no Paley)")
    print(f"  At p=19: degree-6+ terms dominate → Interval wins")

    # Correlation functions
    print(f"\n{'='*75}")
    print("ISING CORRELATION FUNCTIONS: ⟨σ_i σ_j⟩ vs distance")
    print("=" * 75)

    for p in [7, 11]:
        m = (p - 1) // 2
        r = results[p]
        walsh, all_H = full_walsh_expansion(p)

        # Two-point correlations at the ground state
        print(f"\n  p={p}: Two-point correlations at Paley σ_P = {r['sigma_P']}")
        for i in range(m):
            for j in range(i + 1, m):
                corr = r['sigma_P'][i] * r['sigma_P'][j]
                J_ij = r['J'][i][j]
                chord_dist = min(abs(i - j), m - abs(i - j))
                print(f"    ⟨σ_{i+1} σ_{j+1}⟩ = {corr:>+3}, J[{i+1},{j+1}] = {J_ij:>+8.2f}, "
                      f"chord dist = {chord_dist}")

    # The transfer matrix approach
    print(f"\n{'='*75}")
    print("TRANSFER MATRIX APPROACH (sketch)")
    print("=" * 75)

    print("""
  In 1D Ising models, the partition function is computed via transfer matrices:
    Z = Tr(T^N) where T is the 2×2 transfer matrix

  For our model:
    - Spins are chord types 1, ..., m
    - Couplings J[i,j] are NOT nearest-neighbor (they depend on chord geometry)
    - The model is NOT 1D — it's on the complete graph of chords

  However, the QR symmetry gives a SIMPLIFICATION:
    Under QR multiplication, all chords are equivalent (single orbit for p ≡ 3 mod 4)
    So J has structure: J[i,j] depends only on the QR orbit of (i,j)

  For p=7 (m=3):
    J = [[ 0, 3.5, -3.5],
         [ 3.5, 0, -3.5],
         [-3.5, -3.5, 0]]

    The QR group {1,2,4} acts on chords {1,2,3} as:
      1: (1,2,3) → (1,2,3)  (identity)
      2: (1,2,3) → (2,3,1)  (cyclic shift)
      4: (1,2,3) → (3,1,2)  (cyclic shift)

    So J should be circulant on 3 elements...
    J[0,1]=3.5, J[0,2]=-3.5 → NOT circulant (asymmetric)

    But J IS equivariant: J[π(i),π(j)] = J[i,j] for π ∈ QR

    The sign pattern of J encodes which chord pairs "cooperate":
      Chords 1,2 cooperate (J>0): both forward → more HP
      Chords 1,3 and 2,3 anti-cooperate (J<0)

    Paley σ = (1,1,-1) puts chords 1,2 forward and 3 backward.
    This maximizes cooperation: 3.5 + 3.5 + 3.5 = 10.5 (from σ^T J σ / 2)
    """)

    # The susceptibility and specific heat
    print(f"{'='*75}")
    print("SUSCEPTIBILITY AND SPECIFIC HEAT")
    print("=" * 75)

    for p in [7, 11]:
        m = (p - 1) // 2
        r = results[p]
        N = 1 << m

        # Compute "thermal" averages at different "temperatures" β
        # Z(β) = Σ_σ exp(β H(σ))
        # ⟨H⟩_β = (1/Z) Σ H(σ) exp(β H(σ))
        # C_v = β² [⟨H²⟩ - ⟨H⟩²]

        walsh, all_H = full_walsh_expansion(p)
        H_list = list(all_H.values())

        print(f"\n  p={p}: Thermal analysis")
        print(f"  {'β':>8} {'⟨H⟩':>12} {'Var(H)':>12} {'entropy':>10} {'ground σ':>15}")

        for beta in [0, 0.001, 0.01, 0.1, 1.0, 10.0]:
            # Compute Z and averages
            log_weights = [beta * H for H in H_list]
            max_lw = max(log_weights)
            weights = [math.exp(lw - max_lw) for lw in log_weights]
            Z = sum(weights)

            avg_H = sum(w * H for w, H in zip(weights, H_list)) / Z
            avg_H2 = sum(w * H**2 for w, H in zip(weights, H_list)) / Z
            var_H = avg_H2 - avg_H**2

            # Entropy
            probs = [w / Z for w in weights]
            entropy = -sum(p_i * math.log(p_i + 1e-300) for p_i in probs)

            # Ground state at this temperature
            idx_max = max(range(N), key=lambda i: log_weights[i])
            sigmas_list = list(all_H.keys())
            ground = sigmas_list[idx_max]

            print(f"  {beta:>8.3f} {avg_H:>12.2f} {var_H:>12.2f} {entropy:>10.4f} {str(ground):>15}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
