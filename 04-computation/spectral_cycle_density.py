import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
spectral_cycle_density.py
kind-pasteur-2026-03-07-S39b

Investigate the spectral/eigenvalue connection to H-maximization
and cycle density, following INV-055 (Linial-Morgenstern).

Key questions:
1. Does the spectral radius lambda_1(A) correlate with H(T)?
2. Does Paley T_7 have extremal spectral properties?
3. Can eigenvalues predict cycle densities?

Also: Newton inequality margin analysis for root gap investigation (INV-065).
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from tournament_fast import c3_from_score, c4_fast, c5_fast, alpha2_from_trace
from math import comb
import random

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("WARNING: numpy not available. Skipping spectral analysis.")


def cycle_densities(T):
    """Compute normalized cycle densities d_k = c_k / C(n, k)."""
    n = len(T)
    c3 = c3_from_score(T)
    c4 = c4_fast(T)
    c5 = c5_fast(T) if n >= 5 else 0
    return {
        3: c3 / comb(n, 3) if n >= 3 else 0,
        4: c4 / comb(n, 4) if n >= 4 else 0,
        5: c5 / comb(n, 5) if n >= 5 else 0,
    }


if HAS_NUMPY:
    # ============================================================
    # 1. Spectral analysis at n=7
    # ============================================================
    print("=" * 70)
    print("SPECTRAL ANALYSIS: eigenvalues vs H(T) at n=7")
    print("=" * 70)

    n = 7
    m = n * (n - 1) // 2

    # Paley T_7
    T7 = [[0]*7 for _ in range(7)]
    QR = {1, 2, 4}
    for i in range(7):
        for j in range(7):
            if i != j and (j - i) % 7 in QR:
                T7[i][j] = 1

    A7 = np.array(T7, dtype=float)
    evals_paley = sorted(np.linalg.eigvals(A7), key=lambda x: -abs(x))
    h_paley = hamiltonian_path_count(T7)

    print(f"\nPaley T_7:")
    print(f"  H = {h_paley}")
    print(f"  Eigenvalues: {[f'{e:.3f}' for e in evals_paley]}")
    print(f"  Spectral radius |lambda_1| = {abs(evals_paley[0]):.4f}")
    print(f"  c3={c3_from_score(T7)}, c4={c4_fast(T7)}, c5={c5_fast(T7)}")

    # Compare with random tournaments
    random.seed(42)
    sample = [random.randint(0, (1 << m) - 1) for _ in range(500)]

    data = []
    for bits in sample:
        T = tournament_from_bits(n, bits)
        A = np.array(T, dtype=float)
        evals = np.linalg.eigvals(A)
        spec_rad = max(abs(e) for e in evals)
        h = hamiltonian_path_count(T)
        c3 = c3_from_score(T)
        c5 = c5_fast(T)
        alpha2 = alpha2_from_trace(T)
        data.append({
            'bits': bits, 'H': h, 'spec_rad': spec_rad,
            'c3': c3, 'c5': c5, 'alpha2': alpha2,
            'evals': sorted(evals, key=lambda x: -abs(x))
        })

    # Correlation: spectral radius vs H
    spec_rads = [d['spec_rad'] for d in data]
    hs = [d['H'] for d in data]
    corr = np.corrcoef(spec_rads, hs)[0, 1]
    print(f"\n  Correlation(spectral radius, H): {corr:.4f}")

    # Does Paley have max spectral radius?
    max_spec_rad = max(spec_rads)
    paley_spec_rad = abs(evals_paley[0])
    print(f"  Paley spectral radius: {paley_spec_rad:.4f}")
    print(f"  Max spectral radius in sample: {max_spec_rad:.4f}")
    print(f"  Paley is max: {abs(paley_spec_rad - max_spec_rad) < 0.01}")

    # Top 10 by H
    top_h = sorted(data, key=lambda d: -d['H'])[:10]
    print(f"\n  Top 10 by H:")
    print(f"  {'H':>5} {'spec_rad':>10} {'c3':>4} {'c5':>4} {'alpha2':>7}")
    for d in top_h:
        print(f"  {d['H']:>5} {d['spec_rad']:>10.4f} {d['c3']:>4} {d['c5']:>4} {d['alpha2']:>7}")

    # Bottom 10 by H
    bot_h = sorted(data, key=lambda d: d['H'])[:10]
    print(f"\n  Bottom 10 by H:")
    for d in bot_h:
        print(f"  {d['H']:>5} {d['spec_rad']:>10.4f} {d['c3']:>4} {d['c5']:>4} {d['alpha2']:>7}")

    # ============================================================
    # 2. Cycle density comparison: Paley vs others
    # ============================================================
    print("\n" + "=" * 70)
    print("CYCLE DENSITY: Paley vs random vs transitive")
    print("=" * 70)

    # Transitive tournament
    T7_trans = [[0]*7 for _ in range(7)]
    for i in range(7):
        for j in range(i+1, 7):
            T7_trans[i][j] = 1

    paley_dens = cycle_densities(T7)
    trans_dens = cycle_densities(T7_trans)

    print(f"\n  Paley T_7 densities: d3={paley_dens[3]:.4f}, d4={paley_dens[4]:.4f}, d5={paley_dens[5]:.4f}")
    print(f"  Transitive densities: d3={trans_dens[3]:.4f}, d4={trans_dens[4]:.4f}, d5={trans_dens[5]:.4f}")

    # Average random tournament
    avg_d3 = np.mean([cycle_densities(tournament_from_bits(n, d['bits']))[3] for d in data])
    avg_d5 = np.mean([cycle_densities(tournament_from_bits(n, d['bits']))[5] for d in data])
    print(f"  Random avg densities: d3={avg_d3:.4f}, d5={avg_d5:.4f}")
    print(f"  Expected d3 for random: {1/4:.4f} (1/4 = each triple is a 3-cycle with prob 1/4)")

    # Does Paley maximize ALL cycle densities?
    max_d3 = max(cycle_densities(tournament_from_bits(n, d['bits']))[3] for d in data)
    max_d5 = max(cycle_densities(tournament_from_bits(n, d['bits']))[5] for d in data)
    print(f"\n  Paley d3={paley_dens[3]:.4f}, max d3 in sample={max_d3:.4f}")
    print(f"  Paley d5={paley_dens[5]:.4f}, max d5 in sample={max_d5:.4f}")
    print(f"  Paley maximizes d3: {abs(paley_dens[3] - max_d3) < 0.001}")
    print(f"  Paley maximizes d5: {abs(paley_dens[5] - max_d5) < 0.001}")

    # ============================================================
    # 3. Skew-adjacency matrix eigenvalues
    # ============================================================
    print("\n" + "=" * 70)
    print("SKEW-ADJACENCY EIGENVALUES: S = A - A^T")
    print("=" * 70)

    # S is skew-symmetric, so eigenvalues are purely imaginary
    S7 = A7 - A7.T
    evals_S = np.linalg.eigvals(S7)
    imag_parts = sorted([e.imag for e in evals_S], reverse=True)
    print(f"\n  Paley T_7 skew eigenvalues (imaginary parts): {[f'{x:.3f}' for x in imag_parts]}")
    print(f"  Max imaginary part = {max(imag_parts):.4f}")

    # For regular tournaments: S eigenvalues relate to Paley character sums
    # For Paley T_7: eigenvalues should be related to Gauss sums

    # Compare with random
    max_imag_sample = max(
        max(e.imag for e in np.linalg.eigvals(np.array(tournament_from_bits(n, d['bits']), dtype=float) -
            np.array(tournament_from_bits(n, d['bits']), dtype=float).T))
        for d in data[:100]
    )
    print(f"  Max imaginary eigenvalue in sample: {max_imag_sample:.4f}")

    # ============================================================
    # 4. Newton inequality margin analysis (root gap, INV-065)
    # ============================================================
    print("\n" + "=" * 70)
    print("NEWTON INEQUALITY MARGINS for I(Omega(T), x)")
    print("=" * 70)

    for n in [7, 8]:
        m = n * (n - 1) // 2
        random.seed(42)
        if n <= 7:
            sample_bits = [random.randint(0, (1 << m) - 1) for _ in range(200)]
        else:
            sample_bits = [random.randint(0, (1 << m) - 1) for _ in range(100)]

        print(f"\n  n={n}:")

        margins = []
        for bits in sample_bits:
            T = tournament_from_bits(n, bits)
            nn = len(T)

            # Compute independence polynomial coefficients of Omega(T)
            # For n<=8, alpha_max = floor(n/3) = 2 (n=7,8)
            c3 = c3_from_score(T)
            c5 = c5_fast(T) if nn >= 5 else 0
            alpha_1 = c3 + c5
            alpha_2 = alpha2_from_trace(T)

            # I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2
            # (at n=7,8: alpha_3 = 0 for 3-cycles; need to check for n=8 with 5+3 pairs)
            # For now, assume alpha_3 = 0 at n <= 8.

            coeffs = [1, alpha_1, alpha_2]

            # Newton inequality at k=1: alpha_1^2 >= 2 * alpha_0 * alpha_2
            # = alpha_1^2 >= 2 * alpha_2
            if alpha_2 > 0:
                margin = alpha_1 ** 2 - 2 * alpha_2
                ratio = alpha_1 ** 2 / (2 * alpha_2) if alpha_2 > 0 else float('inf')
                margins.append({'margin': margin, 'ratio': ratio, 'a1': alpha_1, 'a2': alpha_2, 'bits': bits})

        if margins:
            min_margin = min(m['margin'] for m in margins)
            min_ratio = min(m['ratio'] for m in margins)
            worst = min(margins, key=lambda m: m['ratio'])
            print(f"    Newton inequality: alpha_1^2 >= 2*alpha_2")
            print(f"    Min margin: {min_margin} (always >= 0 for real roots)")
            print(f"    Min ratio alpha_1^2/(2*alpha_2): {min_ratio:.4f}")
            print(f"    Worst case: alpha_1={worst['a1']}, alpha_2={worst['a2']}, "
                  f"margin={worst['margin']}, ratio={worst['ratio']:.4f}")

            # At what alpha_1/alpha_2 ratio does real-rootedness barely hold?
            # For quadratic ax^2 + bx + c, discriminant = b^2 - 4ac.
            # I(x) = alpha_2*x^2 + alpha_1*x + 1
            # Disc = alpha_1^2 - 4*alpha_2 > 0 for real roots
            disc_margins = []
            for m in margins:
                disc = m['a1']**2 - 4*m['a2']
                disc_margins.append(disc)

            min_disc = min(disc_margins)
            print(f"    Discriminant margin: min = {min_disc} "
                  f"(need > 0 for real roots)")
            print(f"    All discriminants > 0: {all(d > 0 for d in disc_margins)}")

    # ============================================================
    # 5. Root gap at n=9 (approach to failure)
    # ============================================================
    print("\n" + "=" * 70)
    print("ROOT GAP AT n=9: approaching complex roots")
    print("=" * 70)

    n = 9
    m = n * (n - 1) // 2
    random.seed(42)
    sample_bits = [random.randint(0, (1 << m) - 1) for _ in range(200)]

    real_root_count = 0
    complex_root_count = 0
    min_disc_real = float('inf')

    for bits in sample_bits:
        T = tournament_from_bits(n, bits)
        nn = len(T)

        # Compute alpha_1, alpha_2, alpha_3
        c3 = c3_from_score(T)
        c5 = c5_fast(T)
        alpha_1 = c3 + c5  # No 7- or 9-cycles contribute to alpha_1?
        # At n=9: 7-cycles exist. alpha_1 = c3 + c5 + c7 + c9.
        # But computing c7 is expensive. For the Newton inequality,
        # the degree of I(Omega, x) is alpha(Omega) = max independent set.
        # At n=9: max = floor(9/3) = 3 for 3-cycles.
        # But we also have 5-cycles. A 5-cycle + nothing disjoint at n=9.
        # Actually 5+3 = 8 < 9, so a 5-cycle and a disjoint 3-cycle IS possible at n=9!
        # And 3+3+3 = 9, so three disjoint 3-cycles are possible at n=9.
        # So alpha_max could be 3.

        # For now, approximate with alpha_2 from 3-cycles only (OK at n<=7)
        alpha_2 = alpha2_from_trace(T)

        # Need c7 for alpha_1... skip for now
        # Just check if alpha_1, alpha_2 satisfy Newton at degree 2
        if alpha_2 > 0:
            disc = alpha_1**2 - 4*alpha_2
            if disc > 0:
                real_root_count += 1
                if disc < min_disc_real:
                    min_disc_real = disc
            else:
                complex_root_count += 1

    print(f"  n=9 (200 sampled, using c3+c5 only for alpha_1, 3-cycle pairs for alpha_2):")
    print(f"    Real roots (disc > 0): {real_root_count}")
    print(f"    Complex roots (disc <= 0): {complex_root_count}")
    print(f"    Min discriminant among real: {min_disc_real}")

else:
    print("Skipping spectral analysis (numpy not available)")
    print("Install numpy: pip install numpy")


print("\nDone.")
