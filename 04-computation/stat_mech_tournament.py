#!/usr/bin/env python3
"""
stat_mech_tournament.py — kind-pasteur-2026-03-14

Statistical mechanics of tournaments via the partition function H(T) = I(Omega(T), 2).

KEY FRAMEWORK:
  H(T) = I(Omega(T), 2) is the partition function of the hard-core lattice gas
  on the conflict graph Omega(T) at fugacity lambda=2.

EXPLORATIONS:
  1. Boltzmann distribution P(T) = exp(beta * H(T)) / Z(beta) on tournaments
  2. Free energy F(beta) = -log(Z(beta)) / beta
  3. Specific heat C(beta) = beta^2 * Var(H)
  4. Order parameter: score variance
  5. Susceptibility chi = d<H>/d(beta)

Exhaustive for n=5 (1024 tournaments), sampled for n=6 (10000 random).
"""

import numpy as np
from itertools import combinations
from math import comb, factorial
from collections import defaultdict
import random
import sys

random.seed(42)
np.random.seed(42)


# ─────────────────────────────────────────────
#  Core tournament utilities
# ─────────────────────────────────────────────

def all_tournaments(n):
    """Enumerate all 2^C(n,2) tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield mask, A


def random_tournament(n):
    """Generate a uniformly random tournament on n vertices."""
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def scores(A, n):
    """Out-degree sequence of tournament."""
    return sorted([sum(A[i]) for i in range(n)])


def score_variance(A, n):
    """Variance of the score sequence."""
    s = [sum(A[i]) for i in range(n)]
    mean = sum(s) / n
    return sum((x - mean) ** 2 for x in s) / n


def count_directed_3cycles(A, n):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                # Check both orientations of the 3-cycle on {i, j, k}
                if A[i][j] and A[j][k] and A[k][i]:
                    count += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    count += 1
    return count


def compute_H(A, n):
    """
    Compute H(T) = number of directed Hamiltonian paths in T.
    Uses DP (Held-Karp) for efficiency: O(n^2 * 2^n).

    H(T) = I(Omega(T), 2) by the OCF theorem (THM-002),
    where Omega has ALL directed odd cycles as vertices.
    Direct path counting is simpler and faster for small n.
    """
    # dp[mask][v] = number of Hamiltonian paths ending at vertex v
    # using exactly the vertices in mask
    dp = [[0] * n for _ in range(1 << n)]

    # Base: paths of length 1
    for v in range(n):
        dp[1 << v][v] = 1

    # Fill DP
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

    full_mask = (1 << n) - 1
    return sum(dp[full_mask][v] for v in range(n))


# ─────────────────────────────────────────────
#  Statistical mechanics computations
# ─────────────────────────────────────────────

def stat_mech_analysis(n, H_values, sv_values, c3_values, beta_range, label=""):
    """
    Given arrays of H(T) and score_variance(T) for a collection of tournaments,
    compute statistical mechanics quantities as function of beta.
    """
    H_arr = np.array(H_values, dtype=np.float64)
    sv_arr = np.array(sv_values, dtype=np.float64)
    c3_arr = np.array(c3_values, dtype=np.float64)

    results = []
    print(f"\n{'='*80}")
    print(f"  STATISTICAL MECHANICS OF TOURNAMENTS — n={n} {label}")
    print(f"{'='*80}")
    print(f"  Number of tournaments: {len(H_arr)}")
    print(f"  H range: [{int(H_arr.min())}, {int(H_arr.max())}]")
    print(f"  Mean H: {H_arr.mean():.4f}")
    print(f"  H std: {H_arr.std():.4f}")
    print(f"  Score variance range: [{sv_arr.min():.4f}, {sv_arr.max():.4f}]")
    print(f"  3-cycle count range: [{int(c3_arr.min())}, {int(c3_arr.max())}]")

    # Distribution of H values
    H_counts = defaultdict(int)
    for h in H_arr:
        H_counts[int(h)] += 1
    print(f"\n  H value distribution:")
    for h_val in sorted(H_counts.keys()):
        count = H_counts[h_val]
        frac = count / len(H_arr)
        print(f"    H={h_val:5d}: {count:6d} tournaments ({frac*100:6.2f}%)")

    print(f"\n  {'beta':>8s} {'log_Z':>12s} {'<H>':>10s} {'<H^2>':>12s} "
          f"{'C(beta)':>10s} {'<sv>':>10s} {'<c3>':>10s} {'chi':>10s}")
    print(f"  {'-'*8} {'-'*12} {'-'*10} {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")

    for beta in beta_range:
        # Use log-sum-exp for numerical stability
        log_weights = beta * H_arr
        log_Z = np.max(log_weights) + np.log(np.sum(np.exp(log_weights - np.max(log_weights))))

        # Boltzmann weights (normalized)
        log_p = log_weights - log_Z
        p = np.exp(log_p)

        # Expectation values
        mean_H = np.sum(p * H_arr)
        mean_H2 = np.sum(p * H_arr ** 2)
        var_H = mean_H2 - mean_H ** 2
        specific_heat = beta ** 2 * var_H if beta != 0 else var_H

        mean_sv = np.sum(p * sv_arr)
        mean_c3 = np.sum(p * c3_arr)

        # Susceptibility: chi = d<H>/d(beta) = Var(H)
        chi = var_H

        results.append({
            'beta': beta,
            'log_Z': log_Z,
            'mean_H': mean_H,
            'mean_H2': mean_H2,
            'var_H': var_H,
            'specific_heat': specific_heat,
            'mean_sv': mean_sv,
            'mean_c3': mean_c3,
            'chi': chi,
        })

        print(f"  {beta:8.2f} {log_Z:12.4f} {mean_H:10.4f} {mean_H2:12.2f} "
              f"{specific_heat:10.4f} {mean_sv:10.6f} {mean_c3:10.4f} {chi:10.4f}")

    return results


def detect_phase_transition(results):
    """
    Detect phase transitions by finding peaks in specific heat
    and rapid changes in the order parameter.
    """
    betas = [r['beta'] for r in results]
    C_vals = [r['specific_heat'] for r in results]
    sv_vals = [r['mean_sv'] for r in results]
    chi_vals = [r['chi'] for r in results]

    # Find peak of specific heat
    C_max_idx = np.argmax(C_vals)
    C_max_beta = betas[C_max_idx]
    C_max_val = C_vals[C_max_idx]

    # Find peak of susceptibility
    chi_max_idx = np.argmax(chi_vals)
    chi_max_beta = betas[chi_max_idx]
    chi_max_val = chi_vals[chi_max_idx]

    # Compute derivative of <score_var> w.r.t. beta (numerical)
    dsv_dbeta = np.gradient(sv_vals, [r['beta'] for r in results])
    dsv_max_idx = np.argmax(np.abs(dsv_dbeta))
    dsv_max_beta = betas[dsv_max_idx]

    print(f"\n  PHASE TRANSITION ANALYSIS:")
    print(f"  {'='*60}")
    print(f"  Specific heat peak: C_max = {C_max_val:.4f} at beta_c = {C_max_beta:.2f}")
    print(f"  Susceptibility peak: chi_max = {chi_max_val:.4f} at beta = {chi_max_beta:.2f}")
    print(f"  Max |d<sv>/d(beta)| at beta = {dsv_max_beta:.2f}")

    # Check if specific heat peak is sharp (ratio of peak to average)
    C_avg = np.mean(C_vals)
    sharpness = C_max_val / C_avg if C_avg > 0 else 0
    print(f"  Specific heat sharpness (peak/avg): {sharpness:.4f}")

    if sharpness > 3:
        print(f"  ** STRONG phase transition signal (sharpness > 3) **")
    elif sharpness > 2:
        print(f"  ** Moderate phase transition signal (sharpness > 2) **")
    else:
        print(f"  ** Weak or no phase transition (sharpness < 2) **")

    # Report order parameter at various regimes
    print(f"\n  ORDER PARAMETER <score_var> at key betas:")
    for r in results:
        if r['beta'] in [-2, -1, 0, 1, 2, 3, 4, 5]:
            print(f"    beta = {r['beta']:5.1f}: <sv> = {r['mean_sv']:.6f}, "
                  f"<H> = {r['mean_H']:.4f}, <c3> = {r['mean_c3']:.4f}")

    return C_max_beta, C_max_val, chi_max_beta, chi_max_val


def boltzmann_analysis(n, H_values, sv_values, c3_values, scores_list):
    """
    Additional analysis of the Boltzmann distribution.
    """
    H_arr = np.array(H_values, dtype=np.float64)
    sv_arr = np.array(sv_values, dtype=np.float64)
    c3_arr = np.array(c3_values, dtype=np.float64)

    print(f"\n  BOLTZMANN DISTRIBUTION ANALYSIS (n={n}):")
    print(f"  {'='*60}")

    # Natural Boltzmann weight: P(T) = H(T) / sum H(T)
    # This corresponds to beta=1 with the log-partition function
    Z_natural = np.sum(H_arr)
    P_natural = H_arr / Z_natural

    print(f"  Natural partition function Z_1 = sum H(T) = {Z_natural:.0f}")
    print(f"  Expected value: n! * 2^C(n,2) / 2^(n-1) = "
          f"{factorial(n) * 2**comb(n,2) / 2**(n-1):.0f}")
    print(f"    (This is sum_T H(T) = sum_T |ham paths| from double counting)")

    # Boltzmann-weighted properties at beta=1
    mean_H_boltz = np.sum(P_natural * H_arr)
    mean_sv_boltz = np.sum(P_natural * sv_arr)
    mean_c3_boltz = np.sum(P_natural * c3_arr)

    # Uniform properties
    mean_H_uniform = np.mean(H_arr)
    mean_sv_uniform = np.mean(sv_arr)
    mean_c3_uniform = np.mean(c3_arr)

    print(f"\n  Comparison: Uniform vs Boltzmann-weighted (beta=1):")
    print(f"  {'Property':>20s} {'Uniform':>12s} {'Boltzmann':>12s} {'Ratio':>10s}")
    print(f"  {'-'*20} {'-'*12} {'-'*12} {'-'*10}")
    print(f"  {'<H>':>20s} {mean_H_uniform:12.4f} {mean_H_boltz:12.4f} "
          f"{mean_H_boltz/mean_H_uniform:10.4f}")
    print(f"  {'<score_var>':>20s} {mean_sv_uniform:12.6f} {mean_sv_boltz:12.6f} "
          f"{mean_sv_boltz/mean_sv_uniform:10.4f}")
    print(f"  {'<c3>':>20s} {mean_c3_uniform:12.4f} {mean_c3_boltz:12.4f} "
          f"{mean_c3_boltz/mean_c3_uniform:10.4f}")

    # Entropy of the Boltzmann distribution
    S_boltz = -np.sum(P_natural[P_natural > 0] * np.log(P_natural[P_natural > 0]))
    S_uniform = np.log(len(H_arr))
    print(f"\n  Entropy: S_Boltzmann = {S_boltz:.4f}, S_uniform = {S_uniform:.4f}")
    print(f"  Effective number of states: exp(S) = {np.exp(S_boltz):.1f} "
          f"(out of {len(H_arr)})")
    print(f"  Information content: S_uniform - S_Boltzmann = {S_uniform - S_boltz:.4f} nats")

    # Most probable tournament type
    # Group by score sequence
    score_groups = defaultdict(list)
    for i, s in enumerate(scores_list):
        score_groups[tuple(s)].append(i)

    print(f"\n  Score sequence analysis under Boltzmann distribution:")
    score_probs = {}
    for s_seq, indices in sorted(score_groups.items()):
        p_total = sum(P_natural[i] for i in indices)
        avg_H = np.mean([H_arr[i] for i in indices])
        score_probs[s_seq] = p_total

    # Sort by probability
    for s_seq, prob in sorted(score_probs.items(), key=lambda x: -x[1])[:8]:
        indices = score_groups[s_seq]
        avg_H = np.mean([H_arr[i] for i in indices])
        avg_sv = np.mean([sv_arr[i] for i in indices])
        print(f"    scores={list(s_seq)}: P={prob:.4f}, count={len(indices)}, "
              f"<H>={avg_H:.1f}, <sv>={avg_sv:.4f}")


def correlation_analysis(H_values, sv_values, c3_values, n):
    """Analyze correlations between H, score_variance, and cycle counts."""
    H_arr = np.array(H_values, dtype=np.float64)
    sv_arr = np.array(sv_values, dtype=np.float64)
    c3_arr = np.array(c3_values, dtype=np.float64)

    print(f"\n  CORRELATION ANALYSIS (n={n}):")
    print(f"  {'='*60}")

    # Pearson correlations
    corr_H_sv = np.corrcoef(H_arr, sv_arr)[0, 1]
    corr_H_c3 = np.corrcoef(H_arr, c3_arr)[0, 1]
    corr_sv_c3 = np.corrcoef(sv_arr, c3_arr)[0, 1]

    print(f"  Pearson correlations:")
    print(f"    corr(H, score_var) = {corr_H_sv:+.6f}")
    print(f"    corr(H, c3)        = {corr_H_c3:+.6f}")
    print(f"    corr(score_var, c3) = {corr_sv_c3:+.6f}")

    # H vs score_var: key relationship
    # Regular tournaments (sv=0 at odd n) should have high H
    if n % 2 == 1:
        regular_mask = sv_arr == 0
        if np.any(regular_mask):
            H_regular = H_arr[regular_mask]
            H_non_regular = H_arr[~regular_mask]
            print(f"\n  Regular vs non-regular tournaments:")
            print(f"    Regular: count={np.sum(regular_mask)}, "
                  f"<H>={H_regular.mean():.4f}, max H={int(H_regular.max())}")
            print(f"    Non-regular: count={np.sum(~regular_mask)}, "
                  f"<H>={H_non_regular.mean():.4f}, max H={int(H_non_regular.max())}")


# ─────────────────────────────────────────────
#  Main execution
# ─────────────────────────────────────────────

def run_exhaustive(n):
    """Run exhaustive analysis for small n."""
    print(f"\n{'#'*80}")
    print(f"#  EXHAUSTIVE ANALYSIS: n={n} ({2**comb(n,2)} tournaments)")
    print(f"{'#'*80}")

    H_values = []
    sv_values = []
    c3_values = []
    scores_list = []

    total = 2 ** comb(n, 2)
    for idx, (mask, A) in enumerate(all_tournaments(n)):
        if idx % 200 == 0:
            print(f"  Processing tournament {idx}/{total}...", file=sys.stderr)

        H = compute_H(A, n)
        sv = score_variance(A, n)
        c3 = count_directed_3cycles(A, n)
        s = scores(A, n)

        H_values.append(H)
        sv_values.append(sv)
        c3_values.append(c3)
        scores_list.append(s)

    # Verify sum of H
    sum_H = sum(H_values)
    expected_sum = factorial(n) * 2 ** comb(n, 2) // 2 ** (n - 1)
    print(f"\n  VERIFICATION: sum H(T) = {sum_H}")
    print(f"  Expected: n! * 2^C(n,2) / 2^(n-1) = {expected_sum}")
    print(f"  Match: {sum_H == expected_sum}")

    # Beta range
    beta_range = np.arange(-2.0, 5.05, 0.1)

    # Stat mech analysis
    results = stat_mech_analysis(n, H_values, sv_values, c3_values,
                                  beta_range, label="(EXHAUSTIVE)")

    # Phase transition detection
    detect_phase_transition(results)

    # Boltzmann analysis
    boltzmann_analysis(n, H_values, sv_values, c3_values, scores_list)

    # Correlation analysis
    correlation_analysis(H_values, sv_values, c3_values, n)

    return results


def run_sampled(n, num_samples=10000):
    """Run sampled analysis for larger n."""
    print(f"\n{'#'*80}")
    print(f"#  SAMPLED ANALYSIS: n={n} ({num_samples} random tournaments)")
    print(f"{'#'*80}")

    H_values = []
    sv_values = []
    c3_values = []
    scores_list = []

    for idx in range(num_samples):
        if idx % 1000 == 0:
            print(f"  Sampling tournament {idx}/{num_samples}...", file=sys.stderr)

        A = random_tournament(n)
        H = compute_H(A, n)
        sv = score_variance(A, n)
        c3 = count_directed_3cycles(A, n)
        s = scores(A, n)

        H_values.append(H)
        sv_values.append(sv)
        c3_values.append(c3)
        scores_list.append(s)

    # Beta range
    beta_range = np.arange(-2.0, 5.05, 0.1)

    # Stat mech analysis
    results = stat_mech_analysis(n, H_values, sv_values, c3_values,
                                  beta_range, label=f"(SAMPLED, {num_samples} tournaments)")

    # Phase transition detection
    detect_phase_transition(results)

    # Boltzmann analysis
    boltzmann_analysis(n, H_values, sv_values, c3_values, scores_list)

    # Correlation analysis
    correlation_analysis(H_values, sv_values, c3_values, n)

    return results


def summary_comparison(results_5, results_6):
    """Compare phase transition behavior across n=5 and n=6."""
    print(f"\n{'#'*80}")
    print(f"#  CROSS-n COMPARISON")
    print(f"{'#'*80}")

    for n_label, results in [("n=5", results_5), ("n=6", results_6)]:
        C_vals = [r['specific_heat'] for r in results]
        betas = [r['beta'] for r in results]
        C_max_idx = np.argmax(C_vals)
        print(f"\n  {n_label}:")
        print(f"    C_max = {C_vals[C_max_idx]:.4f} at beta = {betas[C_max_idx]:.2f}")
        print(f"    <H> range: [{min(r['mean_H'] for r in results):.2f}, "
              f"{max(r['mean_H'] for r in results):.2f}]")
        print(f"    <sv> range: [{min(r['mean_sv'] for r in results):.6f}, "
              f"{max(r['mean_sv'] for r in results):.6f}]")

    # Check if beta_c shifts with n
    beta_c_5 = [r['beta'] for r in results_5][np.argmax([r['specific_heat'] for r in results_5])]
    beta_c_6 = [r['beta'] for r in results_6][np.argmax([r['specific_heat'] for r in results_6])]
    print(f"\n  Specific heat peak locations:")
    print(f"    n=5: beta_c = {beta_c_5:.2f}")
    print(f"    n=6: beta_c = {beta_c_6:.2f}")
    if abs(beta_c_5 - beta_c_6) < 0.3:
        print(f"    Peak locations are CLOSE — suggests universal beta_c")
    else:
        print(f"    Peak locations DIFFER — beta_c may shift with n")

    # Free energy comparison at selected betas
    print(f"\n  Free energy F(beta) = -log Z(beta) / beta:")
    print(f"  {'beta':>8s} {'F(n=5)':>12s} {'F(n=6)':>12s}")
    print(f"  {'-'*8} {'-'*12} {'-'*12}")
    for r5, r6 in zip(results_5, results_6):
        if abs(r5['beta'] - round(r5['beta'])) < 0.01 and r5['beta'] != 0:
            F5 = -r5['log_Z'] / r5['beta']
            F6 = -r6['log_Z'] / r6['beta']
            print(f"  {r5['beta']:8.1f} {F5:12.4f} {F6:12.4f}")


def finite_size_scaling(all_results):
    """
    Finite-size scaling analysis: how does the specific heat peak scale with n?
    """
    print(f"\n  FINITE-SIZE SCALING:")
    print(f"  {'='*60}")
    print(f"  If C_max ~ n^alpha, a first-order transition has alpha = d (dimension).")
    print(f"  For tournaments, the 'volume' scales as C(n,2) ~ n^2.")

    ns = sorted(all_results.keys())
    if len(ns) >= 2:
        C_maxes = []
        for n in ns:
            C_vals = [r['specific_heat'] for r in all_results[n]]
            C_maxes.append(max(C_vals))

        print(f"\n  {'n':>4s} {'C(n,2)':>8s} {'C_max':>12s} {'log(C_max)':>12s}")
        for n, cm in zip(ns, C_maxes):
            print(f"  {n:4d} {comb(n,2):8d} {cm:12.4f} {np.log(cm):12.4f}")

        if len(ns) >= 2:
            log_ns = [np.log(comb(n, 2)) for n in ns]
            log_Cs = [np.log(cm) for cm in C_maxes]
            # Linear fit in log-log
            if len(ns) == 2:
                alpha = (log_Cs[1] - log_Cs[0]) / (log_ns[1] - log_ns[0])
            else:
                alpha = np.polyfit(log_ns, log_Cs, 1)[0]
            print(f"\n  Effective exponent alpha (C_max ~ V^alpha, V=C(n,2)): {alpha:.4f}")
            print(f"  (alpha=1 for first-order, alpha<1 for continuous/crossover)")


if __name__ == "__main__":
    print("=" * 80)
    print("  STATISTICAL MECHANICS OF TOURNAMENTS")
    print("  H(T) = I(Omega(T), 2) as hard-core lattice gas partition function")
    print("=" * 80)

    # Also do n=4 exhaustive (64 tournaments) for finite-size scaling
    print("\n\n" + "=" * 80)
    print("  n=4 EXHAUSTIVE (for finite-size scaling)")
    print("=" * 80)

    H4 = []
    sv4 = []
    c34 = []
    sc4 = []
    for mask, A in all_tournaments(4):
        H4.append(compute_H(A, 4))
        sv4.append(score_variance(A, 4))
        c34.append(count_directed_3cycles(A, 4))
        sc4.append(scores(A, 4))

    beta_range = np.arange(-2.0, 5.05, 0.1)
    results_4 = stat_mech_analysis(4, H4, sv4, c34, beta_range, label="(EXHAUSTIVE)")
    detect_phase_transition(results_4)
    boltzmann_analysis(4, H4, sv4, c34, sc4)
    correlation_analysis(H4, sv4, c34, 4)

    # n=5 exhaustive
    results_5 = run_exhaustive(5)

    # n=6 sampled
    results_6 = run_sampled(6, num_samples=10000)

    # Cross-n comparison
    summary_comparison(results_5, results_6)

    # Finite-size scaling
    all_results = {4: results_4, 5: results_5, 6: results_6}
    finite_size_scaling(all_results)

    print(f"\n\n{'='*80}")
    print(f"  SUMMARY OF KEY FINDINGS")
    print(f"{'='*80}")

    for n, label in [(4, "n=4"), (5, "n=5"), (6, "n=6")]:
        res = all_results[n]
        C_vals = [r['specific_heat'] for r in res]
        betas = [r['beta'] for r in res]
        idx = np.argmax(C_vals)
        sv_at_0 = [r['mean_sv'] for r in res if abs(r['beta']) < 0.01][0]
        sv_at_5 = [r['mean_sv'] for r in res if abs(r['beta'] - 5.0) < 0.01][0]
        H_at_0 = [r['mean_H'] for r in res if abs(r['beta']) < 0.01][0]
        H_at_5 = [r['mean_H'] for r in res if abs(r['beta'] - 5.0) < 0.01][0]
        print(f"\n  {label}:")
        print(f"    Specific heat peak: beta_c = {betas[idx]:.2f}, C_max = {C_vals[idx]:.4f}")
        print(f"    <H> at beta=0 (uniform): {H_at_0:.4f}")
        print(f"    <H> at beta=5 (near max): {H_at_5:.4f}")
        print(f"    <score_var> at beta=0: {sv_at_0:.6f}")
        print(f"    <score_var> at beta=5: {sv_at_5:.6f}")
        print(f"    Score var change: {sv_at_5 - sv_at_0:+.6f} "
              f"({'DECREASES' if sv_at_5 < sv_at_0 else 'INCREASES'} toward H-maximizers)")

    print(f"\n  INTERPRETATION:")
    print(f"  The H(T) = I(Omega(T), 2) partition function defines a natural")
    print(f"  statistical ensemble on tournaments. Key physics:")
    print(f"  - beta > 0 favors H-maximizers (regular/Paley tournaments)")
    print(f"  - beta < 0 favors H-minimizers (near-transitive tournaments)")
    print(f"  - The specific heat peak indicates a crossover (or phase transition)")
    print(f"    between the 'disordered' (uniform) and 'ordered' (H-maximizing) regimes")
    print(f"  - Score variance serves as an order parameter:")
    print(f"    low sv = regular = high H, high sv = transitive = low H")

    print(f"\n{'='*80}")
    print(f"  END OF STATISTICAL MECHANICS ANALYSIS")
    print(f"{'='*80}")
