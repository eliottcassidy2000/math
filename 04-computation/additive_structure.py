#!/usr/bin/env python3
"""
additive_structure.py — Additive combinatorics of QR vs interval sets

The trace alternation theorem reduces to: QR has more k-sum-zero solutions
than {1,...,m} when k = 1 mod 4, fewer when k = 3 mod 4.

This script investigates the additive structure to understand WHY.

Key questions:
1. What is the structure of k-sum-zero solutions for QR vs Interval?
2. Can we decompose M_k - N_k into understandable terms?
3. Is there a Fourier-analytic proof?

Author: kind-pasteur-2026-03-12-S56c
"""

import cmath
import math
from collections import Counter


def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p - 1) // 2, p) == 1


def count_k_sums(S_list, k, p):
    """Count k-tuples from S summing to 0 mod p, using DP."""
    S = list(S_list)
    # dp[j][r] = number of j-tuples from S summing to r mod p
    dp = [{} for _ in range(k + 1)]
    dp[0][0] = 1
    for j in range(k):
        for r, cnt in dp[j].items():
            for s in S:
                new_r = (r + s) % p
                dp[j + 1][new_r] = dp[j + 1].get(new_r, 0) + cnt
    return dp[k].get(0, 0)


def main():
    print("=" * 70)
    print("ADDITIVE STRUCTURE: QR vs INTERVAL k-SUM SOLUTIONS")
    print("=" * 70)

    # ================================================================
    # Section 1: M_k and N_k for all odd k
    # ================================================================
    print("\n--- Section 1: M_k - N_k for all primes and odd k ---")
    print("  M_k = #{k-tuples from QR summing to 0 mod p}")
    print("  N_k = #{k-tuples from {1,...,m} summing to 0 mod p}")
    print("  Delta_k = p*(M_k - N_k)")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        QR = [j for j in range(1, p) if is_qr(j, p)]
        INT = list(range(1, m + 1))

        print(f"\n  p={p}, QR={QR}, INT={INT}:")
        print(f"    {'k':>4} {'M_k':>12} {'N_k':>12} {'M_k-N_k':>12} {'Delta':>12} {'k%4':>4} {'winner':>8}")

        for k in range(3, p + 1, 2):
            M_k = count_k_sums(QR, k, p)
            N_k = count_k_sums(INT, k, p)
            delta = p * (M_k - N_k)
            winner = "QR" if M_k > N_k else "INT" if N_k > M_k else "TIE"

            print(f"    {k:>4} {M_k:>12} {N_k:>12} {M_k - N_k:>12} {delta:>12} {k%4:>4} {winner:>8}")

    # ================================================================
    # Section 2: Fourier analysis of M_k - N_k
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 2: Fourier decomposition of M_k - N_k")
    print("=" * 70)
    print("""
    M_k = (1/p) * sum_{t=0}^{p-1} [sum_{s in QR} omega^{ts}]^k
        = (1/p) * [m^k + S_P(k)]

    N_k = (1/p) * sum_{t=0}^{p-1} [sum_{s in INT} omega^{ts}]^k
        = (1/p) * [m^k + S_I(k)]

    M_k - N_k = (S_P(k) - S_I(k)) / p = Delta_k / p

    Now decompose S_P - S_I by Fourier mode:
    S_P(k) = sum_{t=1}^{p-1} lambda_t^k  (Paley eigenvalues)
    S_I(k) = sum_{t=1}^{p-1} mu_t^k      (Interval eigenvalues)

    Per-mode difference: lambda_t^k - mu_t^k for t = 1,...,p-1

    Since lambda_{p-t} = conj(lambda_t) and mu_{p-t} = conj(mu_t):
    S_P - S_I = 2*sum_{t=1}^m Re(lambda_t^k - mu_t^k)
    """)

    for p in [7, 11]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        QR = [j for j in range(1, p) if is_qr(j, p)]
        INT = list(range(1, m + 1))

        # Paley eigenvalues
        lam = [sum(omega ** (t * s) for s in QR) for t in range(p)]
        # Interval eigenvalues
        mu = [sum(omega ** (t * s) for s in INT) for t in range(p)]

        print(f"\n  p={p}:")
        print(f"    Per-mode decomposition:")
        print(f"    {'t':>4} {'|lam|':>8} {'|mu|':>8} {'arg_lam':>10} {'arg_mu':>10}")
        for t in range(1, m + 1):
            print(f"    {t:>4} {abs(lam[t]):>8.4f} {abs(mu[t]):>8.4f} "
                  f"{cmath.phase(lam[t]):>10.4f} {cmath.phase(mu[t]):>10.4f}")

        print(f"\n    Re(lam_t^k - mu_t^k) per mode:")
        print(f"    {'t':>4}", end="")
        for k in range(3, min(p + 1, 14), 2):
            print(f" {'k=' + str(k):>12}", end="")
        print()

        for t in range(1, m + 1):
            print(f"    {t:>4}", end="")
            for k in range(3, min(p + 1, 14), 2):
                diff = (lam[t] ** k - mu[t] ** k).real
                print(f" {diff:>12.2f}", end="")
            print()

        # Sum over all t (should give Delta/2 since we pair conjugates)
        print(f"\n    Sum (= Delta/2):", end="")
        for k in range(3, min(p + 1, 14), 2):
            total = sum((lam[t] ** k - mu[t] ** k).real for t in range(1, m + 1))
            print(f" {total:>12.2f}", end="")
        print()

    # ================================================================
    # Section 3: Phase difference analysis
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 3: Phase difference lam vs mu per mode")
    print("=" * 70)
    print("  Key: lam_t and mu_t have DIFFERENT magnitudes and phases.")
    print("  The alternation comes from the PHASE STRUCTURE.\n")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        QR = [j for j in range(1, p) if is_qr(j, p)]
        INT = list(range(1, m + 1))

        lam = [sum(omega ** (t * s) for s in QR) for t in range(p)]
        mu = [sum(omega ** (t * s) for s in INT) for t in range(p)]

        theta = math.atan(math.sqrt(p))
        phi_pred = math.pi * (m + 1) / p

        print(f"  p={p}: Paley phase theta={theta:.4f}, Interval phase phi={phi_pred:.4f}")
        print(f"    theta - pi/2 = {theta - math.pi/2:.4f}")
        print(f"    phi   - pi/2 = {phi_pred - math.pi/2:.4f}")

        print(f"    {'t':>4} {'arg(lam)':>10} {'arg(mu)':>10} {'phase_diff':>12} {'chi(t)':>8}")
        for t in range(1, m + 1):
            phase_l = cmath.phase(lam[t])
            phase_m = cmath.phase(mu[t])
            chi_t = "QR" if is_qr(t, p) else "NQ"

            # Paley: arg(lam_t) = pi - theta if chi(t)=1 (QR), -theta if chi(t)=-1 (NQR)
            # Wait: lam_t = (chi(t)*g - 1)/2 where g = i*sqrt(p)
            # For QR: (g-1)/2 = (i*sqrt(p)-1)/2. arg = pi - arctan(sqrt(p)) = pi - theta
            # For NQR: (-g-1)/2 = (-i*sqrt(p)-1)/2. arg = -(pi - arctan(sqrt(p))) = -(pi-theta)
            expected_l = (math.pi - theta) if is_qr(t, p) else -(math.pi - theta)
            # Normalize to (-pi, pi]
            while expected_l > math.pi:
                expected_l -= 2 * math.pi
            while expected_l <= -math.pi:
                expected_l += 2 * math.pi

            # Interval: arg(mu_t) = pi*t*(m+1)/p
            expected_m = math.pi * t * (m + 1) / p
            while expected_m > math.pi:
                expected_m -= 2 * math.pi
            while expected_m <= -math.pi:
                expected_m += 2 * math.pi

            print(f"    {t:>4} {phase_l:>10.4f} ({expected_l:>7.4f}) "
                  f"{phase_m:>10.4f} ({expected_m:>7.4f}) "
                  f"{(phase_l - phase_m):>12.4f} {chi_t:>8}")

    # ================================================================
    # Section 4: H crossover analysis — exact OCF comparison
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 4: Exact H from trace sums (approximation via tr/k)")
    print("=" * 70)
    print("  H_approx = 1 + sum_{k=3,5,...,p} 2^{(k-1)/2} * tr(A^k)/k")
    print("  This is approximate because tr(A^k)/k != c_k for k >= 7.")
    print("  But it captures the TREND correctly.\n")

    for p in [7, 11, 19, 23]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        QR = [j for j in range(1, p) if is_qr(j, p)]
        INT = list(range(1, m + 1))

        eigs_p = [sum(omega ** (t * s) for s in QR) for t in range(p)]
        eigs_i = [sum(omega ** (t * s) for s in INT) for t in range(p)]

        H_P = 1
        H_I = 1
        for k in range(3, p + 1, 2):
            tr_p = sum(e ** k for e in eigs_p).real
            tr_i = sum(e ** k for e in eigs_i).real
            weight = 2 ** ((k - 1) // 2)
            H_P += weight * tr_p / k
            H_I += weight * tr_i / k

        margin = (H_P - H_I) / H_I * 100
        winner = "PALEY" if H_P > H_I else "INTERVAL"
        print(f"  p={p}: H_P_approx={H_P:.0f}, H_I_approx={H_I:.0f}, "
              f"margin={margin:+.2f}%, winner={winner}")

    # ================================================================
    # Section 5: Additive energy comparison
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 5: Additive energy E(S) = #{(a,b,c,d) in S^4 : a+b=c+d}")
    print("=" * 70)
    print("  Additive energy measures how 'structured' a set is additively.")
    print("  Higher energy = more additive structure = more sum-zero solutions.\n")

    for p in [7, 11, 19, 23, 31]:
        m = (p - 1) // 2
        QR = [j for j in range(1, p) if is_qr(j, p)]
        INT = list(range(1, m + 1))

        # Additive energy
        def additive_energy(S, p):
            """Compute E(S) = #{(a,b,c,d): a+b = c+d mod p}"""
            sum_counts = Counter()
            for a in S:
                for b in S:
                    sum_counts[(a + b) % p] += 1
            return sum(c * c for c in sum_counts.values())

        E_QR = additive_energy(QR, p)
        E_INT = additive_energy(INT, p)

        # Normalized: E/|S|^3 (for random set, E ~ |S|^3/p + |S|^2)
        print(f"  p={p}: E(QR)={E_QR}, E(INT)={E_INT}, "
              f"E(QR)/m^3={E_QR/m**3:.4f}, E(INT)/m^3={E_INT/m**3:.4f}")

    # ================================================================
    # Section 6: Sumset structure comparison
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 6: Sumset structure — k-fold sumset density")
    print("=" * 70)

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        QR = frozenset(j for j in range(1, p) if is_qr(j, p))
        INT = frozenset(range(1, m + 1))

        print(f"\n  p={p}:")
        # kS = {s1+...+sk mod p : si in S}
        for label, S in [("QR", QR), ("INT", INT)]:
            S_list = list(S)
            current = set(s % p for s in S_list)
            print(f"    {label}: |S|={len(S)}, ", end="")
            for k in range(2, 8):
                new_set = set()
                for a in current:
                    for s in S_list:
                        new_set.add((a + s) % p)
                current = new_set
                if len(current) == p:
                    print(f"|{k}S|={p}=Z_p", end=" ")
                    break
                print(f"|{k}S|={len(current)}", end=" ")
            print()

    # ================================================================
    # Section 7: The key asymmetry — QR is multiplicatively closed
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 7: QR multiplicative closure and its additive consequences")
    print("=" * 70)
    print("""
    QR_p is CLOSED under multiplication: a*b in QR if a,b in QR.
    Interval {1,...,m} is NOT multiplicatively closed.

    Key property: For t in QR, the map s -> t*s permutes QR.
    So the k-sum sum_{i} s_i = 0 mod p is equivalent to
    sum_{i} (t*s_i) = t * sum_i s_i = 0 mod p (since t != 0).
    The set of solutions is preserved by QR-multiplication.

    This means M_k(QR) has a large symmetry group (QR itself),
    while N_k(INT) has only trivial symmetry.

    The multiplicative structure of QR creates UNIFORM distribution
    of k-sums over residue classes, while INT creates BIASED distribution
    (sums cluster near k*m/2 = k*(p-1)/4).

    M_k = (1/p) * sum_{t=0}^{p-1} |hat{1_QR}(t)|^{2k} ... no wait,
    M_k counts k-tuples summing to 0, not 2k.

    Actually: p*M_k = sum_t (sum_{s in S} omega^{ts})^k = sum_t hat_S(t)^k.
    This is the k-th moment of the DFT of 1_S.

    For QR: |hat_QR(t)| = sqrt((p+1)/4) = const for all t != 0.
    For INT: |hat_INT(t)| varies (Dirichlet kernel).

    So p*M_k = m^k + (p-1)*((p+1)/4)^{k/2} * cos(k*Theta)  [some aggregate phase]
    p*N_k = m^k + sum_{t=1}^{p-1} |hat_INT(t)|^k * cos(k*arg(hat_INT(t)))
    """)

    # Verify: |hat_QR(t)| constant
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        QR = [j for j in range(1, p) if is_qr(j, p)]
        INT = list(range(1, m + 1))

        hat_QR = [abs(sum(omega ** (t * s) for s in QR)) for t in range(1, p)]
        hat_INT = [abs(sum(omega ** (t * s) for s in INT)) for t in range(1, p)]

        print(f"  p={p}:")
        print(f"    |hat_QR|: min={min(hat_QR):.4f}, max={max(hat_QR):.4f}, "
              f"all equal: {max(hat_QR) - min(hat_QR) < 1e-10}")
        print(f"    |hat_INT|: min={min(hat_INT):.4f}, max={max(hat_INT):.4f}, "
              f"ratio max/min={max(hat_INT)/min(hat_INT):.2f}")

    # ================================================================
    # Section 8: Exact H comparison using actual cycle counts
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 8: OCF-based H comparison (using N_k solutions)")
    print("=" * 70)
    print("  The OCF gives: H = 1 + sum_{k odd} 2^{(k-1)/2} * c_k")
    print("  For circulant tournaments: c_k = tr(A^k)/k exactly only for k<=5.")
    print("  But tr(A^k) = sum_{j=0}^{p-1} lambda_j^k, so")
    print("  Delta_H = sum_{k odd} 2^{(k-1)/2} * (tr_P^k - tr_I^k)/k")
    print("          = sum_{k odd} 2^{(k-1)/2} * Delta_k / k")
    print("  (This is APPROXIMATE for k>=7 due to non-simple walks.)\n")

    # Let's compute the cumulative contribution and see where crossover happens
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        QR = [j for j in range(1, p) if is_qr(j, p)]
        INT = list(range(1, m + 1))

        eigs_p = [sum(omega ** (t * s) for s in QR) for t in range(p)]
        eigs_i = [sum(omega ** (t * s) for s in INT) for t in range(p)]

        print(f"  p={p}: Cumulative Delta_H (partial sums up to k):")
        cum = 0
        for k in range(3, p + 1, 2):
            tr_p = sum(e ** k for e in eigs_p).real
            tr_i = sum(e ** k for e in eigs_i).real
            weight = 2 ** ((k - 1) // 2)
            cum += weight * (tr_p - tr_i) / k
            sign = "+" if cum > 0 else "-"
            print(f"    k<={k}: cum={cum:>16.1f} ({sign}, k%4={k%4})")

    # ================================================================
    # Section 9: Exact crossover prime prediction
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 9: Crossover prime prediction")
    print("=" * 70)
    print("  At each p = 3 mod 4, check if the approximate H_P > H_I.\n")

    for p in [3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        QR = [j for j in range(1, p) if is_qr(j, p)]
        INT = list(range(1, m + 1))

        eigs_p = [sum(omega ** (t * s) for s in QR) for t in range(p)]
        eigs_i = [sum(omega ** (t * s) for s in INT) for t in range(p)]

        # Approximate H via trace sums
        H_P = 1
        H_I = 1
        for k in range(3, p + 1, 2):
            tr_p = sum(e ** k for e in eigs_p).real
            tr_i = sum(e ** k for e in eigs_i).real
            weight = 2 ** ((k - 1) // 2)
            H_P += weight * tr_p / k
            H_I += weight * tr_i / k

        margin_pct = (H_P - H_I) / H_I * 100
        winner = "PALEY" if H_P > H_I else "INTERVAL"
        print(f"  p={p:>3}: margin = {margin_pct:>+10.4f}% => {winner}")


if __name__ == '__main__':
    main()
