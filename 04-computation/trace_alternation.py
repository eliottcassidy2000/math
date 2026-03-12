"""
trace_alternation.py — Why do Paley and Interval alternate trace dominance?

PATTERN: For Paley primes p = 3 mod 4:
  Paley wins tr(A^k)/k at k = 5, 9, 13, ... (k = 1 mod 4)
  Interval wins tr(A^k)/k at k = 7, 11, 15, ... (k = 3 mod 4)
  c_3 = constant (k = 3: tie)
  Split is EXACTLY even (p=7: 1-1, p=11: 2-2, p=19: 4-4)

HYPOTHESIS: This is related to the quadratic residue structure.
  Paley has eigenvalue phases arg(lam_k) = +/- theta
  where theta depends on chi(k) (Legendre symbol).
  The sign pattern of chi(k) for k=1,...,(p-1)/2 determines
  which trace powers are positive.

Author: kind-pasteur-2026-03-12-S56c
"""

import math
import cmath
import sys

sys.path.insert(0, '04-computation')


def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p - 1) // 2, p) == 1


def eigenvalues(S, p):
    omega = cmath.exp(2j * cmath.pi / p)
    return [sum(omega ** (k * s) for s in S) for k in range(p)]


def trace_Ak(eigs, k):
    return sum(e ** k for e in eigs).real


def main():
    print("=" * 70)
    print("TRACE ALTERNATION PATTERN ANALYSIS")
    print("=" * 70)

    # Detailed trace comparison for each prime
    for p in [7, 11, 19, 23, 31, 43]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2

        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        S_interval = frozenset(range(1, m + 1))

        eigs_p = eigenvalues(S_paley, p)
        eigs_i = eigenvalues(S_interval, p)

        print(f"\n  p={p} (m={m}):")
        print(f"  {'k':>4} {'tr_P/k':>14} {'tr_I/k':>14} {'P>I?':>6} {'k mod 4':>7}")

        paley_wins_mod = {1: 0, 3: 0}
        interval_wins_mod = {1: 0, 3: 0}

        for k in range(3, p + 1, 2):
            tr_p = trace_Ak(eigs_p, k) / k
            tr_i = trace_Ak(eigs_i, k) / k
            k_mod4 = k % 4
            if abs(tr_p - tr_i) < 0.5:
                winner = "TIE"
            elif tr_p > tr_i:
                winner = "P"
                if k_mod4 in paley_wins_mod:
                    paley_wins_mod[k_mod4] += 1
            else:
                winner = "I"
                if k_mod4 in interval_wins_mod:
                    interval_wins_mod[k_mod4] += 1

            if k <= 19 or k == p:
                print(f"  {k:>4} {tr_p:>14.1f} {tr_i:>14.1f} {winner:>6} {k_mod4:>7}")

        print(f"    Paley wins at k mod 4: {paley_wins_mod}")
        print(f"    Interval wins at k mod 4: {interval_wins_mod}")

    # ================================================================
    # Analytical explanation attempt
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"ANALYTICAL EXPLANATION")
    print(f"{'=' * 60}")

    print("""
    For Paley tournament on Z_p (p = 3 mod 4):
      S = QR_p = {a^2 mod p : a = 1,...,(p-1)/2}
      lambda_j = sum_{s in QR} omega^{js} = (chi(j)*g - 1)/2

    where chi = Legendre symbol, g = sum_{a=1}^{p-1} chi(a)*omega^a (Gauss sum).
    |g|^2 = p, arg(g) = ? (depends on p).

    Key: chi(j) = +1 if j is QR, -1 if j is NQR.

    So: lambda_j = (g - 1)/2 if j is QR
        lambda_j = (-g - 1)/2 if j is NQR

    For the trace:
      tr(A^k) = m^k + sum_{j=1}^{p-1} lambda_j^k
      = m^k + sum_{j in QR} ((g-1)/2)^k + sum_{j in NQR} ((-g-1)/2)^k
      = m^k + |QR| * ((g-1)/2)^k + |NQR| * ((-g-1)/2)^k
      = m^k + m * [(g-1)^k + (-g-1)^k] / 2^k

    Now (g-1)^k + (-g-1)^k:
    - If k is even: (-g-1)^k = (g+1)^k, so we get (g-1)^k + (g+1)^k
    - If k is odd: (-g-1)^k = -(g+1)^k, so we get (g-1)^k - (g+1)^k

    For k odd: (g-1)^k - (g+1)^k
    = sum_{j=0}^k C(k,j) g^j [(-1)^{k-j} - 1]
    = sum_{j=0}^k C(k,j) g^j * {0 if k-j even, -2 if k-j odd}
    = -2 * sum_{j even, j<k} C(k,j) g^j

    Actually let me redo this. Since k is odd:
    (g-1)^k - (g+1)^k = sum C(k,j)g^j(-1)^{k-j} - sum C(k,j)g^j
    = sum C(k,j)g^j [(-1)^{k-j} - 1]

    When k-j is even (j has same parity as k, i.e., j is odd): (-1)^{k-j} = 1, term = 0
    When k-j is odd (j is even): (-1)^{k-j} = -1, term = -2*C(k,j)*g^j

    So: (g-1)^k - (g+1)^k = -2 * sum_{j=0,2,4,...} C(k,j) g^j
    """)

    # Verify for specific primes
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        eigs_p = eigenvalues(S_paley, p)

        # Compute Gauss sum
        omega = cmath.exp(2j * cmath.pi / p)
        g = sum((1 if is_qr(j, p) else -1) * omega ** j for j in range(1, p))

        print(f"\n  p={p}: g = {g.real:.6f} + {g.imag:.6f}i, |g| = {abs(g):.6f}")
        print(f"    Expected |g| = sqrt({p}) = {math.sqrt(p):.6f}")

        # Verify eigenvalue formula
        for j in [1, 2, 3]:
            lam_formula = ((1 if is_qr(j, p) else -1) * g - 1) / 2
            lam_actual = eigs_p[j]
            print(f"    j={j}: chi(j)={'QR' if is_qr(j,p) else 'NQ'}, "
                  f"formula={lam_formula.real:.4f}+{lam_formula.imag:.4f}i, "
                  f"actual={lam_actual.real:.4f}+{lam_actual.imag:.4f}i, "
                  f"match={abs(lam_formula - lam_actual) < 1e-10}")

    # ================================================================
    # Now compare with Interval
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"INTERVAL TRACE FORMULA")
    print(f"{'=' * 60}")

    print("""
    For Interval S = {1,...,m} with m = (p-1)/2:
      lambda_j = sum_{s=1}^{m} omega^{js} = omega^j * (1 - omega^{jm}) / (1 - omega^j)

    The DIFFERENCE in traces:
      Delta_k = tr(A^k)_Paley - tr(A^k)_Interval
      = sum_{j=1}^{p-1} [lambda_j(Paley)^k - lambda_j(Interval)^k]

    For Paley: lambda_j = (chi(j)*g - 1)/2
    For Interval: lambda_j = sum_{s=1}^m omega^{js} (Dirichlet kernel)

    The sign of Delta_k alternates with k mod 4.
    """)

    # Compute Delta_k and verify alternation
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        S_interval = frozenset(range(1, m + 1))

        eigs_p = eigenvalues(S_paley, p)
        eigs_i = eigenvalues(S_interval, p)

        print(f"\n  p={p}:")
        print(f"  {'k':>4} {'Delta_k':>14} {'sign':>6} {'k mod 4':>7}")

        for k in range(3, p + 1, 2):
            delta = trace_Ak(eigs_p, k) - trace_Ak(eigs_i, k)
            sign = "+" if delta > 0 else "-" if delta < 0 else "0"
            expected = "+" if k % 4 == 1 else "-"
            match = "ok" if sign == expected else "BREAK"
            if k <= 19 or k == p:
                print(f"  {k:>4} {delta:>14.2f} {sign:>6} {k%4:>7}  {match}")

    # ================================================================
    # Phase coherence analysis at k=5
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"PHASE COHERENCE AT k=5 (Paley advantage)")
    print(f"{'=' * 60}")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        S_interval = frozenset(range(1, m + 1))

        eigs_p = eigenvalues(S_paley, p)
        eigs_i = eigenvalues(S_interval, p)

        k = 5
        # For Paley: all |lam| = sqrt(p)/2, so tr = m^5 + (sqrt(p)/2)^5 * sum exp(5i*theta_j)
        sum_phase_p = sum(cmath.exp(1j * k * cmath.phase(e)) for e in eigs_p[1:])
        sum_phase_i = sum(cmath.exp(1j * k * cmath.phase(e)) for e in eigs_i[1:])

        # For Paley: all magnitudes equal
        mag_p = abs(eigs_p[1])
        phase_factor_p = sum_phase_p.real  # only real part matters for trace

        # For Interval: magnitudes vary
        mag_weighted_i = sum(abs(e) ** k * cmath.exp(1j * k * cmath.phase(e))
                            for e in eigs_i[1:])

        tr5_p = trace_Ak(eigs_p, 5)
        tr5_i = trace_Ak(eigs_i, 5)

        print(f"\n  p={p}:")
        print(f"    Paley: mag^5 = {mag_p**k:.4f}, phase factor = {phase_factor_p:.4f}")
        print(f"    tr5_Paley = m^5 + mag^5 * phase_factor = "
              f"{m**k} + {mag_p**k * phase_factor_p:.2f} = {tr5_p:.2f}")
        print(f"    Interval: weighted phase sum = {mag_weighted_i.real:.4f}")
        print(f"    tr5_Interval = m^5 + {mag_weighted_i.real:.2f} = {tr5_i:.2f}")
        print(f"    Delta = {tr5_p - tr5_i:.2f}")


if __name__ == '__main__':
    main()
