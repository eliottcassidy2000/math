#!/usr/bin/env python3
"""
trace_alternation_proof.py — Analytical proof of the trace alternation theorem

THEOREM (HYP-481): For primes p = 3 mod 4, the trace difference
  Delta_k = tr(A^k)_Paley - tr(A^k)_Interval
has sign(Delta_k) = +1 if k = 1 mod 4, -1 if k = 3 mod 4, for odd k >= 5.

PROOF STRATEGY:
  1. Paley trace: tr(A^k) = m^k + m * Re[(g-1)^k - (g+1)^k] / 2^k
     where g = i*sqrt(p) (Gauss sum for p = 3 mod 4)
  2. Interval trace: tr(A^k) = m^k + sum_{j=1}^{p-1} lambda_j^k
     where lambda_j = sum_{s=1}^{m} omega^{js} (Dirichlet kernel)
  3. The m^k terms cancel, so Delta_k depends on eigenvalue k-th powers only.
  4. For Paley: all |lambda| = sqrt(p)/2, phases determined by chi(j).
  5. For Interval: |lambda_j| = |sin(pi*j*m/p) / sin(pi*j/p)|, one dominant.

KEY INSIGHT: The sign alternation comes from the Gauss sum being purely imaginary.
  g = i*sqrt(p) for p = 3 mod 4.
  (g-1)^k - (g+1)^k for k odd:
  = (i*sqrt(p) - 1)^k - (i*sqrt(p) + 1)^k
  Write alpha = i*sqrt(p), then:
  = (alpha - 1)^k - (alpha + 1)^k = -2 * sum_{j even} C(k,j) * alpha^j

  Since alpha = i*sqrt(p): alpha^j = (i*sqrt(p))^j = i^j * p^{j/2}
  For j even: alpha^j = (-1)^{j/2} * p^{j/2} (REAL)

  So (g-1)^k - (g+1)^k = -2 * sum_{j=0,2,...,k-1} C(k,j) * (-1)^{j/2} * p^{j/2}

  This is REAL. The Paley excess trace is:
  m * Re[(g-1)^k - (g+1)^k] / 2^k = m * [(g-1)^k - (g+1)^k] / 2^k

  Now factor out the leading term: j = k-1 (largest even index < k):
  Main term: -2 * C(k,k-1) * (-1)^{(k-1)/2} * p^{(k-1)/2}
           = -2k * (-1)^{(k-1)/2} * p^{(k-1)/2}

  For k = 1 mod 4: (k-1)/2 = even, so (-1)^{(k-1)/2} = +1
    Main term = -2k * p^{(k-1)/2} < 0
    But we have [(g-1)^k - (g+1)^k], and the PALEY trace has
    m * this / 2^k in the NON-zero eigenvalue sum. Need sign carefully.

  Actually let's be more careful about the whole expression.

Author: kind-pasteur-2026-03-12-S56c
"""

import cmath
import math
import sys


def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p - 1) // 2, p) == 1


def main():
    print("=" * 70)
    print("TRACE ALTERNATION: ANALYTICAL PROOF")
    print("=" * 70)

    # ================================================================
    # SECTION 1: Verify the Paley trace formula
    # ================================================================
    print("\n--- Section 1: Paley trace formula verification ---")
    print("Formula: tr_Paley(k) = m^k + m * Re[(g-1)^k - (g+1)^k] / 2^k")
    print("where g = Gauss sum, |g|^2 = p, g = i*sqrt(p) for p = 3 mod 4\n")

    for p in [7, 11, 19, 23]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        # Gauss sum
        g = sum((1 if is_qr(j, p) else -1) * omega ** j for j in range(1, p))

        # Paley eigenvalues (direct)
        S_paley = [j for j in range(1, p) if is_qr(j, p)]
        eigs_paley = []
        for k in range(p):
            eigs_paley.append(sum(omega ** (k * s) for s in S_paley))

        # Interval eigenvalues (direct)
        S_interval = list(range(1, m + 1))
        eigs_interval = []
        for k in range(p):
            eigs_interval.append(sum(omega ** (k * s) for s in S_interval))

        print(f"  p={p}: g = {g.real:.6f} + {g.imag:.6f}i")
        print(f"    i*sqrt(p) = {(1j * math.sqrt(p)).imag:.6f}i")
        print(f"    |g - i*sqrt(p)| = {abs(g - 1j * math.sqrt(p)):.2e}")

        # Verify formula for each odd k
        print(f"    {'k':>4} {'tr_formula':>14} {'tr_direct':>14} {'match':>6} {'Delta':>14} {'sgn':>4} {'k%4':>4}")
        for k in range(3, min(p + 1, 26), 2):
            # Formula: tr = m^k + m * Re[(g-1)^k - (g+1)^k] / 2^k
            paley_excess = m * ((g - 1) ** k - (g + 1) ** k).real / (2 ** k)
            tr_formula = m ** k + paley_excess

            # Direct
            tr_direct_p = sum(e ** k for e in eigs_paley).real
            tr_direct_i = sum(e ** k for e in eigs_interval).real

            delta = tr_direct_p - tr_direct_i
            sgn = "+" if delta > 0.5 else "-" if delta < -0.5 else "0"

            match = abs(tr_formula - tr_direct_p) < 0.01

            if k <= 15 or k == p:
                print(f"    {k:>4} {tr_formula:>14.2f} {tr_direct_p:>14.2f} {'ok' if match else 'FAIL':>6} {delta:>14.2f} {sgn:>4} {k%4:>4}")

    # ================================================================
    # SECTION 2: Decompose the Gauss sum expression
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 2: Gauss sum decomposition for p = 3 mod 4")
    print("=" * 70)
    print("""
    For p = 3 mod 4: g = i*sqrt(p)  (since g^2 = (-1)^{(p-1)/2} * p = -p)

    Paley eigenvalue sum (non-zero terms):
      S_P = sum_{j=1}^{p-1} lambda_j^k = m * [(g-1)^k + (-g-1)^k] / 2^k

    For k ODD: (-g-1)^k = -(g+1)^k
      S_P = m * [(g-1)^k - (g+1)^k] / 2^k

    With g = i*sqrt(p), write alpha = sqrt(p):
      g - 1 = i*alpha - 1,  |g-1|^2 = alpha^2 + 1 = p + 1
      g + 1 = i*alpha + 1,  |g+1|^2 = alpha^2 + 1 = p + 1

    So |lambda_j| = sqrt(p+1)/2 for ALL j != 0  (confirmed!)

    Phase: arg(g-1) = pi - arctan(alpha) (second quadrant)
           arg(g+1) = arctan(alpha) (first quadrant)

    (g-1)^k = (p+1)^{k/2} * e^{ik*[pi - arctan(alpha)]}
    (g+1)^k = (p+1)^{k/2} * e^{ik*arctan(alpha)}

    Difference (k odd):
    (g-1)^k - (g+1)^k = (p+1)^{k/2} * [e^{ik(pi-theta)} - e^{ik*theta}]
      where theta = arctan(sqrt(p))

    e^{ik(pi-theta)} = e^{ikpi} * e^{-ik*theta} = (-1)^k * e^{-ik*theta}
    For k odd: = -e^{-ik*theta}

    So: -e^{-ik*theta} - e^{ik*theta} = -2*cos(k*theta)

    (g-1)^k - (g+1)^k = -2*(p+1)^{k/2} * cos(k*theta)

    S_P = m * [-2*(p+1)^{k/2} * cos(k*theta)] / 2^k
        = -m * (p+1)^{k/2} * cos(k*theta) / 2^{k-1}

    where theta = arctan(sqrt(p)).
    """)

    # Verify this formula
    print("  Verification of S_P = -m * (p+1)^{k/2} * cos(k*theta) / 2^{k-1}:")
    for p in [7, 11, 19, 23, 31, 43]:
        m = (p - 1) // 2
        theta = math.atan(math.sqrt(p))
        omega = cmath.exp(2j * cmath.pi / p)

        S_paley = [j for j in range(1, p) if is_qr(j, p)]
        eigs = [sum(omega ** (k * s) for s in S_paley) for k in range(p)]

        print(f"\n  p={p}, theta = {theta:.6f} rad = {math.degrees(theta):.2f} deg:")
        all_ok = True
        for k in range(3, min(p + 1, 20), 2):
            S_P_direct = sum(e ** k for e in eigs[1:]).real
            S_P_formula = -m * (p + 1) ** (k / 2) * math.cos(k * theta) / 2 ** (k - 1)
            ok = abs(S_P_direct - S_P_formula) < 0.01
            if not ok:
                all_ok = False
            if k <= 13:
                print(f"    k={k}: direct={S_P_direct:.4f}, formula={S_P_formula:.4f}, ok={ok}")
        print(f"    All match: {all_ok}")

    # ================================================================
    # SECTION 3: Interval eigenvalue sum
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 3: Interval eigenvalue sum decomposition")
    print("=" * 70)
    print("""
    For Interval S = {1,...,m}, lambda_j = sum_{s=1}^m omega^{js}.

    This is a geometric series: lambda_j = omega^j * (1 - omega^{jm}) / (1 - omega^j)

    |lambda_j| = |sin(pi*j*m/p) / sin(pi*j/p)|

    S_I = sum_{j=1}^{p-1} lambda_j^k

    The interval eigenvalues are NOT uniform in magnitude.
    But by symmetry: lambda_{p-j} = conj(lambda_j) - 1 ... no wait.

    Actually: lambda_j = sum_{s=1}^m omega^{js}
    lambda_{p-j} = sum_{s=1}^m omega^{(p-j)s} = sum omega^{-js} = conj(lambda_j)

    So eigenvalues come in conjugate pairs: {lambda_j, conj(lambda_j)} for j=1,...,m.

    S_I = sum_{j=1}^m [lambda_j^k + conj(lambda_j)^k] = 2 * sum_{j=1}^m Re(lambda_j^k)

    Write lambda_j = r_j * e^{i*phi_j}, then:
    S_I = 2 * sum_{j=1}^m r_j^k * cos(k*phi_j)

    COMPARISON:
    S_P = -m * (p+1)^{k/2} * cos(k*theta) / 2^{k-1}
    S_I = 2 * sum_{j=1}^m r_j^k * cos(k*phi_j)

    Delta_k = S_P - S_I
    """)

    # ================================================================
    # SECTION 4: Why does sign(Delta_k) depend on k mod 4?
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 4: The k mod 4 mechanism")
    print("=" * 70)

    # For Paley: S_P = -m * (p+1)^{k/2} * cos(k*theta) / 2^{k-1}
    # The sign of S_P depends on cos(k*theta).
    # theta = arctan(sqrt(p)) is close to pi/2 for large p.
    # In fact theta = pi/2 - arctan(1/sqrt(p)) ~ pi/2 - 1/sqrt(p).
    #
    # So k*theta ~ k*pi/2 - k/sqrt(p).
    # cos(k*theta) ~ cos(k*pi/2 - k/sqrt(p))
    #             = cos(k*pi/2)*cos(k/sqrt(p)) + sin(k*pi/2)*sin(k/sqrt(p))
    #
    # For k odd: cos(k*pi/2) = 0, sin(k*pi/2) = (-1)^{(k-1)/2}
    # So cos(k*theta) ~ (-1)^{(k-1)/2} * sin(k/sqrt(p))
    #                 ~ (-1)^{(k-1)/2} * k/sqrt(p)  (for small k/sqrt(p))
    #
    # Therefore S_P ~ -m * (p+1)^{k/2} / 2^{k-1} * (-1)^{(k-1)/2} * k/sqrt(p)
    #
    # For k = 1 mod 4: (k-1)/2 even, (-1)^{(k-1)/2} = +1
    #   S_P ~ -m * positive * (+1) * k/sqrt(p) < 0
    #   (Paley eigenvalue sum is NEGATIVE)
    #
    # For k = 3 mod 4: (k-1)/2 odd, (-1)^{(k-1)/2} = -1
    #   S_P ~ -m * positive * (-1) * k/sqrt(p) > 0
    #   (Paley eigenvalue sum is POSITIVE)
    #
    # The INTERVAL eigenvalue sum S_I is dominated by lambda_1 (Dirichlet peak).
    # lambda_1 has large |lambda_1| and phase phi_1 close to 0 (since S={1,...,m}
    # makes lambda_1 the peak of the kernel).
    #
    # So S_I ~ 2 * r_1^k * cos(k*phi_1) where r_1 >> other r_j.
    # Since phi_1 ~ 0, cos(k*phi_1) ~ 1 for all k.
    # So S_I is ALWAYS POSITIVE and grows like r_1^k.
    #
    # CONCLUSION:
    # Delta_k = S_P - S_I
    # k = 1 mod 4: S_P < 0, S_I > 0 => could go either way (need magnitude)
    # Wait, this gives Delta < 0, meaning INTERVAL wins. But we observed Paley wins!
    #
    # Let me re-examine. The TRACE is m^k + S, so:
    # tr_P = m^k + S_P,  tr_I = m^k + S_I
    # Delta_k = tr_P - tr_I = S_P - S_I
    #
    # At p=7, k=5: tr_P = 14301, tr_I = 13671, Delta = 630 > 0 (Paley wins)
    # So S_P > S_I at k=5. Let me check.

    print("\n  Checking S_P and S_I signs and magnitudes:")
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        theta = math.atan(math.sqrt(p))
        omega = cmath.exp(2j * cmath.pi / p)

        S_paley = [j for j in range(1, p) if is_qr(j, p)]
        eigs_p = [sum(omega ** (k * s) for s in S_paley) for k in range(p)]

        S_interval = list(range(1, m + 1))
        eigs_i = [sum(omega ** (k * s) for s in S_interval) for k in range(p)]

        print(f"\n  p={p}:")
        print(f"    theta = {theta:.6f} ({math.degrees(theta):.2f} deg), pi/2 = {math.pi/2:.6f}")
        print(f"    {'k':>4} {'S_P':>14} {'S_I':>14} {'Delta':>14} {'cos(k*th)':>12} {'(-1)^(k-1)/2':>14}")

        for k in range(3, min(p + 1, 22), 2):
            S_P = sum(e ** k for e in eigs_p[1:]).real
            S_I = sum(e ** k for e in eigs_i[1:]).real
            delta = S_P - S_I
            cos_kth = math.cos(k * theta)
            sign_factor = (-1) ** ((k - 1) // 2)

            print(f"    {k:>4} {S_P:>14.2f} {S_I:>14.2f} {delta:>14.2f} {cos_kth:>12.6f} {sign_factor:>14}")

    # ================================================================
    # SECTION 5: The correct decomposition
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 5: Correct sign analysis")
    print("=" * 70)

    # S_P = -m * (p+1)^{k/2} * cos(k*theta) / 2^{k-1}
    # cos(k*theta) for theta near pi/2:
    # k=3: cos(3*theta) ~ cos(3*pi/2 - 3/sqrt(p)) ~ sin(3/sqrt(p)) > 0
    #       => S_P < 0
    # k=5: cos(5*theta) ~ cos(5*pi/2 - 5/sqrt(p)) ~ -sin(5/sqrt(p)) < 0
    #       => S_P > 0 ... wait, let me be more careful.
    #
    # theta = pi/2 - epsilon, where epsilon = arctan(1/sqrt(p)) > 0
    # k*theta = k*pi/2 - k*epsilon
    # cos(k*pi/2 - k*eps):
    #   k=3: cos(3pi/2 - 3eps) = cos(3pi/2)cos(3eps) + sin(3pi/2)sin(3eps)
    #       = 0*cos(3eps) + (-1)*sin(3eps) = -sin(3eps)
    #       So cos(3*theta) < 0, and S_P = -m*...*(-) = positive
    #   k=5: cos(5pi/2 - 5eps) = cos(pi/2 - 5eps) = sin(5eps)
    #       Wait, 5pi/2 = 2pi + pi/2, so cos(5pi/2 - 5eps) = cos(pi/2 - 5eps) = sin(5eps) > 0
    #       So S_P = -m*...*(+) = negative
    #   k=7: cos(7pi/2 - 7eps) = cos(3pi/2 - 7eps) = -sin(7eps) < 0
    #       S_P = -m*...*(−) = positive
    #   k=9: cos(9pi/2 - 9eps) = cos(pi/2 - 9eps) = sin(9eps) > 0
    #       S_P = -m*...*(+) = negative
    #
    # Pattern: S_P > 0 for k = 3, 7, 11, ... (k = 3 mod 4)
    #          S_P < 0 for k = 5, 9, 13, ... (k = 1 mod 4)
    #
    # Meanwhile S_I > 0 always (dominated by r_1^k * cos(k*phi_1) ~ r_1^k > 0)
    #
    # So Delta_k = S_P - S_I:
    #   k = 3 mod 4: S_P > 0, S_I > 0 => sign depends on magnitudes
    #   k = 1 mod 4: S_P < 0, S_I > 0 => Delta < 0 (Interval wins)
    #
    # But the observation is OPPOSITE: Paley wins at k = 1 mod 4!
    # This means I have a sign error somewhere. Let me recheck directly.

    print("\n  Direct check of S_P sign vs cos(k*theta) sign:")
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        theta = math.atan(math.sqrt(p))
        eps = math.atan(1.0 / math.sqrt(p))

        print(f"\n  p={p}, eps = {eps:.6f} ({math.degrees(eps):.2f} deg):")

        omega = cmath.exp(2j * cmath.pi / p)
        S_paley = [j for j in range(1, p) if is_qr(j, p)]
        eigs_p = [sum(omega ** (k * s) for s in S_paley) for k in range(p)]
        eigs_i = [sum(omega ** (k * s) for s in range(1, m + 1)) for k in range(p)]

        for k in range(3, min(p + 1, 20), 2):
            S_P = sum(e ** k for e in eigs_p[1:]).real
            S_I = sum(e ** k for e in eigs_i[1:]).real
            cos_kth = math.cos(k * theta)
            S_P_formula = -m * (p + 1) ** (k / 2.0) * cos_kth / 2 ** (k - 1)

            # Check which quadrant k*theta is in
            kth_mod_2pi = (k * theta) % (2 * math.pi)
            quadrant = int(kth_mod_2pi / (math.pi / 2)) + 1

            print(f"    k={k}: cos(k*th)={cos_kth:>10.6f}, S_P_formula={S_P_formula:>12.2f}, "
                  f"S_P_direct={S_P:>12.2f}, S_I={S_I:>12.2f}, "
                  f"Delta={S_P - S_I:>12.2f}, k%4={k%4}")

    # ================================================================
    # SECTION 6: The division by k matters for H!
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 6: H comparison via OCF cycle contributions")
    print("=" * 70)
    print("""
    H = 1 + sum_{k odd} 2^{(k-1)/2} * c_k  (OCF where c_k = directed k-cycles)

    For k <= 5: c_k = tr(A^k)/k exactly (HYP-476).
    For k >= 7: c_k = tr(A^k)/k - correction_k.

    The alternation in tr(A^k) does NOT directly translate to H because:
    1. The weights 2^{(k-1)/2} GROW with k
    2. c_k != tr(A^k)/k for k >= 7
    3. The OCF uses ALL odd k simultaneously

    BUT: the TRACE ALTERNATION tells us about the SPECTRAL structure.
    It means Paley and Interval have complementary spectral advantages:
    - At k = 1 mod 4: Paley's flat spectrum wins (higher tr(A^k))
    - At k = 3 mod 4: Interval's concentrated spectrum wins

    The question is: which weight class dominates for H?
    """)

    # Check: what fraction of H comes from each k-class?
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        S_paley = [j for j in range(1, p) if is_qr(j, p)]
        eigs_p = [sum(omega ** (k * s) for s in S_paley) for k in range(p)]
        eigs_i = [sum(omega ** (k * s) for s in range(1, m + 1)) for k in range(p)]

        print(f"\n  p={p}:")
        print(f"    {'k':>4} {'tr_P/k':>14} {'tr_I/k':>14} {'2^((k-1)/2)':>12} {'weighted_diff':>14} {'k%4':>4}")

        total_weighted = 0
        for k in range(3, p + 1, 2):
            tr_p = sum(e ** k for e in eigs_p).real
            tr_i = sum(e ** k for e in eigs_i).real
            ck_p = tr_p / k  # approximate for k >= 7
            ck_i = tr_i / k
            weight = 2 ** ((k - 1) // 2)
            wd = weight * (ck_p - ck_i)
            total_weighted += wd

            if k <= 15 or k == p:
                print(f"    {k:>4} {ck_p:>14.1f} {ck_i:>14.1f} {weight:>12} {wd:>14.1f} {k%4:>4}")

        print(f"    Total weighted difference: {total_weighted:.1f}")
        # This should approximate H(Paley) - H(Interval) when trace = cycle count

    # ================================================================
    # SECTION 7: The crossover as a function of p
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 7: Asymptotic analysis of the crossover")
    print("=" * 70)
    print("""
    S_P = -m * (p+1)^{k/2} * cos(k*theta) / 2^{k-1}

    As p -> infinity: theta -> pi/2, and for fixed k:
      cos(k*theta) -> cos(k*pi/2) = 0 (for k odd)

    More precisely: cos(k*theta) ~ (-1)^{(k+1)/2} * k/sqrt(p)

    So: |S_P| ~ m * (p+1)^{k/2} * k/sqrt(p) / 2^{k-1}
            ~ (p/2) * p^{k/2} * k/sqrt(p) / 2^{k-1}
            = k * p^{(k+1)/2-1} / 2^k
            = k * p^{(k-1)/2} / 2^k

    Meanwhile: |S_I| ~ 2 * r_1^k where r_1 ~ p/pi
    So |S_I| ~ 2 * (p/pi)^k

    Ratio: |S_I| / |S_P| ~ 2 * (p/pi)^k / (k * p^{(k-1)/2} / 2^k)
                         = 2^{k+1} * p^k / (pi^k * k * p^{(k-1)/2})
                         = 2^{k+1} * p^{(k+1)/2} / (pi^k * k)

    For k >= 3 and large p: this ratio -> infinity.

    So S_I ALWAYS dominates S_P eventually (as p grows).

    The crossover prime for each k is where |S_I(k)| ~ |S_P(k)|.

    For H: the OCF sums contributions from ALL k = 3,5,...,p.
    Small k contributions favor Paley (k=5: Paley wins trace).
    Large k contributions favor Interval (Dirichlet peak dominates).

    As p grows: more large-k terms, tipping balance toward Interval.
    This is WHY there's a crossover prime.
    """)

    # Numerical verification of the asymptotic ratio
    print("  |S_I(k)| / |S_P(k)| ratio for k=5:")
    print(f"    {'p':>4} {'|S_P(5)|':>14} {'|S_I(5)|':>14} {'ratio':>10} {'asymp_ratio':>14}")
    for p in [7, 11, 19, 23, 31, 43, 67, 83, 107]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        theta = math.atan(math.sqrt(p))
        k = 5

        S_P_mag = abs(m * (p + 1) ** (k / 2.0) * math.cos(k * theta) / 2 ** (k - 1))

        # Interval: lambda_1 = sum omega^s for s=1..m
        # |lambda_1| = sin(pi*m/p) / sin(pi/p) ~ (p/pi)*sin(pi*m/p) ~ p/pi for m=floor(p/2)
        r1 = abs(math.sin(math.pi * m / p) / math.sin(math.pi / p))
        S_I_approx = 2 * r1 ** k  # dominant eigenvalue contribution

        ratio = S_I_approx / S_P_mag if S_P_mag > 0 else float('inf')
        asymp = 2 ** (k + 1) * p ** ((k + 1) / 2.0) / (math.pi ** k * k)

        print(f"    {p:>4} {S_P_mag:>14.2f} {S_I_approx:>14.2f} {ratio:>10.4f} {asymp:>14.2f}")

    # ================================================================
    # SECTION 8: Key k=p contribution (Hamiltonian cycles)
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 8: k=p contribution (Hamiltonian cycles)")
    print("=" * 70)
    print("""
    At k = p: tr(A^p)/p = c_p + correction (non-simple walks).

    For Paley: S_P(p) = -m * (p+1)^{p/2} * cos(p*theta) / 2^{p-1}
    For Interval: S_I(p) dominated by r_1^p where r_1 ~ p/pi

    The Paley eigenvalue magnitude is ((p+1)/4)^{1/2} ~ sqrt(p)/2.
    So |lambda|^p ~ (sqrt(p)/2)^p = p^{p/2} / 2^p.

    The interval dominant eigenvalue: |lambda_1|^p ~ (p/pi)^p.

    Ratio: (p/pi)^p / (p^{p/2}/2^p) = (2/pi)^p * p^{p/2} -> infinity.

    So at k=p, the interval eigenvalue dominates EXPONENTIALLY.
    This is the fundamental reason the interval wins for large p.
    """)

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        theta = math.atan(math.sqrt(p))

        r1 = abs(math.sin(math.pi * m / p) / math.sin(math.pi / p))
        paley_mag = math.sqrt((p + 1) / 4)

        ratio_per_eig = (r1 / paley_mag) ** p
        print(f"  p={p}: r_1={r1:.4f}, |lam_Paley|={paley_mag:.4f}, "
              f"(r_1/|lam_P|)^p = {ratio_per_eig:.2e}")

    # ================================================================
    # SECTION 9: Summary of the proof structure
    # ================================================================
    print(f"\n{'=' * 70}")
    print("SUMMARY: Proof structure for trace alternation + crossover")
    print("=" * 70)
    print("""
    TRACE ALTERNATION (HYP-481):
    1. Paley eigenvalue sum: S_P(k) = -m*(p+1)^{k/2}*cos(k*theta)/2^{k-1}
       where theta = arctan(sqrt(p)), m = (p-1)/2.
    2. theta = pi/2 - epsilon with epsilon = arctan(1/sqrt(p)) > 0.
    3. cos(k*theta) for odd k has sign pattern determined by k mod 4:
       k = 1 mod 4 (k=5,9,13,...): cos(k*theta) > 0 => S_P < 0
       k = 3 mod 4 (k=3,7,11,...): cos(k*theta) < 0 => S_P > 0
    4. Interval S_I ~ 2*r_1^k > 0 always (Dirichlet peak dominates).
    5. Since S_I > 0 always:
       k = 3 mod 4: S_P > 0, S_I > 0, and |S_P| > |S_I| at small p => Paley wins
       k = 1 mod 4: S_P < 0, S_I > 0 => Interval wins (ALWAYS)

    Wait — this says Interval wins at k=1 mod 4 (k=5,9,...), but observation
    is that PALEY wins at k=1 mod 4!

    LET ME RECHECK THE SIGN...
    """)

    # Direct sign check
    print("  CRITICAL: Direct S_P sign check at k=5")
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        theta = math.atan(math.sqrt(p))
        k = 5

        cos_val = math.cos(k * theta)
        S_P_formula = -m * (p + 1) ** (k / 2.0) * cos_val / 2 ** (k - 1)

        # Direct computation
        omega = cmath.exp(2j * cmath.pi / p)
        S_pal = [j for j in range(1, p) if is_qr(j, p)]
        eigs = [sum(omega ** (kk * s) for s in S_pal) for kk in range(p)]
        S_P_direct = sum(e ** k for e in eigs[1:]).real
        S_I_direct = sum(e ** k for e in [sum(omega ** (kk * s) for s in range(1, m + 1)) for kk in range(1, p)]).real

        tr_P = m ** k + S_P_direct
        tr_I = m ** k + S_I_direct

        print(f"  p={p}, k={k}:")
        print(f"    cos(5*theta) = {cos_val:.6f}")
        print(f"    S_P_formula = {S_P_formula:.2f}")
        print(f"    S_P_direct  = {S_P_direct:.2f}")
        print(f"    S_I_direct  = {S_I_direct:.2f}")
        print(f"    tr_P = {tr_P:.2f}, tr_I = {tr_I:.2f}")
        print(f"    Delta = tr_P - tr_I = {tr_P - tr_I:.2f}")
        print(f"    S_P > S_I? {S_P_direct > S_I_direct}")

    # ================================================================
    # SECTION 10: Interval eigenvalue sum decomposition
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 10: Full eigenvalue-by-eigenvalue comparison")
    print("=" * 70)

    for p in [7]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        S_pal = [j for j in range(1, p) if is_qr(j, p)]
        S_int = list(range(1, m + 1))

        print(f"\n  p={p}, eigenvalue comparison:")
        print(f"    {'j':>4} {'|lam_P|':>10} {'phase_P':>10} {'|lam_I|':>10} {'phase_I':>10} {'QR':>4}")

        for j in range(1, p):
            lam_p = sum(omega ** (j * s) for s in S_pal)
            lam_i = sum(omega ** (j * s) for s in S_int)
            print(f"    {j:>4} {abs(lam_p):>10.4f} {cmath.phase(lam_p):>10.4f} "
                  f"{abs(lam_i):>10.4f} {cmath.phase(lam_i):>10.4f} "
                  f"{'QR' if is_qr(j, p) else 'NQ':>4}")

        print(f"\n    k-th power sums (j=1..{p-1}):")
        for k in [3, 5, 7]:
            terms_p = [(sum(omega ** (j * s) for s in S_pal)) ** k for j in range(1, p)]
            terms_i = [(sum(omega ** (j * s) for s in S_int)) ** k for j in range(1, p)]

            sum_p = sum(t.real for t in terms_p)
            sum_i = sum(t.real for t in terms_i)

            print(f"    k={k}: sum_P = {sum_p:.2f}, sum_I = {sum_i:.2f}, Delta = {sum_p - sum_i:.2f}")

            # Show individual contributions
            print(f"      Per-eigenvalue (Re part):")
            for j in range(1, p):
                print(f"        j={j}: P={terms_p[j-1].real:>10.2f}, I={terms_i[j-1].real:>10.2f}, "
                      f"diff={terms_p[j-1].real - terms_i[j-1].real:>10.2f}")


if __name__ == '__main__':
    main()
