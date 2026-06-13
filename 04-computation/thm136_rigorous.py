#!/usr/bin/env python3
"""
thm136_rigorous.py -- Toward a rigorous proof of THM-136 (Trace Alternation)

THEOREM: For primes p = 3 mod 4 and odd k in [5, p]:
  sign(tr(A_P^k) - tr(A_I^k)) = (-1)^{(k-3)/2}

PROOF STRATEGY:
1. S_P = -m * P_k(p) / 2^{k-1} where P_k(p) = Re[(1+i*sqrt(p))^k] is EXACT
2. S_I = 2*sum r_j^k cos(k*phi_j) where r_j, phi_j are Dirichlet kernel eigenvalues
3. For all k in [5, p]: |S_I| >> |S_P|, so sign(Delta_k) = -sign(S_I)
4. S_I has the correct sign from dominant eigenvalue analysis

RIGOROUS COMPONENTS:
A. sign(P_k(p)) via polynomial evaluation (exact integer arithmetic)
B. |S_I|/|S_P| > 1 via eigenvalue magnitude bounds
C. sign(S_I) via dominant term + error bound

Author: kind-pasteur-2026-03-12-S57
"""

import math
import cmath
from fractions import Fraction


def P_k_exact(k, p):
    """Compute P_k(p) = Re[(1 + i*sqrt(p))^k] exactly as integer.

    P_k(p) = sum_{l=0}^{(k-1)/2} C(k, 2l) * (-p)^l

    Returns exact integer value.
    """
    result = 0
    for l in range((k - 1) // 2 + 1):
        coeff = math.comb(k, 2 * l)
        result += coeff * ((-p) ** l)
    return result


def interval_eigenvalue_magnitudes(p):
    """Compute exact eigenvalue magnitudes r_j for interval tournament.

    mu_j = sum_{s=1}^{m} omega^{js} where omega = e^{2pi*i/p}
    |mu_j| = sin(j*pi*m/p) / sin(j*pi/p)

    Returns list of r_j for j=1,...,m.
    """
    m = (p - 1) // 2
    r = []
    for j in range(1, m + 1):
        num = abs(math.sin(j * math.pi * m / p))
        den = abs(math.sin(j * math.pi / p))
        r.append(num / den)
    return r


def verify_dominance(p, verbose=True):
    """Verify |S_I| > |S_P| and sign(S_I) correct for all odd k in [5, p].

    Returns True if all conditions verified.
    """
    m = (p - 1) // 2
    r = interval_eigenvalue_magnitudes(p)
    r1 = r[0]
    lam_P = math.sqrt((p + 1) / 4)  # Paley eigenvalue magnitude

    if verbose:
        print(f"\n  p = {p}: m = {m}")
        print(f"  Paley |lambda| = {lam_P:.4f}")
        print(f"  Interval |mu_1| = {r1:.4f}")
        print(f"  Ratio r_1/|lam_P| = {r1 / lam_P:.4f}")
        print(f"  Interval magnitudes: {[f'{x:.4f}' for x in r]}")

    # Phase parameters
    delta = math.pi / (2 * p)  # phi_1 - pi/2
    eps = math.atan(1 / math.sqrt(p))  # pi/2 - theta

    all_ok = True
    for k in range(5, p + 1, 2):
        # === PART A: sign(P_k(p)) ===
        Pk = P_k_exact(k, p)
        theta = math.atan(math.sqrt(p))
        cos_k_theta = math.cos(k * theta)
        expected_Pk_sign = 1 if cos_k_theta > 0 else -1

        # S_P = -m * Pk / 2^{k-1}
        S_P_sign = -1 if Pk > 0 else 1  # sign(-m * Pk / 2^{k-1})

        # Expected: k=1 mod 4 => S_P < 0, k=3 mod 4 => S_P > 0
        if k % 4 == 1:
            paley_ok = S_P_sign == -1
        else:
            paley_ok = S_P_sign == 1

        # === PART B: Dominant term of S_I ===
        # S_I = 2 * sum_{j=1}^{m} r_j^k * cos(k * phi_j)
        # Dominant: 2 * r_1^k * cos(k * phi_1)
        # phi_1 = pi/2 + delta, so cos(k*phi_1) = (-1)^{(k+1)/2} * sin(k*delta)
        sin_k_delta = math.sin(k * delta)
        dominant_sign = ((-1) ** ((k + 1) // 2)) * (1 if sin_k_delta > 0 else -1)
        # dominant_sign gives sign of cos(k*phi_1)
        # S_I ~ 2*r_1^k * cos(k*phi_1), so sign(S_I) ~ dominant_sign

        # Expected: k=1 mod 4 => S_I < 0, k=3 mod 4 => S_I > 0
        if k % 4 == 1:
            interval_expected_sign = -1
        else:
            interval_expected_sign = 1

        # === PART C: Error bound ===
        # |S_I - 2*r_1^k*cos(k*phi_1)| <= 2*sum_{j>=2} r_j^k
        dominant_magnitude = 2 * r1 ** k * abs(math.sin(k * delta))
        error_bound = 2 * sum(rj ** k for rj in r[1:])
        dominance_ratio = dominant_magnitude / error_bound if error_bound > 0 else float('inf')

        # Also: |S_I| > |S_P|?
        S_P_magnitude = m * abs(Pk) / (2 ** (k - 1))
        SI_lower = dominant_magnitude - error_bound  # lower bound on |S_I|
        magnitude_ratio = SI_lower / S_P_magnitude if S_P_magnitude > 0 else float('inf')

        # === PART D: Actual computation for verification ===
        omega = cmath.exp(2j * cmath.pi / p)
        S_int = set(range(1, m + 1))
        eigs_I = [sum(omega ** (j * s) for s in S_int) for j in range(p)]
        actual_S_I = sum(e ** k for e in eigs_I).real - m ** k

        S_qr = set(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
        eigs_P = [sum(omega ** (j * s) for s in S_qr) for j in range(p)]
        actual_S_P = sum(e ** k for e in eigs_P).real - m ** k

        actual_Delta = actual_S_P - actual_S_I
        actual_sign = 1 if actual_Delta > 0 else -1
        expected_sign = 1 if k % 4 == 1 else -1

        sign_ok = actual_sign == expected_sign
        if not sign_ok:
            all_ok = False

        # Can we PROVE the sign from the bounds?
        proved = (dominance_ratio > 1) and (magnitude_ratio > 1) and (
            k * delta < math.pi)  # ensures sin(k*delta) > 0

        if verbose:
            status = "OK" if sign_ok else "FAIL"
            prov = "PROVED" if proved else "numeric"
            print(f"    k={k:3d} (mod4={k % 4}): P_k={Pk:>15d}, "
                  f"S_P_sign={'<0' if S_P_sign < 0 else '>0'}, "
                  f"dom/err={dominance_ratio:>8.1f}, "
                  f"|S_I|/|S_P|={magnitude_ratio:>8.1f}, "
                  f"Delta_sign={'>' if actual_Delta > 0 else '<'}0 "
                  f"[{status}] [{prov}]")

    return all_ok


def prove_P_k_sign(verbose=True):
    """Prove sign(P_k(p)) for all relevant (k, p) pairs.

    P_k(p) = sum_{l=0}^{(k-1)/2} C(k,2l)*(-p)^l

    For k=1 mod 4: P_k(p) > 0 (when k*eps < pi)
    For k=3 mod 4: P_k(p) < 0 (when k*eps < pi)
    """
    if verbose:
        print("\n" + "=" * 70)
        print("SIGN OF P_k(p) = Re[(1+i*sqrt(p))^k]")
        print("=" * 70)

    # Check: does P_k have constant sign for p >= 7?
    for k in range(5, 84, 2):
        signs = {}
        for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
            if p % 4 != 3 or k > p:
                continue
            Pk = P_k_exact(k, p)
            signs[p] = 1 if Pk > 0 else -1

        if not signs:
            continue

        all_same = len(set(signs.values())) == 1
        expected = 1 if k % 4 == 1 else -1
        matches = all(v == expected for v in signs.values())

        if verbose and (not matches or k <= 15):
            vals_str = ", ".join(f"p={p}:{'+'if s>0 else '-'}" for p, s in sorted(signs.items()))
            print(f"  k={k:2d} (mod4={k % 4}): expected {'+'if expected>0 else '-'}, "
                  f"{'ALL MATCH' if matches else 'MISMATCH'}: {vals_str}")


def main():
    print("=" * 70)
    print("RIGOROUS ANALYSIS OF THM-136 (TRACE ALTERNATION)")
    print("=" * 70)

    # Part 1: P_k sign analysis
    prove_P_k_sign()

    # Part 2: Dominance analysis at each prime
    print("\n" + "=" * 70)
    print("DOMINANCE ANALYSIS: |S_I| > |S_P| AND SIGN PROOF")
    print("=" * 70)

    all_proved = True
    for p in [7, 11, 19, 23, 31, 43]:
        ok = verify_dominance(p)
        if not ok:
            all_proved = False
            print(f"  *** FAILURE at p={p} ***")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    if all_proved:
        print("  All signs verified numerically for p = 7, 11, 19, 23, 31, 43")
    else:
        print("  Some verifications failed!")

    # Part 3: Asymptotic analysis
    print("\n" + "=" * 70)
    print("ASYMPTOTIC REGIME: When does dominant-term proof work?")
    print("=" * 70)
    print("  The dominant-term proof requires:")
    print("  1. sin(k*pi/(2p)) > 0  (holds for k <= p since k*pi/(2p) <= pi/2)")
    print("  2. dominant/error > 1   (fails at small k for large p)")
    print("  3. |S_I_lower| > |S_P|  (requires r_1/|lam_P| ratio large enough)")

    for p in [7, 11, 19, 23, 31, 43, 59, 83]:
        m = (p - 1) // 2
        r = interval_eigenvalue_magnitudes(p)
        r1 = r[0]
        lam_P = math.sqrt((p + 1) / 4)
        delta = math.pi / (2 * p)

        # Find smallest k where proof works
        min_proved_k = None
        for k in range(5, p + 1, 2):
            sin_kd = math.sin(k * delta)
            dom = 2 * r1 ** k * sin_kd
            err = 2 * sum(rj ** k for rj in r[1:])
            sp_mag = m * abs(P_k_exact(k, p)) / (2 ** (k - 1))
            if dom > err and (dom - err) > sp_mag:
                min_proved_k = k
                break

        total_k = (p - 3) // 2  # number of odd k in [5, p]
        if min_proved_k:
            proved_k = (p - min_proved_k) // 2 + 1
        else:
            proved_k = 0
        print(f"  p={p:3d}: r1/|lam|={r1 / lam_P:.3f}, "
              f"proved for k>={min_proved_k if min_proved_k else 'NONE':>4}, "
              f"covers {proved_k}/{total_k} k-values "
              f"({100 * proved_k / total_k:.0f}%)" if total_k > 0 else "")

    print("\nDONE.")


if __name__ == '__main__':
    main()
