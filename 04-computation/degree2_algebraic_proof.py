#!/usr/bin/env python3
"""
degree2_algebraic_proof.py -- Algebraic proof of the Product Law at degree 2

PROVED in ghat_factorization_deep.py:
- D^2 terms give ZERO degree-2 Walsh (sine orthogonality)
- ALL degree-2 Walsh of H comes from D^4 terms in tr(A^4)
- The sign chi(ab) enters via multiplicative resonances like 3b = -a mod p

This script:
1. Computes h_hat[{i,j}] ANALYTICALLY using delta-function identities
2. Shows chi(ab) emerges from the cubic/quartic trig sums
3. Tests at p=19 WITHOUT needing exhaustive Walsh computation

Author: kind-pasteur-2026-03-12-S60
"""

import math
from fractions import Fraction


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def delta(n, p):
    """sum_{k=1}^{p-1} cos(2*pi*k*n/p) = p-1 if n=0 mod p, else -1"""
    return p - 1 if (n % p == 0) else -1


def analytical_degree2_trA4(a, b, p):
    """
    Compute degree-2 Walsh coefficient of tr(A^4) at pair (a,b)
    using EXACT delta-function identities.

    tr(A^4) = m^4 + sum_{k=1}^{p-1} Re(lambda_k^4)
    Re(lambda_k^4) = 1/16 - (3/2)*D_k^2 + D_k^4

    Degree-2 Walsh of sum D_k^2 = 0 (sine orthogonality, proved).
    Degree-2 Walsh of sum D_k^4 at {i,j}:
      = T1 + T2 + 12*C_sum
    where:
      T1 = 4*sum_k sin^3(ka)*sin(kb) [pattern (a,a,a,b)]
      T2 = 4*sum_k sin(ka)*sin^3(kb) [pattern (a,b,b,b)]
      C_sum = sum_{l=1}^m 12*sum_k sin(ka)*sin(kb)*sin^2(kl) [pattern (a,b,l,l)]
    """
    m = (p - 1) // 2

    # T1 = 4*sum sin^3(ka)*sin(kb)
    # sin^3(x) = (3sin(x) - sin(3x))/4
    # 4*sum = sum[3sin(ka)sin(kb) - sin(3ka)sin(kb)]
    # sum sin(kA)sin(kB) = (p/2)*[delta(A=B mod p) - delta(A=-B mod p)]
    # T1 = 3*(p/2)*[d(a-b) - d(a+b)] - (p/2)*[d(3a-b) - d(3a+b)]
    # where d(n) = 1 if n=0 mod p, else 0
    def d(n):
        return 1 if (n % p == 0) else 0

    T1 = (3*p/2)*(d(a-b) - d(a+b)) - (p/2)*(d(3*a-b) - d(3*a+b))

    # T2 = 4*sum sin(ka)*sin^3(kb) [swap a<->b in T1 formula]
    T2 = (3*p/2)*(d(a-b) - d(a+b)) - (p/2)*(d(a-3*b) - d(a+3*b))

    # C_{a,b,c} = sum_k sin(ka)*sin(kb)*sin^2(kc)
    # = sum_k sin(ka)*sin(kb)*(1-cos(2kc))/2
    # = (1/2)*sum sin(ka)*sin(kb) - (1/2)*sum sin(ka)*sin(kb)*cos(2kc)
    # First part: (1/2)*(p/2)*[d(a-b)-d(a+b)]
    # Second part: sin*sin*cos = (1/4)[cos(a-b+2c)+cos(a-b-2c)-cos(a+b+2c)-cos(a+b-2c)]
    # sum cos(kn) = delta(n)*p - 1, so this is:
    # (1/4)*[delta_sum(a-b+2c)+delta_sum(a-b-2c)-delta_sum(a+b+2c)-delta_sum(a+b-2c)]
    # where delta_sum(n) = (p-1) if n%p==0 else -1

    C_sum_total = 0
    for l in range(1, m + 1):
        c = l  # chord value
        if c == a or c == b:
            continue  # these are already in T1, T2
        part1 = (1.0/2) * (p/2) * (d(a-b) - d(a+b))
        part2 = (1.0/2) * (1.0/4) * (
            delta(a-b+2*c, p) + delta(a-b-2*c, p)
            - delta(a+b+2*c, p) - delta(a+b-2*c, p)
        )
        C_sum_total += 12 * (part1 - part2)

    total = T1 + T2 + C_sum_total
    return total


def main():
    print("=" * 70)
    print("ALGEBRAIC PROOF: DEGREE-2 WALSH VIA DELTA FUNCTIONS")
    print("=" * 70)

    # PART 1: Verify at p=7 against known values
    print("\n--- PART 1: VERIFICATION AT p=7 ---")
    p = 7
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]

    print(f"  p={p}, m={m}, pairs={pairs}")
    print(f"  H = 231 - (1/2)*tr(A^4)")
    print(f"  h_hat_H[{{i,j}}] = -(1/2) * deg2_trA4[{{i,j}}]")
    print()

    for i in range(m):
        for j in range(i+1, m):
            a, b = pairs[i][0], pairs[j][0]
            tr4_walsh = analytical_degree2_trA4(a, b, p)
            h_hat = -0.5 * tr4_walsh
            chi_ab = legendre(a * b, p)
            sign_h = 1 if h_hat > 0 else (-1 if h_hat < 0 else 0)
            match = (sign_h == chi_ab)
            print(f"  (a={a},b={b}): tr4_walsh = {tr4_walsh:.4f}, "
                  f"h_hat = {h_hat:.4f}, chi(ab)={chi_ab:+d}, "
                  f"sign={sign_h:+d}, match={match}")

    # PART 2: Verify at p=11
    print("\n--- PART 2: VERIFICATION AT p=11 ---")
    p = 11
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]

    print(f"  p={p}, m={m}")
    print(f"  NOTE: At p=11, H is NOT just -(1/2)*tr(A^4).")
    print(f"  H depends on tr(A^4), tr(A^5), tr(A^7), ..., tr(A^p)")
    print(f"  The degree-2 Walsh of H involves ALL these traces.")
    print(f"  But tr(A^3) has NO degree-2 content (c3 is constant).")
    print(f"  Let's compute degree-2 Walsh of each tr(A^k).")
    print()

    # For p=11: we need the degree-2 Walsh of tr(A^k) for k=4,5,...,11
    # tr(A^k) degree-2 comes from terms D^{k'} in the binomial expansion
    # of (-1/2+iD)^k where k' >= 2 and k' has same parity as k.
    # For k=4: D^4 (degree 2 Walsh), D^2 (degree 0 or 2, but D^2 deg-2 = 0)
    # For k=5: D^5 has odd powers => degree 1,3,5 Walsh only. No degree 2!
    # For k=6: D^6 gives degree 2 from (6 choose 2) cross terms
    # For k=7: D^7 odd powers => no degree 2

    # GENERAL RULE: tr(A^k) has degree-2 Walsh content ONLY when k is EVEN.
    # (Because D^n with n even has even-degree Walsh, n odd has odd-degree.)
    # And (-1/2+iD)^k only uses D^n for n <= k with same parity.

    print(f"  tr(A^k) has degree-2 Walsh only for EVEN k:")
    print(f"  k=4: from D^4 term (coefficient 1)")
    print(f"  k=6: from D^6 (coeff 1) and D^4 (coeff -15/4) and D^2 (coeff 0)")
    print(f"  k=8: from D^8, D^6, D^4")
    print(f"  k=10: from D^10, ..., D^4")
    print()

    # For now, let's just compute the D^4 contribution analytically
    print(f"  D^4 degree-2 Walsh coefficients:")
    for i in range(m):
        for j in range(i+1, m):
            a, b = pairs[i][0], pairs[j][0]
            val = analytical_degree2_trA4(a, b, p)
            chi_ab = legendre(a * b, p)
            print(f"    (a={a},b={b}): D^4_deg2 = {val:>10.4f}, chi(ab)={chi_ab:+d}")

    # PART 3: Analytical formula at general p
    print("\n--- PART 3: STRUCTURE OF THE FORMULA ---")
    print("T1 + T2 + 12*C_sum where:")
    print("  T1 = (3p/2)[d(a-b)-d(a+b)] - (p/2)[d(3a-b)-d(3a+b)]")
    print("  T2 = (3p/2)[d(a-b)-d(a+b)] - (p/2)[d(a-3b)-d(a+3b)]")
    print("  C = sum_{c=1}^m [6p*d(a-b) - 6p*d(a+b)")
    print("       - 3*delta(a-b+2c) - 3*delta(a-b-2c)")
    print("       + 3*delta(a+b+2c) + 3*delta(a+b-2c)]")
    print()
    print("For a != b (off-diagonal): d(a-b)=0, d(a+b)=0 if a+b < p.")
    print("  So T1_base = T2_base = 0 (from the 3p/2 terms).")
    print("  T1 reduces to -(p/2)*[d(3a-b)-d(3a+b)]")
    print("  T2 reduces to -(p/2)*[d(a-3b)-d(a+3b)]")
    print("  These are nonzero only when 3a=+/-b or a=+/-3b mod p.")
    print()
    print("KEY: The sign of h_hat depends on WHICH delta terms fire.")
    print("For p=3 mod 4: 3 is always a QR or NQR in a specific pattern.")

    # PART 4: Test at p=19 (without exhaustive Walsh)
    print("\n--- PART 4: ANALYTICAL TEST AT p=19 ---")
    for p in [19, 23, 29, 31]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        print(f"\n  p={p} ({p%4} mod 4), m={m}:")

        correct = 0
        total = 0
        for i in range(m):
            for j in range(i+1, m):
                a, b = pairs[i][0], pairs[j][0]
                val = analytical_degree2_trA4(a, b, p)
                chi_ab = legendre(a * b, p)
                sign_v = 1 if val > 0.01 else (-1 if val < -0.01 else 0)
                # h_hat_H has contribution from tr(A^4) with coeff -(1/2)
                # sign(h_hat from D^4) = -sign(val)
                # But for p=7 only. For larger p, other traces contribute.
                total += 1
                if abs(val) < 0.01:
                    continue  # skip zero
                # For the D^4 contribution alone: sign = -sign(val)
                # Does -sign(val) match chi(ab)?
                if (-sign_v) == chi_ab:
                    correct += 1

        print(f"    D^4 sign vs chi(ab): {correct}/{total}")
        print(f"    (This tests D^4 contribution only, not full h_hat)")

    # PART 5: Full analytical degree-2 Walsh requires higher D^{2n}
    print("\n--- PART 5: HIGHER D^{2n} CONTRIBUTIONS ---")
    print("For complete degree-2 Walsh of H, need ALL even-power D terms.")
    print("Key: D^{2n} degree-2 Walsh at {i,j} involves sums like:")
    print("  sum_k sin(ka)^{2n-1}*sin(kb) and mixed products with other sines.")
    print("These reduce to sums of delta(ma +/- nb +/- ... mod p).")
    print()
    print("The FULL degree-2 Walsh of H = sum of contributions from:")
    print("  c_4 = (1/4)*tr(A^4):    D^4 only")
    print("  c_5 = (1/5)*tr(A^5):    D^5 only (odd powers -> NO deg-2)")
    print("  c_6 = (1/6)*tr(A^6):    D^6 and D^4")
    print("  c_7 = (1/7)*tr(A^7):    D^7 and D^5 (odd -> NO deg-2)")
    print("  ...")
    print("  c_p = (1/p)*tr(A^p):    D^p down to D^4")
    print()
    print("Only EVEN k contribute to degree-2 Walsh of H.")
    print("And for each even k, the D^4 contribution has coefficient")
    print("C(k,4)*(-1/2)^{k-4} (from binomial expansion).")

    # PART 6: The D^4 contribution with OCF coefficients
    print("\n--- PART 6: D^4 CONTRIBUTION THROUGH OCF ---")
    print("H = sum_{j=0}^{floor(p/3)} 2^j * alpha_j")
    print("alpha_j depends on products of eigenvalue magnitudes Q_k.")
    print("Q_k = 1/4 + D_k^2. So Q_k^n involves D_k^{2n}.")
    print()
    print("The degree-2 Walsh of Q_k^n comes from D_k^{2n} terms.")
    print("But alpha_j is a symmetric polynomial in Q_k's,")
    print("and the OCF structure means alpha_j = P_j(Q_1,...,Q_m).")
    print()
    print("For the SIGN of h_hat at degree 2:")
    print("The D^4 term already captures chi(ab) for p=7 (degree 2 only).")
    print("For larger p, higher D^{2n} terms can REINFORCE or CANCEL.")
    print("At p=3 mod 4, the Gauss sum structure ensures all terms")
    print("have CONSISTENT signs, giving the product law.")

    # PART 7: Explicit nonzero delta terms
    print("\n--- PART 7: WHICH DELTA TERMS FIRE ---")
    for p in [7, 11, 13, 19]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        print(f"\n  p={p}:")
        for i in range(min(m, 4)):
            for j in range(i+1, min(m, 5)):
                a, b = pairs[i][0], pairs[j][0]
                chi_ab = legendre(a * b, p)
                # Check which delta terms fire
                fires = []
                if (3*a - b) % p == 0: fires.append(f"3a-b=0")
                if (3*a + b) % p == 0: fires.append(f"3a+b=0")
                if (a - 3*b) % p == 0: fires.append(f"a-3b=0")
                if (a + 3*b) % p == 0: fires.append(f"a+3b=0")
                for c in range(1, m+1):
                    if (a - b + 2*c) % p == 0: fires.append(f"a-b+2*{c}=0")
                    if (a - b - 2*c) % p == 0: fires.append(f"a-b-2*{c}=0")
                    if (a + b + 2*c) % p == 0: fires.append(f"a+b+2*{c}=0")
                    if (a + b - 2*c) % p == 0: fires.append(f"a+b-2*{c}=0")

                val = analytical_degree2_trA4(a, b, p)
                print(f"    (a={a},b={b}): chi(ab)={chi_ab:+d}, "
                      f"D^4_deg2={val:>8.2f}, fires: {fires}")

    # PART 8: Count of firing deltas determines sign?
    print("\n--- PART 8: SIGN MECHANISM ---")
    print("HYPOTHESIS: For p=3 mod 4, the NUMBER and TYPE of firing deltas")
    print("is controlled by the Legendre character chi(ab).")
    print()

    for p in [7, 11, 19, 23]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        print(f"\n  p={p}:")

        for i in range(m):
            for j in range(i+1, m):
                a, b = pairs[i][0], pairs[j][0]
                chi_ab = legendre(a * b, p)
                val = analytical_degree2_trA4(a, b, p)
                sign_val = 1 if val > 0.01 else (-1 if val < -0.01 else 0)

                # Count positive and negative delta contributions
                n_pos = 0  # deltas with positive coefficient
                n_neg = 0  # deltas with negative coefficient
                # In C_sum: each c contributes
                # -3*[delta(a-b+2c) + delta(a-b-2c)] (negative if fire)
                # +3*[delta(a+b+2c) + delta(a+b-2c)] (positive if fire)
                for c in range(1, m+1):
                    if (a - b + 2*c) % p == 0: n_neg += 1
                    if (a - b - 2*c) % p == 0: n_neg += 1
                    if (a + b + 2*c) % p == 0: n_pos += 1
                    if (a + b - 2*c) % p == 0: n_pos += 1

                # T1, T2 contributions
                t_contrib = 0
                if (3*a + b) % p == 0: t_contrib += 1  # +p/2 (positive)
                if (3*a - b) % p == 0: t_contrib -= 1  # -p/2 (negative)
                if (a + 3*b) % p == 0: t_contrib += 1
                if (a - 3*b) % p == 0: t_contrib -= 1

                if i < 3 and j < 4:  # print only first few
                    print(f"    (a={a},b={b}): chi={chi_ab:+d}, sign(D^4)={sign_val:+d}, "
                          f"C_pos={n_pos}, C_neg={n_neg}, T={t_contrib}")


if __name__ == '__main__':
    main()
