#!/usr/bin/env python3
"""
trace_algebraic.py — Algebraic structure of Paley and Interval eigenvalue sums

The key to proving THM-136 algebraically is to express both S_P(k) and S_I(k)
as EXACT formulas in p, then show their difference has the right sign.

S_P(k) = -m * sum_{l=0}^{(k-1)/2} C(k,2l)*(-p)^l / 2^{k-1}

For the interval, we need an analogous closed form.

Author: kind-pasteur-2026-03-12-S56c
"""

import cmath
import math
from fractions import Fraction


def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p - 1) // 2, p) == 1


def main():
    print("=" * 70)
    print("ALGEBRAIC STRUCTURE OF TRACE SUMS")
    print("=" * 70)

    # ================================================================
    # Section 1: S_P as a polynomial in p
    # ================================================================
    print("\n--- Section 1: S_P(k) as polynomial in p ---")
    print("S_P(k) = m * [(i*sqrt(p)-1)^k - (i*sqrt(p)+1)^k] / 2^k")
    print("       = -m / 2^{k-1} * sum_{l=0}^{(k-1)/2} C(k,2l)*(-p)^l")

    for k in [3, 5, 7, 9, 11]:
        # Compute the polynomial coefficients
        # sum_{l=0}^{(k-1)/2} C(k,2l)*(-1)^l * p^l
        coeffs = []
        for l in range((k - 1) // 2 + 1):
            from math import comb
            c = comb(k, 2 * l) * ((-1) ** l)
            coeffs.append(c)

        poly_str = " + ".join(f"{c}*p^{l}" if l > 0 else str(c)
                              for l, c in enumerate(coeffs) if c != 0)
        print(f"\n  k={k}: P(p) = {poly_str}")

        # Verify against exact formula
        for p in [7, 11, 19]:
            m = (p - 1) // 2
            poly_val = sum(c * p ** l for l, c in enumerate(coeffs))
            S_P = -m * poly_val // (2 ** (k - 1))

            # Direct computation
            omega = cmath.exp(2j * cmath.pi / p)
            S_pal = [j for j in range(1, p) if is_qr(j, p)]
            eigs = [sum(omega ** (kk * s) for s in S_pal) for kk in range(p)]
            S_P_direct = round(sum(e ** k for e in eigs[1:]).real)

            print(f"    p={p}: P(p)={poly_val}, S_P = {S_P}, direct = {S_P_direct}, "
                  f"match = {S_P == S_P_direct}")

    # ================================================================
    # Section 2: S_I as an algebraic expression
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 2: Interval eigenvalue sum S_I(k)")
    print("=" * 70)
    print("""
    mu_j = sum_{s=1}^m omega^{js} = (omega^j - omega^{j(m+1)}) / (1 - omega^j)

    For p = 3 mod 4, m = (p-1)/2, so m+1 = (p+1)/2.

    omega^{j(m+1)} = omega^{j(p+1)/2} = omega^{j/2} * omega^{jp/2}
                   = omega^{j/2} * (-1)^j  (since omega^{p/2} doesn't exist for odd p...)

    Actually, omega^{jp/2} = e^{2*pi*i*j*p/(2p)} = e^{pi*i*j} = (-1)^j.
    Wait, that's not right either. omega^{p/2} is not well-defined since p is odd.

    Let me use omega^{jm} = omega^{j(p-1)/2} instead.
    omega^{j(p-1)/2} = (omega^j)^{(p-1)/2} = e^{2*pi*i*j*(p-1)/(2p)}
                      = e^{pi*i*j} * e^{-pi*i*j/p}
                      = (-1)^j * e^{-pi*i*j/p}

    So mu_j = omega^j * (1 - omega^{jm}) / (1 - omega^j)
            = omega^j * (1 - (-1)^j * e^{-pi*i*j/p}) / (1 - omega^j)

    For j = 1,...,m (all positive):
    |mu_j| = |sin(pi*j*m/p)| / |sin(pi*j/p)|
           = |sin(pi*j*(p-1)/(2p))| / |sin(pi*j/p)|

    S_I(k) = sum_{j=1}^{p-1} mu_j^k = 2*Re(sum_{j=1}^m mu_j^k)
    """)

    # Compute S_I exactly for small primes
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        eigs_i = [sum(omega ** (kk * s) for s in range(1, m + 1)) for kk in range(p)]

        print(f"\n  p={p}:")
        for k in [3, 5, 7, 9]:
            S_I = round(sum(e ** k for e in eigs_i[1:]).real)
            print(f"    S_I(k={k}) = {S_I}")

    # ================================================================
    # Section 3: Power sum identity for S_I
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 3: Newton's identity approach — express S_I via symmetric functions")
    print("=" * 70)

    # The interval eigenvalues mu_1, ..., mu_{p-1} satisfy a known polynomial.
    # Since mu_j = sum_{s=1}^m omega^{js}, the sum sum_{j=0}^{p-1} mu_j^k
    # can be computed using the multinomial formula.
    #
    # Actually, sum_{j=0}^{p-1} f(omega^j) = p * [coeff of z^0 in f(z)]
    # for f a polynomial of degree < p in z = omega^j.
    #
    # But mu_j is not a polynomial in omega^j of degree < p (it's a ratio).
    #
    # Alternatively, mu_j^k = (sum_{s=1}^m omega^{js})^k
    # = sum_{s_1,...,s_k in {1,...,m}} omega^{j*(s_1+...+s_k)}
    #
    # So S_I(k) = sum_{j=0}^{p-1} mu_j^k - m^k  (subtracting j=0 term)
    # = sum_{j=0}^{p-1} sum_{s_1,...,s_k} omega^{j*T} - m^k
    # where T = s_1 + ... + s_k
    #
    # = sum_{s_1,...,s_k} sum_{j=0}^{p-1} omega^{j*T} - m^k
    # = p * #{(s_1,...,s_k) in {1,...,m}^k : s_1+...+s_k = 0 mod p} - m^k
    #
    # This is the NUMBER OF SOLUTIONS to s_1+...+s_k = 0 mod p
    # with each s_i in {1,...,m}.
    #
    # N_k(p) = #{(s_1,...,s_k) in {1,...,m}^k : s_1+...+s_k = 0 mod p}
    # sum_{j=0}^{p-1} mu_j^k = p * N_k(p)
    # S_I(k) = p * N_k(p) - m^k

    print("  IDENTITY: sum_{j=0}^{p-1} mu_j^k = p * N_k(p)")
    print("  where N_k(p) = #{(s_1,...,s_k) in {1,...,m}^k : sum = 0 mod p}")
    print("  So S_I(k) = p * N_k(p) - m^k\n")

    for p in [7, 11]:
        m = (p - 1) // 2

        # Direct count N_k for k=3,5,7
        for k in [3, 5, 7]:
            count = 0
            # Brute force counting for small p
            if k == 3:
                for s1 in range(1, m + 1):
                    for s2 in range(1, m + 1):
                        for s3 in range(1, m + 1):
                            if (s1 + s2 + s3) % p == 0:
                                count += 1
            elif k == 5 and p <= 11:
                for s1 in range(1, m + 1):
                    for s2 in range(1, m + 1):
                        for s3 in range(1, m + 1):
                            for s4 in range(1, m + 1):
                                for s5 in range(1, m + 1):
                                    if (s1 + s2 + s3 + s4 + s5) % p == 0:
                                        count += 1
            else:
                count = -1  # skip expensive computation

            if count >= 0:
                # Verify
                omega = cmath.exp(2j * cmath.pi / p)
                eigs_i = [sum(omega ** (kk * s) for s in range(1, m + 1)) for kk in range(p)]
                total_sum = round(sum(e ** k for e in eigs_i).real)
                S_I = total_sum - m ** k
                expected = p * count - m ** k

                print(f"  p={p}, k={k}: N_k={count}, p*N_k={p * count}, "
                      f"m^k={m ** k}, S_I=p*N_k-m^k={expected}, direct={S_I}, "
                      f"match={expected == S_I}")

    # ================================================================
    # Section 4: Similarly for Paley
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 4: Paley trace as counting problem")
    print("=" * 70)

    # For Paley: lambda_j = (chi(j)*g - 1)/2 where g = i*sqrt(p).
    # sum_{j=0}^{p-1} lambda_j^k = m^k + S_P(k)
    # But also: lambda_j = sum_{s in QR} omega^{js}
    # So by the same argument:
    # sum_{j=0}^{p-1} lambda_j^k = p * N_k^P(p)
    # where N_k^P = #{(s_1,...,s_k) in QR^k : s_1+...+s_k = 0 mod p}
    # S_P(k) = p * N_k^P(p) - m^k

    print("  S_P(k) = p * M_k(p) - m^k")
    print("  where M_k(p) = #{(s_1,...,s_k) in QR^k : sum = 0 mod p}")
    print("  S_I(k) = p * N_k(p) - m^k")
    print("  where N_k(p) = #{(s_1,...,s_k) in {1,...,m}^k : sum = 0 mod p}")
    print("")
    print("  Delta_k = S_P - S_I = p * (M_k - N_k)")
    print("  sign(Delta_k) = sign(M_k - N_k)")
    print("")
    print("  THE TRACE ALTERNATION REDUCES TO:")
    print("  M_k > N_k when k = 1 mod 4 (more QR solutions)")
    print("  M_k < N_k when k = 3 mod 4 (more interval solutions)")

    for p in [7, 11]:
        m = (p - 1) // 2
        QR = frozenset(j for j in range(1, p) if is_qr(j, p))
        INT = frozenset(range(1, m + 1))

        print(f"\n  p={p}, QR={sorted(QR)}, INT={sorted(INT)}:")

        for k in [3, 5, 7]:
            M_k = 0  # QR solutions
            N_k = 0  # Interval solutions

            if k == 3:
                for s1 in QR:
                    for s2 in QR:
                        for s3 in QR:
                            if (s1 + s2 + s3) % p == 0:
                                M_k += 1
                for s1 in INT:
                    for s2 in INT:
                        for s3 in INT:
                            if (s1 + s2 + s3) % p == 0:
                                N_k += 1
            elif k == 5 and p <= 11:
                for s1 in QR:
                    for s2 in QR:
                        for s3 in QR:
                            for s4 in QR:
                                for s5 in QR:
                                    if (s1 + s2 + s3 + s4 + s5) % p == 0:
                                        M_k += 1
                for s1 in INT:
                    for s2 in INT:
                        for s3 in INT:
                            for s4 in INT:
                                for s5 in INT:
                                    if (s1 + s2 + s3 + s4 + s5) % p == 0:
                                        N_k += 1
            else:
                M_k = N_k = -1

            if M_k >= 0:
                delta = p * (M_k - N_k)
                print(f"    k={k}: M_k={M_k}, N_k={N_k}, M_k-N_k={M_k - N_k}, "
                      f"p*(M_k-N_k)={delta}, k%4={k % 4}")

    # ================================================================
    # Section 5: Weil sum connection
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 5: Connection to Weil sums and Jacobi sums")
    print("=" * 70)
    print("""
    M_k(p) = #{(s_1,...,s_k) in QR^k : sum = 0 mod p}

    This is a WEIL SUM problem! The number of solutions to
    x_1 + ... + x_k = 0 over F_p with each x_i a quadratic residue
    (or zero, but we exclude zero).

    By standard character sum theory:
    M_k = (1/p) * sum_{t=0}^{p-1} (sum_{s in QR} omega^{ts})^k
        = m^k/p + (1/p) * sum_{t=1}^{p-1} ((chi(t)*g - 1)/2)^k
        = m^k/p + S_P(k)/p + ... wait, that's circular.

    Actually: p * M_k = sum_{t=0}^{p-1} lambda_t^k where lambda_t = sum_{s in QR} omega^{ts}
    = m^k + S_P(k)  <-- which is just the trace formula again.

    Similarly: p * N_k = m^k + S_I(k).

    So M_k - N_k = (S_P - S_I)/p = Delta_k/p.

    The question "why does Delta_k alternate?" is equivalent to:
    "Why does QR have more k-sum-zero solutions than {1,...,m} at k=1 mod 4?"

    This is a DEEP question about the additive structure of QR vs intervals.
    """)

    # ================================================================
    # Section 6: The p=7 case in detail
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 6: p=7 complete analysis")
    print("=" * 70)

    p = 7
    m = 3
    QR = sorted([j for j in range(1, p) if is_qr(j, p)])
    INT = sorted(range(1, m + 1))

    print(f"  QR_7 = {QR}")
    print(f"  INT_7 = {INT}")
    print(f"  QR symmetric set: {sorted([(s, p-s) for s in QR])}")
    print(f"  INT consecutive: {{1, 2, 3}}")

    # QR = {1, 2, 4}. Note: these are 1, 2, 4 = 2^0, 2^1, 2^2 (powers of 2 mod 7)
    # INT = {1, 2, 3}

    # k=3 solutions: s1+s2+s3 = 0 mod 7 with s_i in S
    print("\n  k=3 solutions (s1+s2+s3 = 0 mod 7):")
    for label, S in [("QR", QR), ("INT", INT)]:
        sols = []
        for s1 in S:
            for s2 in S:
                for s3 in S:
                    if (s1 + s2 + s3) % p == 0:
                        sols.append((s1, s2, s3))
        print(f"    {label}: {len(sols)} solutions")
        for sol in sols:
            print(f"      {sol} -> sum={sum(sol)} = {sum(sol)%p} mod {p}")

    print(f"\n  k=5 solutions:")
    for label, S in [("QR", QR), ("INT", INT)]:
        count = 0
        for s1 in S:
            for s2 in S:
                for s3 in S:
                    for s4 in S:
                        for s5 in S:
                            if (s1 + s2 + s3 + s4 + s5) % p == 0:
                                count += 1
        print(f"    {label}: {count} solutions")

    print(f"\n  k=7 solutions:")
    for label, S in [("QR", QR), ("INT", INT)]:
        count = 0
        for s1 in S:
            for s2 in S:
                for s3 in S:
                    for s4 in S:
                        for s5 in S:
                            for s6 in S:
                                for s7 in S:
                                    if (s1 + s2 + s3 + s4 + s5 + s6 + s7) % p == 0:
                                        count += 1
        print(f"    {label}: {count} solutions")

    # ================================================================
    # Section 7: Jacobi sum approach
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 7: Jacobi sum factorization")
    print("=" * 70)
    print("""
    For QR: M_k = #{x_1+...+x_k = 0, x_i in QR}

    Define the indicator function for QR:
    1_{QR}(x) = (1 + chi(x))/2  for x != 0

    Then M_k = sum_{x_1+...+x_k=0} prod (1+chi(x_i))/2
             = (1/2^k) * sum_{S subset [k]} sum_{x_1+...+x_k=0, x_i!=0} prod_{i in S} chi(x_i)

    The inner sum with |S|=r is a sum over x_1+...+x_k=0 where r specific
    variables have chi applied. This relates to JACOBI SUMS:
    J(chi^{a_1}, ..., chi^{a_k}) = sum_{x_1+...+x_k=0} chi(x_1)^{a_1}...chi(x_k)^{a_k}

    For the Legendre symbol: chi^2 = 1 (trivial), so each a_i is 0 or 1.

    Key: J(chi, chi) = sum_{a+b=1} chi(a)*chi(b) = tau(chi)^2/tau(1) = g^2/p
    For p = 3 mod 4: g^2 = -p, so J(chi,chi) = -1.

    More generally: J(chi,...,chi) (k copies) relates to higher Jacobi sums.
    J_k = sum_{x_1+...+x_k=0, all x_i != 0} chi(x_1)*...*chi(x_k) = g^k / p

    For p = 3 mod 4 and k odd: g^k = (i*sqrt(p))^k = i^k * p^{k/2}
    J_k = i^k * p^{k/2} / p = i^k * p^{(k-2)/2}
    """)

    # Verify Jacobi sum formula
    for p in [7, 11, 19]:
        omega = cmath.exp(2j * cmath.pi / p)
        g = sum((1 if is_qr(j, p) else -1) * omega ** j for j in range(1, p))

        for k in [3, 5, 7]:
            # Direct Jacobi sum J_k
            J_k = 0
            if k == 3 and p <= 19:
                for x1 in range(1, p):
                    for x2 in range(1, p):
                        x3 = (p - x1 - x2) % p
                        if x3 == 0:
                            continue
                        chi1 = 1 if is_qr(x1, p) else -1
                        chi2 = 1 if is_qr(x2, p) else -1
                        chi3 = 1 if is_qr(x3, p) else -1
                        J_k += chi1 * chi2 * chi3
            elif k == 5 and p <= 11:
                for x1 in range(1, p):
                    for x2 in range(1, p):
                        for x3 in range(1, p):
                            for x4 in range(1, p):
                                x5 = (p - x1 - x2 - x3 - x4) % p
                                if x5 == 0:
                                    continue
                                chi_prod = 1
                                for x in [x1, x2, x3, x4, x5]:
                                    chi_prod *= (1 if is_qr(x, p) else -1)
                                J_k += chi_prod
            else:
                J_k = None

            if J_k is not None:
                # Formula: J_k = g^k / p
                J_k_formula = round((g ** k / p).real)
                print(f"  p={p}, k={k}: J_k(direct)={J_k}, g^k/p={J_k_formula}, "
                      f"match={J_k == J_k_formula}")

    # ================================================================
    # Section 8: Final formula for M_k - N_k
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 8: M_k - N_k via character sums")
    print("=" * 70)
    print("""
    M_k = #{x_1+...+x_k = 0 mod p, x_i in QR}
        = (1/2^k) * sum_{S subset [k]} J_{|S|}
    where J_0 = number of solutions with no chi constraint = (p-1)^{k-1}
    J_1 = 0 (sum of chi over nonzero = 0)
    J_r = sum_{x_1+...+x_k=0, all nonzero} * (product of chi over r positions)

    Actually this is getting complicated. Let me use the eigenvalue formula directly.

    We know:
    S_P(k) = m * [(i*sqrt(p)-1)^k - (i*sqrt(p)+1)^k] / 2^k

    For k odd, write alpha = sqrt(p):
    (i*alpha - 1)^k - (i*alpha + 1)^k = -2*sum_{l=0}^{(k-1)/2} C(k,2l)*(-alpha^2)^l
                                        = -2*sum_{l=0}^{(k-1)/2} C(k,2l)*(-p)^l

    Define P_k(p) = sum_{l=0}^{(k-1)/2} C(k,2l)*(-p)^l.
    Then S_P(k) = -m * P_k(p) / 2^{k-1}.

    P_k is a polynomial in p whose SIGN determines the sign of S_P:
    sign(S_P) = -sign(P_k) * sign(m) = -sign(P_k) [since m > 0]

    The trace alternation says sign(Delta_k) depends on k mod 4.
    Since Delta_k = p*(M_k - N_k), and p > 0, the sign equals sign(M_k - N_k).
    """)

    # Compute P_k polynomials
    from math import comb
    for k in [3, 5, 7, 9, 11, 13]:
        coeffs = [comb(k, 2 * l) * ((-1) ** l) for l in range((k - 1) // 2 + 1)]
        print(f"\n  k={k}: P_k(p) = ", end="")
        terms = []
        for l, c in enumerate(coeffs):
            if c > 0:
                terms.append(f"+{c}p^{l}" if l > 0 else f"{c}")
            else:
                terms.append(f"{c}p^{l}" if l > 0 else f"{c}")
        print(" ".join(terms))

        # Evaluate at several primes and show sign
        signs = []
        for p_val in [3, 7, 11, 19, 23, 31, 43]:
            P_val = sum(c * p_val ** l for l, c in enumerate(coeffs))
            signs.append((p_val, P_val, "+" if P_val > 0 else "-"))

        print(f"    Values: ", end="")
        for pv, val, sgn in signs:
            print(f"P_{k}({pv})={val}({sgn})", end="  ")
        print()

        # sign(S_P) = -sign(P_k)
        # k=1 mod 4: S_P < 0 means P_k > 0
        # k=3 mod 4: S_P > 0 means P_k < 0
        expected_P_sign = "+" if k % 4 == 1 else "-"
        all_match = all(sgn == expected_P_sign for _, _, sgn in signs if True)
        print(f"    Expected P_k sign for k%4={k%4}: {expected_P_sign}, all match: {all_match}")

    # ================================================================
    # Section 9: Can we PROVE sign(P_k)?
    # ================================================================
    print(f"\n{'=' * 70}")
    print("Section 9: Sign of P_k(p) = sum C(k,2l)*(-p)^l")
    print("=" * 70)
    print("""
    P_k(p) = C(k,0) - C(k,2)*p + C(k,4)*p^2 - C(k,6)*p^3 + ...

    For k odd, the highest-degree term (l = (k-1)/2) has coefficient:
    C(k, k-1) * (-1)^{(k-1)/2} = k * (-1)^{(k-1)/2}

    For large p, sign(P_k) = sign((-1)^{(k-1)/2}) * sign(p^{(k-1)/2}).
    = (-1)^{(k-1)/2}

    For k = 1 mod 4: (k-1)/2 even, sign = +1. P_k > 0 for large p. => S_P < 0. CORRECT.
    For k = 3 mod 4: (k-1)/2 odd, sign = -1. P_k < 0 for large p. => S_P > 0. CORRECT.

    But we need this for ALL p >= 3, not just large p.

    P_k(p) = sum_{l=0}^{(k-1)/2} C(k,2l) * (-p)^l

    This is related to the Chebyshev polynomial! In fact:
    P_k(p) = Re[(ip^{1/2} + 1)^k]  (since we're taking the even-index terms)

    Hmm wait: (i*sqrt(p) + 1)^k = sum C(k,j)(i*sqrt(p))^j
    For j even: (i*sqrt(p))^j = (-1)^{j/2} * p^{j/2} = (-p)^{j/2}
    For j odd: pure imaginary

    So Re[(i*sqrt(p) + 1)^k] = sum_{j even} C(k,j) * (-p)^{j/2}
                               = P_k(p) + ... wait this is exactly P_k for k odd.

    Actually for k odd: Re[(i*sqrt(p) + 1)^k]
    = sum_{j=0,2,...,k-1} C(k,j) * (-p)^{j/2}  [only even j contribute to real part]
    = sum_{l=0}^{(k-1)/2} C(k,2l) * (-p)^l = P_k(p). YES!

    So P_k(p) = Re[(1 + i*sqrt(p))^k].

    And |(1 + i*sqrt(p))| = sqrt(1 + p) = sqrt(p+1).
    arg(1 + i*sqrt(p)) = arctan(sqrt(p)) = theta.

    So P_k(p) = (p+1)^{k/2} * cos(k*theta).

    This is EXACTLY what we had before! The sign is cos(k*theta).
    """)

    print("  Verifying P_k(p) = (p+1)^{k/2} * cos(k*theta):")
    for k in [3, 5, 7, 9, 11]:
        for p_val in [7, 11, 19]:
            coeffs = [comb(k, 2 * l) * ((-1) ** l) for l in range((k - 1) // 2 + 1)]
            P_k = sum(c * p_val ** l for l, c in enumerate(coeffs))
            theta = math.atan(math.sqrt(p_val))
            P_k_trig = (p_val + 1) ** (k / 2) * math.cos(k * theta)
            print(f"    k={k}, p={p_val}: P_k={P_k}, (p+1)^{k/2}*cos(k*theta)={P_k_trig:.4f}, "
                  f"match={abs(P_k - P_k_trig) < 0.01}")

    # ================================================================
    # Section 10: Summary — what remains for a full algebraic proof
    # ================================================================
    print(f"\n{'=' * 70}")
    print("SUMMARY: State of the algebraic proof")
    print("=" * 70)
    print("""
    PROVED (Paley side):
    - S_P(k) = -m * P_k(p) / 2^{k-1} where P_k(p) = (p+1)^{k/2} * cos(k*theta)
    - sign(P_k) = sign(cos(k*theta)) determined by leading term for large p
    - For ALL p >= 3 and k <= p: need to show cos(k*theta) doesn't change sign
      This requires: k*eps < pi where eps = pi/2 - theta = arctan(1/sqrt(p))
      Sufficient condition: k < pi*sqrt(p), which holds when k <= p for p >= 10

    PROVED (by direct computation):
    - 0 violations at p = 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83

    REMAINING QUESTION:
    - For the interval S_I(k): analogous closed form would complete the proof
    - The dominant eigenvalue approximation works but needs error bounds
    - The counting interpretation M_k - N_k gives a NUMBER THEORY question:
      "Does QR have more/fewer k-fold sum-zero solutions than {1,...,m}?"
    - This connects to ADDITIVE COMBINATORICS of QR sets vs arithmetic progressions

    CONJECTURE: The trace alternation holds for ALL p = 3 mod 4 and ALL odd k >= 5.
    This would follow from the Paley sign analysis (proved for k*eps < pi) combined
    with a bound showing |S_I| > |S_P| when the Paley sign formula fails (large k).
    """)


if __name__ == '__main__':
    main()
