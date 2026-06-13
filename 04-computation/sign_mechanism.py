#!/usr/bin/env python3
"""
sign_mechanism.py -- WHY does sign(h_hat[{i,j}]) = chi(a_i * a_j)?

From degree2_algebraic_proof.py:
- D^4 gives nonzero degree-2 Walsh ONLY when 3a = +/-b mod p
- The sign of this D^4 term is controlled by which of 3a=b or 3a=-b fires
- For most pairs, D^4 = 0, so sign comes from D^6 or higher

This script:
1. Classifies pairs by their 3-multiplicative relation
2. Shows the T contribution (from 3a=+/-b) determines sign when nonzero
3. Shows the C_sum contribution sign is predictable via Legendre
4. Computes D^6 degree-2 Walsh to see if it fills the gaps

Author: kind-pasteur-2026-03-12-S60
"""

import math
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def delta_sum(n, p):
    """sum_{k=1}^{p-1} cos(2*pi*k*n/p) = p-1 if n=0 mod p, else -1"""
    return p - 1 if (n % p == 0) else -1


def d(n, p):
    return 1 if (n % p == 0) else 0


def D4_deg2(a, b, p):
    """Analytical D^4 degree-2 Walsh at (a,b) summed over k=1..p-1.
    Excludes c=a and c=b from the C_sum (those are in T1/T2)."""
    m = (p - 1) // 2
    T1 = (3*p/2)*(d(a-b,p) - d(a+b,p)) - (p/2)*(d(3*a-b,p) - d(3*a+b,p))
    T2 = (3*p/2)*(d(a-b,p) - d(a+b,p)) - (p/2)*(d(a-3*b,p) - d(a+3*b,p))
    C_sum = 0
    for c in range(1, m + 1):
        if c == a or c == b:
            continue
        C_sum += -(3.0/2) * (
            delta_sum(a-b+2*c, p) + delta_sum(a-b-2*c, p)
            - delta_sum(a+b+2*c, p) - delta_sum(a+b-2*c, p)
        )
    return T1 + T2 + C_sum


def D6_deg2_numerical(a, b, p):
    """Numerical D^6 degree-2 Walsh at (a,b).
    D^6 = (sum sigma_l S_l)^6. Degree-2 Walsh at {i,j} requires
    products of 6 sines where indices i,j each appear odd times."""
    m = (p - 1) // 2
    chords = list(range(1, m + 1))  # chord values

    # For each k, the degree-2 contribution from D_k^6 at {a,b}:
    # Patterns (count_a, count_b, rest) with count_a odd, count_b odd, sum=6:
    # (1,1,4), (1,3,2), (1,5,0), (3,1,2), (3,3,0), (5,1,0)
    # For each pattern, enumerate l-choices for the "rest" part

    total = 0
    for k in range(1, p):
        S = {}
        for c in chords:
            S[c] = math.sin(2 * math.pi * k * c / p)

        # Brute force: compute the degree-2 Walsh of D_k^6
        # (1/2^m) * sum_sigma D_k^6 * sigma_a_idx * sigma_b_idx
        # Since D_k = sum sigma_l * S[l], we need degree-2 component

        # Alternative: enumerate multinomial contributions
        # Each term is C(6; n_1,...,n_m) * prod S[l]^{n_l}
        # where n_a is odd, n_b is odd, all others even

        contrib = 0

        # Pattern (5,1,0): 6!/(5!1!) = 6
        contrib += 6 * S[a]**5 * S[b]
        # Pattern (1,5,0): 6!/(1!5!) = 6
        contrib += 6 * S[a] * S[b]**5
        # Pattern (3,3,0): 6!/(3!3!) = 20
        contrib += 20 * S[a]**3 * S[b]**3
        # Pattern (3,1,2) with c != a,b: 6!/(3!1!2!) = 60
        for c in chords:
            if c == a or c == b:
                continue
            contrib += 60 * S[a]**3 * S[b] * S[c]**2
        # Pattern (1,3,2) with c != a,b: 60
        for c in chords:
            if c == a or c == b:
                continue
            contrib += 60 * S[a] * S[b]**3 * S[c]**2
        # Pattern (1,1,4) with single c: 6!/(1!1!4!) = 30
        for c in chords:
            if c == a or c == b:
                continue
            contrib += 30 * S[a] * S[b] * S[c]**4
        # Pattern (1,1,2,2) with two distinct c,d != a,b: 6!/(1!1!2!2!) = 180
        other = [c for c in chords if c != a and c != b]
        for ci in range(len(other)):
            for cj in range(ci+1, len(other)):
                c1, c2 = other[ci], other[cj]
                contrib += 180 * S[a] * S[b] * S[c1]**2 * S[c2]**2

        total += contrib
    return total


def main():
    print("=" * 70)
    print("SIGN MECHANISM: D^4 and D^6 degree-2 Walsh")
    print("=" * 70)

    # PART 1: D^4 sign pattern analysis
    print("\n--- PART 1: D^4 SIGN vs CHI(ab) ---")
    print("D^4_deg2 nonzero iff 3a = +/-b mod p (T1/T2 fire)")
    print("Sign determined by WHICH fires and by C_sum")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")

        nonzero = 0
        zero = 0
        match_d4 = 0
        for a in range(1, m+1):
            for b in range(a+1, m+1):
                val = D4_deg2(a, b, p)
                chi_ab = legendre(a*b, p)
                if abs(val) > 0.01:
                    nonzero += 1
                    sign_d4 = 1 if val > 0 else -1
                    # D^4 sign vs -chi(ab): at p=7, h_hat = -(1/2)*D^4_deg2
                    # sign(h_hat) = -sign(D^4_deg2)
                    if -sign_d4 == chi_ab:
                        match_d4 += 1
                else:
                    zero += 1

        total = nonzero + zero
        print(f"    Nonzero D^4: {nonzero}/{total}, zero: {zero}/{total}")
        print(f"    -sign(D^4) matches chi(ab): {match_d4}/{nonzero}")
        print(f"    Pattern: 3a=b or a=3b mod p gives {nonzero} pairs")

        # What fraction are handled by D^4 alone?
        # At p=7: 3/3 = 100%. At p=11: 4/10 = 40%. Growing sparsity.

    # PART 2: D^6 fills the gaps?
    print("\n--- PART 2: D^6 DEGREE-2 WALSH ---")
    for p in [7, 11]:
        m = (p - 1) // 2
        print(f"\n  p={p}:")

        for a in range(1, m+1):
            for b in range(a+1, m+1):
                d4 = D4_deg2(a, b, p)
                d6 = D6_deg2_numerical(a, b, p)
                chi_ab = legendre(a*b, p)
                print(f"    (a={a},b={b}): D^4={d4:>10.2f}, D^6={d6:>10.2f}, "
                      f"chi(ab)={chi_ab:+d}")

    # PART 3: Full degree-2 Walsh from tr(A^k) for k=4,6,8,10
    print("\n--- PART 3: COMBINED TRACE CONTRIBUTION AT p=11 ---")
    p = 11
    m = (p - 1) // 2

    # (-1/2 + iD)^k expansion coefficients for D^{2n} term in Re part:
    # Re((-1/2+iD)^k) = sum_{j even} C(k,j)*(-1/2)^{k-j}*(iD)^j
    # (iD)^j = i^j * D^j. For j even: i^j = (-1)^{j/2}. Real.
    # So Re = sum_{n=0}^{k/2} C(k,2n)*(-1/2)^{k-2n}*(-1)^n * D^{2n}
    # Coefficient of D^{2n}: C(k,2n)*(-1/2)^{k-2n}*(-1)^n

    def coeff_D2n_in_trk(k, n):
        """Coefficient of D^{2n} in Re((-1/2+iD)^k)"""
        from math import comb
        if 2*n > k:
            return 0
        return comb(k, 2*n) * (-0.5)**(k - 2*n) * (-1)**n

    # tr(A^k) = m^k + sum_{freq} Re(lambda^k)
    # Degree-2 Walsh of sum_freq Re(lambda^k) at {a,b}:
    # = sum_n coeff_D2n * (degree-2 Walsh of sum_freq D^{2n})
    # For n=1 (D^2): degree-2 Walsh is 0 (orthogonality)
    # For n=2 (D^4): computed above
    # For n>=3 (D^6, D^8, ...): computed numerically

    print(f"\n  p={p}, m={m}")
    print(f"  Coefficients of D^{{2n}} in Re((-1/2+iD)^k):")
    for k in [4, 6, 8, 10]:
        coeffs = [coeff_D2n_in_trk(k, n) for n in range(k//2 + 1)]
        print(f"    k={k}: {[f'{c:.4f}' for c in coeffs]}")

    # Now compute the combined degree-2 Walsh
    # H = sum_j 2^j * alpha_j where alpha_j depends on cycle counts
    # H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
    # alpha_1 = c3 + c5 + c7 + c9 + c11 (all odd cycles)
    # But alpha_j counts INDEPENDENT SETS in Omega(T), not simple cycles
    # The key: c_k = (1/k)*tr(A^k), so tr Walsh is cycle-count Walsh

    # For degree-2: only even-k traces contribute. The degree-2 Walsh of H
    # = sum_k (dH/dc_k) * (degree-2 Walsh of c_k)
    # = sum_{even k} (dH/dc_k) * (1/k) * (degree-2 Walsh of tr(A^k))

    # Since we need the actual dH/dc_k weights (OCF-dependent), let me
    # instead compute the full degree-2 Walsh numerically and compare

    print(f"\n  Full numerical degree-2 Walsh at p={p}:")
    pairs = [(s, p - s) for s in range(1, m + 1)]

    # Compute H for all 2^m orientations
    from itertools import combinations

    def build_adj(pp, S):
        S_set = set(S)
        A = [[0]*pp for _ in range(pp)]
        for i in range(pp):
            for s in S_set:
                A[i][(i + s) % pp] = 1
        return A

    def count_ham_paths(A, n):
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if dp[mask][v] == 0:
                    continue
                for w in range(n):
                    if mask & (1 << w):
                        continue
                    if A[v][w]:
                        dp[mask | (1 << w)][w] += dp[mask][v]
        full = (1 << n) - 1
        return sum(dp[full][v] for v in range(n))

    H_vals = {}
    for bits in range(1 << m):
        sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_vals[bits] = H

    # Walsh decomposition: degree-2 only
    for i in range(m):
        for j in range(i+1, m):
            a, b = pairs[i][0], pairs[j][0]
            S_set = frozenset([i, j])

            val = 0
            for bits in range(1 << m):
                chi = 1
                for idx in S_set:
                    if not (bits & (1 << idx)):
                        chi *= -1
                val += H_vals[bits] * chi
            h_hat = val / (1 << m)

            d4 = D4_deg2(a, b, p)
            d6 = D6_deg2_numerical(a, b, p)
            chi_ab = legendre(a*b, p)
            sign_h = 1 if h_hat > 0 else (-1 if h_hat < 0 else 0)

            print(f"    (a={a},b={b}): h_hat={h_hat:>10.2f}, D^4={d4:>8.2f}, "
                  f"D^6={d6:>8.2f}, chi={chi_ab:+d}, sign={sign_h:+d}, "
                  f"match={sign_h==chi_ab}")

    # PART 4: The multiplicative resonance pattern
    print("\n--- PART 4: MULTIPLICATIVE RESONANCE ---")
    print("For pairs (a,b), D^4_deg2 != 0 iff 3a = +/-b mod p.")
    print("The VALUE is always +/- p (up to scaling).")
    print("The SIGN is: +p if 3a = -b (i.e., a+3b = 0), -p if 3a = b.")
    print()
    print("chi(ab) at these pairs:")
    for p in [7, 11, 19, 23]:
        m = (p - 1) // 2
        chi_3 = legendre(3, p)
        print(f"\n  p={p}, chi(3)={chi_3:+d}:")
        for a in range(1, m+1):
            for b in range(a+1, m+1):
                if (3*a - b) % p == 0:
                    chi_ab = legendre(a*b, p)
                    print(f"    3a=b: a={a},b={b}, ab={a*b}, "
                          f"chi(ab)={chi_ab:+d}, chi(3a^2)={legendre(3*a*a,p):+d}")
                if (3*a + b) % p == 0:
                    chi_ab = legendre(a*b, p)
                    print(f"    3a=-b: a={a},b={b}, ab={a*b}, "
                          f"chi(ab)={chi_ab:+d}, chi(-3a^2)={legendre((-3*a*a)%p,p):+d}")
                if (a - 3*b) % p == 0:
                    chi_ab = legendre(a*b, p)
                    print(f"    a=3b: a={a},b={b}, ab={a*b}, "
                          f"chi(ab)={chi_ab:+d}, chi(3b^2)={legendre(3*b*b,p):+d}")
                if (a + 3*b) % p == 0:
                    chi_ab = legendre(a*b, p)
                    print(f"    a=-3b: a={a},b={b}, ab={a*b}, "
                          f"chi(ab)={chi_ab:+d}, chi(-3b^2)={legendre((-3*b*b)%p,p):+d}")


if __name__ == '__main__':
    main()
