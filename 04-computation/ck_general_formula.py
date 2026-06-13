#!/usr/bin/env python3
"""
ck_general_formula.py -- Verify and extend the general Walsh formula for c_k

THEOREM (j=2 contribution): For circulant tournaments on Z_p, for k <= 5:
  h_hat_{c_k}[{a,b}] = F(k) * eps_3(a,b) * p / 8    if q(a,b) = 3
                      = 0                               otherwise
where F(k) = (-1)^{k+1} * (k-1)(k-2)(k-3) / (3 * 2^{k-4}).

This gives: F(3)=0, F(4)=-2, F(5)=4.

Magnitudes: |h_hat_{c_k}| = |F(k)| * p / 8.
  k=3: 0
  k=4: p/4
  k=5: p/2

EXTENSION: For k >= 6 (if c_k = tr/k), the j=3 term contributes q=5 content:
  h_hat_{c_k}^{(j=3)} involves sum sin^5 * sin, which resolves at q=5.

This script:
1. Verifies the j=2 formula at all primes 7 <= p <= 23
2. Computes the j=3 contribution analytically
3. Checks whether c_6 = tr/6 (are all closed 6-walks simple?)

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def resonance_level(a, b, p):
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def resonance_sign_q3(ga, gb, p):
    """Sign eps for q=3 resonance."""
    if (3 * ga - gb) % p == 0:
        return -1
    elif (3 * ga + gb) % p == 0:
        return +1
    elif (ga - 3 * gb) % p == 0:
        return -1
    elif (ga + 3 * gb) % p == 0:
        return +1
    return 0


def resonance_sign_q5(ga, gb, p):
    """Sign eps for q=5 resonance."""
    if (5 * ga - gb) % p == 0:
        return -1
    elif (5 * ga + gb) % p == 0:
        return +1
    elif (ga - 5 * gb) % p == 0:
        return -1
    elif (ga + 5 * gb) % p == 0:
        return +1
    return 0


def F_j2(k):
    """j=2 coefficient: (-1)^{k+1} * (k-1)(k-2)(k-3) / (3 * 2^{k-4})"""
    if k < 4:
        return 0  # C(k,4)=0 for k<4
    return (-1)**(k+1) * (k-1)*(k-2)*(k-3) / (3 * 2**(k-4))


def G_j3(k):
    """j=3 coefficient for q=5 contribution.

    From Re(lambda^k): D^6 term has coefficient C(k,6)*C^{k-6}*(-1)^3.
    Walsh_2[D^6][{a,b}] has a term involving sum sin^5(ga*t)*sin(gb*t),
    which resolves at q=5.

    sin^5(x) = (10*sin(x) - 5*sin(3x) + sin(5x)) / 16

    So sum sin^5(ga*t)*sin(gb*t) = (1/16) * [10*sum sin(ga)sin(gb) - 5*sum sin(3ga)sin(gb) + sum sin(5ga)sin(gb)]

    The first two sums are zero (q=1 resolved to 0, q=3 resolved separately).
    The q=5 part: sum sin(5ga*t)*sin(gb*t) = +-p/2 if 5ga = +-gb mod p.

    So the j=3 q=5 contribution is:
    (1/k) * C(k,6) * (-1/2)^{k-6} * (-1)^3 * [term from D^6 Walsh involving s_a^5 s_b]

    Need to work out Walsh_2[D^6] more carefully.
    """
    if k < 6:
        return 0
    # D^6 = (sum sigma_j s_j)^6
    # Walsh_2[D^6][{a,b}] involves terms where the product of 6 sigmas * sigma_a * sigma_b
    # has each index appearing an even number of times.
    #
    # The q=5 contribution comes from terms like sigma_a^5 * sigma_b * sigma_c^0... = sigma_a * sigma_b
    # (since sigma^2 = 1 for all sigma in {+-1}).
    # So we need: n_a + 1 and n_b + 1 both even, and all other n_c even.
    # n_a odd, n_b odd, sum n_j = 6.
    #
    # Partitions:
    # (5,1,0,...): n_a=5, n_b=1. 6!/(5!*1!) = 6 terms per such choice. => s_a^5 * s_b
    # (1,5,0,...): n_a=1, n_b=5. 6 terms. => s_a * s_b^5
    # (3,3,0,...): n_a=3, n_b=3. 6!/(3!3!) = 20 terms. => s_a^3 * s_b^3
    # (3,1,2,...): n_a=3, n_b=1, n_c=2. 6!/(3!1!2!) = 60 per c. Sum over c: 60*s_a^3*s_b*sum_c s_c^2
    # (1,3,2,...): similarly => 60*s_a*s_b^3*sum_c s_c^2
    # (1,1,4,...): n_a=1, n_b=1, n_c=4. 6!/(1!1!4!) = 30 per c. => 30*s_a*s_b*s_c^4
    # (1,1,2,2,...): n_a=1, n_b=1, n_c=2, n_d=2. 6!/(1!1!2!2!) = 180 per (c,d). => 180*s_a*s_b*s_c^2*s_d^2
    #
    # This is getting complex. Let me just compute the s_a^5*s_b + s_a*s_b^5 terms.
    # Coefficient of s_a^5*s_b in Walsh_2[D^6]: 6 (from 6!/5!1!)
    # Coefficient of s_a*s_b^5: 6

    # The full Walsh_2[D^6] also has terms from s_a^3*s_b^3 (gives q=3 contribution)
    # and lower-order terms (give q=1 contribution = 0).

    # So the q=5-ONLY part of Walsh_2[D^6] comes from:
    # 6*(s_a^5 s_b + s_a s_b^5) + contributions from cross terms...
    # Actually, other terms like s_a^3 s_b s_c^2 also contribute to q=5 if
    # 3*ga = +-5*gb mod p, but those are higher-resonance interactions.
    # Let me focus on the DIRECT q=5 terms.

    # For the q=5 contribution from the s_a^5 s_b + s_a s_b^5 terms:
    # sum_t sin^5(ga*t)*sin(gb*t) = (1/16)[10*0 - 5*(±p/2 if q=3) + (±p/2 if q=5)]
    # For a pure q=5 pair (q≠3): the q=3 part is 0, so:
    # sum sin^5(ga*t)*sin(gb*t) = (1/16)(±p/2) = ±p/32

    # Similarly: sum sin(ga*t)*sin^5(gb*t) = ±p/32

    # Total from 6*(s_a^5 s_b + s_a s_b^5):
    # 6 * 2 * (±p/32) = ±6p/16 = ±3p/8

    # Wait, I need to account for both terms:
    # 6 * sum s_a^5 s_b + 6 * sum s_a s_b^5

    # For a q=5 pair where 5ga ≡ gb but ga ≢ 5gb:
    # sum s_a^5 s_b = (1/16)(sin(5ga*t)sin(gb*t) part) = (1/16)(p/2) = p/32
    # sum s_a s_b^5: need 5gb ≡ ±ga. If 5ga ≡ gb, then ga ≡ gb/5. 5gb ≡ 5*(5ga) = 25ga.
    # 25 mod p: varies by prime. At p=11: 25 ≡ 3. So 5gb ≡ 3ga. Need 3ga ≡ ±ga → only if 3≡±1, false.
    # So for most q=5 pairs, sum s_a s_b^5 = 0.
    # Exception: if BOTH 5ga ≡ ±gb AND ga ≡ ±5gb (same pair).

    # Actually, if q=5: one of 5ga≡gb, 5ga≡-gb, ga≡5gb, ga≡-5gb holds.
    # The fifth power: sin^5 has sin(5x) component. sum sin(5ga*t)sin(gb*t):
    # 5ga and gb: if 5ga≡gb (the resonance), sum = p/2.
    # For s_a s_b^5: sin(ga*t)sin^5(gb*t) has sin(ga*t)sin(5gb*t) part.
    # 5gb and ga: if 5ga≡gb, then ga≡gb/5. 5gb vs ga: 5gb vs gb/5. 25gb vs gb (mod p). Need 24≡0 mod p. Only for p|24, not for p≥7.
    # So generically, sum s_a s_b^5 contributes 0 for q=5 pairs.

    # Therefore: from s_a^5 s_b + s_a s_b^5, only ONE of the two sums is nonzero for q=5.
    # The contribution is 6 * (±p/32) = ±6p/32 = ±3p/16.

    # Total j=3 coefficient:
    # (1/k) * C(k,6) * (-1/2)^{k-6} * (-1)^3 * 6 * (±p/32)

    # Hmm, this is getting complicated. Let me just compute numerically.
    return None


def compute_c_k_trace(S, p, k):
    """c_k via trace: tr(A^k)/k."""
    omega = cmath.exp(2j * cmath.pi / p)
    total = 0
    for t in range(p):
        lam = sum(omega**(s * t) for s in S)
        total += lam**k
    return total.real / k


def main():
    print("=" * 70)
    print("GENERAL WALSH FORMULA VERIFICATION")
    print("=" * 70)

    print("\nF(k) coefficients (j=2 contribution):")
    for k in range(3, 12):
        print(f"  F({k}) = {F_j2(k):.4f}")

    # Verify at each prime: h_{c_k} = F(k) * eps * p/8 for k=3,4,5
    for p in [7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        if m > 11:
            continue
        n_orient = 1 << m
        pairs = [(j, p - j) for j in range(1, m + 1)]
        chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]

        print(f"\n{'='*60}")
        print(f"p={p}, m={m}")
        print("=" * 60)

        # Compute c_k = tr/k for k=3,4,5,6,7 for all orientations
        all_ck = defaultdict(list)  # k -> [values]
        for bits in range(n_orient):
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                       for i in range(m))
            for k in [3, 4, 5, 6]:
                val = compute_c_k_trace(S, p, k)
                all_ck[k].append(round(val))

        sigma = {}
        for bits in range(n_orient):
            sigma[bits] = tuple(1 if bits & (1 << i) else -1 for i in range(m))

        # For each k and each chord pair, check formula
        print(f"\n  k=3,4,5 verification (j=2 only):")
        max_err = {3: 0, 4: 0, 5: 0}
        for k in [3, 4, 5]:
            for a, b in chord_pairs:
                ga, gb = a + 1, b + 1
                q = resonance_level(ga, gb, p)

                # Numerical
                h_hat = sum(all_ck[k][bits] * sigma[bits][a] * sigma[bits][b]
                            for bits in range(n_orient)) / n_orient

                # Predicted (j=2 only)
                if q == 3:
                    eps = resonance_sign_q3(ga, gb, p)
                    predicted = F_j2(k) * eps * p / 8
                else:
                    predicted = 0

                err = abs(h_hat - predicted)
                max_err[k] = max(max_err[k], err)

            status = "PASS" if max_err[k] < 0.01 else "FAIL"
            print(f"    k={k}: max_error={max_err[k]:.6f} [{status}]")

        # For k=6: check if c_6 = tr/6 by comparing with exact count
        # At k=6, non-simple walks = 3+3 decompositions exist.
        # But c_6 counts even-length simple cycles. The trace counts even-length closed walks.
        # Non-simple 6-walks: decompose into two 3-cycles sharing a vertex.
        # For each vertex v, count pairs of 3-cycles through v.
        print(f"\n  k=6 analysis (are 6-walks simple?):")
        # c_6 via exact cycle enumeration
        if p <= 13:
            S_test = list(range(1, m + 1))  # base orientation
            A = [[0]*p for _ in range(p)]
            for v in range(p):
                for s in S_test:
                    A[v][(v + s) % p] = 1

            from itertools import permutations
            def held_karp_general(A, verts):
                k = len(verts)
                start = 0
                dp = {}
                dp[(1 << start, start)] = 1
                for mask in range(1, 1 << k):
                    if not (mask & (1 << start)):
                        continue
                    for v in range(k):
                        if not (mask & (1 << v)):
                            continue
                        key = (mask, v)
                        if key not in dp or dp[key] == 0:
                            continue
                        cnt = dp[key]
                        for w in range(k):
                            if mask & (1 << w):
                                continue
                            if A[verts[v]][verts[w]]:
                                nkey = (mask | (1 << w), w)
                                dp[nkey] = dp.get(nkey, 0) + cnt
                full = (1 << k) - 1
                total = 0
                for v in range(k):
                    if v == start:
                        continue
                    key = (full, v)
                    if key in dp and dp[key] > 0:
                        if A[verts[v]][verts[start]]:
                            total += dp[key]
                return total

            c6_exact = sum(held_karp_general(A, list(sub))
                           for sub in combinations(range(p), 6))
            c6_trace = compute_c_k_trace(S_test, p, 6)
            N6 = round(c6_trace * 6 - 6 * c6_exact)
            print(f"    c6_exact={c6_exact}, tr/6={c6_trace:.1f}, N6={N6}")
            print(f"    6-walks are simple? {'YES' if N6 == 0 else 'NO (N6=' + str(N6) + ')'}")

        # Degree-2 Walsh of c_6 (from trace, which = exact at k=6 iff N6=0)
        print(f"\n  k=6 Walsh (from trace):")
        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)

            h_hat = sum(all_ck[6][bits] * sigma[bits][a] * sigma[bits][b]
                        for bits in range(n_orient)) / n_orient

            if abs(h_hat) > 0.01:
                eps3 = resonance_sign_q3(ga, gb, p) if q == 3 else 0
                eps5 = resonance_sign_q5(ga, gb, p) if q == 5 else 0
                pred_j2 = F_j2(6) * eps3 * p / 8 if q == 3 else 0
                diff = h_hat - pred_j2
                print(f"    ({a},{b}) q={q}: h_hat={h_hat:>10.4f}, "
                      f"j2_pred={pred_j2:>10.4f}, residual={diff:>10.4f}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
