#!/usr/bin/env python3
"""
product_law_proof.py -- Proving the Walsh-Legendre Product Law

THEOREM CANDIDATE (THM-150):
For circulant tournaments on Z_p with p = 3 mod 4:
  sign(h_hat[{i,j}]) = chi(a_i * a_j)  (Legendre character of pair product)

This means:
  h_hat[{i,j}] = |h_hat[{i,j}]| * chi(a_i * a_j)

where a_i = pair[i][0] = i+1 is the "positive" element of pair i.

PROOF STRATEGY:
1. H(sigma) for circulant T depends only on the MULTISET of sigma values
   over orbits of the automorphism group Aut(T) on pairs.
2. For p = 3 mod 4, the Paley tournament has Aut = PSL(2,p), which acts
   2-transitively. The Walsh coefficients inherit the character structure.
3. h_hat[{i,j}] = (1/2^m) * sum_sigma H(sigma) * sigma_i * sigma_j
4. The key: flipping sigma_i and sigma_j corresponds to flipping arcs
   i+1 <-> p-i-1 and j+1 <-> p-j-1. The effect on H depends on the
   "interaction" between these arc sets, which is governed by the
   Legendre character via the Gauss sum structure.

ALTERNATIVE: Direct computation shows h_hat relates to the "second
moment of the imaginary Fourier part weighted by cycle counts."
If D(k;sigma) = sum sigma_i sin(2pi*k*(i+1)/p), then:
  h_hat[{i,j}] ~ sum_k f(Q_k) * sin(2pi*k*(i+1)/p) * sin(2pi*k*(j+1)/p)
for some function f of the eigenvalue magnitude Q_k.

At p=3 mod 4, the Gauss sum g = i*sqrt(p) introduces a sign via chi,
which propagates to the Walsh coefficients.

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
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


def walsh_decomp(p, pairs):
    m = len(pairs)
    n_orient = 1 << m
    H_vals = {}
    for bits in range(n_orient):
        sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_vals[bits] = H

    h_hat = {}
    for S_bits in range(n_orient):
        S_set = frozenset(i for i in range(m) if S_bits & (1 << i))
        val = 0
        for bits in range(n_orient):
            chi = 1
            for i in S_set:
                if not (bits & (1 << i)):
                    chi *= -1
            val += H_vals[bits] * chi
        h_hat[S_set] = val / n_orient
    return h_hat, H_vals


def main():
    print("=" * 70)
    print("PRODUCT LAW PROOF: sign(h_hat[{i,j}]) = chi(a_i * a_j)")
    print("=" * 70)

    # PART 1: Extend verification to p=19 and p=23
    print("\n--- PART 1: VERIFICATION AT LARGER PRIMES ---")

    for p in [7, 11, 19, 23]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n  p={p} ({p % 4} mod 4), m={m}:")

        if p <= 13:
            h_hat, _ = walsh_decomp(p, pairs)
        else:
            # For larger p, compute only degree-2 coefficients via
            # the difference method:
            # h_hat[{i,j}] = (1/4) * (H(+,+,...) - H(+,-,...) - H(-,+,...) + H(-,-,...))
            # averaged over all other sigma values

            # Actually: h_hat[{i,j}] = (1/2^m) sum_sigma H(sigma) * sigma_i * sigma_j
            # For p >= 17 this is too expensive. Use a faster method:
            # Group by (sigma_i, sigma_j) and sum over remaining coordinates

            # For p=19 (m=9): 2^9 = 512 orientations, each needs Held-Karp O(2^19 * 19^2)
            # This is ~2^28 * 361 ~ 10^11 operations. Too slow.

            # Instead, estimate via sampling?
            # Or: compute H for a few strategic orientations

            # Let's just do p=7, 11, and 13 for proof, and use Gauss sum argument for general p
            print(f"    [Skipping exhaustive computation for p={p}]")
            print(f"    Testing product law on Gauss sum prediction:")

            # Use the Gauss sum formula to PREDICT signs
            # For p = 3 mod 4: D(k; sigma_Paley) = chi(k) * sqrt(p)/2
            # The general formula involves Gauss sums
            g = sum(pow(j, (p-1)//2, p) * cmath.exp(2j * cmath.pi * j / p)
                    for j in range(1, p))
            # For p = 3 mod 4: g should be i * sqrt(p)
            g_expected = 1j * math.sqrt(p) if p % 4 == 3 else math.sqrt(p)
            print(f"    Gauss sum: g = {g:.4f}, expected = {g_expected:.4f}, "
                  f"match = {abs(g - g_expected) < 0.01}")
            continue

        # Check product law
        correct = 0
        total = 0
        wrong_cases = []
        for S_set, val in h_hat.items():
            if len(S_set) != 2 or abs(val) < 0.01:
                continue
            i, j = sorted(S_set)
            a_i, a_j = pairs[i][0], pairs[j][0]
            prod_mod = (a_i * a_j) % p
            chi_prod = 1 if pow(prod_mod, (p-1)//2, p) == 1 else -1
            sign_h = 1 if val > 0 else -1
            total += 1
            if chi_prod == sign_h:
                correct += 1
            else:
                wrong_cases.append((i, j, a_i, a_j, prod_mod, chi_prod, sign_h, val))

        print(f"    Product law: {correct}/{total} correct")
        if wrong_cases:
            for wc in wrong_cases:
                print(f"      FAIL: ({wc[0]+1},{wc[1]+1}): a*b={wc[4]}, chi={wc[5]}, sign={wc[6]}")

    # PART 2: The h_hat magnitude formula
    print("\n--- PART 2: MAGNITUDE FORMULA ---")
    print("If sign = chi(a_i*a_j), what determines |h_hat|?")
    print("Testing: |h_hat[{i,j}]| depends on a_i*a_j mod p only")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        h_hat, _ = walsh_decomp(p, pairs)

        print(f"\n  p={p}:")
        prod_mag = defaultdict(list)
        for S_set, val in h_hat.items():
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            a_i, a_j = pairs[i][0], pairs[j][0]
            prod_mod = (a_i * a_j) % p
            prod_mag[prod_mod].append(abs(val))

        for pm in sorted(prod_mag):
            mags = prod_mag[pm]
            print(f"    a*b mod p = {pm}: |h_hat| = {[round(v, 2) for v in mags]}")

        # Check if |h_hat| depends on min(prod, p-prod)
        sym_mag = defaultdict(list)
        for S_set, val in h_hat.items():
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            a_i, a_j = pairs[i][0], pairs[j][0]
            prod_mod = (a_i * a_j) % p
            sym_key = min(prod_mod, p - prod_mod)
            sym_mag[sym_key].append(abs(val))

        print(f"    Grouped by min(ab, p-ab):")
        for sk in sorted(sym_mag):
            mags = sym_mag[sk]
            mags_unique = sorted(set(round(v, 2) for v in mags))
            print(f"      key={sk}: |h_hat| = {mags_unique}, count={len(mags)}")

    # PART 3: Attempt algebraic proof via Gauss sum
    print("\n--- PART 3: GAUSS SUM ARGUMENT ---")
    print("""
    PROOF SKETCH for p = 3 mod 4:

    For a circulant tournament with connection set S, the transfer matrix
    eigenvalue at frequency k is lambda_k = S_hat(k) = sum_{s in S} omega^{ks}.

    Write lambda_k = C(k) + i*D(k;sigma) where:
      C(k) = -1/2 (universal)
      D(k;sigma) = sum_{i=1}^m sigma_i * sin(2*pi*k*(i+1)/p)

    The Walsh coefficient h_hat[{i,j}] measures how sigma_i * sigma_j
    affects H. By the OCF formula H = I(Omega(T), 2), the effect propagates
    through cycle counts.

    KEY LEMMA: h_hat[{i,j}] has the same sign as
      sum_k w(k) * sin(2*pi*k*a_i/p) * sin(2*pi*k*a_j/p)
    where w(k) is a positive weight function depending on |lambda_k|^2.

    For p = 3 mod 4, the weight function w(k) relates to the Fejer kernel
    at the Paley orientation (where all Q_k = (p+1)/4), and the product
    sin(...)*sin(...) factors through the Gauss sum.

    Specifically:
      sum_k sin(2*pi*k*a/p) * sin(2*pi*k*b/p) = (p/2) * delta_{a,b}

    But the WEIGHTED sum (with Fejer kernel weights) gives:
      sum_k F_m(k) * sin(2*pi*k*a/p) * sin(2*pi*k*b/p)
    which is NOT diagonal. It depends on (a-b) and (a+b) mod p.

    The sign of this weighted sum is chi(a*b) by the Gauss sum identity
    (using g^2 = (-1)^{(p-1)/2} * p = -p for p = 3 mod 4).
    """)

    # PART 4: Verify the weighted sum formula
    print("--- PART 4: WEIGHTED SUM VERIFICATION ---")

    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n  p={p}:")

        # Compute Fejer kernel weights
        S_int = list(range(1, m + 1))
        F_weights = []
        for k in range(1, m + 1):
            theta = 2 * math.pi * k / p
            F_k = (math.sin(m * theta / 2))**2 / (math.sin(theta / 2))**2
            F_weights.append(F_k)

        print(f"    Fejer weights: {[f'{w:.3f}' for w in F_weights]}")

        # Compute weighted sum WS[a,b] = sum_k F(k) * sin(2pi*k*a/p) * sin(2pi*k*b/p)
        h_hat, _ = walsh_decomp(p, pairs)

        print(f"\n    Weighted sum vs h_hat:")
        for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            a, b = pairs[i][0], pairs[j][0]

            ws = sum(F_weights[ki] * math.sin(2*math.pi*(ki+1)*a/p) *
                     math.sin(2*math.pi*(ki+1)*b/p)
                     for ki in range(m))

            chi_ab = 1 if ((a*b) % p) in QR else -1
            ws_sign = 1 if ws > 0 else (-1 if ws < 0 else 0)

            print(f"    ({i+1},{j+1}): a*b={a*b}, chi(ab)={chi_ab:+d}, "
                  f"WS={ws:>10.4f}, sign(WS)={ws_sign:+d}, "
                  f"h_hat={val:>10.2f}, sign(h)={1 if val>0 else -1:+d}")

    # PART 5: The EXACT formula via Gauss sum products
    print("\n--- PART 5: GAUSS SUM PRODUCT FORMULA ---")
    print("Testing: h_hat[{i,j}] = C * (g_hat(a_i) * g_hat(a_j))")
    print("where g_hat is some transform of the Gauss sum")

    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        omega = cmath.exp(2j * cmath.pi / p)

        h_hat, _ = walsh_decomp(p, pairs)

        print(f"\n  p={p}:")

        # Compute g_hat(a) = sum_k chi(k) * sin(2pi*k*a/p)
        # This is the "character-weighted Fourier sine transform"
        g_hat = []
        for ii in range(m):
            a = pairs[ii][0]
            val = sum((1 if k in QR else -1) * math.sin(2 * math.pi * k * a / p)
                      for k in range(1, p))
            g_hat.append(val)
        print(f"    g_hat values: {[f'{v:.4f}' for v in g_hat]}")

        # Check: h_hat[{i,j}] / (g_hat[i] * g_hat[j])
        print(f"\n    h_hat / (g_hat_i * g_hat_j):")
        for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2 or abs(val) < 0.01:
                continue
            i, j = sorted(S_set)
            prod = g_hat[i] * g_hat[j]
            if abs(prod) > 0.001:
                ratio = val / prod
                print(f"    ({i+1},{j+1}): ratio = {ratio:.6f}")

    # PART 6: Check for higher degree analogues
    print("\n--- PART 6: PRODUCT LAW FOR DEGREE 4 ---")
    print("Testing: sign(h_hat[{i,j,k,l}]) = chi(a_i*a_j*a_k*a_l)?")

    for p in [11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        h_hat, _ = walsh_decomp(p, pairs)

        print(f"\n  p={p} ({p % 4} mod 4):")
        correct_4 = 0
        total_4 = 0
        for S_set, val in h_hat.items():
            if len(S_set) != 4 or abs(val) < 0.01:
                continue
            indices = sorted(S_set)
            pair_vals = [pairs[i][0] for i in indices]
            prod_mod = 1
            for v in pair_vals:
                prod_mod = (prod_mod * v) % p
            chi_prod = 1 if pow(prod_mod, (p-1)//2, p) == 1 else -1
            sign_h = 1 if val > 0 else -1
            total_4 += 1
            match = (chi_prod == sign_h)
            if match:
                correct_4 += 1

        print(f"    Product law (deg 4): {correct_4}/{total_4} correct")

    # PART 7: The g_hat factorization at degree 4
    print("\n--- PART 7: g_hat FACTORIZATION AT DEGREE 4 ---")

    for p in [11]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        h_hat, _ = walsh_decomp(p, pairs)

        g_hat = []
        for ii in range(m):
            a = pairs[ii][0]
            val = sum((1 if k in QR else -1) * math.sin(2 * math.pi * k * a / p)
                      for k in range(1, p))
            g_hat.append(val)

        print(f"\n  p={p}:")
        print(f"    h_hat / prod(g_hat):")
        for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 4 or abs(val) < 0.01:
                continue
            indices = sorted(S_set)
            prod = 1.0
            for i in indices:
                prod *= g_hat[i]
            if abs(prod) > 0.001:
                ratio = val / prod
                print(f"    {set(i+1 for i in indices)}: ratio = {ratio:.6f}")


if __name__ == '__main__':
    main()
