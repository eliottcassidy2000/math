#!/usr/bin/env python3
"""
gauss_sum_sign_law.py -- Gauss sum approach to Walsh-Legendre Sign Law

THE SIGN LAW (p = 3 mod 4):
  sign(h_hat[S]) = chi(prod_{i in S} (i+1))

where chi is the Legendre symbol mod p, S is a subset of {0,...,m-1},
and h_hat[S] is the Walsh-Hadamard coefficient of H(sigma).

APPROACH: Connect Walsh coefficients to Gauss sums via the Fourier
decomposition S_hat(k;sigma) = -1/2 + i * D(k;sigma).

The key insight is that D(k;sigma) = sum_i sigma_i * sin(2*pi*k*(i+1)/p).
At the Paley point, sigma_i = chi(i+1), so:
  D(k; sigma_Paley) = sum_i chi(i+1) * sin(2*pi*k*(i+1)/p)
                     = Im(sum_{j=1}^m chi(j) * omega^{kj})
which is the imaginary part of a partial Gauss sum.

GAUSS SUM: g = sum_{j=1}^{p-1} chi(j) * omega^j
  For p = 3 mod 4: g = i * sqrt(p)  (pure imaginary)
  For p = 1 mod 4: g = sqrt(p)  (real)

Can we express h_hat[S] in terms of products of Gauss-like sums?

Author: kind-pasteur-2026-03-12-S59c
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


def legendre(a, p):
    """Legendre symbol (a/p)"""
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def gauss_sum(p):
    """Compute the quadratic Gauss sum g = sum_{j=0}^{p-1} chi(j) * omega^j"""
    omega = cmath.exp(2j * cmath.pi / p)
    g = sum(legendre(j, p) * omega**j for j in range(p))
    return g


def main():
    print("=" * 70)
    print("GAUSS SUM STRUCTURE AND WALSH-LEGENDRE SIGN LAW")
    print("=" * 70)

    # Part 1: Gauss sum values
    print("\n--- PART 1: GAUSS SUMS ---")
    print(f"{'p':>4} {'mod4':>5} {'Re(g)':>12} {'Im(g)':>12} {'|g|^2':>10} {'sqrt(p)':>10} {'g/i*sqrt(p)':>12}")
    print("-" * 70)

    for p in [5, 7, 11, 13, 17, 19, 23, 29]:
        g = gauss_sum(p)
        sqp = math.sqrt(p)
        if p % 4 == 3:
            ratio = g / (1j * sqp)
        else:
            ratio = g / sqp
        print(f"{p:>4} {p%4:>5} {g.real:>12.6f} {g.imag:>12.6f} {abs(g)**2:>10.4f} {sqp:>10.4f} {ratio.real:>+12.6f}")

    # Part 2: Partial Gauss sums and D(k; sigma_Paley)
    print("\n--- PART 2: D(k; sigma_Paley) = Im(partial Gauss sum) ---")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n  p={p}, m={m}, QR={sorted(QR)}")

        # sigma_Paley: sigma_i = chi(i+1) for i=0,...,m-1
        sigma_paley = tuple(1 if (i+1) in QR else -1 for i in range(m))
        print(f"  sigma_Paley = {sigma_paley}")

        # D(k; sigma_Paley)
        for k in range(1, m+1):
            D_k = sum(sigma_paley[i] * math.sin(2 * math.pi * k * (i+1) / p)
                      for i in range(m))

            # Partial Gauss sum: sum_{j=1}^m chi(j) * omega^{kj}
            partial_g = sum(legendre(j, p) * omega**(k*j) for j in range(1, m+1))

            # Full Gauss sum rotated: g_k = sum_{j=1}^{p-1} chi(j) * omega^{kj} = chi(k) * g
            # (Gauss sum multiplication property)
            g_k = sum(legendre(j, p) * omega**(k*j) for j in range(1, p))
            g = gauss_sum(p)
            chi_k = legendre(k, p)

            print(f"    k={k}: D={D_k:>+10.6f}, Im(partial_g)={partial_g.imag:>+10.6f}, "
                  f"g_k/g={g_k/g:.6f}, chi(k)={chi_k:>+2}")

    # Part 3: Connection between D and Gauss sums
    print("\n--- PART 3: D(k;sigma_Paley) FROM GAUSS SUM ---")
    print("Testing: D(k;sigma_Paley) = Im(chi(k)*g)/2 = chi(k)*sqrt(p)/2 for p=3 mod 4")

    for p in [7, 11, 19, 23]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        omega = cmath.exp(2j * cmath.pi / p)
        g = gauss_sum(p)

        sigma_paley = tuple(1 if (i+1) in QR else -1 for i in range(m))

        print(f"\n  p={p} ({p%4} mod 4):")
        match = True
        for k in range(1, m+1):
            D_k = sum(sigma_paley[i] * math.sin(2*math.pi*k*(i+1)/p) for i in range(m))

            # Full Gauss sum: sum_{j=1}^{p-1} chi(j) omega^{kj} = chi(k) * g
            chi_k = legendre(k, p)
            predicted_D = chi_k * g.imag / 2  # But is this right?

            # Actually, D = sum_{j=1}^m chi(j) sin(2pi*k*j/p)
            # = Im(sum_{j=1}^m chi(j) omega^{kj})
            # = Im(partial Gauss sum with j only up to m)

            # The full sum: sum_{j=1}^{p-1} chi(j) sin(2pi*kj/p) = Im(g_k) = Im(chi(k)*g)
            # For p=3 mod 4: g = i*sqrt(p), so chi(k)*g = chi(k)*i*sqrt(p)
            # Im(chi(k)*i*sqrt(p)) = chi(k)*sqrt(p) (since chi(k) is real)
            full_D = sum(legendre(j,p) * math.sin(2*math.pi*k*j/p) for j in range(1,p))

            # But our D only sums j=1..m, not j=1..p-1
            # sum_{j=m+1}^{p-1} chi(j) sin(2pi*kj/p) = ?
            # j -> p-j: sin(2pi*k*(p-j)/p) = -sin(2pi*kj/p)
            # chi(p-j) = chi(-1)*chi(j) = (-1)^{(p-1)/2} * chi(j)
            # For p=3 mod 4: chi(-1)=-1, so chi(p-j)=-chi(j)
            # => sum_{j=m+1}^{p-1} = sum_{j=1}^{m} chi(p-j)*sin(2pi*k*(p-j)/p)
            #                       = sum_{j=1}^m (-chi(j))*(-sin(2pi*kj/p))
            #                       = sum_{j=1}^m chi(j)*sin(2pi*kj/p) = D_k
            # So full_D = 2*D_k for p=3 mod 4!

            predicted = full_D / 2
            err = abs(D_k - predicted)

            sqp = math.sqrt(p)
            predicted2 = chi_k * sqp / 2  # Im(chi(k)*g) / 2 = chi(k)*sqrt(p)/2
            err2 = abs(D_k - predicted2)

            if err > 1e-8 or err2 > 1e-8:
                match = False

            print(f"    k={k}: D={D_k:>+10.6f}, full_D/2={predicted:>+10.6f} (err={err:.2e}), "
                  f"chi(k)*sqrt(p)/2={predicted2:>+10.6f} (err={err2:.2e})")

        if match:
            print(f"  ==> CONFIRMED: D(k;sigma_Paley) = chi(k)*sqrt(p)/2 for all k")

    # Part 4: Consequence for Q_k at Paley
    print("\n--- PART 4: Q_k(Paley) = 1/4 + p/4 = (p+1)/4 ---")
    print("Since D(k) = chi(k)*sqrt(p)/2, D^2 = p/4 for all k")
    print("=> Q_k = 1/4 + p/4 = (p+1)/4 for all k")
    print("This PROVES the flat Fourier spectrum from Gauss sum theory!")

    # Part 5: What about Walsh coefficients?
    print("\n--- PART 5: WALSH-TO-GAUSS CONNECTION ---")
    print("h_hat[S] = (1/2^m) * sum_sigma H(sigma) * prod_{i in S} sigma_i")
    print("At sigma_Paley: prod_{i in S} sigma_i = prod_{i in S} chi(i+1) = chi(prod(i+1))")
    print("=> H(sigma_Paley) = sum_S h_hat[S] * chi(prod gaps in S)")
    print("")
    print("The Sign Law says: h_hat[S] * chi(prod gaps in S) = |h_hat[S]|")
    print("i.e., every Walsh term is POSITIVE at sigma_Paley")
    print("=> H(sigma_Paley) = sum |h_hat[S]| >= sum h_hat[S]*prod(sigma) for any sigma")

    for p in [7, 11]:
        if p > 11:
            continue  # too slow
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        pairs = [(s, p - s) for s in range(1, m + 1)]

        # Compute all H
        H_dict = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            H_dict[sigma] = H

        # Walsh transform
        n = 1 << m
        h_hat = {}
        for bits in range(n):
            S_idx = tuple(i for i in range(m) if bits & (1 << i))
            coeff = 0
            for sigma_bits in range(n):
                sigma = tuple((1 if sigma_bits & (1 << i) else -1) for i in range(m))
                prod = 1
                for i in S_idx:
                    prod *= sigma[i]
                coeff += H_dict[sigma] * prod
            h_hat[S_idx] = coeff / n

        print(f"\n  p={p}, m={m}:")

        # Check sign law
        sigma_paley = tuple(1 if (i+1) in QR else -1 for i in range(m))
        all_match = True
        for S_idx, coeff in sorted(h_hat.items(), key=lambda x: len(x[0])):
            if len(S_idx) == 0 or len(S_idx) % 2 == 1:
                continue
            gaps = [i+1 for i in S_idx]
            prod_gap = 1
            for g in gaps:
                prod_gap *= g
            chi_val = legendre(prod_gap, p)
            sign_coeff = 1 if coeff > 0 else (-1 if coeff < 0 else 0)
            match = (sign_coeff == chi_val)
            if not match:
                all_match = False

            # Check: is |h_hat[S]| related to a Gauss sum product?
            print(f"    S={S_idx}, deg={len(S_idx)}, gaps={gaps}, "
                  f"prod={prod_gap}, chi={chi_val:>+2}, "
                  f"h_hat={coeff:>+12.4f}, sign={sign_coeff:>+2}, match={'Y' if match else 'N'}")

        print(f"  Sign Law holds: {all_match}")

    # Part 6: prod Q_k primality investigation
    print("\n" + "=" * 70)
    print("PART 6: PROD Q_k VALUES AND PRIMALITY")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        pairs = [(s, p - s) for s in range(1, m + 1)]

        prod_vals = set()
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

            Q_vals = []
            for k in range(1, m + 1):
                val = sum(omega ** (k * s) for s in S)
                Q_vals.append(abs(val)**2)

            prod_Q = round(math.prod(Q_vals))
            prod_vals.add(prod_Q)

        print(f"\n  p={p}, m={m}:")
        for val in sorted(prod_vals):
            # Check primality
            is_prime = val > 1
            if val > 1:
                for d in range(2, int(math.sqrt(val)) + 1):
                    if val % d == 0:
                        is_prime = False
                        break

            # Factor
            if not is_prime and val > 1:
                factors = []
                n = val
                for d in range(2, int(math.sqrt(n)) + 1):
                    while n % d == 0:
                        factors.append(d)
                        n //= d
                if n > 1:
                    factors.append(n)
                print(f"    prod Q_k = {val} = {'*'.join(map(str,factors))} (COMPOSITE)")
            elif is_prime:
                print(f"    prod Q_k = {val} (PRIME)")
            else:
                print(f"    prod Q_k = {val}")

    # Part 7: The resultant / discriminant connection
    print("\n" + "=" * 70)
    print("PART 7: CHARACTERISTIC POLYNOMIAL AND DISCRIMINANT")
    print("=" * 70)
    print("Q_1,...,Q_m are roots of a degree-m polynomial with integer coefficients.")
    print("prod Q_k = (-1)^m * constant term = e_m (the last ESF).")
    print("The DISCRIMINANT of this polynomial might have number-theoretic meaning.")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        pairs = [(s, p - s) for s in range(1, m + 1)]

        # Just compute for Interval
        sigma = tuple([1]*m)
        S = sorted(pairs[i][0] for i in range(m))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        # Compute characteristic polynomial: prod (t - Q_k)
        # Expand via Newton's identities
        e_vals = [1.0]
        for j in range(1, m + 1):
            ej = sum(math.prod(Q_vals[i] for i in combo)
                    for combo in combinations(range(m), j))
            e_vals.append(round(ej))

        # Discriminant = prod_{i<j} (Q_i - Q_j)^2
        disc = 1.0
        for i in range(m):
            for j in range(i+1, m):
                disc *= (Q_vals[i] - Q_vals[j])**2

        print(f"\n  p={p} (Interval), m={m}:")
        print(f"    Q_vals = {[f'{q:.6f}' for q in Q_vals]}")
        print(f"    e_j = {[int(e) for e in e_vals[1:]]}")
        print(f"    prod Q_k = {round(math.prod(Q_vals))}")
        print(f"    Discriminant = {disc:.4f} (rounded: {round(disc)})")

        # Check if disc is a perfect square times p^something
        disc_int = round(disc)
        if disc_int > 0:
            # Factor out powers of p
            n = disc_int
            p_power = 0
            while n % p == 0:
                n //= p
                p_power += 1
            print(f"    disc / p^{p_power} = {n}")

    # Part 8: Deep structure - can we express h_hat[S] via Jacobi sums?
    print("\n" + "=" * 70)
    print("PART 8: WALSH COEFFICIENTS — MAGNITUDE STRUCTURE")
    print("=" * 70)
    print("If h_hat[S] ~ product of Gauss/Jacobi sums, |h_hat[S]| ~ p^{|S|/2}")

    for p in [7, 11]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        pairs = [(s, p - s) for s in range(1, m + 1)]

        H_dict = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            H_dict[sigma] = H

        n = 1 << m
        h_hat = {}
        for bits in range(n):
            S_idx = tuple(i for i in range(m) if bits & (1 << i))
            coeff = 0
            for sigma_bits in range(n):
                sigma = tuple((1 if sigma_bits & (1 << i) else -1) for i in range(m))
                prod = 1
                for i in S_idx:
                    prod *= sigma[i]
                coeff += H_dict[sigma] * prod
            h_hat[S_idx] = coeff / n

        print(f"\n  p={p}, m={m}:")
        print(f"    {'deg':>4} {'|S|':>4} {'|h_hat|':>12} {'|h_hat|/p^(d/2)':>16} {'chi*h_hat':>12}")
        for S_idx, coeff in sorted(h_hat.items(), key=lambda x: (len(x[0]), x[0])):
            if abs(coeff) < 1e-10:
                continue
            d = len(S_idx)
            gaps = [i+1 for i in S_idx]
            prod_gap = 1
            for g in gaps:
                prod_gap *= g
            chi_val = legendre(prod_gap, p)

            normalized = abs(coeff) / (p ** (d/2)) if d > 0 else abs(coeff)
            chi_h = chi_val * coeff if d % 2 == 0 else 0
            print(f"    {d:>4} {d:>4} {abs(coeff):>12.4f} {normalized:>16.8f} {chi_h:>+12.4f}")


if __name__ == '__main__':
    main()
