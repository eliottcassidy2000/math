#!/usr/bin/env python3
"""
walsh_magnitude_formula.py -- Finding the formula for Walsh coefficient magnitudes

OBSERVED:
  p=7, deg-2: |h_hat| = 3.5 (one magnitude)
  p=11, deg-2: |h_hat| in {272.25, 8.25}  (two magnitudes)
  p=11, deg-4: |h_hat| = 118.25 (one magnitude)
  p=13, deg-2: |h_hat| in {2632.5, 778.375} (two magnitudes)
  p=13, deg-4: |h_hat| in {5209.75, 3508.375, 219.375} (three magnitudes)

CONJECTURE: The number of distinct magnitudes at degree d equals the
number of orbits of C(m,d)-element subsets of {1,...,m} under the
multiplier group action.

More precisely, the magnitude depends on the "type" of the subset S,
where the type is determined by the multiset of gap differences
{g_i - g_j mod p : i,j in S}.

This script:
1. Computes the magnitudes and their algebraic expressions
2. Tests the orbit structure under the multiplier action
3. Seeks a closed-form formula for |h_hat|

Author: kind-pasteur-2026-03-12-S59c
"""

from math import comb
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


def walsh_decomposition(p):
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    H_dict = {}
    for bits in range(1 << m):
        sigma = []
        S = []
        for i, (a, b) in enumerate(pairs):
            if bits & (1 << i):
                S.append(a)
                sigma.append(+1)
            else:
                S.append(b)
                sigma.append(-1)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_dict[tuple(sigma)] = H

    n = 1 << m
    h_hat = {}
    for bits in range(n):
        S_idx = [i for i in range(m) if bits & (1 << i)]
        coeff = 0
        for sigma_bits in range(n):
            sigma = [(1 if sigma_bits & (1 << i) else -1) for i in range(m)]
            H = H_dict[tuple(sigma)]
            prod = 1
            for i in S_idx:
                prod *= sigma[i]
            coeff += H * prod
        h_hat[tuple(S_idx)] = coeff / n

    return h_hat, H_dict


def multiplier_orbit(S, m, p):
    """Compute the orbit of subset S of {0,...,m-1} under multiplier action.

    The multiplier a maps chord index i (gap g=i+1) to chord index j
    where j+1 = a*(i+1) mod p, with canonicalization: if j+1 > m, use p-(j+1).
    """
    orbit = set()
    S_tuple = tuple(sorted(S))
    orbit.add(S_tuple)

    for a in range(2, p):
        # Map each element of S
        new_S = []
        for i in S:
            g = i + 1
            ag = (a * g) % p
            if ag > m:
                j = p - ag - 1  # canonical index
            else:
                j = ag - 1
            new_S.append(j)
        new_S_tuple = tuple(sorted(new_S))
        orbit.add(new_S_tuple)

    return frozenset(orbit)


def main():
    print("=" * 70)
    print("WALSH COEFFICIENT MAGNITUDE FORMULA")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        h_hat, H_dict = walsh_decomposition(p)

        # Group by degree and magnitude
        for deg in range(0, m+1):
            items = [(S, c) for S, c in h_hat.items()
                     if len(S) == deg and abs(c) > 0.01]
            if not items:
                continue

            print(f"\n  Degree {deg}: {len(items)} non-zero coefficients")

            # Group by magnitude
            mag_groups = defaultdict(list)
            for S, c in items:
                mag = round(abs(c), 6)
                mag_groups[mag].append((S, c))

            for mag in sorted(mag_groups, reverse=True):
                group = mag_groups[mag]
                n_pos = sum(1 for _, c in group if c > 0)
                n_neg = sum(1 for _, c in group if c < 0)

                # Check if all are in the same multiplier orbit
                orbits = set()
                for S, c in group:
                    orb = multiplier_orbit(S, m, p)
                    orbits.add(orb)

                print(f"\n    |h| = {mag} ({n_pos}+ {n_neg}-)")
                print(f"    {len(orbits)} orbit(s)")

                # Product of gaps for each coefficient
                for S, c in group:
                    gaps = tuple(i + 1 for i in S)
                    prod_gaps = 1
                    for g in gaps:
                        prod_gaps = (prod_gaps * g) % p
                    chi = 1 if prod_gaps in QR else -1
                    print(f"      S={set(S)}, gaps={gaps}, "
                          f"prod={prod_gaps}, chi={chi:>+d}, "
                          f"h={c:>+.2f}")

        # ====== RATIONAL RECONSTRUCTION ======
        print(f"\n  --- RATIONAL RECONSTRUCTION ---")
        # Multiply by 2^m to get integers
        scale = 1 << m
        print(f"  Scale factor: 2^{m} = {scale}")

        for deg in range(0, m+1):
            items = [(S, c) for S, c in h_hat.items()
                     if len(S) == deg and abs(c) > 0.01]
            if not items:
                continue

            mags = sorted(set(round(abs(c) * scale) for _, c in items))
            if len(mags) > 0:
                print(f"\n  Degree {deg}: scaled magnitudes = {mags}")
                for mag_s in mags:
                    # Try to factor
                    factors = []
                    rem = mag_s
                    for f in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
                        while rem % f == 0:
                            factors.append(f)
                            rem //= f
                    if rem > 1:
                        factors.append(rem)
                    print(f"    {mag_s} = {'*'.join(str(f) for f in factors) if factors else '1'}")

                # Ratios between magnitudes
                if len(mags) > 1:
                    for i in range(len(mags)):
                        for j in range(i+1, len(mags)):
                            ratio = mags[j] / mags[i] if mags[i] > 0 else float('inf')
                            print(f"    ratio {mags[j]}/{mags[i]} = {ratio:.6f}")

        # ====== FORMULA HUNT ======
        print(f"\n  --- FORMULA HUNT ---")
        mean_H = h_hat[()]
        print(f"  h_hat[{{}}] = mean H = {mean_H:.2f}")
        print(f"  mean H * 2^m = {mean_H * scale:.0f}")

        # Try expressing magnitudes in terms of p, m, mean_H
        for deg in [2, 4]:
            items = [(S, c) for S, c in h_hat.items()
                     if len(S) == deg and abs(c) > 0.01]
            if not items:
                continue

            print(f"\n  Degree {deg}:")
            mags = sorted(set(abs(c) for _, c in items))
            for mag in mags:
                # Normalized
                ratio_to_mean = mag / mean_H if mean_H > 0 else 0
                ratio_to_p = mag * p
                print(f"    |h| = {mag:.4f}")
                print(f"    |h|/mean = {ratio_to_mean:.8f}")
                print(f"    |h|*p = {ratio_to_p:.4f}")
                print(f"    |h|*p/mean = {ratio_to_p/mean_H:.8f}" if mean_H > 0 else "")
                print(f"    |h|*p^2/mean = {mag*p*p/mean_H:.8f}" if mean_H > 0 else "")

    # ====== CROSS-PRIME COMPARISON ======
    print(f"\n{'='*70}")
    print("CROSS-PRIME COMPARISON OF MAGNITUDES")
    print("=" * 70)

    for deg in [2]:
        print(f"\n  Degree {deg}:")
        for p in [7, 11, 13]:
            m = (p - 1) // 2
            h_hat, _ = walsh_decomposition(p)
            items = [(S, c) for S, c in h_hat.items()
                     if len(S) == deg and abs(c) > 0.01]
            mags = sorted(set(abs(c) for _, c in items), reverse=True)
            mean_H = h_hat[()]

            print(f"\n    p={p}, m={m}:")
            for i, mag in enumerate(mags):
                scaled = mag * (1 << m)
                print(f"      mag_{i} = {mag:.4f}, scaled = {scaled:.0f}, "
                      f"|h|/mean = {mag/mean_H:.8f}")

            # Check: does |h|_large / |h|_small = some function of p?
            if len(mags) >= 2:
                ratio = mags[0] / mags[1]
                print(f"      ratio large/small = {ratio:.4f}")
                # Check against p, m, etc
                print(f"      p = {p}, m = {m}")
                print(f"      (m+1) = {m+1}, p-2 = {p-2}")
                print(f"      ratio vs m = {ratio/m:.4f}")
                print(f"      ratio vs (m-1) = {ratio/(m-1):.4f}")
                print(f"      ratio vs (p-2) = {ratio/(p-2):.4f}")


if __name__ == '__main__':
    main()
