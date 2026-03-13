#!/usr/bin/env python3
"""
eigenvalue_walsh_layers.py -- Layer-by-layer eigenvalue contribution to Walsh coefficients

KEY DISCOVERY: The degree-2 Walsh of tr(A^k)/k decomposes as:
  h_hat = sum_{j=2,4,6,...,k} C(k,j)*(-1/2)^{k-j} * [D^j]_{a,b} / k

where D = sum_i sigma_i * i*sin(2*pi*g_i*t/p).

The [D^j]_{a,b} (degree-2 Walsh of D^j) involves j-th order Gauss sums.
Each even j activates a new resonance level:
  j=2: only q=1 (orthogonal, gives 0 for a!=b)
  j=4: activates q=3 (via sin^3 -> sin(3x) -> tripled frequency)
  j=6: activates q=5 (via sin^5 -> sin(5x) -> quintupled frequency)
  j=2r: activates q=2r-1

This script verifies and extracts exact coefficients for all j-layers.

THEOREM CANDIDATES:
  j=4, q=3: coefficient = -eps * C(k,4) * (1/2)^{k-4} * p / k
  j=6, q=3: coefficient = eps * m^2 * p / 8 (at k=6)
  j=6, q=5: coefficient = eps5 * 4 * p / 8 = eps5 * p / 2 (at k=6)

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def resonance_level(a, b, p):
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def resonance_sign(ga, gb, p, q):
    if (q * ga - gb) % p == 0:
        return -1
    elif (q * ga + gb) % p == 0:
        return +1
    elif (ga - q * gb) % p == 0:
        return -1
    elif (ga + q * gb) % p == 0:
        return +1
    return 0


def compute_D_j_walsh2(j, a, b, m, p):
    """Compute sum_t [D^j]_{a,b} numerically.

    D = sum_c sigma_c * d_c where d_c = i*sin(2*pi*g_c*t/p).
    [D^j]_{a,b} = degree-2 Walsh component (sigma_a * sigma_b) of D^j.

    We compute this for each t=1,...,p-1 and sum.
    """
    ga = a + 1
    gb = b + 1
    total = 0

    for t in range(1, p):
        # Build d_c values for all chords
        d = []
        for c in range(m):
            gc = c + 1
            alpha = 2 * math.pi * gc * t / p
            d.append(1j * math.sin(alpha))

        da = d[a]
        db = d[b]
        d_other = [d[c] for c in range(m) if c != a and c != b]

        # Compute [D^j]_{a,b} at this t using the partition formula
        # We need all partitions (n_a, n_b, {n_c}) with:
        #   n_a odd, n_b odd, all n_c even, sum = j
        # The coefficient for partition P = j! / prod(n_i!)
        # and contribution = coefficient * d_a^{n_a} * d_b^{n_b} * prod d_c^{n_c}

        # For efficiency, compute using the identity:
        # [D^j]_{a,b} = (1/4) * [(da+db+S)^j + (da+db-S)^j - (da-db+S)^j - (da-db-S)^j]
        # where S = sum d_c for c not a,b... NO this doesn't work for Walsh.

        # Actually, the correct formula uses:
        # [D^j]_{a,b} = (1/4) * sum over (eps_a, eps_b) in {+1,-1}^2:
        #   eps_a * eps_b * (eps_a*da + eps_b*db + sum d_other)^j
        # This extracts the coefficient of sigma_a * sigma_b.

        # Must average over all 2^{m-2} settings of sigma_c for c != a,b
        n_other = len(d_other)
        val = 0
        for ea in [1, -1]:
            for eb in [1, -1]:
                # Average over all sigma' for other chords
                avg = 0
                for bits_c in range(1 << n_other):
                    S_other = sum(d_other[c] * (1 if bits_c & (1 << c) else -1)
                                 for c in range(n_other))
                    x = ea * da + eb * db + S_other
                    avg += x**j
                avg /= (1 << n_other)
                val += ea * eb * avg
        val /= 4

        total += val

    return total.real  # Should be real after sum over t


def main():
    print("=" * 70)
    print("EIGENVALUE WALSH LAYERS: j=2,4,6,... decomposition")
    print("=" * 70)

    for p in [7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        if m > 11:
            continue
        print(f"\n{'='*60}")
        print(f"p={p}, m={m}")
        print("=" * 60)

        chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]

        for k in [4, 5, 6, 7, 8]:
            if k > p:
                continue
            print(f"\n  k={k}:")
            max_j = min(k, 2 * m)  # j can't exceed 2m (each chord used at most... well it can)
            # Actually j goes up to k

            # For each pair, compute j-layer contributions
            print(f"  {'pair':>8} {'q':>3} {'eps':>4}", end="")
            for j in range(2, k + 1, 2):
                print(f" {'j='+str(j):>10}", end="")
            print(f" {'total':>10} {'coeff_j4':>10} {'coeff_j6':>10}")

            for a, b in chord_pairs:
                ga, gb = a + 1, b + 1
                q = resonance_level(ga, gb, p)
                eps_q = resonance_sign(ga, gb, p, q)

                # Compute each j-layer
                layers = {}
                total_pred = 0
                for j in range(2, k + 1, 2):
                    Dj = compute_D_j_walsh2(j, a, b, m, p)
                    coeff = math.comb(k, j) * (-0.5) ** (k - j)
                    layer = coeff * Dj / k  # divided by k for tr/k
                    layers[j] = layer
                    total_pred += layer

                # Only print if significant
                if abs(total_pred) > 0.01:
                    print(f"  ({a},{b}) q={q:>3} eps={eps_q:>+2d}", end="")
                    for j in range(2, k + 1, 2):
                        print(f" {layers[j]:>10.4f}", end="")

                    # Coefficient: layer / (eps * p / 8)
                    c_j4 = layers.get(4, 0) / (eps_q * p / 8) if eps_q != 0 and abs(layers.get(4, 0)) > 0.01 else 0
                    c_j6 = layers.get(6, 0) / (eps_q * p / 8) if eps_q != 0 and abs(layers.get(6, 0)) > 0.01 else 0

                    # For q=5, use eps5
                    if q == 5 and abs(layers.get(6, 0)) > 0.01:
                        eps5 = resonance_sign(ga, gb, p, 5)
                        c_j6 = layers.get(6, 0) / (eps5 * p / 8) if eps5 != 0 else 0

                    print(f" {total_pred:>10.4f} {c_j4:>10.4f} {c_j6:>10.4f}")

        # Summary: universal j-layer coefficient patterns
        print(f"\n  SUMMARY of j-layer coefficients (normalized by eps*p/8):")
        for k in [4, 5, 6, 7, 8]:
            if k > p:
                continue
            print(f"\n    k={k}:")
            for q_target in [3, 5, 7]:
                coeffs_by_j = defaultdict(list)
                for a, b in chord_pairs:
                    ga, gb = a + 1, b + 1
                    q = resonance_level(ga, gb, p)
                    if q != q_target:
                        continue
                    eps_q = resonance_sign(ga, gb, p, q)

                    for j in range(2, k + 1, 2):
                        Dj = compute_D_j_walsh2(j, a, b, m, p)
                        coeff = math.comb(k, j) * (-0.5) ** (k - j)
                        layer = coeff * Dj / k
                        if abs(layer) > 0.001:
                            normalized = layer / (eps_q * p / 8)
                            coeffs_by_j[j].append(round(normalized, 4))

                for j in sorted(coeffs_by_j.keys()):
                    vals = coeffs_by_j[j]
                    if len(set(round(v, 2) for v in vals)) == 1:
                        print(f"      q={q_target}, j={j}: coeff = {vals[0]:.4f} (UNIVERSAL)")
                    else:
                        print(f"      q={q_target}, j={j}: coeffs = {sorted(set(round(v, 2) for v in vals))} (VARIES)")

    print("\nDONE.")


if __name__ == '__main__':
    main()
