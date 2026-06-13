#!/usr/bin/env python3
"""
j_eigenvalue_transition.py -- Track the Paley J-eigenvalue across primes

The interaction matrix J[i,j] = Walsh degree-2 coefficient for {i,j}
encodes the 2-body Ising interaction on the orientation cube.

DISCOVERY (HYP-512): The Paley eigenvalue of J flips sign!
- p=7: lambda_P = +7.0 (POSITIVE, maximum eigenvalue)
- p=19: lambda_P = -544M (NEGATIVE, 3rd smallest)

Question: at what p does the flip occur? Need J eigenvalues at p=11.

For circulant tournaments, J can be computed from the full orientation cube.
p=11: 2^5 = 32 evaluations (fast!)

Author: kind-pasteur-2026-03-12-S58
"""

import time
import numpy as np


def held_karp_H(p, S_set):
    adj = [0] * p
    for i in range(p):
        for s in S_set:
            adj[i] |= (1 << ((i + s) % p))
    full_mask = (1 << p) - 1
    dp = [0] * ((1 << p) * p)
    dp[1 * p + 0] = 1
    for mask in range(1, 1 << p):
        if not (mask & 1):
            continue
        base = mask * p
        for v in range(p):
            cnt = dp[base + v]
            if cnt == 0:
                continue
            expand = adj[v] & (~mask) & full_mask
            while expand:
                w = expand & (-expand)
                w_idx = w.bit_length() - 1
                dp[(mask | w) * p + w_idx] += cnt
                expand ^= w
    return sum(dp[full_mask * p + v] for v in range(p))


def compute_walsh_J(p):
    """Compute the full orientation cube and extract J matrix."""
    m = (p - 1) // 2
    pairs = [(k+1, p-(k+1)) for k in range(m)]

    # Compute H for all 2^m orientations
    results = {}
    for bits in range(2**m):
        S = []
        sigma = []
        for k in range(m):
            if bits & (1 << k):
                S.append(pairs[k][1])
                sigma.append(1)
            else:
                S.append(pairs[k][0])
                sigma.append(0)
        S.sort()
        h0 = held_karp_H(p, set(S))
        H = p * h0
        results[bits] = (tuple(sigma), S, H)

    # Walsh-Fourier decomposition: compute all degree-2 coefficients
    walsh = {}
    for subset_bits in range(2**m):
        deg = bin(subset_bits).count('1')
        if deg > 2:
            continue
        coeff = 0
        for bits, (sigma, S, H) in results.items():
            parity = 0
            for k in range(m):
                if subset_bits & (1 << k):
                    parity ^= sigma[k]
            sign = 1 - 2 * parity
            coeff += H * sign
        walsh[subset_bits] = coeff / (2**m)

    # Build J matrix
    J = np.zeros((m, m))
    for i in range(m):
        for j in range(i+1, m):
            bits_ij = (1 << i) | (1 << j)
            J[i][j] = walsh[bits_ij]
            J[j][i] = walsh[bits_ij]

    H_0 = walsh[0]
    return J, H_0, results


def main():
    print("=" * 70)
    print("J-EIGENVALUE TRANSITION: TRACKING THE PALEY EIGENVALUE SIGN FLIP")
    print("=" * 70)

    for p in [5, 7, 11, 13]:
        m = (p - 1) // 2
        paley_exists = (p % 4 == 3)

        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m+1))
        pairs = [(k+1, p-(k+1)) for k in range(m)]

        print(f"\n{'='*60}")
        print(f"p={p}, m={m}, {2**m} orientations, Paley {'EXISTS' if paley_exists else 'N/A'}")
        print(f"{'='*60}")

        t0 = time.time()
        J, H_0, results = compute_walsh_J(p)
        elapsed = time.time() - t0
        print(f"  Computation time: {elapsed:.1f}s")
        print(f"  H_0 (mean) = {H_0:.2f}")

        # J eigenvalues
        eigvals = np.linalg.eigvalsh(J)
        print(f"  J eigenvalues: {sorted(eigvals)}")
        print(f"  Max eigenvalue: {max(eigvals):.2f}")
        print(f"  Min eigenvalue: {min(eigvals):.2f}")

        # Paley quadratic form
        if paley_exists:
            S_qr_set = set(S_qr)
            sigma_qr = tuple(0 if pairs[k][0] in S_qr_set else 1 for k in range(m))
            s_qr = np.array([1 - 2*s for s in sigma_qr], dtype=float)
            Q_qr = float(s_qr @ J @ s_qr)
            lambda_P = Q_qr / m
            print(f"\n  Paley sigma = {sigma_qr}")
            print(f"  Q(Paley) = sigma^T J sigma = {Q_qr:.2f}")
            print(f"  lambda_P = Q/m = {lambda_P:.2f}")

            # Check if it's an eigenvector
            Js = J @ s_qr
            if np.linalg.norm(Js) > 0:
                ratio = Js / s_qr
                if np.std(ratio[s_qr != 0]) < abs(np.mean(ratio[s_qr != 0])) * 0.001:
                    print(f"  Paley IS an eigenvector of J with eigenvalue {np.mean(ratio[s_qr != 0]):.2f}")
                else:
                    print(f"  Paley is NOT an eigenvector of J")
                    print(f"    Js/s ratios: {ratio}")

        # Interval quadratic form
        sigma_int = tuple(0 for _ in range(m))
        s_int = np.array([1 - 2*s for s in sigma_int], dtype=float)  # all +1
        Q_int = float(s_int @ J @ s_int)
        print(f"\n  Interval sigma = {sigma_int}")
        print(f"  Q(Interval) = {Q_int:.2f}")
        print(f"  sum(J) = {np.sum(J):.2f}")

        # SDP bound
        max_eig = max(eigvals)
        sdp_bound = H_0 + m * max_eig
        H_max = max(r[2] for r in results.values())
        sdp_gap = sdp_bound - H_max
        print(f"\n  SDP bound = H_0 + m*max_eig = {sdp_bound:.2f}")
        print(f"  Actual max H = {H_max}")
        print(f"  SDP gap = {sdp_gap:.2f} ({'POSITIVE' if sdp_gap > 0 else 'NEGATIVE'})")

        if paley_exists:
            print(f"\n  KEY: lambda_P = {lambda_P:.2f}, sign = {'POSITIVE' if lambda_P > 0 else 'NEGATIVE'}")

    # Now do p=19 using precomputed data
    print(f"\n{'='*60}")
    print(f"p=19, m=9, 512 orientations (from precomputed data)")
    print(f"{'='*60}")
    print(f"  J eigenvalues: [-604M(x2), -544M(x1), +163M(x2), +253M(x2), +460M(x2)]")
    print(f"  Paley eigenvalue = -544M (NEGATIVE)")
    print(f"  Q(Paley) = -4,900,099,958")
    print(f"  Q(Interval) = +1,117,343,600")

    # Summary table
    print(f"\n{'='*70}")
    print("PALEY EIGENVALUE TRANSITION TABLE")
    print("=" * 70)
    print(f"  {'p':>3s}  {'m':>3s}  {'lambda_P':>15s}  {'sign':>8s}  {'max_eig':>15s}  {'SDP_gap':>12s}  {'H winner':>10s}")


if __name__ == '__main__':
    main()
