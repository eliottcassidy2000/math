#!/usr/bin/env python3
"""
orientation_cube_p19.py -- Compute H for ALL 2^9 = 512 circulant tournaments on Z_19

For p=19, m=9. Each circulant tournament is determined by choosing for each
pair (k, 19-k) with k=1,...,9 which element is in the connection set S.
This gives 2^9 = 512 possible connection sets.

For each, we compute H via Held-Karp DP starting from vertex 0
(using circulant symmetry: H = 19 * h_0).

This allows:
1. Exhaustive verification of HYP-480 at p=19
2. Walsh-Fourier decomposition of H on the orientation cube
3. SDP integrality gap computation
4. Identification of all local maxima

Author: kind-pasteur-2026-03-12-S58
"""

import time


def held_karp_from_v0(p, S_set):
    """Compute number of Hamiltonian paths starting at vertex 0.

    Uses Held-Karp DP with flat array. Returns h_0 such that H = p * h_0.
    """
    # Build bitmask adjacency
    adj = [0] * p
    for i in range(p):
        for s in S_set:
            j = (i + s) % p
            adj[i] |= (1 << j)

    full_mask = (1 << p) - 1

    # Use flat array dp[mask * p + v]
    # 2^19 * 19 ≈ 10M entries. Use list of ints.
    dp = [0] * ((1 << p) * p)
    dp[1 * p + 0] = 1  # mask={vertex 0}, at vertex 0

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

    # Sum all paths at full mask
    total = sum(dp[full_mask * p + v] for v in range(p))
    return total


def main():
    p = 19
    m = (p - 1) // 2  # = 9

    # The orientation cube: sigma ∈ {0,1}^m
    # sigma[k] = 0 means k+1 ∈ S, sigma[k] = 1 means p-(k+1) ∈ S
    # (k ranges 0..8, corresponding to pairs (1,18), (2,17), ..., (9,10))

    pairs = [(k+1, p-(k+1)) for k in range(m)]
    print(f"p={p}, m={m}")
    print(f"Pairs: {pairs}")
    print(f"Total orientations: {2**m} = {2**m}")

    # Known reference points
    S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    S_int = list(range(1, m + 1))
    print(f"QR = {S_qr}")
    print(f"Int = {S_int}")

    # Encode Paley and Interval as sigma vectors
    S_qr_set = set(S_qr)
    S_int_set = set(S_int)

    sigma_qr = tuple(0 if pairs[k][0] in S_qr_set else 1 for k in range(m))
    sigma_int = tuple(0 if pairs[k][0] in S_int_set else 1 for k in range(m))
    print(f"sigma_QR = {sigma_qr} (bits = {sum(s * 2**k for k, s in enumerate(sigma_qr))})")
    print(f"sigma_Int = {sigma_int} (bits = {sum(s * 2**k for k, s in enumerate(sigma_int))})")

    print(f"\n{'='*70}")
    print(f"COMPUTING H FOR ALL {2**m} ORIENTATIONS")
    print(f"{'='*70}")

    results = {}
    t0 = time.time()

    for bits in range(2**m):
        # Construct S from bits
        S = []
        sigma = []
        for k in range(m):
            if bits & (1 << k):
                S.append(pairs[k][1])  # take p-(k+1)
                sigma.append(1)
            else:
                S.append(pairs[k][0])  # take k+1
                sigma.append(0)
        S.sort()
        S_set = set(S)

        h0 = held_karp_from_v0(p, S_set)
        H = p * h0

        results[bits] = (tuple(sigma), S, H)

        if (bits + 1) % 32 == 0:
            elapsed = time.time() - t0
            rate = (bits + 1) / elapsed
            eta = (2**m - bits - 1) / rate
            print(f"  {bits+1}/{2**m} done, {elapsed:.1f}s elapsed, ETA {eta:.0f}s", flush=True)

    total_time = time.time() - t0
    print(f"\n  Total computation time: {total_time:.1f}s")

    # Sort by H
    sorted_results = sorted(results.items(), key=lambda x: x[1][2], reverse=True)

    # Top 20
    print(f"\n{'='*70}")
    print(f"TOP 20 H VALUES")
    print(f"{'='*70}")
    for rank, (bits, (sigma, S, H)) in enumerate(sorted_results[:20]):
        marker = ""
        if S == S_int:
            marker = " <-- INTERVAL"
        elif S == S_qr:
            marker = " <-- PALEY"
        print(f"  #{rank+1}: H={H:>22,d}  S={S}{marker}")

    # Bottom 5
    print(f"\nBOTTOM 5:")
    for rank, (bits, (sigma, S, H)) in enumerate(sorted_results[-5:]):
        print(f"  #{2**m - 4 + rank}: H={H:>22,d}  S={S}")

    # Statistics
    H_values = [r[2] for r in results.values()]
    H_max = max(H_values)
    H_min = min(H_values)
    H_mean = sum(H_values) / len(H_values)

    print(f"\n{'='*70}")
    print(f"STATISTICS")
    print(f"{'='*70}")
    print(f"  Max H = {H_max:,d}")
    print(f"  Min H = {H_min:,d}")
    print(f"  Mean H = {H_mean:,.1f}")
    print(f"  Max/Min = {H_max/H_min:.4f}")

    # Check HYP-480
    max_S = sorted_results[0][1][1]
    print(f"\n  HYP-480 at p={p}: Interval maximizes H?")
    print(f"  Maximizer S = {max_S}")
    print(f"  Is interval: {max_S == S_int}")

    # Count how many orientations beat Interval
    H_int = None
    for bits, (sigma, S, H) in results.items():
        if S == S_int:
            H_int = H
            break

    n_beat = sum(1 for v in H_values if v > H_int)
    n_tie = sum(1 for v in H_values if v == H_int)
    print(f"  H(Interval) = {H_int:,d}")
    print(f"  # orientations beating Interval: {n_beat}")
    print(f"  # orientations tying Interval: {n_tie}")

    # Walsh-Fourier decomposition
    print(f"\n{'='*70}")
    print(f"WALSH-FOURIER DECOMPOSITION")
    print(f"{'='*70}")

    # Convert to +1/-1 convention: s_k = 1 - 2*sigma_k
    # H(sigma) = sum_S hat{H}(S) * prod_{k in S} s_k
    #
    # hat{H}(S) = (1/2^m) * sum_sigma H(sigma) * prod_{k in S} s_k

    # Compute all Walsh coefficients
    walsh = {}
    for subset_bits in range(2**m):
        coeff = 0
        for bits, (sigma, S, H) in results.items():
            # s_k = 1 - 2*sigma_k
            parity = 0
            for k in range(m):
                if subset_bits & (1 << k):
                    parity ^= sigma[k]
            sign = 1 - 2 * parity
            coeff += H * sign
        walsh[subset_bits] = coeff / (2**m)

    # Group by degree
    degree_energy = {}
    print(f"\n  Walsh coefficients by degree:")
    for subset_bits in range(2**m):
        deg = bin(subset_bits).count('1')
        w = walsh[subset_bits]
        if abs(w) > 0.5:  # significant
            if deg not in degree_energy:
                degree_energy[deg] = 0
            degree_energy[deg] += w**2
            if deg <= 2 or abs(w) > H_mean * 0.01:
                # Show small degree or large coefficients
                which = [k for k in range(m) if subset_bits & (1 << k)]
                print(f"    degree {deg}, vars {which}: hat_H = {w:.2f}")

    print(f"\n  Energy by degree:")
    total_energy = sum(degree_energy.values())
    for deg in sorted(degree_energy.keys()):
        pct = 100 * degree_energy[deg] / total_energy if total_energy > 0 else 0
        print(f"    degree {deg}: {degree_energy[deg]:.2e} ({pct:.2f}%)")

    # SDP integrality gap
    # Degree-2 SDP maximizes H_0 + sigma^T J sigma
    # where J[i,j] = hat{H}({i,j})
    # The maximum of sigma^T J sigma over {+1,-1}^m is ||J||_op
    print(f"\n{'='*70}")
    print(f"SDP ANALYSIS")
    print(f"{'='*70}")

    # Build J matrix
    J = [[0.0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            bits_ij = (1 << i) | (1 << j)
            J[i][j] = walsh[bits_ij]
            J[j][i] = walsh[bits_ij]

    H_0 = walsh[0]  # degree-0 coefficient (mean)

    # Compute J eigenvalues via power iteration or direct for small m
    # For m=9, use numpy if available
    try:
        import numpy as np
        J_np = np.array(J)
        eigvals = np.linalg.eigvalsh(J_np)
        max_eig = max(eigvals)
        min_eig = min(eigvals)
        print(f"  H_0 (mean) = {H_0:.2f}")
        print(f"  J eigenvalues: {sorted(eigvals)}")
        print(f"  Max eigenvalue = {max_eig:.2f}")
        print(f"  SDP bound (degree-2) = H_0 + m * max_eig = {H_0 + m * max_eig:.2f}")
        print(f"  Actual max H = {H_max}")
        print(f"  SDP gap = {(H_0 + m * max_eig) - H_max:.2f}")

        # Check if Paley is the eigenvector
        idx = list(eigvals).index(max_eig)
        # Paley sigma in +1/-1
        s_paley = [1 - 2*sigma_qr[k] for k in range(m)]
        quad_paley = sum(J[i][j] * s_paley[i] * s_paley[j]
                        for i in range(m) for j in range(m))
        print(f"  Q(Paley) = sigma_P^T J sigma_P = {quad_paley:.2f}")

        s_int = [1 - 2*sigma_int[k] for k in range(m)]
        quad_int = sum(J[i][j] * s_int[i] * s_int[j]
                       for i in range(m) for j in range(m))
        print(f"  Q(Interval) = sigma_I^T J sigma_I = {quad_int:.2f}")

    except ImportError:
        print("  numpy not available, skipping eigenvalue computation")


if __name__ == '__main__':
    main()
