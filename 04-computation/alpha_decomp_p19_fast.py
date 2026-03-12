#!/usr/bin/env python3
"""
alpha_decomp_p19_fast.py -- Cycle counts at p=19 using circulant symmetry

Key optimization: For circulant tournaments, c_k = p * c_k^{(0)} / k
where c_k^{(0)} = number of directed k-cycles through vertex 0.
This reduces to enumerating C(p-1, k-1) subsets instead of C(p, k).

We compute alpha_1 = total directed odd cycles = sum c_k.
Then higher = H - 1 - 2*alpha_1 gives the disjoint-packing contribution.

Author: kind-pasteur-2026-03-12-S58
"""

import time
from itertools import combinations


def build_adj(p, S):
    """Build adjacency matrix for circulant tournament on Z_p."""
    A = [[0]*p for _ in range(p)]
    S_set = set(S)
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_cycles_through_v0(A, p, other_verts):
    """Count directed Hamiltonian cycles through vertex 0 on {0} + other_verts.

    Uses Held-Karp DP with vertex 0 as fixed source.
    Returns number of distinct directed Hamiltonian cycles.
    """
    k = len(other_verts) + 1  # total vertices including 0
    verts = [0] + list(other_verts)

    # Build sub-adjacency indexed 0..k-1
    # Use flat list for speed
    adj = [0] * k  # bitmask adjacency
    for i in range(k):
        for j in range(k):
            if i != j and A[verts[i]][verts[j]]:
                adj[i] |= (1 << j)

    full_mask = (1 << k) - 1

    # DP: dp[mask * k + v] = # Hamiltonian paths from 0 to v using mask
    size = (1 << k) * k
    dp = [0] * size
    dp[1 * k + 0] = 1  # start: mask={0}, at vertex 0

    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        base = mask * k
        for v in range(k):
            cnt = dp[base + v]
            if cnt == 0:
                continue
            # Expand to unvisited neighbors (skip 0 = source during path)
            expand = adj[v] & (~mask) & (full_mask ^ 1)  # exclude bit 0
            while expand:
                w = expand & (-expand)
                w_idx = w.bit_length() - 1
                dp[(mask | w) * k + w_idx] += cnt
                expand ^= w

    # Count: paths using all vertices, ending at v where v->0
    total = 0
    base = full_mask * k
    for v in range(1, k):
        if dp[base + v] and (adj[v] & 1):
            total += dp[base + v]

    return total


def main():
    p = 19
    m = (p - 1) // 2  # = 9

    S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    S_int = list(range(1, m + 1))

    known_H = {
        "Paley": 1172695746915,
        "Interval": 1184212824763
    }

    print("=" * 70)
    print(f"CYCLE COUNTS AND ALPHA DECOMPOSITION AT p={p}")
    print("=" * 70)
    print(f"QR set:       {S_qr}")
    print(f"Interval set: {S_int}")

    other_verts = list(range(1, p))

    for name, S in [("Paley", S_qr), ("Interval", S_int)]:
        A = build_adj(p, S)
        print(f"\n{'='*70}")
        print(f"{name} tournament on Z_{p}")
        print(f"{'='*70}")

        cycle_counts = {}
        total_alpha1 = 0
        grand_t0 = time.time()

        for k in range(3, p + 1, 2):
            n_subsets = 1
            for i in range(k - 1):
                n_subsets = n_subsets * (p - 1 - i) // (i + 1)

            print(f"\n  k={k}: enumerating cycles through vertex 0 ({n_subsets} subsets)...", flush=True)
            t0 = time.time()

            # Count cycles through vertex 0
            c_k_through_0 = 0
            done = 0

            for subset in combinations(other_verts, k - 1):
                hc = count_ham_cycles_through_v0(A, p, subset)
                c_k_through_0 += hc
                done += 1
                if done % 5000 == 0:
                    elapsed = time.time() - t0
                    rate = done / elapsed if elapsed > 0 else float('inf')
                    eta = (n_subsets - done) / rate if rate > 0 else 0
                    print(f"    {done}/{n_subsets} subsets, {c_k_through_0} cycles, "
                          f"{elapsed:.1f}s elapsed, ETA {eta:.0f}s", flush=True)

            elapsed = time.time() - t0

            # By circulant symmetry: c_k = p * c_k^{(0)} / k
            c_k = p * c_k_through_0 // k
            cycle_counts[k] = c_k
            total_alpha1 += c_k

            # Verify integrality
            if p * c_k_through_0 % k != 0:
                print(f"  WARNING: p * c_k_through_0 / k = {p}*{c_k_through_0}/{k} is NOT integer!", flush=True)

            print(f"  k={k}: c_k^(0) = {c_k_through_0}, c_k = {c_k} ({elapsed:.1f}s)", flush=True)

        grand_elapsed = time.time() - grand_t0
        print(f"\n  Total enumeration time: {grand_elapsed:.1f}s")

        print(f"\n  CYCLE COUNTS:")
        for k in sorted(cycle_counts.keys()):
            c = cycle_counts[k]
            nC = 1
            for i in range(k):
                nC = nC * (p - i) // (i + 1)
            print(f"    c_{k:2d} = {c:>20,d}   (C({p},{k})={nC}, ratio={c/nC:.4f})")

        print(f"\n  alpha_1 = {total_alpha1:,d}")
        H = known_H[name]
        higher = H - 1 - 2 * total_alpha1
        print(f"  H (known) = {H:,d}")
        print(f"  2*alpha_1 = {2*total_alpha1:,d}")
        print(f"  higher (sum_{{j>=2}} 2^j alpha_j) = {higher:,d}")
        print(f"  higher / H = {higher/H:.6f} = {100*higher/H:.2f}%")

    # Comparison
    print(f"\n{'='*70}")
    print("ISING DECOMPOSITION COMPARISON")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
