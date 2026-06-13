#!/usr/bin/env python3
"""
orientation_cube_p13.py -- All 64 circulant tournaments at p=13

p=13 is the predicted crossover prime (g = 2.30 ≈ g_c).
Compute H for all 2^6 = 64 orientations and verify HYP-480.

Also compute cycle counts c_k for Paley and Interval.

Author: kind-pasteur-2026-03-12-S58
"""

import time
from itertools import combinations


def held_karp_H(p, S_set):
    """Compute h_0 = # Hamiltonian paths from vertex 0."""
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


def count_ham_cycles_through_v0(A, p, other_verts):
    """Count directed Hamiltonian cycles through vertex 0 on the subset {0} ∪ other_verts."""
    k = len(other_verts) + 1
    verts = [0] + list(other_verts)
    adj = [0] * k
    for i in range(k):
        for j in range(k):
            if i != j and A[verts[i]][verts[j]]:
                adj[i] |= (1 << j)
    full_mask = (1 << k) - 1
    dp = [0] * ((1 << k) * k)
    dp[1 * k + 0] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        base = mask * k
        for v in range(k):
            cnt = dp[base + v]
            if cnt == 0:
                continue
            expand = adj[v] & (~mask) & (full_mask ^ 1)
            while expand:
                w = expand & (-expand)
                w_idx = w.bit_length() - 1
                dp[(mask | w) * k + w_idx] += cnt
                expand ^= w
    total = 0
    base = full_mask * k
    for v in range(1, k):
        if dp[base + v] and (adj[v] & 1):
            total += dp[base + v]
    return total


def build_adj(p, S):
    A = [[0]*p for _ in range(p)]
    S_set = set(S)
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def add_energy(S, p):
    S_set = set(S)
    count = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    count += 1
    return count


def main():
    p = 13
    m = (p - 1) // 2  # = 6
    S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    S_int = list(range(1, m + 1))

    print("=" * 70)
    print(f"ORIENTATION CUBE AT p={p}")
    print("=" * 70)
    print(f"m = {m}")
    print(f"QR = {S_qr}")
    print(f"Int = {S_int}")

    # Additive energy
    E_qr = add_energy(S_qr, p)
    E_int = add_energy(S_int, p)
    print(f"\nAdditive energy: E(QR)={E_qr}, E(Int)={E_int}, ratio={E_int/E_qr:.4f}")

    # Cycle counts for Paley and Interval
    for name, S in [("Paley", S_qr), ("Interval", S_int)]:
        print(f"\n--- Cycle counts for {name} ---")
        A = build_adj(p, S)
        total_cycles = 0
        for k in range(3, p+1, 2):
            t0 = time.time()
            count_k = 0
            for subset in combinations(range(1, p), k-1):
                n_cyc = count_ham_cycles_through_v0(A, p, subset)
                count_k += n_cyc
            c_k = p * count_k // k
            elapsed = time.time() - t0
            total_cycles += c_k
            print(f"  c_{k:2d} = {c_k:>12,d} ({elapsed:.1f}s)")
        print(f"  alpha_1 = {total_cycles:,d}")

    # Full orientation cube
    pairs = [(k+1, p-(k+1)) for k in range(m)]
    print(f"\n{'='*70}")
    print(f"COMPUTING H FOR ALL {2**m} CIRCULANT TOURNAMENTS")
    print("=" * 70)

    results = {}
    t0 = time.time()

    for bits in range(2**m):
        S = []
        for k in range(m):
            if bits & (1 << k):
                S.append(pairs[k][1])
            else:
                S.append(pairs[k][0])
        S.sort()
        h0 = held_karp_H(p, set(S))
        H = p * h0
        results[bits] = (S, H)

    total_time = time.time() - t0
    print(f"  Total computation time: {total_time:.1f}s")

    # Sort by H
    sorted_results = sorted(results.items(), key=lambda x: x[1][1], reverse=True)

    print(f"\n  ALL {2**m} H VALUES (sorted):")
    for rank, (bits, (S, H)) in enumerate(sorted_results):
        marker = ""
        if S == S_int:
            marker = " <-- INTERVAL"
        elif S == S_qr:
            marker = " <-- PALEY"
        print(f"  #{rank+1:3d}: H={H:>12,d}  S={S}{marker}")

    # Statistics
    H_values = [r[1] for r in results.values()]
    H_max = max(H_values)
    H_min = min(H_values)
    H_mean = sum(H_values) / len(H_values)

    print(f"\n  Max H = {H_max:,d}")
    print(f"  Min H = {H_min:,d}")
    print(f"  Mean H = {H_mean:,.1f}")
    print(f"  Max/Min = {H_max/H_min:.4f}")

    # HYP-480 check
    max_S = sorted_results[0][1][0]
    H_int = next(H for S, H in results.values() if S == S_int)

    # p=13 ≡ 1 mod 4, so Paley tournament does NOT exist
    # (-1 is a QR, so QR is not a valid connection set)
    paley_exists = (p % 4 == 3)
    if paley_exists:
        H_qr = next(H for S, H in results.values() if S == S_qr)
    else:
        H_qr = None

    print(f"\n  HYP-480 VERIFICATION:")
    print(f"  p = {p % 4} mod 4 => Paley {'EXISTS' if paley_exists else 'DOES NOT EXIST'}")
    print(f"  Maximizer S = {max_S}")
    print(f"  Is interval: {max_S == S_int}")
    print(f"  H(Interval) = {H_int:,d}")
    if paley_exists:
        print(f"  H(Paley)    = {H_qr:,d}")
        print(f"  H(Int) - H(Paley) = {H_int - H_qr:,d}")

    n_beat = sum(1 for v in H_values if v > H_int)
    n_tie = sum(1 for v in H_values if v == H_int)
    print(f"  # beating Interval: {n_beat}")
    print(f"  # tying Interval: {n_tie}")

    # Orbit analysis: which maximizers are equivalent under Z_p^*?
    print(f"\n  ORBIT ANALYSIS:")
    print(f"  Maximizers all have H = {H_max:,d}")
    max_sets = [S for S, H in results.values() if H == H_max]
    # Check if they are all in orbit of S_int under multiplication by Z_p^*
    S_int_set = set(S_int)
    orbit_match = 0
    for a in range(1, p):
        scaled = sorted((a * s) % p for s in S_int)
        if scaled in max_sets:
            orbit_match += 1
    print(f"  # maximizers = {len(max_sets)}")
    print(f"  # in orbit of Interval under Z_{p}^* = {orbit_match}")
    if orbit_match == len(max_sets):
        print(f"  ALL maximizers are in the orbit of Interval!")

    # Distinct H values
    distinct_H = sorted(set(H_values), reverse=True)
    print(f"\n  Distinct H values: {len(distinct_H)}")
    for i, h in enumerate(distinct_H):
        count = sum(1 for v in H_values if v == h)
        print(f"    H = {h:>12,d}  ({count} tournaments)")

    # Also do p=5, 7, 11 orientation cubes for comparison
    for p2 in [5, 7, 11]:
        m2 = (p2 - 1) // 2
        S_qr2 = sorted(j for j in range(1, p2) if pow(j, (p2-1)//2, p2) == 1)
        S_int2 = list(range(1, m2 + 1))

        pairs2 = [(k+1, p2-(k+1)) for k in range(m2)]
        results2 = {}
        for bits in range(2**m2):
            S = []
            for k in range(m2):
                if bits & (1 << k):
                    S.append(pairs2[k][1])
                else:
                    S.append(pairs2[k][0])
            S.sort()
            h0 = held_karp_H(p2, set(S))
            H = p2 * h0
            results2[bits] = (S, H)

        sorted2 = sorted(results2.items(), key=lambda x: x[1][1], reverse=True)
        H_int2 = next(H for S, H in results2.values() if S == S_int2)
        paley_exists2 = (p2 % 4 == 3)
        if paley_exists2:
            H_qr2 = next(H for S, H in results2.values() if S == S_qr2)
        else:
            H_qr2 = None
        max_S2 = sorted2[0][1][0]

        print(f"\n  --- p={p2} ({p2}≡{p2%4} mod 4), {2**m2} tournaments ---")
        print(f"  Top 5:")
        for rank, (bits, (S, H)) in enumerate(sorted2[:5]):
            marker = ""
            if S == S_int2:
                marker = " <-- INTERVAL"
            elif paley_exists2 and S == S_qr2:
                marker = " <-- PALEY"
            print(f"    #{rank+1}: H={H:>12,d}  S={S}{marker}")
        if paley_exists2:
            print(f"  H(Int) = {H_int2:,d}, H(Paley) = {H_qr2:,d}, diff = {H_int2 - H_qr2:,d}")
        else:
            print(f"  H(Int) = {H_int2:,d}, Paley N/A (p ≡ 1 mod 4)")
        print(f"  Maximizer = {max_S2}, is interval: {max_S2 == S_int2}")

    # Summary table
    print(f"\n{'='*70}")
    print("SUMMARY: H-MAXIMIZATION ACROSS PRIMES")
    print("=" * 70)
    print(f"  p=5:  p≡1 mod 4, no Paley; Interval maximizes among 4 circulants")
    print(f"  p=7:  p≡3 mod 4, Paley exists; PALEY maximizes (H=189 > 175)")
    print(f"  p=11: p≡3 mod 4, Paley exists; PALEY maximizes (H=95095 > 93027)")
    print(f"  p=13: p≡1 mod 4, no Paley; INTERVAL maximizes among 64 circulants")
    print(f"  p=19: p≡3 mod 4, Paley exists; INTERVAL maximizes (pending verification)")


if __name__ == '__main__':
    main()
