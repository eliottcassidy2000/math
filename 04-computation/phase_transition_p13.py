#!/usr/bin/env python3
"""
phase_transition_p13.py -- Verify the phase transition at p=13

p=13 is predicted to be the crossover prime where Interval first beats Paley
(g(13) = 2*sqrt(13)/pi = 2.30 ≈ g_c).

Compute:
1. Full alpha_j decomposition for BOTH Paley and Interval at p=13
2. H via OCF = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
3. Compare with Held-Karp H for verification
4. Complete cycle count c_k for all odd k
5. Additive energy E(S)

p=13: m=6, so 2^6=64 circulant tournaments exist.
Can also do full orientation cube (much faster than p=19).

Author: kind-pasteur-2026-03-12-S58
"""

import time
from itertools import combinations


def build_adj(p, S):
    A = [[0]*p for _ in range(p)]
    S_set = set(S)
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def count_ham_cycles_k3(A, p, subset):
    a, b, c = subset
    if A[a][b] and A[b][c] and A[c][a]:
        return 1
    if A[a][c] and A[c][b] and A[b][a]:
        return 1
    return 0


def count_ham_cycles_on_subset(A, subset):
    k = len(subset)
    if k < 3:
        return 0
    verts = list(subset)
    adj = [0]*k
    for i in range(k):
        for j in range(k):
            if i != j and A[verts[i]][verts[j]]:
                adj[i] |= (1 << j)
    full_mask = (1 << k) - 1
    dp = [0] * ((1 << k) * k)
    dp[1*k + 0] = 1
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


def held_karp_H(p, S_set):
    """Compute H = p * h_0 via Held-Karp from vertex 0."""
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
    total = sum(dp[full_mask * p + v] for v in range(p))
    return total  # this is h_0; H = p * h_0


def enumerate_all_cycles(A, p):
    """Enumerate all directed odd cycles, returning list of frozensets with multiplicities."""
    cycles = []
    for k in range(3, p+1, 2):
        t0 = time.time()
        count = 0
        for subset in combinations(range(p), k):
            if k == 3:
                n_cyc = count_ham_cycles_k3(A, p, subset)
            else:
                n_cyc = count_ham_cycles_on_subset(A, subset)
            if n_cyc > 0:
                fs = frozenset(subset)
                for _ in range(n_cyc):
                    cycles.append(fs)
                count += n_cyc
        elapsed = time.time() - t0
        print(f"    k={k}: {count} directed cycles ({elapsed:.1f}s)")
    return cycles


def compute_alpha(cycles_list):
    """Compute alpha_j for j=0,1,2,3,..."""
    n = len(cycles_list)
    alphas = [1, n]  # alpha_0=1, alpha_1=n

    if n == 0:
        return alphas

    # alpha_2
    alpha_2 = 0
    for i in range(n):
        for j in range(i+1, n):
            if not (cycles_list[i] & cycles_list[j]):
                alpha_2 += 1
    alphas.append(alpha_2)

    # alpha_3
    compatible = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if not (cycles_list[i] & cycles_list[j]):
                compatible[i].append(j)
                compatible[j].append(i)

    alpha_3 = 0
    for i in range(n):
        for j in compatible[i]:
            if j <= i:
                continue
            for k_idx in compatible[j]:
                if k_idx <= j:
                    continue
                if not (cycles_list[i] & cycles_list[k_idx]):
                    alpha_3 += 1
    alphas.append(alpha_3)

    # alpha_4
    alpha_4 = 0
    for i in range(n):
        comp_i = set(j for j in compatible[i] if j > i)
        for j in sorted(comp_i):
            comp_ij = set(k for k in compatible[j] if k > j and k in comp_i
                          if not (cycles_list[i] & cycles_list[k]))
            # Wait, need to also check i vs k already done via comp_i
            comp_ij_real = []
            for k_idx in compatible[j]:
                if k_idx <= j:
                    continue
                if not (cycles_list[i] & cycles_list[k_idx]):
                    comp_ij_real.append(k_idx)
            for k_idx in comp_ij_real:
                for l_idx in compatible[k_idx]:
                    if l_idx <= k_idx:
                        continue
                    if not (cycles_list[i] & cycles_list[l_idx]):
                        if not (cycles_list[j] & cycles_list[l_idx]):
                            alpha_4 += 1
    alphas.append(alpha_4)

    return alphas


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
    print(f"PHASE TRANSITION ANALYSIS AT p={p}")
    print("=" * 70)
    print(f"m = {m}")
    print(f"QR = {S_qr}")
    print(f"Int = {S_int}")
    print(f"Number of circulant tournaments: 2^{m} = {2**m}")

    # Additive energy
    E_qr = add_energy(S_qr, p)
    E_int = add_energy(S_int, p)
    print(f"\nAdditive energy: E(QR)={E_qr}, E(Int)={E_int}, ratio={E_int/E_qr:.4f}")

    # Full alpha decomposition for Paley and Interval
    for name, S in [("Paley", S_qr), ("Interval", S_int)]:
        print(f"\n{'='*60}")
        print(f"{name} tournament: S = {S}")
        print(f"{'='*60}")

        A = build_adj(p, S)

        # 1. Enumerate all cycles
        print("\n  Enumerating all directed odd cycles...")
        cycles = enumerate_all_cycles(A, p)
        print(f"  Total cycles (alpha_1): {len(cycles)}")

        # 2. Compute H via Held-Karp for verification
        print("\n  Computing H via Held-Karp...", flush=True)
        t0 = time.time()
        h0 = held_karp_H(p, set(S))
        H = p * h0
        elapsed = time.time() - t0
        print(f"  H = {H:,d} ({elapsed:.1f}s)")

        # 3. Alpha decomposition
        print("\n  Computing alpha decomposition...")
        alphas = compute_alpha(cycles)
        print(f"  alphas = {alphas}")

        # 4. Verify OCF
        H_ocf = sum(2**j * alphas[j] for j in range(len(alphas)))
        print(f"  H(OCF) = {H_ocf:,d}")
        print(f"  H(HK)  = {H:,d}")
        if H_ocf == H:
            print(f"  ✓ OCF VERIFIED")
        else:
            print(f"  ✗ OCF MISMATCH (diff = {H_ocf - H})")
            # Check if we missed higher alphas
            print(f"  Missing = {H - H_ocf}")

    # Full orientation cube at p=13
    print(f"\n{'='*70}")
    print(f"FULL ORIENTATION CUBE: ALL {2**m} CIRCULANT TOURNAMENTS AT p={p}")
    print(f"{'='*70}")

    pairs = [(k+1, p-(k+1)) for k in range(m)]
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
        S_set = set(S)
        h0 = held_karp_H(p, S_set)
        H = p * h0
        results[bits] = (S, H)
        if (bits + 1) % 16 == 0:
            elapsed = time.time() - t0
            rate = (bits + 1) / elapsed
            eta = (2**m - bits - 1) / rate
            print(f"  {bits+1}/{2**m} done, {elapsed:.1f}s elapsed, ETA {eta:.0f}s", flush=True)

    total_time = time.time() - t0
    print(f"\n  Total computation time: {total_time:.1f}s")

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

    # Identify maximizer
    max_S = sorted_results[0][1][0]
    print(f"\n  HYP-480 at p={p}: Interval maximizes H?")
    print(f"  Maximizer S = {max_S}")
    print(f"  Is interval: {max_S == S_int}")

    H_int = next(H for S, H in results.values() if S == S_int)
    H_qr = next(H for S, H in results.values() if S == S_qr)

    n_beat = sum(1 for v in H_values if v > H_int)
    print(f"\n  H(Interval) = {H_int:,d}")
    print(f"  H(Paley)    = {H_qr:,d}")
    print(f"  H(Int) - H(Paley) = {H_int - H_qr:,d}")
    print(f"  # orientations beating Interval: {n_beat}")

    # Phase transition check
    g = 2 * (p ** 0.5) / 3.14159265
    print(f"\n  Dimensionless coupling g = 2*sqrt(p)/pi = {g:.4f}")
    if H_int > H_qr:
        print(f"  INTERVAL WINS at p={p} (g > g_c)")
    elif H_int < H_qr:
        print(f"  PALEY WINS at p={p} (g < g_c)")
    else:
        print(f"  TIE at p={p}")


if __name__ == '__main__':
    main()
