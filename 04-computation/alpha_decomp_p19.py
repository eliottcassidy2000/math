#!/usr/bin/env python3
"""
alpha_decomp_p19.py -- Full alpha_j decomposition at p=19 for Paley and Interval

Computes:
  1. All directed odd cycles (vertices of Omega) for both tournaments
  2. alpha_1 = total cycle count, broken down by cycle length
  3. alpha_2 = number of vertex-disjoint cycle pairs
  4. Higher alpha_j as feasible
  5. H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ... (verification)

Uses Held-Karp DP on each k-vertex subset to count directed Hamiltonian cycles.

Author: kind-pasteur-2026-03-12-S58
"""

import time
import sys
from itertools import combinations


def build_adj(p, S):
    """Build adjacency matrix for circulant tournament on Z_p with connection set S."""
    A = [[0]*p for _ in range(p)]
    S_set = set(S)
    for i in range(p):
        for s in S_set:
            j = (i + s) % p
            A[i][j] = 1
    return A


def count_ham_cycles_on_subset(A, subset):
    """Count directed Hamiltonian cycles on the given vertex subset using Held-Karp.

    Returns the number of distinct directed Hamiltonian cycles (up to cyclic rotation).
    """
    k = len(subset)
    if k < 3:
        return 0

    # Map subset vertices to indices 0..k-1
    verts = list(subset)
    # Build sub-adjacency
    adj = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                adj[i][j] = A[verts[i]][verts[j]]

    # Held-Karp: fix source = vertex 0
    # dp[mask][v] = number of Hamiltonian paths from 0 to v using vertices in mask
    full_mask = (1 << k) - 1

    # Use dict for sparse DP at small k, flat array for larger k
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1  # start at vertex 0, mask = {0}

    for mask in range(1, 1 << k):
        if not (mask & 1):  # source vertex 0 must be in mask
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            cnt = dp[mask][v]
            for w in range(1, k):  # don't go back to source mid-path
                if mask & (1 << w):
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += cnt

    # Count cycles: paths ending at v where v -> 0
    total_cycles = 0
    for v in range(1, k):
        if dp[full_mask][v] and adj[v][0]:
            total_cycles += dp[full_mask][v]

    return total_cycles


def count_ham_cycles_k3(A, p, subset):
    """Optimized 3-cycle check."""
    a, b, c = subset
    # Check a->b->c->a
    if A[a][b] and A[b][c] and A[c][a]:
        return 1
    # Check a->c->b->a
    if A[a][c] and A[c][b] and A[b][a]:
        return 1
    return 0


def count_ham_cycles_k5(A, subset):
    """Count directed Hamiltonian cycles on 5-vertex subset using bitwise Held-Karp."""
    verts = list(subset)
    k = 5

    # Build sub-adjacency as bitmask
    adj = [0]*k
    for i in range(k):
        for j in range(k):
            if i != j and A[verts[i]][verts[j]]:
                adj[i] |= (1 << j)

    # Held-Karp from source=0
    dp = [0] * (32 * k)  # 2^5 * 5 flat array
    dp[1*k + 0] = 1

    for mask in range(1, 32):
        if not (mask & 1):
            continue
        base = mask * k
        for v in range(k):
            cnt = dp[base + v]
            if cnt == 0:
                continue
            # Expand to neighbors not in mask (skip vertex 0)
            expand = adj[v] & (~mask) & 0x1E  # bits 1-4 only
            while expand:
                w = expand & (-expand)  # lowest set bit
                w_idx = w.bit_length() - 1
                dp[(mask | w) * k + w_idx] += cnt
                expand ^= w

    # Count cycles: full_mask = 31
    total = 0
    base = 31 * k
    for v in range(1, k):
        if dp[base + v] and (adj[v] & 1):  # v -> 0
            total += dp[base + v]
    return total


def enumerate_cycles(A, p, k):
    """Enumerate all directed k-cycles in tournament, returning list of frozensets."""
    cycles = []
    verts = list(range(p))
    count = 0
    total_subsets = 0

    t0 = time.time()
    for subset in combinations(verts, k):
        total_subsets += 1
        if k == 3:
            n_cyc = count_ham_cycles_k3(A, p, subset)
        elif k == 5:
            n_cyc = count_ham_cycles_k5(A, subset)
        else:
            n_cyc = count_ham_cycles_on_subset(A, subset)

        if n_cyc > 0:
            fs = frozenset(subset)
            # Each Hamiltonian cycle on k vertices gives one vertex of Omega
            # If multiple cycles exist on the same vertex set, they are DIFFERENT
            # vertices of Omega (same vertex set, different arc pattern)
            for _ in range(n_cyc):
                cycles.append(fs)
            count += n_cyc

        if total_subsets % 10000 == 0:
            elapsed = time.time() - t0
            print(f"    k={k}: {total_subsets} subsets, {count} cycles found, {elapsed:.1f}s", flush=True)

    elapsed = time.time() - t0
    print(f"  k={k}: DONE. {total_subsets} subsets, {count} directed {k}-cycles, {elapsed:.1f}s", flush=True)
    return cycles


def compute_alpha_decomposition(cycles_list):
    """Given list of cycles (as frozensets), compute alpha_0, alpha_1, alpha_2, ...

    alpha_j = number of collections of j mutually vertex-disjoint cycles.
    """
    n = len(cycles_list)
    alpha_1 = n

    if n == 0:
        return [1, 0]

    # Compute alpha_2: count pairs of vertex-disjoint cycles
    print(f"  Computing alpha_2 from {n} cycles...", flush=True)
    t0 = time.time()
    alpha_2 = 0
    for i in range(n):
        for j in range(i+1, n):
            if not (cycles_list[i] & cycles_list[j]):  # disjoint
                alpha_2 += 1
    elapsed = time.time() - t0
    print(f"  alpha_2 = {alpha_2} ({elapsed:.1f}s)", flush=True)

    # For alpha_3: need triples of mutually disjoint cycles
    # Build adjacency list for "compatible" (disjoint) pairs
    if n <= 5000:
        print(f"  Computing alpha_3...", flush=True)
        t0 = time.time()
        compatible = [[] for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if not (cycles_list[i] & cycles_list[j]):
                    compatible[i].append(j)
                    compatible[j].append(i)

        alpha_3 = 0
        for i in range(n):
            # Find pairs (j, k) with j < k, both > i, disjoint from i and from each other
            comp_i = [j for j in compatible[i] if j > i]
            for ji, j in enumerate(comp_i):
                for k in comp_i[ji+1:]:
                    if not (cycles_list[j] & cycles_list[k]):
                        alpha_3 += 1
        elapsed = time.time() - t0
        print(f"  alpha_3 = {alpha_3} ({elapsed:.1f}s)", flush=True)

        # alpha_4 if feasible
        if n <= 2000 and alpha_3 > 0:
            print(f"  Computing alpha_4...", flush=True)
            t0 = time.time()
            alpha_4 = 0
            for i in range(n):
                comp_i = sorted([j for j in compatible[i] if j > i])
                for ji, j in enumerate(comp_i):
                    comp_ij = [k for k in comp_i[ji+1:] if not (cycles_list[j] & cycles_list[k])]
                    for ki, k in enumerate(comp_ij):
                        for l in comp_ij[ki+1:]:
                            if not (cycles_list[k] & cycles_list[l]):
                                if not (cycles_list[i] & cycles_list[l]):
                                    # Wait, l is in comp_i so disjoint from i
                                    # l in comp_ij so disjoint from j
                                    # Need: disjoint from k (checked)
                                    alpha_4 += 1
            elapsed = time.time() - t0
            print(f"  alpha_4 = {alpha_4} ({elapsed:.1f}s)", flush=True)
            return [1, alpha_1, alpha_2, alpha_3, alpha_4]

        return [1, alpha_1, alpha_2, alpha_3]
    else:
        return [1, alpha_1, alpha_2]


def main():
    p = 19
    m = (p - 1) // 2  # = 9

    # Connection sets
    S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    S_int = list(range(1, m + 1))

    print("=" * 70)
    print(f"ALPHA DECOMPOSITION AT p={p}")
    print("=" * 70)
    print(f"QR set: {S_qr}")
    print(f"Interval set: {S_int}")

    A_P = build_adj(p, S_qr)
    A_I = build_adj(p, S_int)

    for name, A in [("Paley", A_P), ("Interval", A_I)]:
        print(f"\n{'='*70}")
        print(f"{name} tournament on Z_{p}")
        print(f"{'='*70}")

        all_cycles = []
        cycle_counts = {}

        for k in range(3, p + 1, 2):
            n_subsets = 1
            for i in range(k):
                n_subsets = n_subsets * (p - i) // (i + 1)
            print(f"\n  Processing k={k} ({n_subsets} subsets)...", flush=True)

            t0 = time.time()
            cycles_k = enumerate_cycles(A, p, k)
            cycle_counts[k] = len(cycles_k)
            all_cycles.extend(cycles_k)

            # Check if remaining k values will take too long
            elapsed = time.time() - t0
            if elapsed > 300 and k < p - 2:
                print(f"  WARNING: k={k} took {elapsed:.0f}s. Estimating remaining time...", flush=True)
                # Rough estimate for next k
                next_k = k + 2
                next_subsets = 1
                for i in range(next_k):
                    next_subsets = next_subsets * (p - i) // (i + 1)
                ratio = next_subsets / n_subsets * (2**(next_k - k))
                print(f"  Next k={next_k}: ~{n_subsets}->{next_subsets} subsets, est. {elapsed*ratio:.0f}s", flush=True)

        print(f"\n  CYCLE COUNTS BY LENGTH:")
        total_alpha1 = 0
        for k in sorted(cycle_counts.keys()):
            c = cycle_counts[k]
            total_alpha1 += c
            print(f"    c_{k} = {c}")
        print(f"    alpha_1 = {total_alpha1}")

        print(f"\n  ALPHA DECOMPOSITION:")
        alphas = compute_alpha_decomposition(all_cycles)
        for j, a in enumerate(alphas):
            print(f"    alpha_{j} = {a}")

        # Compute H from alphas
        H_computed = sum(a * (2**j) for j, a in enumerate(alphas))
        higher = H_computed - 1 - 2 * alphas[1]
        print(f"\n  H (from alphas, partial) = {H_computed}")
        print(f"  Higher-order contribution = {higher}")

        # Known H values for comparison
        known_H = {"Paley": 1172695746915, "Interval": 1184212824763}
        if name in known_H:
            H_known = known_H[name]
            remaining = H_known - H_computed
            print(f"  H (known) = {H_known}")
            print(f"  Remaining (from alpha_{len(alphas)}+) = {remaining}")

    print(f"\n{'='*70}")
    print("COMPARISON")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
