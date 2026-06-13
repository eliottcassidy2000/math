#!/usr/bin/env python3
"""
Correct full Omega(T) construction and component factorization verification.

Omega(T) vertices = individual directed odd cycles (up to cyclic rotation).
A directed k-cycle on vertex set S is a cyclic ordering (v_1, v_2, ..., v_k)
such that v_i -> v_{i+1} for all i (mod k). Two orderings are the same cycle
iff one is a cyclic rotation of the other.

Two cycles are adjacent in Omega iff they share a tournament vertex.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations, permutations
from collections import defaultdict


def held_karp(n, adj_bits):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            c = dp[S][v]
            if c == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj_bits[v] & (1 << u):
                    dp[S | (1 << u)][u] += c
    return sum(dp[(1 << n) - 1])


def count_directed_cycles_on_set(verts, A):
    """Count directed Hamiltonian cycles on vertex set verts in tournament A.
    Returns count of distinct directed cycles (up to cyclic rotation)."""
    k = len(verts)
    if k < 3 or k % 2 == 0:
        return 0

    # Fix first vertex, permute rest, check if directed cycle
    v0 = verts[0]
    rest = verts[1:]
    count = 0
    for p in permutations(rest):
        seq = (v0,) + p
        is_cycle = True
        for i in range(k):
            if not A[seq[i]][seq[(i+1) % k]]:
                is_cycle = False
                break
        if is_cycle:
            count += 1
    # Each directed cycle is counted (k-1)!/k * k = (k-1)! times when we
    # fix v0 and permute rest. Actually, fixing v0 and permuting rest,
    # each directed cycle through v0 is counted exactly once (since v0 is fixed
    # and the rest follow in the unique order of the cycle).
    # Wait no: a directed cycle (v0, v1, ..., v_{k-1}) with v0 fixed is
    # uniquely determined by the sequence (v1, ..., v_{k-1}). Different
    # cyclic rotations put different vertices first, but we fixed v0, so
    # each cycle is counted exactly once.
    return count


def build_full_omega(n, A):
    """Build the full Omega(T) graph with ALL odd directed cycles as vertices.
    Returns list of (vertex_set, cycle_id) and adjacency."""
    # For each odd-size subset, count directed cycles
    cycles = []  # Each entry: (frozenset of vertices, index)

    for size in range(3, n+1, 2):
        for verts in combinations(range(n), size):
            num_cycles = count_directed_cycles_on_set(verts, A)
            vset = frozenset(verts)
            for _ in range(num_cycles):
                cycles.append(vset)

    # Build adjacency: two cycles adjacent iff they share a vertex
    nc = len(cycles)
    omega_adj = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                omega_adj[i].add(j)
                omega_adj[j].add(i)

    return cycles, omega_adj


def compute_I_at_2(nodes, omega_adj):
    """Compute I(subgraph, 2) for a set of Omega nodes."""
    sz = len(nodes)
    I_val = 0
    for mask in range(1 << sz):
        selected = [nodes[i] for i in range(sz) if mask & (1 << i)]
        independent = True
        for i in range(len(selected)):
            for j in range(i+1, len(selected)):
                if selected[j] in omega_adj[selected[i]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            I_val += 2 ** bin(mask).count('1')
    return I_val


def verify_full_omega_factorization(n):
    """Verify H(T) = I(Omega(T), 2) with full Omega."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    mismatches = 0
    total = 0

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1

        adj_bits = [0]*n
        for i in range(n):
            for j in range(n):
                if i != j and A[i][j]:
                    adj_bits[i] |= (1 << j)
        H = held_karp(n, adj_bits)

        cycles, omega_adj = build_full_omega(n, A)

        # Compute I(Omega, 2)
        I_total = compute_I_at_2(list(range(len(cycles))), omega_adj)

        if H != I_total:
            mismatches += 1
            if mismatches <= 3:
                print(f"MISMATCH: bits={bits}, H={H}, I={I_total}, #cycles={len(cycles)}")
                for i, c in enumerate(cycles):
                    print(f"  cycle {i}: {set(c)}")

        total += 1

    print(f"n={n}: {total} tournaments, {mismatches} mismatches")
    return mismatches


if __name__ == '__main__':
    # Quick test at n=4 first
    print("=== n=4 (quick test) ===")
    verify_full_omega_factorization(4)

    print("\n=== n=5 ===")
    verify_full_omega_factorization(5)
