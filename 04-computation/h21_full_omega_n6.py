#!/usr/bin/env python3
"""
Full Omega component analysis at n=6. Optimized.

At n=6: odd cycle sizes are 3 and 5.
For 3-cycles: each triple has at most 1 directed 3-cycle.
For 5-cycles: each 5-set can have multiple directed 5-cycles.
No 7-cycles since n=6.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations, permutations
from collections import defaultdict
import time


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


def build_omega_n6(A):
    """Build full Omega for n=6 tournament. Cycles: 3-cycles and 5-cycles."""
    n = 6
    cycles = []  # list of frozensets (vertex sets)

    # 3-cycles: exactly 1 per vertex set if cyclic
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    cycles.append(frozenset({a, b, c}))

    # 5-cycles: count directed Hamiltonian cycles per 5-set
    for verts in combinations(range(n), 5):
        v0 = verts[0]
        rest = verts[1:]
        count = 0
        for p in permutations(rest):
            seq = (v0,) + p
            is_cycle = True
            for i in range(5):
                if not A[seq[i]][seq[(i+1) % 5]]:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
        vset = frozenset(verts)
        for _ in range(count):
            cycles.append(vset)

    # Build adjacency
    nc = len(cycles)
    omega_adj = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                omega_adj[i].add(j)
                omega_adj[j].add(i)

    return cycles, omega_adj


def compute_I_at_2(nodes, omega_adj):
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


def main():
    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    total = 1 << m  # 32768

    mismatches = 0
    all_comp_vals = set()
    H_to_comp = defaultdict(set)  # H -> set of component tuples

    t0 = time.time()

    for bits in range(total):
        if bits % 5000 == 0 and bits > 0:
            elapsed = time.time() - t0
            print(f"  progress: {bits}/{total} ({elapsed:.1f}s)")

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

        cycles, omega_adj = build_omega_n6(A)

        if not cycles:
            continue

        nc = len(cycles)
        visited = [False]*nc
        components = []
        for start in range(nc):
            if visited[start]:
                continue
            comp = []
            queue = [start]
            visited[start] = True
            while queue:
                v = queue.pop()
                comp.append(v)
                for u in omega_adj[v]:
                    if not visited[u]:
                        visited[u] = True
                        queue.append(u)
            components.append(comp)

        H_product = 1
        comp_vals = []
        for comp in components:
            I_val = compute_I_at_2(comp, omega_adj)
            comp_vals.append(I_val)
            H_product *= I_val
            all_comp_vals.add(I_val)

        if H != H_product:
            mismatches += 1
            if mismatches <= 3:
                print(f"MISMATCH: bits={bits}, H={H}, product={H_product}")

        H_to_comp[H].add(tuple(sorted(comp_vals)))

    elapsed = time.time() - t0
    print(f"\nn=6: {total} tournaments in {elapsed:.1f}s, {mismatches} mismatches")
    print(f"Component I(C,2) values: {sorted(all_comp_vals)}")
    print(f"I(C,2)=7 achievable: {7 in all_comp_vals}")
    print(f"I(C,2)=21 achievable: {21 in all_comp_vals}")

    # Show component structure for H near 21
    print("\nComponent structures for H in [15..27]:")
    for h in range(15, 28, 2):
        if h in H_to_comp:
            print(f"  H={h}: {sorted(H_to_comp[h])}")
        else:
            print(f"  H={h}: NOT ACHIEVABLE")


if __name__ == '__main__':
    main()
