#!/usr/bin/env python3
"""
Debug: verify component factorization H = product of I(C_i, 2).

The claim was verified exhaustively at n=5,6 with 0 mismatches.
But we just found tournaments with single component I(C,2)=21 and H != 21.

Something is wrong. Let me re-examine.

WAIT: The OCF says H(T) = I(Omega(T), 2) where Omega(T) has vertices = directed 3-cycles
AND 5-cycles AND 7-cycles AND all odd-length directed cycles.

I was only counting 3-cycles! At n=6 there are also 5-cycles.

The component factorization applies to the FULL Omega(T) with ALL odd cycles,
not just the "3-cycle conflict graph."

Let me redo with ALL odd cycles.

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


def find_all_directed_cycles(n, A):
    """Find ALL directed cycles (as vertex sets) of odd length in tournament A."""
    edge_set = set()
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]:
                edge_set.add((i, j))

    cycles = set()

    # For each subset of vertices of odd size >= 3, check if it supports a directed cycle
    for size in range(3, n+1, 2):  # odd sizes only
        for verts in combinations(range(n), size):
            # Check all cyclic permutations for directed cycles
            # A directed cycle on k vertices: need to find a Hamiltonian cycle
            # in the subtournament on these vertices.
            # For small sizes, enumerate permutations.
            vset = frozenset(verts)

            if size == 3:
                a, b, c = verts
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    cycles.add(vset)
            elif size <= 7:
                # Check if any cyclic permutation gives a directed cycle
                found = False
                v0 = verts[0]
                rest = verts[1:]
                for p in permutations(rest):
                    seq = (v0,) + p
                    is_cycle = True
                    for i in range(size):
                        if (seq[i], seq[(i+1) % size]) not in edge_set:
                            is_cycle = False
                            break
                    if is_cycle:
                        found = True
                        break
                if found:
                    cycles.add(vset)

    return cycles


def verify_factorization(n):
    """Verify H = product of I(C_i, 2) using FULL Omega (all odd cycles)."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    mismatches = 0
    total = 0
    component_21_examples = []

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1

        # Compute H via Held-Karp
        adj_bits = [0]*n
        for i in range(n):
            for j in range(n):
                if i != j and A[i][j]:
                    adj_bits[i] |= (1 << j)
        H = held_karp(n, adj_bits)

        # Find ALL odd directed cycles
        all_cycles = find_all_directed_cycles(n, A)

        if not all_cycles:
            if H != 1:
                print(f"ERROR: no cycles but H={H}")
                mismatches += 1
            total += 1
            continue

        # Build Omega graph
        cycle_list = list(all_cycles)
        nc = len(cycle_list)
        omega_adj = [set() for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if cycle_list[i] & cycle_list[j]:
                    omega_adj[i].add(j)
                    omega_adj[j].add(i)

        # Find components
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

        # Compute I(C_i, 2) for each component
        H_product = 1
        comp_vals = []
        for comp in components:
            sz = len(comp)
            I_val = 0
            for mask in range(1 << sz):
                nodes = [comp[i] for i in range(sz) if mask & (1 << i)]
                independent = True
                for i in range(len(nodes)):
                    for j in range(i+1, len(nodes)):
                        if nodes[j] in omega_adj[nodes[i]]:
                            independent = False
                            break
                    if not independent:
                        break
                if independent:
                    I_val += 2**bin(mask).count('1')
            comp_vals.append(I_val)
            H_product *= I_val

        if H != H_product:
            mismatches += 1
            if mismatches <= 3:
                print(f"MISMATCH: bits={bits}, H={H}, product={H_product}, comps={comp_vals}")
                print(f"  cycles: {[set(c) for c in cycle_list]}")

        # Check for I=21 component
        for cv in comp_vals:
            if cv == 21:
                component_21_examples.append((bits, H, comp_vals))

        total += 1

    print(f"\nn={n}: {total} tournaments, {mismatches} mismatches")
    if component_21_examples:
        print(f"  {len(component_21_examples)} tournaments with I(C,2)=21 component")
        H_vals = set(h for _, h, _ in component_21_examples)
        print(f"  Their H values: {sorted(H_vals)}")
        for b, h, cv in component_21_examples[:5]:
            print(f"    bits={b}: H={h}, components={cv}")
    else:
        print(f"  NO tournaments with I(C,2)=21 component")

    return mismatches


if __name__ == '__main__':
    for n in [5, 6]:
        print(f"=== n={n} ===")
        verify_factorization(n)
        print()
