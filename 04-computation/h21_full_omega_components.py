#!/usr/bin/env python3
"""
Component factorization with full Omega(T) and check for I(C,2)=21.

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
    k = len(verts)
    if k < 3 or k % 2 == 0:
        return 0
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
    return count


def build_full_omega(n, A):
    cycles = []
    for size in range(3, n+1, 2):
        for verts in combinations(range(n), size):
            num = count_directed_cycles_on_set(verts, A)
            vset = frozenset(verts)
            for _ in range(num):
                cycles.append(vset)
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


def analyze_components(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    mismatches = 0
    total = 0
    all_component_vals = set()
    comp_21_count = 0
    comp_7_count = 0

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

        if not cycles:
            if H != 1:
                mismatches += 1
            total += 1
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
            all_component_vals.add(I_val)

        if H != H_product:
            mismatches += 1
            if mismatches <= 3:
                print(f"MISMATCH: bits={bits}, H={H}, product={H_product}, comps={comp_vals}")

        if 21 in comp_vals:
            comp_21_count += 1
        if 7 in comp_vals:
            comp_7_count += 1

        total += 1

    print(f"n={n}: {total} tournaments, {mismatches} factorization mismatches")
    print(f"  All component I(C,2) values: {sorted(all_component_vals)}")
    print(f"  Tournaments with I(C,2)=7 component: {comp_7_count}")
    print(f"  Tournaments with I(C,2)=21 component: {comp_21_count}")


if __name__ == '__main__':
    for n in [4, 5]:
        print(f"=== n={n} ===")
        analyze_components(n)
        print()
