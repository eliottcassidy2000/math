#!/usr/bin/env python3
"""
Find tournaments at n=6 where a component of Omega has I(C,2)=21.
Then check what H(T) is for those tournaments.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations
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


def analyze_n6_component_21():
    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    examples = []

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1

        # Find 3-cycles
        cycles = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if (A[a][b] and A[b][c] and A[c][a]) or \
                       (A[a][c] and A[c][b] and A[b][a]):
                        cycles.append(frozenset({a, b, c}))

        if not cycles:
            continue

        nc = len(cycles)
        omega_adj = [set() for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if cycles[i] & cycles[j]:
                    omega_adj[i].add(j)
                    omega_adj[j].add(i)

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

        comp_vals = []
        has_21 = False
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
            if I_val == 21:
                has_21 = True

        if has_21:
            # Compute H
            adj_bits = [0]*n
            for i in range(n):
                for j in range(n):
                    if i != j and A[i][j]:
                        adj_bits[i] |= (1 << j)
            H = held_karp(n, adj_bits)

            # Product of components
            H_product = 1
            for v in comp_vals:
                H_product *= v

            examples.append((bits, H, comp_vals, len(cycles), cycles))

    print(f"Found {len(examples)} tournaments at n=6 with I(C,2)=21 component")
    print()

    # Show first few
    for bits, H, comp_vals, ncycles, cycles in examples[:10]:
        print(f"bits={bits}: H={H}, components={comp_vals}, product={1}")
        product = 1
        for v in comp_vals:
            product *= v
        print(f"  H={H}, product of components={product}, match={H==product}")
        print(f"  {ncycles} 3-cycles: {[set(c) for c in cycles]}")

        # Show the component with I=21
        for comp_idx, cv in enumerate(comp_vals):
            if cv == 21:
                print(f"  Component with I=21 has {len([c for c in cycles])} total cycles")
                break
        print()

    # Summarize H values
    H_vals = set(h for _, h, _, _, _ in examples)
    print(f"H values of tournaments with I(C,2)=21 component: {sorted(H_vals)}")
    print(f"Note: H is product of ALL component values.")
    print(f"If I(C,2)=21 is the ONLY component, then H=21.")
    print(f"If there are other components, H = 21 * product_of_others.")

    single_comp = [(b, h, cv) for b, h, cv, _, _ in examples if len(cv) == 1]
    multi_comp = [(b, h, cv) for b, h, cv, _, _ in examples if len(cv) > 1]
    print(f"\nSingle-component (H=21 directly): {len(single_comp)}")
    print(f"Multi-component (H=21*k): {len(multi_comp)}")

    if single_comp:
        print("FOUND tournaments with H=21!")
        for b, h, cv in single_comp[:5]:
            print(f"  bits={b}: H={h}, components={cv}")
    else:
        print("NO tournaments with H=21 (I=21 always comes with other components)")
        print("Multi-component examples:")
        for b, h, cv in multi_comp[:10]:
            print(f"  bits={b}: H={h}, components={cv}")


if __name__ == '__main__':
    analyze_n6_component_21()
