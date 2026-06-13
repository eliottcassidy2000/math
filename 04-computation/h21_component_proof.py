#!/usr/bin/env python3
"""
Prove I(C,2)=21 impossible for any connected component of Omega(T).

Connected graphs G with I(G,2)=21:
I(G,2) = 1 + 2*|V| + 4*i_2 + 8*i_3 + ...
=> |V| + 2*i_2 + 4*i_3 + ... = 10

Complete list (connected only):
1. K_10: 10 verts, 0 non-edges
2. K_8 - e: 8 verts, 1 non-edge
3. K_6 - M: 6 verts, 2 disjoint non-edges (matching)
4. K_6 - P3: 6 verts, 2 adjacent non-edges (sharing endpoint)
5. complement(P_4): 4 verts, 3 non-edges (complement of path)

Each must be shown impossible as a connected component of Omega(T).

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations
from collections import defaultdict


def three_cycles_through_vertex(n, d):
    """Number of 3-cycles through a vertex of out-degree d in n-tournament."""
    return (n-1)*(n-2)//2 - d*(d-1)//2 - (n-1-d)*(n-2-d)//2


def max_three_cycles_per_vertex(n):
    """Max 3-cycles through any single vertex."""
    return max(three_cycles_through_vertex(n, d) for d in range(n))


def exhaustive_component_check(n):
    """Check if I(C,2)=21 ever appears as a component value."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    component_vals = set()

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

            component_vals.add(I_val)

    return sorted(component_vals)


def analyze_vertex_bounds():
    """Analyze whether K_10 component is possible based on vertex counts."""
    print("=== 3-cycles per vertex bounds ===")
    for n in range(5, 12):
        mx = max_three_cycles_per_vertex(n)
        print(f"  n={n}: max 3-cycles through one vertex = {mx}")

    print()
    print("K_10 requires 10 pairwise-intersecting 3-cycles.")
    print("If they share a common vertex v, need >= 10 through v.")
    print("At n=7: max is 9. IMPOSSIBLE at n=7.")
    print("At n=8: max is 12. Possible in principle, but component isolation needed.")
    print()

    print("K_8-e requires 8 3-cycles, 7 pairwise-intersecting + 1 disjoint pair.")
    print("At n=7: max through one vertex is 9. Could have 7 through common vertex")
    print("  plus 1 disjoint. But disjoint cycle uses 3 OTHER vertices (only 6 available).")


def check_4cycle_exactly_constraint(n):
    """
    Can a tournament on n vertices have EXACTLY 4 directed 3-cycles?
    complement(P4) as component requires exactly 4 3-cycles forming that pattern.
    """
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    count = 0
    patterns = defaultdict(int)

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1

        cycles = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if (A[a][b] and A[b][c] and A[c][a]) or \
                       (A[a][c] and A[c][b] and A[b][a]):
                        cycles.append(frozenset({a, b, c}))

        if len(cycles) != 4:
            continue

        count += 1

        # Check adjacency pattern
        adj_edges = set()
        for i in range(4):
            for j in range(i+1, 4):
                if cycles[i] & cycles[j]:
                    adj_edges.add((i, j))

        num_edges = len(adj_edges)
        # For complement(P4): 3 edges, 3 non-edges
        if num_edges == 3:
            # Check if it is complement(P4)
            # complement(P4) has degree sequence [2,1,1,2] (sorted: [1,1,2,2])
            degs = [0]*4
            for (i, j) in adj_edges:
                degs[i] += 1
                degs[j] += 1
            sorted_degs = sorted(degs)
            if sorted_degs == [1, 1, 2, 2]:
                patterns['complement_P4'] += 1
            else:
                patterns[f'3edges_deg{sorted_degs}'] += 1
        else:
            patterns[f'{num_edges}edges'] += 1

    return count, patterns


if __name__ == '__main__':
    analyze_vertex_bounds()

    print("\n" + "="*60)
    print("\n=== Exhaustive component I(C,2) values ===")
    for n in [5, 6]:
        vals = exhaustive_component_check(n)
        has_21 = 21 in vals
        print(f"n={n}: I(C,2) values = {vals}")
        print(f"  I(C,2)=21 found: {has_21}")

    print("\n" + "="*60)
    print("\n=== Tournaments with exactly 4 directed 3-cycles ===")
    for n in [5, 6, 7]:
        cnt, pats = check_4cycle_exactly_constraint(n)
        print(f"n={n}: {cnt} tournaments with exactly 4 3-cycles")
        for p, c in sorted(pats.items()):
            print(f"  pattern {p}: {c}")
