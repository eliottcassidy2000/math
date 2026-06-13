#!/usr/bin/env python3
"""
Check if 5/7-cycles bridge 3-cycle components.

For multi-component tournaments (3-cycle graph disconnected),
check if adding 5-cycles would merge the components.

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


def get_three_cycles(n, adj_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (adj_bits[i] & (1 << j)):
                A[i][j] = 1
    cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    cycles.append(frozenset({a, b, c}))
    return cycles, A


def count_ham_cycles_on_subset(verts, A):
    """Count directed Hamiltonian cycles on subset (fix first vertex)."""
    k = len(verts)
    v0 = verts[0]
    rest = verts[1:]
    count = 0
    for p in permutations(rest):
        seq = (v0,) + p
        ok = True
        for i in range(k):
            if not A[seq[i]][seq[(i+1) % k]]:
                ok = False
                break
        if ok:
            count += 1
    return count


def get_components(cycles):
    if not cycles:
        return []
    nc = len(cycles)
    adj = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i].add(j)
                adj[j].add(i)
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
            for u in adj[v]:
                if not visited[u]:
                    visited[u] = True
                    queue.append(u)
        components.append(comp)
    return components


def check_bridges(n_max=7):
    """For each multi-component tournament, check if 5/7-cycles bridge components."""
    for n in range(5, n_max+1):
        edges = [(i, j) for i in range(n) for j in range(i+1, n)]
        m = len(edges)
        total = 1 << m

        bridged = 0
        unbridged = 0

        for bits in range(total):
            adj_bits = [0]*n
            for k, (i, j) in enumerate(edges):
                if bits & (1 << k):
                    adj_bits[j] |= (1 << i)
                else:
                    adj_bits[i] |= (1 << j)

            three_cycles, A = get_three_cycles(n, adj_bits)
            if not three_cycles:
                continue

            comps = get_components(three_cycles)
            if len(comps) <= 1:
                continue

            # Multi-component 3-cycle graph. Check for bridging 5/7-cycles.
            comp_verts = []
            for comp in comps:
                verts = set()
                for idx in comp:
                    verts |= three_cycles[idx]
                comp_verts.append(verts)

            # Check 5-cycles (and 7-cycles if n>=7)
            bridges_found = False
            for size in range(5, n+1, 2):
                for verts in combinations(range(n), size):
                    nc = count_ham_cycles_on_subset(verts, A)
                    if nc == 0:
                        continue
                    vset = frozenset(verts)
                    # Check if this cycle touches multiple 3-cycle components
                    touches = set()
                    for ci, cv in enumerate(comp_verts):
                        if vset & cv:
                            touches.add(ci)
                    if len(touches) > 1:
                        bridges_found = True
                        break
                if bridges_found:
                    break

            if bridges_found:
                bridged += 1
                H = held_karp(n, adj_bits)
                H3 = 1  # H from 3-cycles only
                # Actually let me just report
                if bridged <= 3:
                    print(f"  n={n}: BRIDGED tournament bits={bits}, H={held_karp(n, adj_bits)}")
                    print(f"    3-cycle components: {[sorted(cv) for cv in comp_verts]}")
            else:
                unbridged += 1

        print(f"n={n}: {bridged} bridged, {unbridged} truly disconnected")


if __name__ == '__main__':
    check_bridges(7)
