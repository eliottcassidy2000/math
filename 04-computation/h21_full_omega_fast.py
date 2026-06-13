#!/usr/bin/env python3
"""
Fast full Omega component analysis at n=5,6.

Key optimization: instead of enumerating all permutations for 5-cycles,
use Held-Karp DP to count directed Hamiltonian cycles on each 5-subset.

Directed Hamiltonian cycles on k vertices:
  dp[S][v] = number of paths starting at v0 (fixed), visiting exactly S, ending at v.
  Then #cycles = sum_v dp[full_set][v] where v->v0.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations
from collections import defaultdict
import time


def held_karp_full(n, adj_bits):
    """Full n-vertex Held-Karp for H(T)."""
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


def count_ham_cycles(verts, adj_bits_full):
    """Count directed Hamiltonian cycles on vertex subset verts.
    Fix first vertex, count paths through all others that return to first."""
    k = len(verts)
    v0 = verts[0]

    # Local adjacency for subset
    # Map global vertex indices to local [0..k-1]
    local = {v: i for i, v in enumerate(verts)}

    # Local adj bitmasks
    local_adj = [0] * k
    for i, v in enumerate(verts):
        for j, u in enumerate(verts):
            if i != j and (adj_bits_full[v] & (1 << u)):
                local_adj[i] |= (1 << j)

    # DP: dp[S][v] = #paths starting at local vertex 0, visiting set S, ending at v
    full = (1 << k) - 1
    dp = [[0] * k for _ in range(1 << k)]
    dp[1][0] = 1  # start at local vertex 0

    for S in range(1, 1 << k):
        for v in range(k):
            if not (S & (1 << v)):
                continue
            c = dp[S][v]
            if c == 0:
                continue
            for u in range(k):
                if S & (1 << u):
                    continue
                if local_adj[v] & (1 << u):
                    dp[S | (1 << u)][u] += c

    # Count cycles: paths visiting all vertices, ending at v where v->0
    count = 0
    for v in range(1, k):  # v != 0
        if dp[full][v] > 0 and (local_adj[v] & 1):  # v -> 0
            count += dp[full][v]

    return count


def build_full_omega_fast(n, adj_bits):
    """Build full Omega graph efficiently."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (adj_bits[i] & (1 << j)):
                A[i][j] = 1

    cycles = []  # list of frozensets

    # 3-cycles
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    cycles.append(frozenset({a, b, c}))

    # 5-cycles (and 7-cycles if n>=7)
    for size in range(5, n+1, 2):
        for verts in combinations(range(n), size):
            nc = count_ham_cycles(verts, adj_bits)
            vset = frozenset(verts)
            for _ in range(nc):
                cycles.append(vset)

    # Build adjacency
    total_cycles = len(cycles)
    omega_adj = [set() for _ in range(total_cycles)]
    for i in range(total_cycles):
        for j in range(i+1, total_cycles):
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


def analyze(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    total = 1 << m

    mismatches = 0
    all_comp_vals = set()
    H_to_comp = defaultdict(set)
    comp_val_counts = defaultdict(int)

    t0 = time.time()

    for bits in range(total):
        if bits % 5000 == 0 and bits > 0:
            elapsed = time.time() - t0
            rate = bits / elapsed
            eta = (total - bits) / rate
            print(f"  {bits}/{total} ({elapsed:.0f}s, ETA {eta:.0f}s)")

        adj_bits = [0]*n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj_bits[j] |= (1 << i)
            else:
                adj_bits[i] |= (1 << j)

        H = held_karp_full(n, adj_bits)
        cycles, omega_adj = build_full_omega_fast(n, adj_bits)

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
            comp_val_counts[I_val] += 1

        if H != H_product:
            mismatches += 1
            if mismatches <= 5:
                print(f"MISMATCH: bits={bits}, H={H}, product={H_product}, comps={comp_vals}")

        H_to_comp[H].add(tuple(sorted(comp_vals)))

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.1f}s, {mismatches} mismatches")
    print(f"Component I(C,2) values: {sorted(all_comp_vals)}")
    print(f"I(C,2)=7 achievable: {7 in all_comp_vals}")
    print(f"I(C,2)=21 achievable: {21 in all_comp_vals}")

    print("\nComponent value frequencies:")
    for v in sorted(comp_val_counts.keys()):
        print(f"  I={v}: {comp_val_counts[v]} occurrences")

    print("\nComponent structures for H in [15..27]:")
    for h in range(15, 28, 2):
        if h in H_to_comp:
            print(f"  H={h}: {sorted(H_to_comp[h])}")
        else:
            print(f"  H={h}: NOT ACHIEVABLE")

    # Check multi-component factorizations
    print("\nMulti-component factorizations:")
    for h in sorted(H_to_comp.keys()):
        for ct in H_to_comp[h]:
            if len(ct) > 1:
                print(f"  H={h} = {'*'.join(str(v) for v in ct)}")


if __name__ == '__main__':
    print("=== n=5 ===")
    analyze(5)
    print("\n" + "="*60)
    print("\n=== n=6 ===")
    analyze(6)
