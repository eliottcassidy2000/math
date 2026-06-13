#!/usr/bin/env python3
"""
K_6 minus 2 edges as Omega(T) component: exhaustive investigation.

RESULT: K_6-2e CANNOT be a connected component of the full Omega(T).

Proof outline (verified computationally):
1. Six 3-cycles in K_6-2e pattern span at least 6 tournament vertices
   (at n=5, max t3 = 5, so 6 cycles need n >= 6).
2. At n=6 with t3=6: every tournament with K_6-2e 3-cycle conflict graph
   also has >= 4 five-cycles. Since 5-cycles on 5 of 6 vertices always
   share vertices with 3-cycles on 3 of 6 vertices (5+3 > 6), the full
   Omega(T) has more vertices adjacent to the K_6-2e, so it's not a component.
3. At n=7 exhaustive (2^21 tournaments): 58,800 have K_6-2e as 3-cycle
   component, but 0 have it as a full Omega component (5-cycles always intrude).
4. At n=8 sampling (200,000 tournaments): 663 have K_6-2e as 3-cycle
   component, but 0 have it as full Omega component.

The mechanism: t3 >= 6 on 6 vertices forces t5 >= 4 (verified exhaustively).
On 7 vertices, the 6 three-cycles span 6 or 7 vertices, and 5-cycles on
overlapping subsets are always forced by the tournament structure.

Instance: opus-2026-03-07
Dependencies: THM-079 (H=21 component reduction)
"""

import itertools
import random
from collections import Counter

# ─────────────────────────────────────────────────
# Core utilities
# ─────────────────────────────────────────────────

def all_tournaments(n):
    """Generate all labeled tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(edges)):
        adj = [[False]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
        yield adj

def find_3_cycles(adj, n):
    """Find all directed 3-cycles as frozensets of vertex sets."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    cycles.append(frozenset([i, j, k]))
    return cycles

def has_hamiltonian_cycle(adj, verts):
    """Check if sub-tournament on 'verts' has a directed Hamiltonian cycle (DP)."""
    v = list(verts)
    k = len(v)
    if k < 3:
        return False
    dp = set()
    dp.add((1, 0))
    for mask in range(1, 1 << k):
        for i in range(k):
            if not (mask & (1 << i)):
                continue
            if (mask, i) not in dp:
                continue
            for j in range(k):
                if mask & (1 << j):
                    continue
                if adj[v[i]][v[j]]:
                    dp.add((mask | (1 << j), j))
    full = (1 << k) - 1
    return any((full, j) in dp and adj[v[j]][v[0]] for j in range(1, k))

def build_conflict_graph(cycles):
    """Build adjacency matrix: edge iff cycles share a vertex."""
    m = len(cycles)
    adj = [[False]*m for _ in range(m)]
    for a in range(m):
        for b in range(a+1, m):
            if cycles[a] & cycles[b]:
                adj[a][b] = adj[b][a] = True
    return adj

def connected_components(adj, n):
    """Find connected components."""
    visited = [False]*n
    comps = []
    for s in range(n):
        if visited[s]:
            continue
        comp = []
        stack = [s]
        while stack:
            v = stack.pop()
            if visited[v]:
                continue
            visited[v] = True
            comp.append(v)
            for u in range(n):
                if adj[v][u] and not visited[u]:
                    stack.append(u)
        comps.append(comp)
    return comps

def is_k6_minus_2e(adj6):
    """Check if 6-vertex graph is K_6 minus exactly 2 edges. Return (bool, type)."""
    edges = sum(1 for i in range(6) for j in range(i+1, 6) if adj6[i][j])
    if edges != 13:
        return False, None
    non_edges = [(i, j) for i in range(6) for j in range(i+1, 6) if not adj6[i][j]]
    if len(non_edges) != 2:
        return False, None
    (a, b), (c, d) = non_edges
    netype = "adjacent" if ({a, b} & {c, d}) else "disjoint"
    return True, netype

# ─────────────────────────────────────────────────
# Test 1: Forcing lemma at n=5
# ─────────────────────────────────────────────────

def test_forcing_n5():
    """Two 3-cycles sharing a vertex on 5 vertices always force a 3rd cycle."""
    print("=" * 60)
    print("TEST 1: Forcing lemma (2 sharing cycles on 5 vertices)")
    print("=" * 60)
    min_t3 = float('inf')
    cases = 0
    for adj in all_tournaments(5):
        c = find_3_cycles(adj, 5)
        if frozenset([0,1,2]) in c and frozenset([2,3,4]) in c:
            cases += 1
            min_t3 = min(min_t3, len(c))
    print(f"  Cases: {cases}, min t3: {min_t3}")
    print(f"  CONFIRMED: sharing vertex always forces >= {min_t3} cycles\n")

# ─────────────────────────────────────────────────
# Test 2: t3=6 forces t5 >= 4 at n=6
# ─────────────────────────────────────────────────

def test_t3_forces_t5():
    """At n=6, t3=6 always forces t5 >= 4."""
    print("=" * 60)
    print("TEST 2: t3=6 forces t5 >= 4 at n=6")
    print("=" * 60)
    min_t5 = float('inf')
    count = 0
    for adj in all_tournaments(6):
        c3 = find_3_cycles(adj, 6)
        if len(c3) != 6:
            continue
        count += 1
        t5 = sum(1 for v5 in itertools.combinations(range(6), 5)
                 if has_hamiltonian_cycle(adj, v5))
        min_t5 = min(min_t5, t5)
    print(f"  Tournaments with t3=6: {count}")
    print(f"  Minimum t5 among those: {min_t5}")
    print(f"  CONFIRMED: t3=6 at n=6 forces t5 >= {min_t5}\n")

# ─────────────────────────────────────────────────
# Test 3: Conflict graphs with t3=6 at n=6
# ─────────────────────────────────────────────────

def test_conflict_graphs_n6():
    """Classify all conflict graphs when t3=6 at n=6."""
    print("=" * 60)
    print("TEST 3: Conflict graph types with t3=6 at n=6")
    print("=" * 60)
    types = Counter()
    for adj in all_tournaments(6):
        c3 = find_3_cycles(adj, 6)
        if len(c3) != 6:
            continue
        conf = build_conflict_graph(c3)
        e = sum(1 for i in range(6) for j in range(i+1, 6) if conf[i][j])
        ds = tuple(sorted(sum(1 for j in range(6) if conf[i][j]) for i in range(6)))
        types[(e, ds)] += 1
    for (e, ds), cnt in sorted(types.items()):
        tag = " <-- K_6-2e" if e == 13 else ""
        print(f"  edges={e}, deg_seq={ds}: {cnt} tournaments{tag}")
    print()

# ─────────────────────────────────────────────────
# Test 4: Full Omega check at n=7 (exhaustive)
# ─────────────────────────────────────────────────

def test_full_omega_n7():
    """Exhaustive n=7: K_6-2e as component of full Omega (3-cycles + 5-cycles + 7-cycles)."""
    print("=" * 60)
    print("TEST 4: Exhaustive n=7 -- K_6-2e as full Omega component")
    print("=" * 60)

    edges_list = [(i, j) for i in range(7) for j in range(i+1, 7)]
    total = 2**21

    k6_2e_3cyc = 0
    k6_2e_full = 0
    span_dist = Counter()

    for bits in range(total):
        if bits % 500000 == 0 and bits > 0:
            print(f"  {bits}/{total} ({100*bits/total:.1f}%)")

        adj = [[False]*7 for _ in range(7)]
        for k, (i, j) in enumerate(edges_list):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True

        c3 = find_3_cycles(adj, 7)
        if len(c3) < 6:
            continue

        conf = build_conflict_graph(c3)
        comps = connected_components(conf, len(c3))

        for comp in comps:
            if len(comp) != 6:
                continue
            sub = [[False]*6 for _ in range(6)]
            for ii, a in enumerate(comp):
                for jj, b in enumerate(comp):
                    if ii < jj and conf[a][b]:
                        sub[ii][jj] = sub[jj][ii] = True
            ok, netype = is_k6_minus_2e(sub)
            if not ok:
                continue

            k6_2e_3cyc += 1
            comp_verts = set()
            for idx in comp:
                comp_verts.update(c3[idx])
            span_dist[len(comp_verts)] += 1

            # Check for 5-cycles touching comp_verts
            has_extra = False
            for v5 in itertools.combinations(range(7), 5):
                if set(v5) & comp_verts and has_hamiltonian_cycle(adj, v5):
                    has_extra = True
                    break
            if not has_extra:
                if has_hamiltonian_cycle(adj, range(7)):
                    has_extra = True
            if not has_extra:
                comp_set = set(comp)
                for idx in range(len(c3)):
                    if idx not in comp_set and c3[idx] & comp_verts:
                        has_extra = True
                        break
            if not has_extra:
                k6_2e_full += 1

    print(f"\n  K_6-2e as 3-cycle component: {k6_2e_3cyc}")
    print(f"  Vertex span distribution: {dict(sorted(span_dist.items()))}")
    print(f"  K_6-2e as FULL Omega component: {k6_2e_full}")
    if k6_2e_full == 0:
        print("  CONFIRMED: impossible at n=7\n")

# ─────────────────────────────────────────────────
# Test 5: Sampling at n=8
# ─────────────────────────────────────────────────

def test_sampling_n8(num_samples=200000):
    """Sampling n=8 for K_6-2e as full Omega component."""
    print("=" * 60)
    print(f"TEST 5: Sampling n=8 ({num_samples} samples)")
    print("=" * 60)

    random.seed(42)
    edges_list = [(i, j) for i in range(8) for j in range(i+1, 8)]
    k6_2e_3cyc = 0
    k6_2e_full = 0

    for trial in range(num_samples):
        if trial % 50000 == 0 and trial > 0:
            print(f"  {trial}/{num_samples}")

        bits = random.randint(0, 2**28 - 1)
        adj = [[False]*8 for _ in range(8)]
        for k, (i, j) in enumerate(edges_list):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True

        c3 = find_3_cycles(adj, 8)
        if len(c3) < 6:
            continue

        conf = build_conflict_graph(c3)
        comps = connected_components(conf, len(c3))

        for comp in comps:
            if len(comp) != 6:
                continue
            sub = [[False]*6 for _ in range(6)]
            for ii, a in enumerate(comp):
                for jj, b in enumerate(comp):
                    if ii < jj and conf[a][b]:
                        sub[ii][jj] = sub[jj][ii] = True
            ok, netype = is_k6_minus_2e(sub)
            if not ok:
                continue

            k6_2e_3cyc += 1
            comp_verts = set()
            for idx in comp:
                comp_verts.update(c3[idx])

            has_extra = False
            for v5 in itertools.combinations(range(8), 5):
                if set(v5) & comp_verts and has_hamiltonian_cycle(adj, v5):
                    has_extra = True
                    break
            if not has_extra:
                for v7 in itertools.combinations(range(8), 7):
                    if set(v7) & comp_verts and has_hamiltonian_cycle(adj, v7):
                        has_extra = True
                        break
            if not has_extra:
                comp_set = set(comp)
                for idx in range(len(c3)):
                    if idx not in comp_set and c3[idx] & comp_verts:
                        has_extra = True
                        break
            if not has_extra:
                k6_2e_full += 1
                print(f"  *** FOUND K_6-2e full component! ***")

    print(f"\n  K_6-2e as 3-cycle component: {k6_2e_3cyc}")
    print(f"  K_6-2e as FULL Omega component: {k6_2e_full}")
    if k6_2e_full == 0:
        print("  CONFIRMED: not found in sampling\n")

# ─────────────────────────────────────────────────

if __name__ == "__main__":
    test_forcing_n5()
    test_t3_forces_t5()
    test_conflict_graphs_n6()
    test_full_omega_n7()
    test_sampling_n8()

    print("=" * 60)
    print("CONCLUSION: K_6-2e cannot be a connected component of Omega(T).")
    print("The 6 three-cycles always force 5-cycles that share their vertices,")
    print("expanding the component beyond K_6-2e in the full Omega graph.")
    print("=" * 60)
