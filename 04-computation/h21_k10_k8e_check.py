#!/usr/bin/env python3
"""
h21_k10_k8e_check.py
opus-2026-03-07

Verify that (10,0) and (8,1) decompositions of I(Omega(T),2) = 21 are impossible.

Key facts from THM-079:
  - Omega(T) = conflict graph on ALL directed odd cycles (3, 5, 7, ...)
  - H(T) = I(Omega(T), 2)
  - For H=21, need I(Omega(T), 2) = 21
  - Remaining 3-cycle skeleton decompositions: (10,0)=K_10, (8,1)=K_8-e, (6,2)=K_6-2e
  - At n=6, K_6-2e in 3-cycle conflict always forces 5-cycles, pushing I >= 29

This script checks n=7 (exhaustive) and n=8 (sampling):
  - For tournaments whose 3-CYCLE conflict graph has K_10 or K_8-e pattern:
    * How many 5-cycles and 7-cycles are forced?
    * What is I(full Omega(T), 2)?
    * Can I(Omega(T), 2) = 21?

Uses bitmask adjacency for speed.
"""

import random
from collections import defaultdict

def tournament_from_index(n, idx):
    """Create bitmask adjacency from tournament index."""
    adj = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i+1, n):
            if idx & (1 << bit):
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
            bit += 1
    return adj

def random_tournament(n):
    adj = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            if random.randint(0, 1):
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
    return adj

def find_3cycles(adj, n):
    """Find all directed 3-cycles as sorted vertex triples."""
    cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (adj[a] >> b & 1) and (adj[b] >> c & 1) and (adj[c] >> a & 1):
                    cycles.append((a, b, c))
                elif (adj[a] >> c & 1) and (adj[c] >> b & 1) and (adj[b] >> a & 1):
                    cycles.append((a, b, c))
    return cycles

def find_directed_5cycles(adj, n):
    """Find all directed 5-cycles. Return as frozensets of vertices.
    Each vertex SET supporting a 5-cycle appears once (even if multiple
    directed orderings exist), because all directed cycles on the same
    vertex set are the same vertex in our conflict model (they all share
    all vertices, so they're adjacent in Omega and form a clique).

    Actually, each DIRECTED cycle is a separate vertex in Omega.
    We need to count directed cycles, not vertex sets.
    A 5-cycle a->b->c->d->e->a. We enumerate by checking all 5-permutations."""

    # Enumerate directed 5-cycles: sequences (v0,v1,v2,v3,v4) with
    # v0->v1->v2->v3->v4->v0, and we normalize by choosing the rotation
    # starting with the smallest vertex.
    cycles = set()
    verts = list(range(n))

    from itertools import permutations
    for perm in permutations(verts, 5):
        v0, v1, v2, v3, v4 = perm
        if ((adj[v0] >> v1 & 1) and (adj[v1] >> v2 & 1) and
            (adj[v2] >> v3 & 1) and (adj[v3] >> v4 & 1) and (adj[v4] >> v0 & 1)):
            # Normalize: rotate to start with smallest, then pick direction
            # giving the smaller second element
            cycle = [v0, v1, v2, v3, v4]
            # Find rotation starting with min
            min_v = min(cycle)
            idx_min = cycle.index(min_v)
            rotated = tuple(cycle[idx_min:] + cycle[:idx_min])
            # Also consider reverse direction
            rev = tuple(cycle[(idx_min) % 5 - i] if i > 0 else cycle[idx_min] for i in range(5))
            # Simpler: just use frozenset of the directed edges as canonical form
            # Or: use the vertex SET as the canonical form since all directed cycles
            # on the same 5 vertices share all vertices and are thus adjacent in Omega

            # Actually for Omega, each directed cycle is a separate vertex.
            # But cycles on the same vertex set are ALL pairwise adjacent (share all vertices).
            # For the independence polynomial, a clique of size k contributes
            # I(K_k, x) = 1 + kx. So we need to count directed cycles correctly.
            #
            # Canonical form: rotate to min vertex, then pick the direction where
            # the second vertex is smaller.
            rotated_fwd = tuple(cycle[(idx_min + i) % 5] for i in range(5))
            rotated_bwd = tuple(cycle[(idx_min - i) % 5] for i in range(5))
            canonical = min(rotated_fwd, rotated_bwd)
            cycles.add(canonical)

    return list(cycles)

def find_directed_5cycles_fast(adj, n):
    """Faster: enumerate directed 5-cycles by vertex subsets, then check orientations."""
    from itertools import combinations
    cycles = []
    for subset in combinations(range(n), 5):
        # Check all possible directed 5-cycles on this subset
        # A directed 5-cycle is a permutation (v0,v1,v2,v3,v4) with v0->v1->v2->v3->v4->v0
        # On 5 vertices there are 4!/2 = 12 directed 5-cycles (up to rotation and reflection)
        # Actually: number of distinct directed Hamiltonian cycles on 5 vertices = (5-1)!/2 = 12
        # We check each
        from itertools import permutations
        found = []
        for perm in permutations(subset):
            v0, v1, v2, v3, v4 = perm
            if ((adj[v0] >> v1 & 1) and (adj[v1] >> v2 & 1) and
                (adj[v2] >> v3 & 1) and (adj[v3] >> v4 & 1) and (adj[v4] >> v0 & 1)):
                # Canonicalize
                cycle = list(perm)
                min_v = min(cycle)
                idx_min = cycle.index(min_v)
                rotated_fwd = tuple(cycle[(idx_min + i) % 5] for i in range(5))
                rotated_bwd = tuple(cycle[(idx_min - i) % 5] for i in range(5))
                canonical = min(rotated_fwd, rotated_bwd)
                found.append(canonical)
        # Deduplicate within this subset
        for c in set(found):
            cycles.append(c)
    return cycles

def find_directed_7cycles_fast(adj, n):
    """Find directed 7-cycles (only relevant at n=7 where there's exactly one 7-subset)."""
    if n < 7:
        return []
    from itertools import combinations, permutations
    cycles = []
    for subset in combinations(range(n), 7):
        for perm in permutations(subset):
            v = perm
            ok = True
            for i in range(7):
                if not (adj[v[i]] >> v[(i+1) % 7] & 1):
                    ok = False
                    break
            if ok:
                cycle = list(perm)
                min_v = min(cycle)
                idx_min = cycle.index(min_v)
                fwd = tuple(cycle[(idx_min + i) % 7] for i in range(7))
                bwd = tuple(cycle[(idx_min - i) % 7] for i in range(7))
                canonical = min(fwd, bwd)
                cycles.append(canonical)
    return list(set(cycles))

def count_disjoint_pairs_triples(cycles3):
    """Count pairs of 3-cycles with no shared vertex."""
    m = len(cycles3)
    disjoint = 0
    for i in range(m):
        si = set(cycles3[i])
        for j in range(i+1, m):
            if not (si & set(cycles3[j])):
                disjoint += 1
    return disjoint

def build_full_omega_and_compute_I2(all_cycles):
    """Given list of all odd cycles (as tuples of vertices),
    build Omega and compute I(Omega, 2) via recursive DP with memoization."""
    m = len(all_cycles)
    if m == 0:
        return 1

    # Build conflict adjacency using bitmasks
    vsets = [set(c) for c in all_cycles]
    conflict = [0] * m
    for i in range(m):
        for j in range(i+1, m):
            if vsets[i] & vsets[j]:
                conflict[i] |= (1 << j)
                conflict[j] |= (1 << i)

    # Use recursive DP with memoization on vertex subsets
    # I(G[S], 2) for vertex subset S (bitmask)
    # Recurrence: pick lowest-index vertex v in S,
    #   I(G[S], 2) = I(G[S\{v}], 2) + 2 * I(G[S \ N[v]], 2)
    memo = {}

    def dp(subset):
        if subset == 0:
            return 1
        if subset in memo:
            return memo[subset]

        # Pick the lowest vertex in subset
        v = (subset & -subset).bit_length() - 1
        s_minus_v = subset & ~(1 << v)
        s_minus_Nv = subset & ~((1 << v) | conflict[v])

        result = dp(s_minus_v) + 2 * dp(s_minus_Nv)
        memo[subset] = result
        return result

    full = (1 << m) - 1
    return dp(full)

def count_hamiltonian_paths(adj, n):
    """Count Hamiltonian paths using Held-Karp DP."""
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            out = adj[v] & ~mask & full
            tmp = out
            while tmp:
                w = (tmp & -tmp).bit_length() - 1
                dp[mask | (1 << w)][w] += dp[mask][v]
                tmp &= tmp - 1
    return sum(dp[full])

# ========================================================
# Main analysis
# ========================================================

def analyze_tournament(adj, n, compute_H=False):
    """Full analysis: find all odd cycles, check 3-cycle skeleton, compute I(Omega, 2)."""
    c3 = find_3cycles(adj, n)
    t3 = len(c3)

    # 3-cycle conflict pattern
    disjoint3 = count_disjoint_pairs_triples(c3)

    # Check skeleton patterns
    is_k10_skeleton = (t3 == 10 and disjoint3 == 0)
    is_k8e_skeleton = (t3 == 8 and disjoint3 == 1)
    is_k6_2e_skeleton = (t3 == 6 and disjoint3 == 2)

    result = {
        't3': t3, 'disjoint3': disjoint3,
        'is_k10': is_k10_skeleton, 'is_k8e': is_k8e_skeleton,
        'is_k6_2e': is_k6_2e_skeleton,
    }

    # If skeleton matches, compute full Omega
    if is_k10_skeleton or is_k8e_skeleton or is_k6_2e_skeleton:
        c5 = find_directed_5cycles_fast(adj, n)
        t5 = len(c5)
        c7 = find_directed_7cycles_fast(adj, n) if n >= 7 else []
        t7 = len(c7)

        all_cycles = [(c, 3) for c in c3] + [(c, 5) for c in c5] + [(c, 7) for c in c7]
        all_cycle_verts = [c for c, _ in all_cycles]

        total_cycles = len(all_cycle_verts)
        I2 = build_full_omega_and_compute_I2(all_cycle_verts)

        result['t5'] = t5
        result['t7'] = t7
        result['total_cycles'] = total_cycles
        result['I2'] = I2

        if compute_H:
            H = count_hamiltonian_paths(adj, n)
            result['H'] = H

    return result

def full_omega_analysis(adj, n):
    """For a tournament matching a skeleton pattern, compute full Omega and I(Omega,2)."""
    c3 = find_3cycles(adj, n)
    c5 = find_directed_5cycles_fast(adj, n)
    c7 = find_directed_7cycles_fast(adj, n) if n >= 7 else []
    all_cycles = list(c3) + list(c5) + list(c7)
    I2 = build_full_omega_and_compute_I2(all_cycles)
    return len(c3), len(c5), len(c7), len(all_cycles), I2

def run_n7_exhaustive():
    """Exhaustively check all n=7 tournaments."""
    n = 7
    num_edges = n * (n-1) // 2
    total = 1 << num_edges

    print(f"=== n=7 Exhaustive Check ({total} tournaments) ===")
    print()

    k10_count = 0
    k8e_count = 0
    k6_2e_count = 0

    I2_dist_k10 = defaultdict(int)
    I2_dist_k8e = defaultdict(int)
    t5_dist_k10 = defaultdict(int)
    t5_dist_k8e = defaultdict(int)

    for idx in range(total):
        if idx % 500000 == 0:
            print(f"  Progress: {idx}/{total} ({100*idx/total:.1f}%)", flush=True)

        adj = tournament_from_index(n, idx)
        c3 = find_3cycles(adj, n)
        t3 = len(c3)

        # Quick screen: only compute disjoint pairs for m=10 or m=8 or m=6
        if t3 == 10:
            d = count_disjoint_pairs_triples(c3)
            if d == 0:
                # K_10 skeleton! Do full analysis.
                k10_count += 1
                _, t5, t7, total_c, I2 = full_omega_analysis(adj, n)
                I2_dist_k10[I2] += 1
                t5_dist_k10[t5] += 1
                if k10_count <= 3:
                    print(f"  K_10 skeleton: t3=10, t5={t5}, t7={t7}, "
                          f"total={total_c}, I(Omega,2)={I2}")

        elif t3 == 8:
            d = count_disjoint_pairs_triples(c3)
            if d == 1:
                k8e_count += 1
                _, t5, t7, total_c, I2 = full_omega_analysis(adj, n)
                I2_dist_k8e[I2] += 1
                t5_dist_k8e[t5] += 1
                if k8e_count <= 3:
                    print(f"  K_8-e skeleton: t3=8, t5={t5}, t7={t7}, "
                          f"total={total_c}, I(Omega,2)={I2}")

        elif t3 == 6:
            d = count_disjoint_pairs_triples(c3)
            if d == 2:
                k6_2e_count += 1

    print()
    print(f"n=7 Results:")
    print(f"  K_10 skeleton count: {k10_count}")
    if I2_dist_k10:
        print(f"    I(Omega,2) distribution: {dict(sorted(I2_dist_k10.items()))}")
        print(f"    t5 distribution: {dict(sorted(t5_dist_k10.items()))}")
        min_I2 = min(I2_dist_k10.keys())
        print(f"    Min I(Omega,2) = {min_I2}  {'<= 21 POSSIBLE!' if min_I2 <= 21 else '> 21 IMPOSSIBLE'}")

    print(f"  K_8-e skeleton count: {k8e_count}")
    if I2_dist_k8e:
        print(f"    I(Omega,2) distribution: {dict(sorted(I2_dist_k8e.items()))}")
        print(f"    t5 distribution: {dict(sorted(t5_dist_k8e.items()))}")
        min_I2 = min(I2_dist_k8e.keys())
        print(f"    Min I(Omega,2) = {min_I2}  {'<= 21 POSSIBLE!' if min_I2 <= 21 else '> 21 IMPOSSIBLE'}")

    print(f"  K_6-2e skeleton count: {k6_2e_count}")
    print()
    return k10_count, k8e_count

def run_n8_sampling(num_samples=100000):
    """Sample n=8 tournaments, check K_10 and K_8-e skeletons, compute full I(Omega,2)."""
    n = 8
    print(f"=== n=8 Sampling ({num_samples} tournaments) ===")
    print()

    k10_count = 0
    k8e_count = 0

    I2_dist_k10 = defaultdict(int)
    I2_dist_k8e = defaultdict(int)
    t5_dist_k10 = defaultdict(int)
    t5_dist_k8e = defaultdict(int)
    H_dist_k10 = defaultdict(int)
    H_dist_k8e = defaultdict(int)

    for trial in range(num_samples):
        if trial % 25000 == 0 and trial > 0:
            print(f"  Progress: {trial}/{num_samples}", flush=True)

        adj = random_tournament(n)
        c3 = find_3cycles(adj, n)
        t3 = len(c3)

        if t3 == 10:
            d = count_disjoint_pairs_triples(c3)
            if d == 0:
                k10_count += 1
                _, t5, t7, total_c, I2 = full_omega_analysis(adj, n)
                H = count_hamiltonian_paths(adj, n)
                I2_dist_k10[I2] += 1
                t5_dist_k10[t5] += 1
                H_dist_k10[H] += 1
                if k10_count <= 5:
                    print(f"  K_10: t5={t5}, t7={t7}, total={total_c}, "
                          f"I(Omega,2)={I2}, H={H}, match={I2==H}")

        elif t3 == 8:
            d = count_disjoint_pairs_triples(c3)
            if d == 1:
                k8e_count += 1
                _, t5, t7, total_c, I2 = full_omega_analysis(adj, n)
                H = count_hamiltonian_paths(adj, n)
                I2_dist_k8e[I2] += 1
                t5_dist_k8e[t5] += 1
                H_dist_k8e[H] += 1
                if k8e_count <= 5:
                    print(f"  K_8-e: t5={t5}, t7={t7}, total={total_c}, "
                          f"I(Omega,2)={I2}, H={H}, match={I2==H}")

    print()
    print(f"n=8 Results ({num_samples} samples):")
    print(f"  K_10 skeleton count: {k10_count}")
    if I2_dist_k10:
        print(f"    I(Omega,2) distribution: {dict(sorted(I2_dist_k10.items()))}")
        print(f"    t5 distribution: {dict(sorted(t5_dist_k10.items()))}")
        print(f"    H distribution: {dict(sorted(H_dist_k10.items()))}")
        min_I2 = min(I2_dist_k10.keys())
        print(f"    Min I(Omega,2) = {min_I2}  {'<= 21 POSSIBLE!' if min_I2 <= 21 else '> 21 IMPOSSIBLE'}")

    print(f"  K_8-e skeleton count: {k8e_count}")
    if I2_dist_k8e:
        print(f"    I(Omega,2) distribution: {dict(sorted(I2_dist_k8e.items()))}")
        print(f"    t5 distribution: {dict(sorted(t5_dist_k8e.items()))}")
        print(f"    H distribution: {dict(sorted(H_dist_k8e.items()))}")
        min_I2 = min(I2_dist_k8e.keys())
        print(f"    Min I(Omega,2) = {min_I2}  {'<= 21 POSSIBLE!' if min_I2 <= 21 else '> 21 IMPOSSIBLE'}")
    print()
    return k10_count, k8e_count

def run_verification():
    """Quick verification: H = I(Omega_full, 2) on small sample."""
    n = 5
    total = 1 << (n*(n-1)//2)
    print(f"=== Verification: H = I(Omega_full, 2) at n={n} ({total} tournaments) ===")
    mismatches = 0
    for idx in range(total):
        adj = tournament_from_index(n, idx)
        c3 = find_3cycles(adj, n)
        c5 = find_directed_5cycles_fast(adj, n)
        all_cycles = list(c3) + [c for c in c5]
        H = count_hamiltonian_paths(adj, n)
        I2 = build_full_omega_and_compute_I2(all_cycles)
        if H != I2:
            mismatches += 1
            if mismatches <= 3:
                print(f"  MISMATCH idx={idx}: H={H}, I(Omega,2)={I2}, t3={len(c3)}, t5={len(c5)}")
    print(f"  Checked {total}, mismatches: {mismatches}")
    if mismatches == 0:
        print("  PERFECT: H = I(Omega_full, 2) for all n=5 tournaments")
    print()

if __name__ == '__main__':
    # Step 1: Verify H = I(Omega_full, 2) at n=5
    run_verification()

    # Step 2: n=7 exhaustive
    k10_7, k8e_7 = run_n7_exhaustive()

    # Step 3: n=8 sampling
    k10_8, k8e_8 = run_n8_sampling(100000)

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"n=7 exhaustive:")
    print(f"  K_10 3-cycle skeleton: {k10_7} tournaments")
    print(f"  K_8-e 3-cycle skeleton: {k8e_7} tournaments")
    print(f"n=8 sampling (100k):")
    print(f"  K_10 3-cycle skeleton: {k10_8} found")
    print(f"  K_8-e 3-cycle skeleton: {k8e_8} found")
    print()
    print("Key question: Is min I(full Omega, 2) > 21 for all K_10/K_8-e skeleton tournaments?")
    print("If yes, these decompositions are impossible for H=21.")
