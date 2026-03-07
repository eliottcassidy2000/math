#!/usr/bin/env python3
"""
H=21 structural proof: analyze (8,1) and (10,0) cases.

For (10,0): alpha_1=10, need i_2=0 (all pairwise sharing).
  But computationally alpha_1=10 always gives i_2=2.
  WHY? This script investigates the structure.

For (8,1): alpha_1=8, need i_2=1.
  But computationally alpha_1=8 gives i_2 in {0,7} at n=8.
  WHY? This script investigates.

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
import sys
import time

def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1
        yield adj

def find_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

def find_5cycles(adj, n):
    cycles = []
    for verts in combinations(range(n), 5):
        v = list(verts)
        count = 0
        for perm in permutations(range(5)):
            ok = True
            for idx in range(5):
                if not adj[v[perm[idx]]][v[perm[(idx+1)%5]]]:
                    ok = False
                    break
            if ok:
                count += 1
        # Each directed 5-cycle counted 5 times (cyclic rotations)
        num_cycles = count // 5
        for _ in range(num_cycles):
            cycles.append(frozenset(verts))
    return cycles

def find_7cycles(adj, n):
    if n < 7:
        return []
    cycles = []
    for verts in combinations(range(n), 7):
        v = list(verts)
        count = 0
        # Use DP for efficiency
        dp = {}
        dp[(1, 0)] = 1
        for S in range(1, 128):
            for i in range(7):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(7):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        for j in range(1, 7):
            if (127, j) in dp and adj[v[j]][v[0]]:
                count += dp[(127, j)]
        num_cycles = count // 7
        for _ in range(num_cycles):
            cycles.append(frozenset(verts))
    return cycles

def compute_alpha_i2(cycles):
    alpha1 = len(cycles)
    i2 = 0
    for a in range(alpha1):
        for b in range(a+1, alpha1):
            if not (cycles[a] & cycles[b]):
                i2 += 1
    return alpha1, i2

def held_karp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full])

def analyze_alpha10_structure(adj, n, c3, c5, c7):
    """Deep structural analysis when alpha_1=10."""
    all_cycles = c3 + c5 + c7

    # Find the 2 disjoint pairs
    disjoint_pairs = []
    for a in range(len(all_cycles)):
        for b in range(a+1, len(all_cycles)):
            if not (all_cycles[a] & all_cycles[b]):
                la = len(all_cycles[a])
                lb = len(all_cycles[b])
                disjoint_pairs.append((la, lb, all_cycles[a], all_cycles[b]))

    return {
        't3': len(c3),
        't5': len(c5),
        't7': len(c7),
        'i2': len(disjoint_pairs),
        'pair_types': [(a,b) for a,b,_,_ in disjoint_pairs],
        'pair_details': disjoint_pairs,
    }

def analyze_alpha8_structure(adj, n, c3, c5, c7):
    """Deep structural analysis when alpha_1=8."""
    all_cycles = c3 + c5 + c7

    disjoint_pairs = []
    for a in range(len(all_cycles)):
        for b in range(a+1, len(all_cycles)):
            if not (all_cycles[a] & all_cycles[b]):
                la = len(all_cycles[a])
                lb = len(all_cycles[b])
                disjoint_pairs.append((la, lb))

    cycle_types = Counter(len(c) for c in all_cycles)

    return {
        'cycle_types': dict(cycle_types),
        'i2': len(disjoint_pairs),
        'pair_types': disjoint_pairs,
    }


def main():
    for n in [5, 6]:
        print(f"\n{'='*60}")
        print(f"n={n}: Structural analysis of alpha_1=8 and alpha_1=10")
        print(f"{'='*60}")

        alpha10_data = []
        alpha8_data = []
        alpha_i2_summary = defaultdict(lambda: defaultdict(int))

        t0 = time.time()
        count = 0
        for adj in all_tournaments(n):
            c3 = find_3cycles(adj, n)
            c5 = find_5cycles(adj, n)
            c7 = find_7cycles(adj, n) if n >= 7 else []
            all_cycles = c3 + c5 + c7
            alpha1 = len(all_cycles)
            _, i2 = compute_alpha_i2(all_cycles)
            H = held_karp(adj, n)

            alpha_i2_summary[alpha1][i2] += 1

            if alpha1 == 10:
                info = analyze_alpha10_structure(adj, n, c3, c5, c7)
                info['H'] = H
                alpha10_data.append(info)

            if alpha1 == 8:
                info = analyze_alpha8_structure(adj, n, c3, c5, c7)
                info['H'] = H
                alpha8_data.append(info)

            count += 1

        elapsed = time.time() - t0
        print(f"Checked {count} tournaments in {elapsed:.1f}s")

        # Print summary for all alpha_1 values
        print(f"\n--- (alpha_1, i_2) summary ---")
        for a1 in sorted(alpha_i2_summary.keys()):
            if a1 > 14:
                continue
            i2_dict = alpha_i2_summary[a1]
            i2_vals = sorted(i2_dict.keys())
            total = sum(i2_dict.values())
            print(f"  alpha_1={a1:2d}: i_2 in {{{','.join(str(v) for v in i2_vals)}}}, count={total}")

        # Detail for alpha_1=10
        if alpha10_data:
            print(f"\n--- alpha_1=10 detail ({len(alpha10_data)} tournaments) ---")
            type_dist = Counter((d['t3'], d['t5'], d['t7']) for d in alpha10_data)
            i2_dist = Counter(d['i2'] for d in alpha10_data)
            print(f"  i_2 distribution: {dict(sorted(i2_dist.items()))}")
            print(f"  (t3,t5,t7) distribution:")
            for k, v in sorted(type_dist.items()):
                print(f"    {k}: {v}")
            if alpha10_data:
                # Show structure of disjoint pairs
                pair_type_dist = Counter(tuple(sorted(d['pair_types'])) for d in alpha10_data)
                print(f"  Disjoint pair types (cycle lengths):")
                for k, v in sorted(pair_type_dist.items()):
                    print(f"    {k}: {v}")
        else:
            print(f"\n  No tournaments with alpha_1=10 at n={n}")

        # Detail for alpha_1=8
        if alpha8_data:
            print(f"\n--- alpha_1=8 detail ({len(alpha8_data)} tournaments) ---")
            i2_dist = Counter(d['i2'] for d in alpha8_data)
            print(f"  i_2 distribution: {dict(sorted(i2_dist.items()))}")
            type_dist = Counter(tuple(sorted(d['cycle_types'].items())) for d in alpha8_data)
            print(f"  Cycle type distribution:")
            for k, v in sorted(type_dist.items()):
                print(f"    {k}: {v}")
        else:
            print(f"\n  No tournaments with alpha_1=8 at n={n}")

    # At n=7, too slow for Python exhaustive.
    # Instead, sample and look for the structural pattern.
    n = 7
    print(f"\n{'='*60}")
    print(f"n={n}: Sampling alpha_1=8,10 tournaments (10k random)")
    print(f"{'='*60}")

    import random
    random.seed(42)

    alpha10_data = []
    alpha8_data = []
    alpha_i2_summary = defaultdict(lambda: defaultdict(int))

    t0 = time.time()
    for trial in range(10000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        c3 = find_3cycles(adj, n)
        c5 = find_5cycles(adj, n)
        c7 = find_7cycles(adj, n)
        all_cycles = c3 + c5 + c7
        alpha1 = len(all_cycles)
        _, i2 = compute_alpha_i2(all_cycles)

        alpha_i2_summary[alpha1][i2] += 1

        if alpha1 == 10:
            H = held_karp(adj, n)
            info = analyze_alpha10_structure(adj, n, c3, c5, c7)
            info['H'] = H
            alpha10_data.append(info)

        if alpha1 == 8:
            H = held_karp(adj, n)
            info = analyze_alpha8_structure(adj, n, c3, c5, c7)
            info['H'] = H
            alpha8_data.append(info)

    elapsed = time.time() - t0
    print(f"Sampled 10000 in {elapsed:.1f}s")

    for a1 in sorted(alpha_i2_summary.keys()):
        if a1 > 14:
            continue
        i2_dict = alpha_i2_summary[a1]
        i2_vals = sorted(i2_dict.keys())
        total = sum(i2_dict.values())
        print(f"  alpha_1={a1:2d}: i_2 in {{{','.join(str(v) for v in i2_vals)}}}, count={total}")

    if alpha10_data:
        print(f"\nalpha_1=10: {len(alpha10_data)} tournaments")
        i2_dist = Counter(d['i2'] for d in alpha10_data)
        print(f"  i_2 distribution: {dict(sorted(i2_dist.items()))}")
        type_dist = Counter((d['t3'], d['t5'], d['t7']) for d in alpha10_data)
        print(f"  (t3,t5,t7):")
        for k, v in sorted(type_dist.items()):
            print(f"    {k}: {v}")
        pair_type_dist = Counter(tuple(sorted(d['pair_types'])) for d in alpha10_data)
        print(f"  Disjoint pair types:")
        for k, v in sorted(pair_type_dist.items()):
            print(f"    {k}: {v}")

    if alpha8_data:
        print(f"\nalpha_1=8: {len(alpha8_data)} tournaments")
        i2_dist = Counter(d['i2'] for d in alpha8_data)
        print(f"  i_2 distribution: {dict(sorted(i2_dist.items()))}")
        type_dist = Counter(tuple(sorted(d['cycle_types'].items())) for d in alpha8_data)
        print(f"  Cycle type distribution:")
        for k, v in sorted(type_dist.items()):
            print(f"    {k}: {v}")

if __name__ == "__main__":
    main()
