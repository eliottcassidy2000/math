#!/usr/bin/env python3
"""
THEORETICAL PROOF ATTEMPT for H=21 impossibility at ALL n.

Key observations from data:
1. alpha_1=10 forces i_2 >= 2 (proved at n<=7, need general)
2. alpha_1=8 forces i_2 in {0} (not 1) at n<=7

For (10,0):
  alpha_1=10 with t3+t5+t7+... = 10.
  At n<=7: ONLY (3,3) disjoint pairs possible (l1+l2 > n for all other types).
  At n=8: also (3,5) possible.
  At n >= 10: also (3,7) and (5,5) possible.

  But the key constraint: with alpha_1=10, the tournament must have EXACTLY 10
  directed odd cycles. This is a very restrictive condition.

  Moon's formula: t3 = C(n,3) - sum_v C(s_v, 2).
  For t3 to be small enough that t3 + (rest) = 10, the tournament must be
  nearly transitive.

  A nearly-transitive tournament has very few cycles of any kind.

For (8,1):
  alpha_1=8 with i_2=1 means exactly 1 disjoint pair among 8 cycles.

APPROACH: Use the "few cycles" constraint to bound the structure.
If alpha_1=10, the tournament is nearly transitive. In a nearly-transitive
tournament, the cycles are clustered around a small number of "non-transitive"
arcs. This clustering forces many disjoint pairs OR forces all cycles to share.

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations
from collections import Counter, defaultdict
import random
import time


def find_directed_cycles_dp(adj, n, k):
    result = []
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        full = (1 << k) - 1
        for S in range(1, full + 1):
            for i in range(k):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(k):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, k):
            if (full, j) in dp and adj[v[j]][v[0]]:
                count += dp[(full, j)]
        if count > 0:
            result.append((frozenset(verts), count))
    return result


def get_all_cycles(adj, n):
    all_cycles = []
    max_k = n if n % 2 == 1 else n - 1
    for k in range(3, max_k + 1, 2):
        for vs, d in find_directed_cycles_dp(adj, n, k):
            for _ in range(d):
                all_cycles.append(vs)
    return all_cycles


def compute_i2(all_cycles):
    i2 = 0
    for a in range(len(all_cycles)):
        for b in range(a+1, len(all_cycles)):
            if not (all_cycles[a] & all_cycles[b]):
                i2 += 1
    return i2


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


def explore_low_alpha1(n, target_alpha1, samples=50000):
    """Find tournaments with exactly target_alpha1 directed odd cycles."""
    print(f"\nn={n}: Searching for alpha_1={target_alpha1}")

    random.seed(42 + n*100 + target_alpha1)
    results = []

    for trial in range(samples):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        all_cyc = get_all_cycles(adj, n)
        alpha1 = len(all_cyc)

        if alpha1 != target_alpha1:
            continue

        c_by_len = defaultdict(int)
        for c in all_cyc:
            c_by_len[len(c)] += 1

        i2 = compute_i2(all_cyc)
        H = held_karp(adj, n)
        scores = tuple(sorted([sum(adj[i]) for i in range(n)]))

        # Compute backward arcs (arcs going "wrong way" in score ordering)
        sorted_verts = sorted(range(n), key=lambda v: sum(adj[v]))
        back_arcs = 0
        for i in range(n):
            for j in range(i+1, n):
                vi, vj = sorted_verts[i], sorted_verts[j]
                if adj[vj][vi]:  # lower-score beats higher-score
                    back_arcs += 1

        results.append({
            'i2': i2, 'H': H, 'score': scores,
            'comp': dict(c_by_len), 'back_arcs': back_arcs,
        })

    if not results:
        print(f"  Not found in {samples} samples")
        return

    i2_dist = Counter(d['i2'] for d in results)
    H_dist = Counter(d['H'] for d in results)
    comp_dist = Counter(tuple(sorted(d['comp'].items())) for d in results)
    score_dist = Counter(d['score'] for d in results)

    print(f"  Found {len(results)}")
    print(f"  i_2: {dict(sorted(i2_dist.items()))}")
    print(f"  H: {dict(sorted(H_dist.items()))}")
    print(f"  Compositions: {dict(sorted(comp_dist.items()))}")
    if len(score_dist) <= 10:
        print(f"  Scores: {dict(sorted(score_dist.items()))}")

    # Check needed i_2 for H=21
    needed = (10 - target_alpha1) // 2 if (10 - target_alpha1) >= 0 and (10 - target_alpha1) % 2 == 0 else None
    if needed is not None:
        if needed in i2_dist:
            print(f"  ALERT: i_2={needed} exists! Check H values...")
            h_at_needed = [d['H'] for d in results if d['i2'] == needed]
            print(f"    H values at i_2={needed}: {dict(Counter(h_at_needed))}")
        else:
            print(f"  BLOCKED: need i_2={needed}, not achievable")


def main():
    # Check alpha_1=8 and alpha_1=10 at n=8,9
    for n in [8, 9]:
        for a1 in [8, 10]:
            explore_low_alpha1(n, a1, samples=100000)


if __name__ == "__main__":
    main()
