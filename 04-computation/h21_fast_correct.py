#!/usr/bin/env python3
"""
Fast and CORRECT analysis of (alpha_1, i_2) for H=21 investigation.

The DP-based 5-cycle counter had a bug: it used adj[v[i]][v[j]] which
is correct, but the issue was in the PREVIOUS script's find_5cycles_dp
which counted cycles per VERTEX SET rather than per DIRECTED CYCLE.

Actually, let me first verify: does the DP method give the RIGHT count
of directed cycles? Let me compare DP vs brute-force at n=5,6,7.

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
import random
import time


def find_directed_cycles_dp(adj, n, k):
    """Find all directed k-cycles using DP.
    Returns list of (vertex_set, num_directed_cycles) pairs.
    The DP counts Hamiltonian cycles starting at v[0] in each k-subset.
    Each directed cycle is counted k times (once per starting vertex),
    but since we fix the start at v[0], each cycle is counted once.
    Then we need to divide by... nothing actually, since we fix one vertex.
    Wait: a directed cycle visits all k vertices. If we fix v[0] as start,
    each directed cycle passes through v[0] exactly once, so it's counted once.
    But there are k rotations of the same cycle. If we fix v[0] as start vertex,
    we count each directed cycle that goes v[0] -> ... -> ... -> v[0] exactly once.
    But actually, the set of k vertices has k! permutations.
    A directed k-cycle corresponds to (k-1)! / 1 = ... hmm.

    Let me think: a directed k-cycle on {v[0],...,v[k-1]} is a cyclic sequence
    v[pi(0)] -> v[pi(1)] -> ... -> v[pi(k-1)] -> v[pi(0)].
    Two cyclic sequences are the same cycle if one is a rotation of the other.
    So there are k rotations per cycle.
    If we fix v[0] as the START, we get (k-1)! arrangements, each corresponding
    to exactly one directed Hamiltonian cycle. No need to divide.
    BUT: we also need the return edge v[last] -> v[0].

    The DP computes: number of Hamiltonian PATHS from v[0] through all vertices,
    ending at some v[j], where also adj[v[j]][v[0]] = 1.
    This is the number of directed Hamiltonian cycles through v[0].
    Each directed cycle is counted exactly once (with v[0] as start).
    But wait: v[0] is the FIRST vertex in the sorted order. So we're counting
    cycles where v[0] is the starting point. Each cycle has exactly one rotation
    starting at v[0]. So the count is correct: num_directed_cycles.

    CORRECTION: Actually NO. The "starting vertex" is fixed. But there's no
    concept of "start" for a cycle. A directed cycle v0->v1->v2->...->vk-1->v0
    is the SAME cycle as v1->v2->...->vk-1->v0->v1. We want to count
    DISTINCT directed cycles, so we should NOT divide by k because we already
    fixed v[0] as starting point.
    """
    result = []
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1  # start at v[0] (bitmask has bit 0 set)
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
        # Count paths that return to v[0]
        count = 0
        for j in range(1, k):
            if (full, j) in dp and adj[v[j]][v[0]]:
                count += dp[(full, j)]
        # count = number of directed Hamiltonian cycles through v[0]
        # But these are ORDERED cycles starting at v[0].
        # Each DIRECTED cycle has k rotations, but we fixed start = v[0],
        # so each directed cycle is counted exactly once.
        # HOWEVER: the DP counts Hamiltonian PATHS v[0] -> ... -> v[j],
        # and we check adj[v[j]][v[0]]. This gives the number of
        # directed Hamiltonian cycles containing vertex v[0].
        # Since we want TOTAL directed k-cycles (not just through v[0]):
        # each cycle contains v[0] iff v[0] is in the vertex set, which it always is.
        # So count IS the total number of directed k-cycles on this vertex set.
        if count > 0:
            result.append((frozenset(verts), count))
    return result


def find_directed_cycles_brute(adj, n, k):
    """Brute-force: count directed k-cycles on each k-subset."""
    result = []
    for verts in combinations(range(n), k):
        v = list(verts)
        count = 0
        for perm in permutations(range(k)):
            ok = True
            for idx in range(k):
                if not adj[v[perm[idx]]][v[perm[(idx+1)%k]]]:
                    ok = False
                    break
            if ok:
                count += 1
        # Each directed cycle counted k times (k rotations)
        d = count // k
        if d > 0:
            result.append((frozenset(verts), d))
    return result


def verify_dp_vs_brute():
    """Verify DP gives same results as brute-force."""
    print("=== VERIFICATION: DP vs brute-force cycle counting ===")

    for n in [5, 6, 7]:
        random.seed(42)
        mismatches = 0
        total = 0

        for trial in range(min(200, 2**(n*(n-1)//2))):
            adj = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1

            for k in [3, 5]:
                if k > n:
                    continue
                dp_result = dict(find_directed_cycles_dp(adj, n, k))
                bf_result = dict(find_directed_cycles_brute(adj, n, k))

                if dp_result != bf_result:
                    mismatches += 1
                    if mismatches <= 3:
                        print(f"  MISMATCH at n={n}, k={k}, trial={trial}")
                        for vs in set(list(dp_result.keys()) + list(bf_result.keys())):
                            d = dp_result.get(vs, 0)
                            b = bf_result.get(vs, 0)
                            if d != b:
                                print(f"    {sorted(vs)}: DP={d}, BF={b}")
                total += 1

        print(f"  n={n}: {total} checks, {mismatches} mismatches")


def main():
    verify_dp_vs_brute()

    # Now run the corrected analysis
    print("\n=== n=6 EXHAUSTIVE: (alpha_1, i_2) with correct counting ===")
    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]

    alpha_i2 = defaultdict(lambda: defaultdict(int))
    h21_count = 0

    for bits in range(2**15):
        adj = [[0]*n for _ in range(n)]
        for k_idx, (i, j) in enumerate(edges):
            if (bits >> k_idx) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1

        # All odd cycles with correct counting
        all_cycles = []
        for vs, d in find_directed_cycles_dp(adj, n, 3):
            for _ in range(d):
                all_cycles.append(vs)
        for vs, d in find_directed_cycles_dp(adj, n, 5):
            for _ in range(d):
                all_cycles.append(vs)

        alpha1 = len(all_cycles)
        i2 = 0
        for a in range(alpha1):
            for b in range(a+1, alpha1):
                if not (all_cycles[a] & all_cycles[b]):
                    i2 += 1

        alpha_i2[alpha1][i2] += 1

        # Check H=21
        # H = 1 + 2*alpha_1 + 4*i2 + ... (need full IP)
        # Quick: held_karp
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
        H = sum(dp[full])
        if H == 21:
            h21_count += 1

    print(f"  H=21 count at n=6: {h21_count}")
    for a1 in sorted(alpha_i2.keys()):
        if a1 > 20:
            continue
        i2_dict = alpha_i2[a1]
        i2_vals = sorted(i2_dict.keys())
        total = sum(i2_dict.values())
        needed = (10 - a1) // 2 if 0 <= 10 - a1 and (10 - a1) % 2 == 0 else None
        tag = ""
        if needed is not None and needed >= 0:
            if needed in i2_vals:
                tag = f" [i_2={needed} EXISTS]"
            else:
                tag = f" [BLOCKED: need i_2={needed}]"
        print(f"  alpha_1={a1:2d}: i_2 in {{{','.join(str(v) for v in i2_vals)}}}, count={total}{tag}")

    # n=7 sampling
    print("\n=== n=7 SAMPLING (10k) ===")
    n = 7
    random.seed(42)
    alpha_i2 = defaultdict(lambda: defaultdict(int))
    h21_count = 0

    t0 = time.time()
    for trial in range(10000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        all_cycles = []
        for vs, d in find_directed_cycles_dp(adj, n, 3):
            for _ in range(d):
                all_cycles.append(vs)
        for vs, d in find_directed_cycles_dp(adj, n, 5):
            for _ in range(d):
                all_cycles.append(vs)
        for vs, d in find_directed_cycles_dp(adj, n, 7):
            for _ in range(d):
                all_cycles.append(vs)

        alpha1 = len(all_cycles)
        if alpha1 > 25:
            continue

        i2 = 0
        for a in range(alpha1):
            for b in range(a+1, alpha1):
                if not (all_cycles[a] & all_cycles[b]):
                    i2 += 1

        alpha_i2[alpha1][i2] += 1

        # H
        dp_h = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp_h[1 << v][v] = 1
        full = (1 << n) - 1
        for S in range(1, full + 1):
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                if dp_h[S][v] == 0:
                    continue
                for u in range(n):
                    if S & (1 << u):
                        continue
                    if adj[v][u]:
                        dp_h[S | (1 << u)][u] += dp_h[S][v]
        H = sum(dp_h[full])
        if H == 21:
            h21_count += 1

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s, H=21 count: {h21_count}")

    for a1 in sorted(alpha_i2.keys()):
        if a1 < 4 or a1 > 16:
            continue
        i2_dict = alpha_i2[a1]
        i2_vals = sorted(i2_dict.keys())
        total = sum(i2_dict.values())
        needed = (10 - a1) // 2 if 0 <= 10 - a1 and (10 - a1) % 2 == 0 else None
        tag = ""
        if needed is not None and needed >= 0:
            if needed in i2_vals:
                tag = f" [i_2={needed} EXISTS]"
            else:
                tag = f" [BLOCKED: need i_2={needed}]"
        print(f"  alpha_1={a1:2d}: i_2 in {{{','.join(str(v) for v in i2_vals)}}}, count={total}{tag}")


if __name__ == "__main__":
    main()
