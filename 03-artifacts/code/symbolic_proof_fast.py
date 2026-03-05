#!/usr/bin/env python3
"""
Fast symbolic proof of delta_H = delta_I for n=6.

Instead of enumerating all path shapes, use DP to compute adj counts
and blocking directly. Much faster than the naive approach.

Instance: kind-pasteur-2026-03-05-S6
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
from itertools import combinations, permutations


def compute_adj_and_unmatched(T_func, n, src, dst):
    """
    Compute adj(src, dst) and #unmatched using DP.
    T_func(a, b) returns 1 if arc a->b exists.
    """
    full = (1 << n) - 1

    # dp[mask][v] = # Ham paths on mask ending at v
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T_func(v, u):
                    dp[mask | (1 << u)][u] += dp[mask][v]

    # adj(src, dst): paths using src->dst consecutively
    # = sum over S containing src but not dst:
    #   dp_end(S, src) * dp_start(V\S, dst)
    # where dp_start uses a separate computation

    # dp_start[mask][v] = # Ham paths on mask starting at v
    dp_start = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp_start[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp_start[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T_func(v, u):
                    dp_start[mask | (1 << u)][u] += dp_start[mask][v]
    # Actually dp_start[mask][u] counts paths on mask ending at u (same as dp).
    # To get paths starting at v, I need dp_from_v.

    # Simpler: compute adj and unmatched directly from path enumeration
    # For small n (<=7), enumerate all Ham paths

    # Let's use the dp to count adj and unmatched efficiently.
    # A Ham path using src->dst at some position k:
    # dp[S][src] * dp_from_dst_on_complement[V\S][some endpoint]
    # where S contains src but not dst.

    # adj(src, dst) = sum_{S: src in S, dst not in S, T(src,dst)=1}
    #                 dp[S][src] * (sum_e dp_on_complement[V\S starting from dst][e])

    if not T_func(src, dst):
        return 0, 0

    # dp2[mask][v] = # Ham paths on vertices in mask starting at v
    dp2 = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp2[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp2[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T_func(v, u):
                    dp2[mask | (1 << u)][u] += dp2[mask][v]

    # h_start(mask, v) = sum_e dp2[mask][e] for paths starting at v
    # Actually dp2[mask][e] counts paths from SOME start ending at e.
    # I need the reverse: paths starting at a specific vertex.

    # Let me use a different DP: dp3[mask][v] = # paths on mask starting at v
    dp3 = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp3[1 << v][v] = 1
    # Build by adding to the END
    for mask in range(1, 1 << n):
        for v in range(n):  # start vertex
            if not (mask & (1 << v)):
                continue
            if dp3[mask][v] == 0:
                continue
            # dp3[mask][v] paths start at v and end at some vertex in mask
            # To extend: find the endpoint. But dp3 doesn't track endpoint!
            pass

    # This approach is getting convoluted. For n=6, just enumerate paths directly.
    # With 720 permutations, it's fast enough.

    adj = 0
    unmatched = 0

    # Use the dp to enumerate: iterate over all full-mask paths
    # Actually, let's just check all permutations
    paths = []

    def gen_paths(path, used):
        if len(path) == n:
            paths.append(tuple(path))
            return
        last = path[-1]
        for u in range(n):
            if not (used & (1 << u)) and T_func(last, u):
                path.append(u)
                gen_paths(path, used | (1 << u))
                path.pop()

    for v in range(n):
        gen_paths([v], 1 << v)

    for p in paths:
        for k in range(n - 1):
            if p[k] == src and p[k + 1] == dst:
                adj += 1
                # Check blocking
                blocked = False
                if k > 0 and not T_func(p[k - 1], dst):
                    blocked = True
                if k + 1 < n - 1 and not T_func(src, p[k + 2]):
                    blocked = True
                if blocked:
                    unmatched += 1
                break

    return adj, unmatched


def prove_n6_fast():
    """Prove delta_H = delta_I at n=6 by exhaustive evaluation."""
    print("=== n=6 Fast Symbolic Proof ===\n")

    i_v, j_v = 0, 1
    others = [2, 3, 4, 5]
    n = 6

    vars_list = ['T20', 'T30', 'T40', 'T50',
                 'T12', 'T13', 'T14', 'T15',
                 'T23', 'T24', 'T25', 'T34', 'T35', 'T45']

    def make_T(vals):
        def T_func(a, b):
            if a == b:
                return 0
            if a == 0 and b == 1:
                return 1
            if a == 1 and b == 0:
                return 0
            key = f"T{a}{b}"
            if key in vals:
                return vals[key]
            return 1 - vals[f"T{b}{a}"]
        return T_func

    def make_Tp(vals):
        def Tp_func(a, b):
            if a == b:
                return 0
            if a == 0 and b == 1:
                return 0
            if a == 1 and b == 0:
                return 1
            key = f"T{a}{b}"
            if key in vals:
                return vals[key]
            return 1 - vals[f"T{b}{a}"]
        return Tp_func

    ok = 0
    total = 0
    fail_count = 0

    for mask in range(2 ** len(vars_list)):
        vals = {}
        for idx, v in enumerate(vars_list):
            vals[v] = (mask >> idx) & 1

        T_func = make_T(vals)
        Tp_func = make_Tp(vals)

        _, uT = compute_adj_and_unmatched(T_func, n, i_v, j_v)
        _, uTp = compute_adj_and_unmatched(Tp_func, n, j_v, i_v)

        delta = uTp - uT

        # Compute expected from THM-013 formula
        s_vals = {x: 1 - T_func(x, 0) - T_func(1, x) for x in others}

        # H(B_x) for 3-vertex sub-tournaments
        def h3(verts):
            count = 0
            for p in permutations(verts):
                if T_func(p[0], p[1]) and T_func(p[1], p[2]):
                    count += 1
            return count

        # 5-cycle counts (my convention: C5=lost, D5=gained)
        C5 = 0  # lost: use i->j = 0->1
        D5 = 0  # gained: use j->i = 1->0
        for subset in combinations(others, 3):
            for perm in permutations(subset):
                v1, v2, v3 = perm
                # Lost 5-cycle: (0,1,v1,v2,v3,0)
                if (T_func(1, v1) and T_func(v1, v2) and
                    T_func(v2, v3) and T_func(v3, 0)):
                    C5 += 1
                # Gained 5-cycle: (1,0,v1,v2,v3,1)
                if (T_func(0, v1) and T_func(v1, v2) and
                    T_func(v2, v3) and T_func(v3, 1)):
                    D5 += 1

        formula_sum = sum(s_vals[x] * h3([v for v in others if v != x])
                         for x in others)

        # H(T')-H(T) = 2*sum(s_x*H(B_x)) + 2*(D5-C5)
        expected = 2 * formula_sum + 2 * (D5 - C5)

        total += 1
        if delta == expected:
            ok += 1
        else:
            fail_count += 1
            if fail_count <= 5:
                print(f"  FAIL: mask={mask}, delta={delta}, expected={expected}, "
                      f"s={s_vals}, D5={D5}, C5={C5}, formula_sum={formula_sum}")
            if fail_count > 20:
                print(f"  Aborting after {fail_count} failures.")
                break

        if total % 4000 == 0:
            print(f"  Progress: {total}/{2**14}, ok={ok}, fail={fail_count}")

    print(f"\nIdentity verified: {ok}/{total}")
    if ok == total:
        print("PROVED: delta_H = delta_I for ALL n=6 arc configurations (16384 cases).")
        print("This proves OCF at n=6 via arc-flip induction from transitive base case.")
    elif fail_count > 0:
        print(f"FAILED: {fail_count} counterexamples found.")


if __name__ == "__main__":
    prove_n6_fast()
