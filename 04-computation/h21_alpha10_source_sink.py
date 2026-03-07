#!/usr/bin/env python3
"""
Does alpha_1=10 force a source or sink at n >= 8?

From Moon's formula: t3 = C(n,3) - sum_v C(s_v, 2)
For alpha_1 = 10: need t3 + t5 + t7 + ... = 10.

At n=8: C(8,3) = 56. So sum_v C(s_v, 2) = 56 - t3.
Since t3 <= 10, sum_v C(s_v, 2) >= 46.

Moon bound: sum C(s_v, 2) >= n * C((n-1)/2, 2) for n odd.
At n=8: sum >= 46 is very high. Let's check score constraints.

Key hypothesis: alpha_1=10 at n>=8 forces score sequence to have 0 or n-1.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from collections import Counter
from math import comb

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


def main():
    # ANALYTICAL approach: what score sequences allow t3 <= 10?
    print("=== Score sequences with t3 <= 10 ===")

    for n in [7, 8, 9, 10]:
        print(f"\nn={n}:")
        max_t3 = comb(n, 3)
        # t3 = C(n,3) - sum_v C(s_v, 2)
        # For t3 <= 10: sum_v C(s_v, 2) >= C(n,3) - 10
        threshold = max_t3 - 10
        print(f"  C(n,3) = {max_t3}")
        print(f"  Need sum C(s_v, 2) >= {threshold}")

        # Minimum sum C(s_v, 2) over all score sequences
        # Regular tournament (n odd): s_v = (n-1)/2 for all v
        # sum = n * C((n-1)/2, 2)
        if n % 2 == 1:
            reg_sum = n * comb((n-1)//2, 2)
            print(f"  Regular sum C(s_v, 2) = {reg_sum}")
            if reg_sum >= threshold:
                print(f"    -> Regular tournaments have t3 <= {max_t3 - reg_sum}")
            else:
                print(f"    -> Regular CAN have t3 up to {max_t3 - reg_sum}")

        # Near-transitive: scores (0, 1, 2, ..., n-1) gives
        # sum C(k, 2) for k=0,...,n-1 = sum_{k=0}^{n-1} k(k-1)/2
        trans_sum = sum(comb(k, 2) for k in range(n))
        print(f"  Transitive sum C(s_v, 2) = {trans_sum}")
        print(f"    -> Transitive t3 = {max_t3 - trans_sum}")

        # Find ALL score sequences with sum C(s_v, 2) >= threshold
        # Score sequence: s_0 <= s_1 <= ... <= s_{n-1}, sum = C(n,2), s_i in {0,...,n-1}
        # Landau's condition: sum_{i=0}^{k-1} s_i >= C(k,2) for all k

        def gen_scores(n):
            """Generate valid tournament score sequences."""
            target = n * (n-1) // 2

            def backtrack(idx, rem, prev, seq):
                if idx == n:
                    if rem == 0:
                        yield tuple(seq)
                    return
                lo = max(prev, 0)
                hi = min(n - 1, rem - (n - idx - 1) * prev if prev >= 0 else n - 1)
                for s in range(lo, n):
                    if rem - s < 0:
                        break
                    new_seq = seq + [s]
                    # Landau check: sum of first k scores >= C(k,2)
                    k = len(new_seq)
                    if sum(new_seq) < comb(k, 2):
                        continue
                    # Upper bound check
                    remaining = n - k
                    max_possible = sum(new_seq) + remaining * (n - 1)
                    if max_possible < target:
                        continue
                    yield from backtrack(idx + 1, rem - s, s, new_seq)

            yield from backtrack(0, target, 0, [])

        count_total = 0
        count_low_t3 = 0
        low_t3_scores = []

        for scores in gen_scores(n):
            count_total += 1
            s_sum = sum(comb(s, 2) for s in scores)
            t3 = max_t3 - s_sum
            if t3 <= 10:
                count_low_t3 += 1
                has_source_sink = (0 in scores) or ((n-1) in scores)
                low_t3_scores.append((scores, t3, has_source_sink))

        print(f"  Total score sequences: {count_total}")
        print(f"  With t3 <= 10: {count_low_t3}")

        if low_t3_scores:
            all_have_ss = all(ss for _, _, ss in low_t3_scores)
            print(f"  ALL have source/sink: {all_have_ss}")
            for s, t3, ss in low_t3_scores:
                ss_mark = "SS" if ss else "NO-SS"
                print(f"    {s}: t3={t3} [{ss_mark}]")

        # KEY CHECK: does t3 <= 10 force source/sink?
        no_ss = [s for s, t3, ss in low_t3_scores if not ss]
        if no_ss:
            print(f"  WARNING: {len(no_ss)} score(s) with t3<=10 and NO source/sink!")
            for s in no_ss:
                print(f"    {s}")
        else:
            print(f"  CONFIRMED: t3 <= 10 => source/sink at n={n}")


if __name__ == "__main__":
    main()
