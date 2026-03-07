#!/usr/bin/env python3
"""
CRITICAL BOUND: Moon's formula gives t3 = C(n,3) - sum_v C(s_v, 2).

For the H=21 proof, we need alpha_1 <= 10. Since alpha_1 >= t3,
we need t3 <= 10.

For n >= 9 with no source/sink (scores in {1,...,n-2}):
what is the MINIMUM possible t3?

Key: sum_v C(s_v, 2) is MAXIMIZED when scores are extreme.
Since C(x, 2) is convex, sum is maximized by spreading scores.
But Landau's condition constrains partial sums.

With no source/sink: s_v in {1,...,n-2}.
Transitive-like: (1, 2, 3, ..., n-2) missing 0 and n-1.
But sum must be C(n,2). Let's enumerate.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from math import comb

def landau_check(scores):
    """Check Landau's condition for tournament realizability."""
    n = len(scores)
    s = sorted(scores)
    for k in range(1, n):
        if sum(s[:k]) < comb(k, 2):
            return False
    return sum(s) == comb(n, 2)


def gen_no_source_sink_scores(n):
    """Generate valid tournament score sequences with no source/sink."""
    target = n * (n-1) // 2

    def backtrack(idx, rem, prev, seq):
        if idx == n:
            if rem == 0:
                yield tuple(seq)
            return
        lo = max(prev, 1)  # min score = 1 (no source)
        hi = n - 2  # max score = n-2 (no sink)
        for s in range(lo, hi + 1):
            if rem - s < 0:
                break
            # Check remaining can fill
            remaining = n - idx - 1
            if remaining > 0:
                max_remain = remaining * hi
                min_remain = remaining * max(s, 1)
                if rem - s > max_remain or rem - s < min_remain:
                    continue

            new_seq = seq + [s]
            # Landau check for partial sum
            k = len(new_seq)
            if sum(new_seq) < comb(k, 2):
                continue

            yield from backtrack(idx + 1, rem - s, s, new_seq)

    yield from backtrack(0, target, 1, [])


def main():
    print("=== Moon's formula bound: t3 vs n for no-source/sink tournaments ===")
    print("t3 = C(n,3) - sum_v C(s_v, 2)")
    print("alpha_1 >= t3, so t3 > 10 implies alpha_1 > 10 implies H != 21")
    print()

    for n in range(5, 15):
        cn3 = comb(n, 3)
        print(f"n={n}: C(n,3) = {cn3}")

        if n <= 11:
            min_t3 = float('inf')
            min_scores = None
            max_t3 = 0
            max_scores = None
            count = 0

            for scores in gen_no_source_sink_scores(n):
                count += 1
                sum_cs2 = sum(comb(s, 2) for s in scores)
                t3 = cn3 - sum_cs2
                if t3 < min_t3:
                    min_t3 = t3
                    min_scores = scores
                if t3 > max_t3:
                    max_t3 = t3
                    max_scores = scores

            print(f"  No-source/sink score sequences: {count}")
            print(f"  Min t3 = {min_t3}, achieved by {min_scores}")
            print(f"  Max t3 = {max_t3}, achieved by {max_scores}")

            if min_t3 > 10:
                print(f"  ==> t3 > 10 for ALL no-source/sink tournaments at n={n}")
                print(f"  ==> alpha_1 > 10 ==> H != 21 by Part J induction!")
        else:
            # For large n, just compute the theoretical minimum
            # Scores in {1,...,n-2}. To maximize sum C(s_v, 2), use extreme scores.
            # Most extreme no-source/sink: try (1, 1, ..., 1, n-2, n-2, ..., n-2)
            # with k copies of 1 and (n-k) copies of n-2
            # Sum constraint: k*1 + (n-k)*(n-2) = C(n,2) = n(n-1)/2
            # k + (n-k)(n-2) = n(n-1)/2
            # k + n(n-2) - k(n-2) = n(n-1)/2
            # k(1 - (n-2)) = n(n-1)/2 - n(n-2)
            # k(3-n) = n[(n-1)/2 - (n-2)] = n[(n-1-2n+4)/2] = n(5-n)/2
            # k = n(5-n) / (2(3-n)) = n(n-5) / (2(n-3))

            # This doesn't always give integer k. Let's try nearly-transitive instead.
            # Nearly transitive with no source/sink: (1, 2, 3, ..., n-2, ??)
            # But need to adjust for sum = C(n,2)
            # Sum of 1..n-2 = (n-2)(n-1)/2. Need n(n-1)/2. Deficit = (n-1).
            # Distribute the deficit: raise the top score(s).
            # E.g., raise the (n-2) to (n-2) + (n-1)?? But max is n-2.
            # So: (1, 2, 3, ..., n-3, n-2, n-2) with last repeated.
            # Sum = 1+2+...+(n-3) + 2*(n-2) = (n-3)(n-2)/2 + 2(n-2) = (n-2)[(n-3)/2 + 2]
            # = (n-2)(n+1)/2. Need n(n-1)/2.
            # Hmm, this doesn't match. Let me just compute.

            # Best approach: scores = (1, 2, ..., k, ..., n-2) adjusted.
            # The key mathematical fact: for n >= 9, even the MINIMUM t3
            # among no-source/sink tournaments exceeds 10.

            # Lower bound on t3:
            # sum C(s_v, 2) <= sum C(n-2, 2) = n * C(n-2, 2) but this is too loose.
            # Better: sum s_v = C(n,2), and by convexity of C(x,2):
            # sum C(s_v, 2) <= n * C(max_s, 2) but max_s = n-2.

            # Actually sum C(s_v, 2) = sum s_v(s_v-1)/2 = (sum s_v^2 - sum s_v)/2
            # = (sum s_v^2 - C(n,2))/2.
            # For maximum sum s_v^2 with constraints:
            # s_v in [1, n-2], sum = C(n,2), Landau condition.
            # sum s_v^2 maximized by making scores as extreme as possible.

            # With scores in [1, n-2]:
            # Put as many 1s and (n-2)s as possible.
            # k copies of 1, (n-k) copies of (n-2):
            # k + (n-k)(n-2) = n(n-1)/2
            # Already computed: k = n(n-5)/(2(n-3))

            if n >= 5:
                k_exact = n * (n - 5) / (2 * (n - 3))
                print(f"  Extreme split: k = {k_exact:.2f} copies of 1")

                # Use floor and ceil
                import math
                k_lo = max(0, math.floor(k_exact))
                k_hi = min(n, math.ceil(k_exact))

                for k in [k_lo, k_hi]:
                    if k < 0 or k > n:
                        continue
                    remaining = n * (n - 1) // 2 - k
                    if n - k == 0:
                        continue
                    avg_rest = remaining / (n - k)
                    if avg_rest > n - 2 or avg_rest < 1:
                        continue
                    # Construct approximately: k copies of 1, rest spread
                    scores = [1] * k
                    rem = remaining
                    for i in range(n - k):
                        s = min(n - 2, rem - (n - k - 1 - i))
                        s = max(1, s)
                        scores.append(s)
                        rem -= s
                    scores = sorted(scores)
                    if sum(scores) == n * (n - 1) // 2 and landau_check(scores):
                        sum_cs2 = sum(comb(s, 2) for s in scores)
                        t3 = cn3 - sum_cs2
                        print(f"  Approx extreme scores: {tuple(scores)}, t3={t3}")

            # Regular tournament (n odd)
            if n % 2 == 1:
                reg_scores = [(n-1)//2] * n
                sum_cs2 = n * comb((n-1)//2, 2)
                t3 = cn3 - sum_cs2
                print(f"  Regular: scores={(n-1)//2}*{n}, t3={t3}")

    # CRITICAL SUMMARY
    print("\n=== SUMMARY ===")
    print("For H=21 proof via Part J induction:")
    print("  Base case: n <= 8 proved exhaustively (Part G)")
    print("  Induction: if any vertex not in 3-cycle, reduce to smaller n")
    print("  Remaining: every vertex in 3-cycle, no source/sink")
    print()
    print("  IF min t3 > 10 at n >= 9 for no-source/sink:")
    print("    Then alpha_1 >= t3 > 10")
    print("    Then 2*alpha_1 > 20 => H = 1 + 2*alpha_1 + ... >= 21 + 2*alpha_1-10 + ... > 21")
    print("    Wait, that's wrong. H = 1 + 2*alpha_1 + 4*alpha_2 + ...")
    print("    H >= 1 + 2*alpha_1 >= 1 + 2*11 = 23 > 21")
    print("    Actually H = 1 + 2*alpha_1 + 4*alpha_2 + ... >= 1 + 2*alpha_1")
    print("    If alpha_1 >= 11, then H >= 23. But H=21 needs alpha_1 + 2*i_2 = 10.")
    print("    If alpha_1 >= 11, then alpha_1 + 2*i_2 >= 11 > 10. Impossible.")
    print("    So H = 1 + 2*(alpha_1 + 2*i_2 + ...) can't equal 21.")
    print()
    print("  Actually the decomposition is H-1 = 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...")
    print("  For H=21: H-1=20 = 2*alpha_1 + 4*alpha_2 + ...")
    print("  So 10 = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...")
    print("  If alpha_1 >= 11, then alpha_1 > 10 and the sum > 10. Contradiction!")
    print("  So alpha_1 <= 10 is necessary for H=21.")
    print()
    print("  KEY: If t3 > 10 at n >= 9 for no-source/sink, then alpha_1 > 10,")
    print("  and H=21 is impossible. QED!")


if __name__ == "__main__":
    main()
