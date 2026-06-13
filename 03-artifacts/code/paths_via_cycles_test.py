#!/usr/bin/env python3
"""
Test the conjecture:
  #{Ham paths using arc i->j in T} = 2 * sum_{odd C using i->j} H(T[V\V(C)])

with convention H(empty tournament) = 1.

If true, this would be a fundamental identity connecting Ham paths to
cycle-weighted complement counts. Combined with the algebraic formula
for delta_I (verified), it would prove OCF = Claim A by induction.

Instance: kind-pasteur-2026-03-05-S5
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
import random


def H_complement(T, verts_to_remove):
    """H(T[V \ verts_to_remove]) with convention H(empty) = 1."""
    n = len(T)
    keep = [v for v in range(n) if v not in verts_to_remove]
    if len(keep) == 0:
        return 1  # convention: empty tournament has 1 Ham path
    # Build subtournament
    k = len(keep)
    sub = [[0]*k for _ in range(k)]
    for a in range(k):
        for b in range(k):
            sub[a][b] = T[keep[a]][keep[b]]
    return hamiltonian_path_count(sub)


def paths_using_arc(T, i, j):
    """Count Ham paths in T that use arc i->j (i immediately before j)."""
    n = len(T)
    if not T[i][j]:
        return 0

    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    # dp_used tracks paths that have used i->j
    dp_used = [[0]*n for _ in range(1 << n)]

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[mask][v]
            cu = dp_used[mask][v]
            if c == 0 and cu == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    nmask = mask | (1 << u)
                    dp[nmask][u] += c
                    if v == i and u == j:
                        dp_used[nmask][u] += cu + c
                    else:
                        dp_used[nmask][u] += cu

    return sum(dp_used[full][v] for v in range(n))


def cycles_using_arc(T, i, j):
    """Find all odd cycles in T that use arc i->j."""
    all_cyc = find_odd_cycles(T)
    result = []
    for c in all_cyc:
        L = len(c)
        for k in range(L):
            if c[k] == i and c[(k+1) % L] == j:
                result.append(c)
                break
    return result


def test_conjecture(n, exhaustive=True, n_random=200):
    """Test #{paths using i->j} = 2 * sum_{C using i->j} H(complement)."""
    print(f"\n=== n={n} {'(exhaustive)' if exhaustive else f'(random {n_random})'} ===")

    rng = random.Random(42)
    total = 0
    ok = 0
    fail = 0
    fail_examples = []

    if exhaustive:
        tournaments = list(all_tournaments(n))
    else:
        tournaments = [random_tournament(n, rng) for _ in range(n_random)]

    for T in tournaments:
        for i in range(n):
            for j in range(n):
                if i == j or not T[i][j]:
                    continue

                p = paths_using_arc(T, i, j)
                cycs = cycles_using_arc(T, i, j)

                rhs = 0
                for c in cycs:
                    h_comp = H_complement(T, set(c))
                    rhs += h_comp

                rhs *= 2

                total += 1
                if p == rhs:
                    ok += 1
                else:
                    fail += 1
                    if len(fail_examples) < 5:
                        fail_examples.append({
                            'i': i, 'j': j,
                            'paths': p, 'rhs': rhs,
                            'n_cycles': len(cycs),
                            'cycle_lengths': [len(c) for c in cycs],
                        })

    print(f"Total (T, i->j) pairs: {total}")
    print(f"  Match: {ok}/{total}")
    print(f"  Fail: {fail}/{total}")

    if fail_examples:
        print("  Failure examples:")
        for ex in fail_examples:
            print(f"    arc {ex['i']}->{ex['j']}: paths={ex['paths']}, "
                  f"2*sum_H={ex['rhs']}, #cycles={ex['n_cycles']}, "
                  f"lengths={ex['cycle_lengths']}")

    return fail == 0


def main():
    print("=== Testing: #{paths using i->j} = 2 * sum_{C using i->j} H(complement) ===")

    for n in range(3, 6):
        test_conjecture(n, exhaustive=True)

    test_conjecture(6, exhaustive=False, n_random=100)

    print("\n=== Summary ===")
    print("If the conjecture holds, it directly implies OCF/Claim A by induction!")


if __name__ == "__main__":
    main()
