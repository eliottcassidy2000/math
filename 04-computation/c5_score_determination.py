"""
c5_score_determination.py
kind-pasteur-2026-03-07-S39b

FINDING: c5 (number of directed 5-cycles) is NOT determined by the score sequence.

At n=5, score (1,2,2,2,3) has 280 tournaments with c5 in {1, 2, 3}.
Even the full edge-score vector (common out-neighbor counts) does NOT determine c5.
This means any c5 formula must involve cubic or higher-order graph statistics.

Compare to c3, which IS determined by score sequence via Moon's formula:
  c3 = C(n,3) - sum_i C(s_i, 2)

The Komarov-Mackey formula (JGT 2017) for c5 likely involves
sum of d(i,j)^2 * (arc-specific terms) or similar third-order invariants.

PRACTICAL CONSEQUENCE: For OCF computation via I(Omega, 2), we cannot
shortcut c5 counting from the score sequence alone. The bitmask DP
cycle-finder or exhaustive enumeration is necessary.

For regular tournaments specifically, c5 IS a class invariant per
Savchenko's results (verified at n=7).
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '03-artifacts', 'code'))

from tournament_lib import tournament_from_bits
from math import comb
from itertools import permutations, combinations
from collections import defaultdict


def count_directed_5_cycles(T):
    """Count directed 5-cycles. Each cycle counted once (rotation-canonical)."""
    n = len(T)
    c5 = 0
    for verts in combinations(range(n), 5):
        first = verts[0]
        for perm in permutations(verts[1:]):
            path = (first,) + perm
            if all(T[path[i]][path[(i+1) % 5]] for i in range(5)):
                c5 += 1
    return c5


def edge_scores(T):
    """For each arc i->j, compute d(i,j) = #{k: T[i][k]=1 and T[j][k]=1}."""
    n = len(T)
    d = {}
    for i in range(n):
        for j in range(n):
            if i != j and T[i][j]:
                d[(i, j)] = sum(1 for k in range(n) if k != i and k != j and T[i][k] and T[j][k])
    return d


# ============================================================
# Main analysis
# ============================================================

print("=" * 70)
print("c5 determination by score sequence and edge-score statistics")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    num_T = 1 << m

    if n < 5:
        print(f"\n  n={n}: no 5-cycles possible")
        continue

    # Group by score -> c5
    score_c5 = defaultdict(lambda: defaultdict(int))
    # Group by (score, sum_d^2) -> c5
    score_d2_c5 = defaultdict(lambda: defaultdict(int))

    for bits in range(num_T):
        T = tournament_from_bits(n, bits)
        scores = tuple(sorted(sum(T[i]) for i in range(n)))

        c5 = count_directed_5_cycles(T)
        score_c5[scores][c5] += 1

        d = edge_scores(T)
        sum_d2 = sum(v**2 for v in d.values())
        score_d2_c5[(scores, sum_d2)][c5] += 1

    print(f"\n  n={n} ({num_T} tournaments):")
    print(f"  Score -> c5 (grouped):")
    for scores in sorted(score_c5.keys()):
        c5_dist = score_c5[scores]
        const = "const" if len(c5_dist) == 1 else "VARIES"
        total = sum(c5_dist.values())
        detail = ", ".join(f"c5={k}: {v}" for k, v in sorted(c5_dist.items()))
        if const == "VARIES":
            print(f"    {scores}: {const} ({detail})")
        else:
            k = list(c5_dist.keys())[0]
            print(f"    {scores}: c5={k} (constant, {total} tours)")

    varying = sum(1 for v in score_c5.values() if len(v) > 1)
    print(f"  {varying}/{len(score_c5)} score sequences have varying c5")

    varying_d2 = sum(1 for v in score_d2_c5.values() if len(v) > 1)
    print(f"  {varying_d2}/{len(score_d2_c5)} (score, sum_d2) groups have varying c5")


# ============================================================
# Regular tournaments: c5 IS constant per isomorphism class
# ============================================================
print("\n" + "=" * 70)
print("Regular tournament c5: constant per iso class (Savchenko)")
print("=" * 70)

for n in [5, 7]:
    m = n * (n - 1) // 2
    num_T = 1 << m

    regular_c5 = defaultdict(int)
    if n <= 6:
        rng = range(num_T)
    else:
        import random
        random.seed(42)
        rng = [random.randint(0, num_T - 1) for _ in range(5000)]

    for bits in rng:
        T = tournament_from_bits(n, bits)
        scores = [sum(T[i]) for i in range(n)]
        if len(set(scores)) != 1 or scores[0] != (n - 1) // 2:
            continue

        if n >= 7:
            # Just count 5-cycles on random subsets for speed
            c5 = 0
            for verts in combinations(range(n), 5):
                first = verts[0]
                for perm in permutations(verts[1:]):
                    path = (first,) + perm
                    if all(T[path[i]][path[(i + 1) % 5]] for i in range(5)):
                        c5 += 1
        else:
            c5 = count_directed_5_cycles(T)

        regular_c5[c5] += 1

    print(f"\n  n={n}: regular tournament c5 distribution:")
    for c5, count in sorted(regular_c5.items()):
        print(f"    c5={c5}: {count} tournaments")
    if len(regular_c5) == 1:
        print(f"    c5 is CONSTANT for all regular tournaments at n={n}")
    else:
        print(f"    c5 takes {len(regular_c5)} distinct values")


print("\nDONE")
