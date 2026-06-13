#!/usr/bin/env python3
"""
INV-011: Does alpha_1(T) ≡ c_3(T) (mod 2)?

alpha_1 = total number of directed odd cycles in T
c_3 = number of directed 3-cycles in T

At n=3,4: alpha_1 = c_3 (only 3-cycles exist), so trivially equal mod 2.
At n=5: 5-cycles first appear. Question: does c_5 always have even parity?
At n=6: both 5-cycles contribute. More complex.

Moon's formula: c_3 = C(n,3) - sum_v C(d+(v), 2) — determined by score sequence.
If alpha_1 ≡ c_3 (mod 2), then H(T) mod 4 is score-sequence-determined.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, random_tournament
from math import comb
from itertools import permutations
from collections import defaultdict

def count_3cycles(T):
    n = len(T)
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    c3 += 1
                elif T[i][k] and T[k][j] and T[j][i]:
                    c3 += 1
    return c3

def count_directed_5cycles(T):
    """Count directed 5-cycles."""
    n = len(T)
    count = 0
    for perm in permutations(range(n), 5):
        a, b, c, d, e = perm
        if T[a][b] and T[b][c] and T[c][d] and T[d][e] and T[e][a]:
            count += 1
    return count // 5  # Each directed 5-cycle counted 5 times (rotations)

def count_directed_7cycles(T):
    """Count directed 7-cycles (only for n>=7)."""
    n = len(T)
    if n < 7:
        return 0
    count = 0
    # Fix first vertex to remove rotational symmetry
    for perm in permutations(range(1, n), 6):
        path = [0] + list(perm)
        if len(path) == 7 and all(T[path[i]][path[(i+1) % 7]] for i in range(7)):
            count += 1
    return count  # Already counted once per directed cycle (first vertex fixed)

def score_sequence(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)]))

print("=" * 70)
print("INV-011: alpha_1 vs c_3 mod 2")
print("=" * 70)

for n in [3, 4, 5, 6, 7]:
    print(f"\n--- n = {n} ---")
    m = n * (n - 1) // 2

    if n <= 6:
        mode = 'exhaustive'
        total = 1 << m
    else:
        mode = 'random'
        total = 500

    mismatch_count = 0
    total_checked = 0
    c5_parity_stats = defaultdict(int)
    score_determines_alpha1_mod2 = True
    score_alpha1_mod2 = {}  # score_seq -> set of alpha1 mod 2 values

    for idx in range(total):
        if mode == 'exhaustive':
            T = tournament_from_bits(n, idx)
        else:
            T = random_tournament(n)

        c3 = count_3cycles(T)
        c5 = count_directed_5cycles(T)

        if n >= 7:
            c7 = count_directed_7cycles(T)
        else:
            c7 = 0

        alpha_1 = c3 + c5 + c7
        total_checked += 1

        # Check alpha_1 ≡ c_3 (mod 2)
        if (alpha_1 % 2) != (c3 % 2):
            mismatch_count += 1

        # Track c5 parity
        c5_parity_stats[c5 % 2] += 1

        # Track score sequence → alpha_1 mod 2
        ss = score_sequence(T)
        if ss not in score_alpha1_mod2:
            score_alpha1_mod2[ss] = set()
        score_alpha1_mod2[ss].add(alpha_1 % 2)

        if total_checked % 100000 == 0 and n <= 6:
            print(f"  ... {total_checked}/{total}")

    print(f"  Checked: {total_checked}")
    print(f"  alpha_1 ≡ c_3 (mod 2) mismatches: {mismatch_count}")

    if mismatch_count == 0:
        print(f"  CONJECTURE HOLDS: alpha_1 ≡ c_3 (mod 2) for all n={n} tournaments")
        print(f"  => c_5 + c_7 is always even!")
    else:
        pct = 100 * mismatch_count / total_checked
        print(f"  CONJECTURE FAILS: {pct:.1f}% mismatch rate")
        print(f"  c_5 parity: even={c5_parity_stats[0]}, odd={c5_parity_stats[1]}")

    # Check if score sequence determines alpha_1 mod 2
    ambiguous = [(ss, vals) for ss, vals in score_alpha1_mod2.items() if len(vals) > 1]
    if ambiguous:
        print(f"  Score sequence does NOT determine alpha_1 mod 2: {len(ambiguous)} ambiguous classes")
        for ss, vals in ambiguous[:3]:
            print(f"    Score {ss}: alpha_1 mod 2 can be {vals}")
    else:
        print(f"  Score sequence DETERMINES alpha_1 mod 2 ({len(score_alpha1_mod2)} classes)")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
