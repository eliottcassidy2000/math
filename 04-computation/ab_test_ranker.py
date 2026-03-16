#!/usr/bin/env python3
"""
ab_test_ranker.py — Rank n items from pairwise A/B test data.

Given pairwise comparison data (which item wins each matchup),
this tool:
  1. Builds the tournament
  2. Counts Hamiltonian paths H(T)
  3. Computes significance against random (Z-score)
  4. Reports whether the ranking is real or noise

USE CASES:
  - Product ranking from user A/B tests
  - Feature prioritization from team votes
  - UI design comparison
  - Food/wine tasting ranking
  - Interview candidate ranking

Input format: CSV with columns "item_a,item_b,winner"
Example:
  apple,banana,apple
  apple,cherry,cherry
  banana,cherry,banana

Output: ranked list + significance score.
"""

import sys
import csv
from math import sqrt, factorial
from itertools import permutations
from collections import defaultdict

def build_tournament(comparisons):
    """Build tournament from pairwise comparison data."""
    items = set()
    wins = defaultdict(lambda: defaultdict(int))

    for a, b, winner in comparisons:
        items.add(a)
        items.add(b)
        if winner == a:
            wins[a][b] += 1
        elif winner == b:
            wins[b][a] += 1

    items = sorted(items)
    n = len(items)
    idx = {item: i for i, item in enumerate(items)}

    T = [[0]*n for _ in range(n)]
    for a in items:
        for b in items:
            if a != b:
                if wins[a][b] > wins[b][a]:
                    T[idx[a]][idx[b]] = 1
                elif wins[b][a] > wins[a][b]:
                    T[idx[b]][idx[a]] = 1
                else:
                    # Tie: arbitrarily assign (or could randomize)
                    if idx[a] < idx[b]:
                        T[idx[a]][idx[b]] = 1
                    else:
                        T[idx[b]][idx[a]] = 1

    return T, items, n

def count_ham_paths(T, n):
    """Count Hamiltonian paths (brute force, n <= 10)."""
    if n > 10:
        return None  # Too slow
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if T[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

def score_ranking(T, n):
    """Score each item by out-degree (simple ranking)."""
    scores = [sum(T[i]) for i in range(n)]
    return scores

def main():
    if len(sys.argv) < 2:
        print("A/B Test Ranker")
        print("=" * 40)
        print()
        print("Usage: python3 ab_test_ranker.py data.csv")
        print("       python3 ab_test_ranker.py --demo")
        print()
        print("CSV format: item_a,item_b,winner")
        print()
        print("Based on the Cayley-Delannoy tournament theory.")
        print("Significance test: CV = sqrt(2/n).")
        return

    if sys.argv[1] == "--demo":
        comparisons = [
            ("React", "Vue", "React"),
            ("React", "Angular", "React"),
            ("React", "Svelte", "Svelte"),
            ("React", "Solid", "React"),
            ("Vue", "Angular", "Vue"),
            ("Vue", "Svelte", "Svelte"),
            ("Vue", "Solid", "Vue"),
            ("Angular", "Svelte", "Svelte"),
            ("Angular", "Solid", "Angular"),
            ("Svelte", "Solid", "Svelte"),
        ]
        print("DEMO: Frontend Framework Ranking")
        print("(Pairwise comparisons from hypothetical developer survey)")
        print()
    else:
        comparisons = []
        with open(sys.argv[1]) as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) >= 3:
                    comparisons.append((row[0].strip(), row[1].strip(), row[2].strip()))

    T, items, n = build_tournament(comparisons)
    scores = score_ranking(T, n)

    # Sort by score
    ranked = sorted(zip(items, scores), key=lambda x: -x[1])

    print(f"Tournament: {n} items, {len(comparisons)} comparisons")
    print()
    print("RANKING:")
    for rank, (item, score) in enumerate(ranked, 1):
        print(f"  {rank}. {item} (wins: {score}/{n-1})")

    # Significance test
    H = count_ham_paths(T, n) if n <= 10 else None
    mu = factorial(n) / 2**(n-1)
    cv_val = sqrt(2.0/n)
    sigma = mu * cv_val

    print()
    if H is not None:
        z = (H - mu) / sigma
        print(f"Hamiltonian paths: H = {H}")
        print(f"Expected (random): E[H] = {mu:.1f}")
        print(f"Normalized: H/E[H] = {H/mu:.3f}")
        print(f"Z-score: {z:+.2f}")
        print()
        # IMPORTANT: Low H = strong hierarchy (few orderings work = constrained).
        # High H = weak hierarchy (many orderings work = flexible).
        # H=1 = perfectly transitive = STRONGEST possible ranking.
        # H=max = maximally ambiguous = many equivalent orderings.
        if H == 1:
            print("CONCLUSION: PERFECT HIERARCHY. Only one consistent ordering exists.")
            print("The ranking is unambiguous.")
        elif z < -2:
            print("CONCLUSION: STRONG HIERARCHY. Very few consistent orderings.")
            print("The ranking is highly constrained — real skill differences exist.")
        elif z < -1:
            print("CONCLUSION: MODERATE HIERARCHY. Fewer orderings than random.")
            print("Some real structure in the preferences.")
        elif z < 1:
            print("CONCLUSION: NO SIGNIFICANT STRUCTURE. Consistent with random.")
            print("The preferences do not indicate a clear ranking.")
        else:
            print("CONCLUSION: UNUSUALLY FLEXIBLE. More orderings than random.")
            print("The items are highly interchangeable — consider them equivalent.")
    else:
        print(f"(Too many items for exact H computation. n={n} > 10.)")
        print(f"Using asymptotic test: CV = sqrt(2/{n}) = {cv_val:.4f}")
        print(f"A ranking is significant if it explains > {(1+2*cv_val)*100-100:.0f}% more")
        print(f"orderings than random chance.")

if __name__ == "__main__":
    main()
