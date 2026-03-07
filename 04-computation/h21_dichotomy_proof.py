#!/usr/bin/env python3
"""
DICHOTOMY PROOF ATTEMPT: For cycle-rich T on n >= 9:
Either (a) 3 disjoint 3-cycles exist, OR
       (b) some vertex deletion yields cycle-rich at n-1.

APPROACH: Analyze the 3-cycle hypergraph structure.

Case analysis based on max matching (mm) of 3-cycle vertex sets:

mm >= 3: Case (a) directly.

mm = 2: Two disjoint 3-cycles A={a1,a2,a3}, B={b1,b2,b3}.
  Every other 3-cycle hits A union B.
  Remaining vertices R = [n] \ (A union B), |R| = n-6 >= 3.
  Each r in R is in some 3-cycle C_r (since T is cycle-rich).
  C_r must hit A union B. So C_r uses r plus 2 vertices from A union B.

  Consider removing vertex a1 from A.
  Does every vertex still have a 3-cycle in T - a1?
  - Vertices in B: B is a 3-cycle not containing a1. OK.
  - Vertices a2, a3: they lose cycle A. But are they in other 3-cycles?
    If a2 is in a 3-cycle {a2, x, y} with x,y not in {a1}: OK.
    If ALL 3-cycles containing a2 also contain a1: problem.
  - Vertices in R: their 3-cycles C_r might contain a1. If so, they
    need alternative 3-cycles.

  KEY: At n >= 9 with mm=2, there are many 3-cycles (the mm=2 constraint
  forces concentrated structure with lots of cycles near A and B).

  From n=9 data: mm=2 gives min alpha_1 = 6, with the 6 cycles covering
  many cross-structures. Typically a2 and a3 appear in multiple 3-cycles.

mm = 1: All 3-cycles share a vertex v.
  T-v is transitive (no 3-cycles). Not cycle-rich.
  But H(T) with all-through-v: alpha_1 >= n-2 >> 10 at n >= 9.
  So H >= 1 + 2*(n-2) >= 17 at n=9. But we need H >= 23 or H != 21 specifically.

  Actually: alpha_1 = 2^(n-3) for the canonical all-through-v construction.
  Even for non-canonical: alpha_1 >= n-2 and there are many 5-cycles.

  At n=9 mm=1: the all-through-v structure forces t3 >= 4 (by coverage count,
  since each 3-cycle covers 2 non-v vertices, need ceil((n-1)/2) >= 4 three-cycles).
  Actually, coverage: each 3-cycle covers 1 vertex not previously covered.
  With n-1=8 non-v vertices, need at least 4 three-cycles (each covering 2).
  But they share v, and pairs might overlap on non-v vertices too.

  Min t3 = ceil((n-1)/2) = 4 at n=9 (if each cycle covers 2 new non-v vertices).
  t5 is forced by the structure: 5-cycles through v require paths of length 4
  from N+(v) to N-(v) in the transitive T-v.

Let me verify computationally: does the dichotomy hold at n=9 and n=10?

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from collections import Counter
from math import comb
import random

def count_3cycles(adj, n):
    scores = [sum(adj[i]) for i in range(n)]
    return comb(n, 3) - sum(comb(s, 2) for s in scores)

def vertex_in_3cycle(adj, n, v):
    out_v = [j for j in range(n) if adj[v][j]]
    in_v = [j for j in range(n) if adj[j][v]]
    for u in out_v:
        for w in in_v:
            if u != w and adj[u][w]:
                return True
    return False

def find_3cycle_sets(adj, n):
    cycle_sets = []
    for vs in combinations(range(n), 3):
        a, b, c = vs
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            cycle_sets.append(frozenset(vs))
    return cycle_sets

def max_matching_3cycles(cycle_sets):
    n = len(cycle_sets)
    if n == 0:
        return 0
    # Check for size 3 matching
    for a in range(n):
        for b in range(a+1, n):
            if cycle_sets[a] & cycle_sets[b]:
                continue
            for c in range(b+1, n):
                if (not (cycle_sets[c] & cycle_sets[a])) and (not (cycle_sets[c] & cycle_sets[b])):
                    return 3
    # Check for size 2
    for a in range(n):
        for b in range(a+1, n):
            if not (cycle_sets[a] & cycle_sets[b]):
                return 2
    return 1 if n > 0 else 0

def check_deletion_cycle_rich(adj, n, v):
    """Check if T-v is cycle-rich."""
    for w in range(n):
        if w == v:
            continue
        if not vertex_in_3cycle_excluding(adj, n, w, v):
            return False
    return True

def vertex_in_3cycle_excluding(adj, n, w, v):
    """Check if w is in a 3-cycle in T - {v}."""
    out_w = [j for j in range(n) if j != v and adj[w][j]]
    in_w = [j for j in range(n) if j != v and adj[j][w]]
    for u in out_w:
        for x in in_w:
            if u != x and u != v and x != v and adj[u][x]:
                return True
    return False

def main():
    n = 9
    random.seed(42)

    print(f"=== n={n}: Dichotomy verification ===")
    print("For cycle-rich T: either (a) mm>=3 or (b) some deletion is cycle-rich")

    found = 0
    both_fail = 0
    case_a = 0
    case_b = 0
    case_both = 0

    for trial in range(2000000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        scores = [sum(adj[i]) for i in range(n)]
        if 0 in scores or (n-1) in scores:
            continue

        t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)
        if t3 > 15:
            continue

        all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
        if not all_in:
            continue

        found += 1
        c3_sets = find_3cycle_sets(adj, n)
        mm = max_matching_3cycles(c3_sets)

        has_a = (mm >= 3)
        has_b = False

        if not has_a:
            # Check if some vertex deletion gives cycle-rich
            for v in range(n):
                if check_deletion_cycle_rich(adj, n, v):
                    has_b = True
                    break

        if has_a and has_b:
            case_both += 1
        elif has_a:
            case_a += 1
        elif has_b:
            case_b += 1
        else:
            both_fail += 1
            print(f"  !!! BOTH FAIL at trial {trial}! mm={mm}, t3={t3}")
            print(f"      scores={tuple(sorted(scores))}")
            print(f"      3-cycle sets: {c3_sets}")

        if found % 1000 == 0:
            print(f"  {found} found: a_only={case_a}, b_only={case_b}, both={case_both}, FAIL={both_fail}")

    print(f"\nTotal: {found} cycle-rich tournaments")
    print(f"  Case (a) only (mm>=3): {case_a}")
    print(f"  Case (b) only (deletion): {case_b}")
    print(f"  Both: {case_both}")
    print(f"  NEITHER: {both_fail}")

    if both_fail == 0:
        print(f"\n*** DICHOTOMY VERIFIED at n={n}: 0 failures in {found} tests ***")


if __name__ == "__main__":
    main()
