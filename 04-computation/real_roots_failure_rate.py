#!/usr/bin/env python3
"""
What fraction of n=9 tournaments have non-real-rooted I(Omega(T), x)?
What characterizes the failing tournaments?

We use Omega_3 (3-cycles only) as a proxy since:
1. It's much faster to compute
2. The full Omega failure implies Omega_3 failure (in the known counterexample)

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament
from math import comb
import numpy as np

def find_3cycles(T):
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i,j,k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i,j,k))
    return cycles

def count_a2_a3(cycles):
    m = len(cycles)
    cycle_sets = [frozenset(c) for c in cycles]
    disjoint_after = [[] for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if not (cycle_sets[i] & cycle_sets[j]):
                disjoint_after[i].append(j)
    a2 = sum(len(d) for d in disjoint_after)
    a3 = 0
    for i in range(m):
        ui = cycle_sets[i]
        for j in disjoint_after[i]:
            uij = ui | cycle_sets[j]
            for k in disjoint_after[j]:
                if k > j and not (cycle_sets[k] & uij):
                    a3 += 1
    return a2, a3

n = 9
samples = 100000
print(f"Sampling {samples} random n={n} tournaments...")

fail_omega3 = 0
total_deg3 = 0
fail_triples = []
score_seqs = {}

for trial in range(samples):
    T = random_tournament(n)
    c3_list = find_3cycles(T)
    a1 = len(c3_list)
    if a1 < 3:
        continue

    a2, a3 = count_a2_a3(c3_list)
    if a3 == 0:
        continue

    total_deg3 += 1
    disc = 18*a1*a2*a3 - 4*a1**3*a3 + a1**2*a2**2 - 4*a2**3 - 27*a3**2

    if disc < 0:
        fail_omega3 += 1
        fail_triples.append((a1, a2, a3, disc))
        scores = tuple(sorted(sum(row) for row in T))
        score_seqs[scores] = score_seqs.get(scores, 0) + 1

    if (trial + 1) % 25000 == 0:
        print(f"  ... {trial+1}/{samples}: {fail_omega3} failures so far")

print(f"\nResults:")
print(f"  Total tournaments: {samples}")
print(f"  Degree-3 polynomials (a3 > 0): {total_deg3}")
print(f"  Omega_3 real-root failures: {fail_omega3}")
print(f"  Failure rate (among deg-3): {fail_omega3/total_deg3*100:.3f}%")
print(f"  Failure rate (overall): {fail_omega3/samples*100:.4f}%")

if fail_triples:
    print(f"\n  Distinct failing (a1, a2, a3) triples:")
    from collections import Counter
    triple_counts = Counter((a1,a2,a3) for a1,a2,a3,d in fail_triples)
    for (a1,a2,a3), count in triple_counts.most_common(20):
        disc = 18*a1*a2*a3 - 4*a1**3*a3 + a1**2*a2**2 - 4*a2**3 - 27*a3**2
        print(f"    ({a1:2d}, {a2:2d}, {a3:2d}): disc={disc:7d}, count={count}")

    print(f"\n  Score sequences of failing tournaments:")
    for scores, count in sorted(score_seqs.items(), key=lambda x: -x[1])[:15]:
        print(f"    {scores}: {count} times")

    # Characterize: what's special about these?
    # Look at score variance
    print(f"\n  Score variance of failing tournaments:")
    variances = []
    for scores, count in score_seqs.items():
        mean = sum(scores) / len(scores)
        var = sum((s - mean)**2 for s in scores) / len(scores)
        variances.append((var, scores, count))
    variances.sort(reverse=True)
    for var, scores, count in variances[:10]:
        print(f"    var={var:.2f}: {scores} ({count} times)")
