#!/usr/bin/env python3
"""
Analyze the typed independence polynomial at n=7 for all tournaments.

I_typed(y3, y5, y7) = 1 + t3*y3 + t5*y5 + t7*y7 + bc*y3^2

where t3 = #directed 3-cycles, t5 = #5-cycles, t7 = #7-cycles,
bc = #disjoint 3-3 pairs.

At y3=y5=y7=2: I_typed(2,2,2) = 1 + 2*t3 + 2*t5 + 2*t7 + 4*bc = I(Omega, 2) = H(T).

Question: does I_typed carry more information than I(Omega, x)?
I(Omega, x) = 1 + (t3+t5+t7)*x + bc*x^2

The typed polynomial separates alpha_1 into t3, t5, t7.
The standard polynomial loses this separation.

Key test: are there non-isomorphic tournaments with same I(Omega,x)
but different I_typed?

opus-2026-03-07-S37
"""
from itertools import combinations, permutations
from collections import defaultdict
import random

n = 7

def make_tournament(bits):
    """Create tournament from bit vector."""
    A = [[0]*n for _ in range(n)]
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    for k, (i, j) in enumerate(edges):
        if bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A

def count_typed(A):
    """Count t3, t5, t7, bc for tournament A."""
    edge_set = set()
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]:
                edge_set.add((i, j))

    # Count directed 3-cycles
    t3 = 0
    for triple in combinations(range(n), 3):
        a, b, c = triple
        for perm in [(a,b,c), (a,c,b)]:
            x, y, z = perm
            if (x,y) in edge_set and (y,z) in edge_set and (z,x) in edge_set:
                t3 += 1
                break  # one directed cycle per vertex set

    # Count directed 5-cycles
    t5 = 0
    for five in combinations(range(n), 5):
        for p in permutations(five):
            if all((p[i], p[(i+1) % 5]) in edge_set for i in range(5)):
                t5 += 1
                break  # one direction per vertex set... no, need to count ALL directed cycles
    # Actually: for each 5-vertex set, count the number of distinct directed 5-cycles
    # A directed 5-cycle: v0->v1->v2->v3->v4->v0.
    # Two directed cycles are the same iff they're the same cyclic sequence.
    # 5! orderings / 5 (cyclic) = 24. Check each.
    t5 = 0
    for five in combinations(range(n), 5):
        # Count distinct directed 5-cycles on these vertices
        seen = set()
        for p in permutations(five):
            if all((p[i], p[(i+1) % 5]) in edge_set for i in range(5)):
                # Canonical form: start from min vertex
                min_idx = p.index(min(p))
                canonical = tuple(p[min_idx:] + p[:min_idx])
                if canonical not in seen:
                    seen.add(canonical)
                    t5 += 1

    # Count directed 7-cycles (Hamiltonian)
    t7 = 0
    seen7 = set()
    for p in permutations(range(n)):
        if all((p[i], p[(i+1) % n]) in edge_set for i in range(n)):
            min_idx = list(p).index(0)
            canonical = tuple(list(p)[min_idx:] + list(p)[:min_idx])
            if canonical not in seen7:
                seen7.add(canonical)
                t7 += 1

    # Count bc: disjoint 3-cycle pairs
    three_cycles = []
    for triple in combinations(range(n), 3):
        a, b, c = triple
        for perm in [(a,b,c), (a,c,b)]:
            x, y, z = perm
            if (x,y) in edge_set and (y,z) in edge_set and (z,x) in edge_set:
                three_cycles.append(frozenset(triple))
                break

    bc = 0
    for i in range(len(three_cycles)):
        for j in range(i+1, len(three_cycles)):
            if len(three_cycles[i] & three_cycles[j]) == 0:
                bc += 1

    return t3, t5, t7, bc

# Sample random tournaments and collect (I_std, I_typed)
signatures = defaultdict(list)  # I_std -> list of I_typed
typed_to_H = {}

num_samples = 50000
for seed in range(num_samples):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3, t5, t7, bc = count_typed(A)
    alpha1 = t3 + t5 + t7
    H = 1 + 2*alpha1 + 4*bc

    I_std = (alpha1, bc)  # Standard I(Omega, x) = 1 + alpha1*x + bc*x^2
    I_typed = (t3, t5, t7, bc)

    if I_std not in signatures:
        signatures[I_std] = set()
    signatures[I_std].add(I_typed)

    typed_key = I_typed
    if typed_key not in typed_to_H:
        typed_to_H[typed_key] = H

# Find cases where same I_std has different I_typed
print(f"Sampled {num_samples} random tournaments on n=7")
print(f"Distinct I_std signatures: {len(signatures)}")
print(f"Distinct I_typed signatures: {len(typed_to_H)}")

collisions = 0
for I_std, typed_set in signatures.items():
    if len(typed_set) > 1:
        collisions += 1
        if collisions <= 5:
            alpha1, bc = I_std
            H = 1 + 2*alpha1 + 4*bc
            print(f"\nI_std = (alpha1={alpha1}, alpha2={bc}), H={H}:")
            for t in sorted(typed_set):
                print(f"  I_typed = (t3={t[0]}, t5={t[1]}, t7={t[2]}, bc={t[3]})")

print(f"\nTotal I_std values with multiple typed versions: {collisions}")
