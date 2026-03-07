#!/usr/bin/env python3
"""
WHY does W(0) = 0 for all odd-t3 GS tilings at n=5?

W(0) = 1 - t3 + 2*t5 at n=5.
W(0) = 0 iff t3 = 1 + 2*t5.

For ALL GS tilings with odd t3 at n=5, we had W(0) = 0.
So: odd t3 implies t3 = 1 + 2*t5 for GS tilings.

Is this true for ALL n=5 tournaments with odd t3? Or only GS ones?
Is this a GS-specific constraint?

Also: at n=5, the tournament is determined by 6 edges.
t3 + C(5,3) - 2*t3 = C(5,3) - t3 is the count of transitive triples.
Actually: total triples = C(5,3) = 10, each is cyclic or transitive.
For 3 vertices with backbone i->i+1, the third vertex determines the type.

Let me check: for NON-GS tournaments, does t3 odd imply t3 = 1 + 2*t5?

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_5cycles(A, n):
    t5 = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t5 += 1
    return t5 // 5

n = 5
m = num_tiling_bits(n)

print(f"t3 vs t5 for ALL n={n} tournaments (2^{m} = {2**m})")
print("="*60)

# Check all
odd_t3_constraint = True
even_t3_constraint = True
t3_t5_pairs = defaultdict(int)

for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    t5 = count_5cycles(A, n)
    t3_t5_pairs[(t3, t5)] += 1

    if t3 % 2 == 1:
        if t3 != 1 + 2*t5:
            odd_t3_constraint = False

print("(t3, t5) -> count:")
for (t3, t5) in sorted(t3_t5_pairs.keys()):
    w0 = 1 - t3 + 2*t5
    print(f"  t3={t3}, t5={t5}: count={t3_t5_pairs[(t3,t5)]}, t3%2={t3%2}, W(0)={w0}")

print(f"\nOdd t3 always has t3 = 1 + 2*t5? {odd_t3_constraint}")

# The constraint: at n=5, every tournament has exactly one 5-vertex set.
# On that set, the tournament is either:
# - Has a directed Hamiltonian cycle (5-cycle) or not
# But t5 counts directed 5-cycles = directed Hamiltonian cycles!

# Relation between t3 and t5 at n=5:
# Total edges = C(5,2) = 10
# Each vertex has out-degree d_i, sum d_i = 10
# t3 = C(5,3) - sum C(d_i, 2) + ... (Kendall-Babington-Smith formula)
# Actually: t3 = C(n,3) - sum_i C(d_i, 2)  <-- this is for n=5?

# Let me verify: for score sequence (d_0, ..., d_4) with sum = 10
# t3 = C(5,3) - sum C(d_i, 2) = 10 - sum d_i(d_i-1)/2

for bits in [0, 1, 2**m-1]:
    A = tournament_from_tiling(n, bits)
    scores = [sum(A[i][j] for j in range(n)) for i in range(n)]
    t3 = count_3cycles(A, n)
    t3_formula = 10 - sum(s*(s-1)//2 for s in scores)
    print(f"\n  bits={bits}: scores={sorted(scores)}, t3={t3}, formula={t3_formula}, match={t3==t3_formula}")

# So t3 is determined by score sequence!
# And t5 at n=5?
print(f"\n{'='*60}")
print("Score sequence -> (t3, t5) mapping")
print("="*60)

score_to_invariants = defaultdict(set)
for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
    t3 = count_3cycles(A, n)
    t5 = count_5cycles(A, n)
    score_to_invariants[scores].add((t3, t5))

for scores in sorted(score_to_invariants.keys()):
    invariants = sorted(score_to_invariants[scores])
    print(f"  scores={scores}: (t3,t5) = {invariants}")

print(f"\n  => t3 is determined by scores (single value per score seq)")
print(f"  => t5 may or may not be (check if multiple values per score)")
