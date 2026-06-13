#!/usr/bin/env python3
"""
WHY does odd t3 force t3 = 1 + 2*t5 at n=5?

Equivalently: for ALL n=5 tournaments, t3 - 2*t5 is in {0, 1}.
And t3 - 2*t5 = 1 iff t3 is odd, t3 - 2*t5 = 0 iff... wait let me check.

Actually from the data:
  t3=0, t5=0: t3-2*t5 = 0 (even t3)
  t3=1, t5=0: t3-2*t5 = 1 (odd t3)
  t3=2, t5=0: t3-2*t5 = 2 (even t3)  <-- NOT 0!
  t3=3, t5=1: t3-2*t5 = 1 (odd t3)
  t3=4, t5=1: t3-2*t5 = 2 (even t3)
  t3=4, t5=2: t3-2*t5 = 0 (even t3)
  t3=4, t5=3: t3-2*t5 = -2 (even t3)
  t3=5, t5=2: t3-2*t5 = 1 (odd t3)

So W(0) = 1 - (t3-2*t5) takes values {-1, 0, 1, 3}.
W(0) = 0 iff t3 - 2*t5 = 1, which happens exactly when t3 is odd.

This is a STRONG statement: at n=5, odd t3 forces a unique t5 value.

Actually this is related to:
  H(T) = 1 + 2*t3 + 2*t5 (OCF at n=5, since Omega has no bc terms)
  W(0) = 1 - t3 + 2*t5

So H(T) - W(0) = 3*t3, and H(T) + W(0) = 2 + t3 + 4*t5.

For odd t3 (W(0) = 0): H = 1 + 2*t3 + 2*t5 = 1 + 2*t3 + t3 - 1 = 3*t3.
So H = 3*t3 when t3 is odd at n=5!

Let's verify and check what determines the constraint.

kind-pasteur-2026-03-07-S26
"""
from itertools import combinations
from math import comb

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

def count_directed_cycles(A, n, length):
    """Count directed cycles of given length."""
    from itertools import permutations
    count = 0
    for combo in combinations(range(n), length):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[(i+1)%length]] for i in range(length)):
                count += 1
    return count // length

n = 5
m = num_tiling_bits(n)

print("H(T) = 3*t3 when t3 is odd at n=5?")
print("="*50)

for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    t5 = count_directed_cycles(A, n, 5)

    # H via OCF
    H = 1 + 2*t3 + 2*t5

    if t3 % 2 == 1:
        expected_H = 3*t3
        print(f"  bits={bits}: t3={t3} (odd), t5={t5}, H={H}, 3*t3={expected_H}, match={H==expected_H}")

print(f"\n{'='*50}")
print("Why t3 odd => t5 = (t3-1)/2 at n=5?")
print("="*50)
print("""
At n=5, we have C(5,3) = 10 triples.
t3 = number of directed 3-cycles (2 per cyclic triple).
Transitive triples: 10 - t3.

Score formula: t3 = 10 - sum C(d_i, 2) where d_i are out-degrees.
sum d_i = C(5,2) = 10.

For n=5, the independence polynomial of Omega has:
I(Omega, x) = 1 + t3*x + ? * x^2

Actually at n=5: there is only one 5-vertex set.
The t5 count is the number of directed Hamiltonian cycles in T.

Key constraint: Omega(T) at n=5 has:
- t3 vertices (one per 3-cycle vertex set)
- Edges between vertex sets that share 2 vertices
- i_2 = number of independent pairs = t5 coefficient in I(Omega,2)

So H = I(Omega(T), 2) = 1 + 2*t3 + 2*t5 (since i_2 = t5 at n=5??)
Wait, that can't be right. Let me recalculate.
""")

# Actually compute I(Omega, x) at n=5
from itertools import combinations as comb_iter

for bits in [0, 1, 42, 63]:
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    t5 = count_directed_cycles(A, n, 5)

    # Find all 3-cycle vertex sets
    cycle_sets = []
    for triple in comb_iter(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            cycle_sets.append(frozenset(triple))

    # Build Omega: edges between overlapping cycle sets
    edges = []
    for a, b in comb_iter(range(len(cycle_sets)), 2):
        if len(cycle_sets[a] & cycle_sets[b]) > 0:
            edges.append((a, b))

    # Independent sets of Omega
    # i_0 = 1, i_1 = t3, i_2 = ...
    i_2 = 0
    for a, b in comb_iter(range(len(cycle_sets)), 2):
        if len(cycle_sets[a] & cycle_sets[b]) == 0:
            i_2 += 1

    print(f"  bits={bits}: t3={t3}, t5={t5}, cycle_sets={len(cycle_sets)}, i_2={i_2}")
    print(f"    H via OCF = 1 + 2*{t3} + 4*{i_2} = {1 + 2*t3 + 4*i_2}")
    scores = [sum(A[i][j] for j in range(n)) for i in range(n)]
    print(f"    scores = {sorted(scores)}")

    # Wait: OCF uses 2^{parts(I)} * I_S(T). For single cycles: weight 2.
    # For pairs of disjoint cycles: weight 4 (2^2).
    # So H = 1 + 2*alpha_1 + 4*alpha_2 where alpha_1 = #{cycle vertex sets}
    # and alpha_2 = #{pairs of disjoint cycle vertex sets}.
    # But alpha_1 = t3 (number of directed 3-cycles),
    # and alpha_2 = i_2 (independent pairs in Omega).

    # At n=5: alpha_2 = #{disjoint 3-cycle pairs} (disjoint = no shared vertex)
    # Two 3-element subsets of [5] are disjoint iff... impossible!
    # C(5,3) = 10 triples, but any two 3-element subsets of [5] share at least 1 vertex.
    # So i_2 = 0 always at n=5!

print(f"\n  At n=5: any two 3-element subsets share >= 1 vertex.")
print(f"  So i_2 = 0 ALWAYS.")
print(f"  Therefore H = 1 + 2*t3 + 0 = 1 + 2*t3 at n=5!")
print(f"  Wait, but 5-cycles contribute via t5...")
print(f"  Let me recheck: at n=5, the 5-cycle IS a separate independent set type.")
print(f"  OCF: H = 1 + 2*t3 + 2*t5 (t5 = directed 5-cycles, each is a single vertex in Omega_5)")

# But actually OCF involves Omega which includes BOTH 3-cycles and 5-cycles!
# At n=5: Omega has vertices for 3-cycle sets AND 5-cycle sets.
# A 5-cycle vertex set = all 5 vertices (only one such set).
# It conflicts with ALL 3-cycle sets (shares vertices).
# So in the independence polynomial, t5 vertices are always in conflict with t3 vertices.

# OCF: I(Omega, 2) = 1 + 2*(t3 + t5) + 4*i_2
# where i_2 counts independent pairs in the FULL Omega.
# Since 5-cycle vertex sets conflict with everything, i_2 = 0.
# So H = 1 + 2*t3 + 2*t5.

# And W(0) = 1 - t3 + 2*t5 from THM-059.
# H = 1 + 2*t3 + 2*t5
# W(0) = 1 - t3 + 2*t5
# H - W(0) = 3*t3
# H + W(0) = 2 + t3 + 4*t5

# For odd t3: H = 1 + 2*t3 + 2*t5 and W(0) = 0
# => H = 1 + 2*t3 + 2*t5, but we need W(0)=0 => t3 = 1 + 2*t5 => t5 = (t3-1)/2
# => H = 1 + 2*t3 + t3 - 1 = 3*t3

# Why is t5 = (t3-1)/2 forced when t3 is odd?
# t5 = number of directed Hamiltonian cycles on 5 vertices.
# For a tournament on 5 vertices, there's a classical formula:
# t5 = C(5,5) - 5*t3 + sum something...

# Actually: the number of Hamiltonian cycles in a tournament on n vertices
# can be computed from the cycle counts. For n=5:
# Total cycles = t3 + t5 (only 3-cycles and 5-cycles)
# Each vertex is in exactly... hmm.

# Let's just check: what determines t5 given t3?
print(f"\n{'='*50}")
print("t5 as a function of t3 at n=5")
print("="*50)

from collections import defaultdict
t3_to_t5 = defaultdict(set)
for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t3_to_t5[t3].add(t5)

for t3_val in sorted(t3_to_t5.keys()):
    t5_vals = sorted(t3_to_t5[t3_val])
    pred = (t3_val - 1) / 2 if t3_val % 2 == 1 else "N/A"
    print(f"  t3={t3_val}: t5 in {t5_vals}, predicted (odd): {pred}")

# Note: for the REGULAR tournament (t3=5), t5=2.
# t3=5, t5=2: t3 = 1 + 2*2 = 5. ✓

# This is actually related to a classical result.
# For a tournament on 5 vertices with score sequence (d_0,...,d_4):
# t3 = 10 - sum C(d_i, 2)
# And there's a formula for t5 involving higher moments of the score sequence.

# Let me look up: Moon's formula for Hamilton cycles in tournaments.
# t5 = something involving score sequence...

# Actually, at n=5, there's a result: t5 = 12 - sum d_i^2 + (stuff)
# Let me just check directly.

print(f"\n{'='*50}")
print("t5 formula from score sequence?")
print("="*50)

from collections import defaultdict
score_to_t3_t5 = defaultdict(set)
for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
    t3 = count_3cycles(A, n)
    t5 = count_directed_cycles(A, n, 5)
    score_to_t3_t5[scores].add((t3, t5))

for scores in sorted(score_to_t3_t5.keys()):
    pairs = sorted(score_to_t3_t5[scores])
    sum_sq = sum(s**2 for s in scores)
    sum_cube = sum(s**3 for s in scores)
    print(f"  scores={scores}: sum_d^2={sum_sq}, sum_d^3={sum_cube}, (t3,t5) = {pairs}")

# For score (1,2,2,2,3):
# t3=4 is forced, but t5 can be 1, 2, or 3.
# So t5 is NOT determined by score sequence alone!
# But the constraint t3 odd => t5 = (t3-1)/2 still holds.

# This means the constraint is deeper than score sequences.
