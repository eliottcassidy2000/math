#!/usr/bin/env python3
"""
EVEN BLOCK VANISHING THEOREM

THEOREM: In the position-topology decomposition of e_{2k}, a block of ODD length
(involving EVEN number of vertices) contributes ZERO.

Equivalently: only blocks of EVEN length (involving ODD number of vertices)
contribute to the coefficients c_{n-1-2k}.

PROOF:
For a block of length m, the contribution involves c_0(S) of the (m+1)-vertex
subtournament. But c_0(S) is the constant term of tr(M_S(r)).

At even vertex count (m+1 even, so m odd):
tr(M_S(r)) has ONLY ODD powers of r (not even).
Therefore c_0(S) = 0.

Why only odd powers at even n? Because the complement map s_i -> -s_i
sends prod(r+s_i) -> prod(r-s_i) = (-1)^{n-1} prod(-r+s_i).
At even n-1 (odd n): (-1)^{n-1} = 1, so tr(M(r)) = tr(M(-r)) [even powers only].
At odd n-1 (even n): (-1)^{n-1} = -1, so tr(M(r)) = -tr(M(-r)) [odd powers only].

(More precisely: summing over all permutations P of complement T^op,
we get the same multiset of products since P -> P is a bijection.
So tr(M_T(r)) = (-1)^{n-1} * tr(M_T(-r)).)

CONSEQUENCE:
The surviving patterns for e_{2k} are exactly those whose run structure
has ALL runs of EVEN length (2, 4, 6, ...).

EXAMPLES:
At n=7, e_4 has patterns:
  (4,): single block of 4 -> 5 vertices -> contributes
  (2,2): two blocks of 2 -> 3+3 vertices -> contributes

At n=9, e_6 has patterns:
  (6,): contributes (7v)
  (4,2): contributes (5v + 3v)
  (3,3): ZERO (4v + 4v, even blocks!)
  (2,2,2): contributes (3v + 3v + 3v)

At n=11, e_8 has patterns:
  (8,): contributes (9v)
  (6,2): contributes (7v + 3v)
  (4,4): contributes (5v + 5v)
  (4,2,2): contributes (5v + 3v + 3v)
  (2,2,2,2): contributes (3v + 3v + 3v + 3v)
  (5,3): ZERO (6v + 4v)
  (3,3,2): ZERO (4v + 4v + 3v — one even block)

PARTITIONS OF EVEN NUMBERS INTO EVEN PARTS:
These are precisely the partitions counted by the number of even partitions of 2k.

opus-2026-03-06-S11b (continued^6)
"""
from itertools import combinations, permutations
from collections import defaultdict

def run_structure(S):
    S_sorted = sorted(S)
    runs = []
    current = 1
    for i in range(1, len(S_sorted)):
        if S_sorted[i] == S_sorted[i-1] + 1:
            current += 1
        else:
            runs.append(current)
            current = 1
    runs.append(current)
    return tuple(sorted(runs, reverse=True))

def has_isolated(S):
    S_set = set(S)
    return any((s-1) not in S_set and (s+1) not in S_set for s in S)

def has_odd_block(S):
    S_sorted = sorted(S)
    runs = []
    current = 1
    for i in range(1, len(S_sorted)):
        if S_sorted[i] == S_sorted[i-1] + 1:
            current += 1
        else:
            runs.append(current)
            current = 1
    runs.append(current)
    return any(r % 2 == 1 for r in runs)

# Verify computationally at n=9
n = 9
print(f"Pattern analysis at n={n}:")
print(f"Positions 0..{n-2}")
print()

for subset_size in range(2, n, 2):
    patterns = defaultdict(int)
    for S in combinations(range(n-1), subset_size):
        if has_isolated(S):
            continue
        rs = run_structure(S)
        patterns[rs] += 1

    print(f"e_{subset_size} (position {subset_size}-subsets without isolated elements):")
    for rs, count in sorted(patterns.items()):
        odd_block = any(r % 2 == 1 for r in rs)
        n_verts = sum(r + 1 for r in rs)
        status = "ZERO (odd block!)" if odd_block else "CONTRIBUTES"
        print(f"  {rs}: {count} subsets, {n_verts} vertices -> {status}")
    print()

# Verify the ZERO claim by direct computation
print("=" * 60)
print("VERIFICATION: (3,3) pattern gives exactly zero at n=9")
import random
random.seed(42)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

# Find (3,3) patterns
for S in combinations(range(n-1), 6):
    if has_isolated(S):
        continue
    rs = run_structure(S)
    if rs == (3, 3):
        # Direct computation
        contrib = 0.0
        for p in permutations(range(n)):
            prod_val = 1.0
            for pos in S:
                prod_val *= (A[p[pos]][p[pos+1]] - 0.5)
            contrib += prod_val
        print(f"  Pattern {S}: contribution = {contrib:.10f}")

# Also verify c_0 = 0 for all n=4 tournaments
print()
print("Verification: c_0 = 0 for ALL n=4 tournaments")
n4 = 4
for bits in range(64):
    A4 = [[0]*n4 for _ in range(n4)]
    edges4 = [(i,j) for i in range(n4) for j in range(i+1, n4)]
    for idx, (i,j) in enumerate(edges4):
        if bits & (1 << idx):
            A4[i][j] = 1
        else:
            A4[j][i] = 1
    c0 = sum(
        eval("*".join(f"({A4[p[k]][p[k+1]]}-0.5)" for k in range(n4-1)))
        for p in permutations(range(n4))
    )
    if abs(c0) > 1e-10:
        print(f"  COUNTEREXAMPLE at bits={bits}: c_0 = {c0}")
        break
else:
    print(f"  Confirmed: c_0 = 0 for all {64} tournaments on 4 vertices")

# Similarly for n=6
print()
print("Verification: c_0 = 0 for sample n=6 tournaments")
n6 = 6
for trial in range(20):
    random.seed(trial * 17 + 3)
    A6 = [[0]*n6 for _ in range(n6)]
    for i in range(n6):
        for j in range(i+1, n6):
            if random.random() < 0.5:
                A6[i][j] = 1
            else:
                A6[j][i] = 1
    c0 = sum(
        eval("*".join(f"({A6[p[k]][p[k+1]]}-0.5)" for k in range(n6-1)))
        for p in permutations(range(n6))
    )
    if abs(c0) > 1e-10:
        print(f"  COUNTEREXAMPLE at trial={trial}: c_0 = {c0}")
        break
else:
    print(f"  Confirmed: c_0 = 0 for all 20 sample n=6 tournaments")

print("\nDONE")
