#!/usr/bin/env python3
"""
Check for mixed-direction cycle collections at n=9.

At n=7, no mixed collections exist. At n=9 with 3 disjoint 3-cycles,
there might be room for one T-cycle and one T^op-cycle.

opus-2026-03-07-S37
"""
from itertools import permutations, combinations
import random

n = 9

def find_directed_3cycles(n, edge_set):
    """Find all directed 3-cycles."""
    cycles_T = []
    cycles_Top = []
    opp = set((j, i) for (i, j) in edge_set)

    for triple in combinations(range(n), 3):
        a, b, c = triple
        # Check all 2 orientations of 3-cycle
        for perm in [(a,b,c), (a,c,b)]:
            x, y, z = perm
            if (x,y) in edge_set and (y,z) in edge_set and (z,x) in edge_set:
                cycles_T.append(perm)
            if (x,y) in opp and (y,z) in opp and (z,x) in opp:
                cycles_Top.append(perm)
    return cycles_T, cycles_Top

def check_mixed_disjoint_pairs(n, edge_set):
    """Check if there exist a T-3-cycle and T^op-3-cycle that are vertex-disjoint."""
    opp = set((j, i) for (i, j) in edge_set)

    T_cycles = []
    Top_cycles = []

    for triple in combinations(range(n), 3):
        a, b, c = triple
        for perm in [(a,b,c), (a,c,b)]:
            x, y, z = perm
            if (x,y) in edge_set and (y,z) in edge_set and (z,x) in edge_set:
                T_cycles.append(frozenset(perm))
            if (x,y) in opp and (y,z) in opp and (z,x) in opp:
                Top_cycles.append(frozenset(perm))

    # Remove duplicates (same vertex set, different starting point)
    T_vsets = set(T_cycles)
    Top_vsets = set(Top_cycles)

    mixed_pairs = 0
    for tv in T_vsets:
        for topv in Top_vsets:
            if len(tv & topv) == 0:
                mixed_pairs += 1
                if mixed_pairs <= 3:
                    print(f"    Mixed pair: T-cycle on {sorted(tv)}, T^op-cycle on {sorted(topv)}")

    return mixed_pairs, len(T_vsets), len(Top_vsets)

print(f"Checking for mixed-direction disjoint 3-cycle pairs at n={n}...")
print()

found_any = False
for seed in range(100):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    edge_set = set()
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]:
                edge_set.add((i, j))

    mp, nt, ntop = check_mixed_disjoint_pairs(n, edge_set)
    if mp > 0:
        if not found_any:
            print(f"  seed={seed}: {mp} mixed pairs (T-cycles: {nt}, T^op-cycles: {ntop})")
            found_any = True
        elif seed < 10:
            print(f"  seed={seed}: {mp} mixed pairs")

if not found_any:
    print("  NO mixed disjoint pairs found in 100 random tournaments!")
else:
    print(f"\n  Mixed pairs exist! Checking the full Rédei-Berge sum...")
    # For the first tournament with mixed pairs, verify if they contribute to the sum
    # (they should, as individual 3-cycles in T and T^op)

    # But wait: in the Rédei-Berge framework, we need the PERMUTATION
    # to have ALL nontrivial cycles be directed cycles of T or T^op.
    # A permutation with one T-cycle and one T^op-cycle IS a valid permutation
    # in the sum. Let me check this contributes.
    pass

# Actually, let me re-examine the original check.
# The issue might be simpler: at n=5 and n=7, there simply isn't ROOM
# for two disjoint 3-cycles (one T, one T^op) because 3+3=6 <= 7
# and the vertex constraints are too tight.

# At n=9: 3+3=6 <= 9, so there's room. Let me check if there are
# actually disjoint T and T^op 3-cycles.

print(f"\n=== Direct 3-cycle count at n=9 ===")
rng = random.Random(42)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if rng.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

edge_set = set()
opp_set = set()
for i in range(n):
    for j in range(n):
        if i != j:
            if A[i][j]:
                edge_set.add((i, j))
            else:
                opp_set.add((i, j))

# Every directed 3-cycle (a,b,c) with a->b->c->a in T
# corresponds to the cycle (a,c,b) with a->c->b->a in T^op
# (since a->c in T^op means c->a in T, etc.)
# So cycle reversal maps T-3-cycles to T^op-3-cycles.

# The vertex SETS are the same! A 3-cycle on {a,b,c} in T gives
# a reversed 3-cycle on the SAME {a,b,c} in T^op.

# Therefore: for 3-cycles, the T-cycle vertex sets and T^op-cycle
# vertex sets are IDENTICAL. Two disjoint T and T^op cycles means
# two disjoint triples, one with T-direction, one with T^op-direction.
# These are vertex-disjoint, so both are independent in the sense
# of the Rédei-Berge sum.

# But can they COEXIST? Of course! Just take any two vertex-disjoint
# triples {a,b,c} and {d,e,f}. One triple forms a T-3-cycle,
# the other forms a T^op-3-cycle (equivalently, a T-3-cycle in reverse).

# Wait, EVERY triple in a tournament contains exactly one directed 3-cycle
# OR is acyclic (transitive). If it's a 3-cycle, the cycle goes in one
# direction (say T). The reverse goes in T^op.

# So a "T^op 3-cycle on {d,e,f}" is just the reversal of the
# "T 3-cycle on {d,e,f}" (if one exists). But as a DIRECTED CYCLE,
# d->f->e->d (T^op) is different from d->e->f->d (T).

# As a PERMUTATION: the T-cycle (d,e,f) maps d->e->f->d.
# The T^op cycle on the same vertices: d->f->e->d, which as a
# permutation is (d,f,e).

# So the permutation sigma with T-cycle (a,b,c) and T^op-cycle (d,f,e)
# (and fixed points g,h,i) IS a valid permutation in the Rédei-Berge sum.
# Its type is (3,3,1,1,1).

# phi = sum of (len-1) for T-cycles only = 2 (from the T-3-cycle).
# Sign = (-1)^2 = +1.
# Specialization: p_3 * p_3 * p_1^3 -> 2 * 2 * 1 = 4.
# Contribution: +4.

# But there's ALSO the permutation with T^op-cycle (a,c,b) and T-cycle (d,e,f).
# phi = 2 (from the T-cycle). Sign = +1. Same contribution +4.

# And the permutation with T-cycle (a,b,c) and T-cycle (d,e,f).
# phi = 2 + 2 = 4. Sign = +1. Contribution +4.

# And T^op (a,c,b) and T^op (d,f,e).
# phi = 0. Sign = +1. Contribution +4.

# So there are 4 permutations total, each contributing +4.
# But the OCF says this collection of vertex sets {a,b,c} and {d,e,f}
# should contribute 2^2 = 4 (one term for the pair of disjoint cycles).

# Wait, there are TWO disjoint cycle collections here:
# 1. T-cycle {a,b,c} + T-cycle {d,e,f} — this is ONE independent set
#    in Omega(T) of size 2.
# 2. T^op-cycle {a,b,c} + T^op-cycle {d,e,f} — this is ONE independent
#    set in Omega(T^op) of size 2.

# But in the OCF, we only count T-cycle collections. So collection 1
# contributes 2^2 = 4.

# In the FULL Rédei-Berge sum (T and T^op mixed):
# We get 4 permutations, each contributing +4, total = 16.
# But the OCF gives H(T) once, not 4 times.

# This contradicts our earlier finding that mixed = T_only = 189 for Paley!
# Unless the mixed permutations DON'T exist...

# The key: at n=5 and n=7, there were NO mixed-direction permutations.
# Is this a general theorem for tournaments?

# CLAIM: In a tournament, if C is a directed 3-cycle on {a,b,c}
# and D is a directed 3-cycle of T^op on disjoint {d,e,f}, then
# this creates a VALID permutation in S(T, T^op) of type (3,3,...).

# But our computation said there are 0 such permutations at n=7!
# How is this possible?

# Let me think... Wait, maybe my check was wrong. Let me re-examine.

# At n=7, for Paley: every triple has a 3-cycle (since score is regular).
# So both T-cycles and T^op-cycles exist on every triple.
# Take any two disjoint triples. Each has a T-cycle and a T^op-cycle.
# The permutation with T-cycle on first triple and T^op-cycle on second
# SHOULD be valid.

# Unless my code had a bug. Let me recheck with explicit examples.

print(f"\n=== Manual check at n=7 (Paley) ===")
n7 = 7
QR = {1, 2, 4}
paley7 = set()
for i in range(n7):
    for j in range(n7):
        if i != j and (j - i) % n7 in QR:
            paley7.add((i, j))

# Find a T-3-cycle
for a, b, c in combinations(range(n7), 3):
    if (a,b) in paley7 and (b,c) in paley7 and (c,a) in paley7:
        print(f"T-cycle: ({a},{b},{c})")
        # Find a T^op-cycle on disjoint vertices
        remaining = [v for v in range(n7) if v not in {a,b,c}]
        for d, e, f in combinations(remaining, 3):
            opp = set((j, i) for (i, j) in paley7)
            if (d,e) in opp and (e,f) in opp and (f,d) in opp:
                print(f"  T^op-cycle: ({d},{e},{f}) (i.e. T has {d}<-{e}<-{f}<-{d})")
                # Check: is this a valid permutation in the Rédei-Berge sum?
                # sigma: a->b, b->c, c->a, d->e, e->f, f->d, remaining fixed
                # Wait: T^op cycle (d,e,f) means d->e in T^op, so e->d in T.
                # As a permutation: sigma(d)=e, sigma(e)=f, sigma(f)=d.
                # Check if (d,e,f) is a directed cycle of T^op:
                print(f"    (d,e) in T^op: {(d,e) in opp}, (e,f) in T^op: {(e,f) in opp}, (f,d) in T^op: {(f,d) in opp}")
                break
        break
