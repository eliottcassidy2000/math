#!/usr/bin/env python3
"""
gammoid_matroid_test.py — Do vertex-disjoint odd cycles in tournaments form a matroid?

BACKGROUND:
- Omega(T) = conflict graph of odd directed cycles in T
- I(Omega(T), 2) = sum over independent sets S of 2^{|S|} = H(T) (OCF)
- Independent sets of Omega(T) correspond to collections of
  vertex-disjoint odd directed cycles in T

QUESTION: Does the collection of vertex-disjoint odd cycle sets form a MATROID?
If so, which kind? A gammoid (linkage matroid on a digraph)?

A GAMMOID on digraph D with sources S:
  Independent sets = subsets of V that can be linked to S by vertex-disjoint paths.

For our case:
  "Elements" = odd directed cycles of T
  "Independent" = pairwise vertex-disjoint collection

This is an INDEPENDENCE SYSTEM but is it a matroid?
A matroid must satisfy the EXCHANGE AXIOM:
  If I, J are independent with |I| < |J|, then exists j in J\I with I+j independent.

ALSO: Is this related to the TRANSVERSAL MATROID or PARTITION MATROID?

Author: opus-2026-03-07-S45
"""
from itertools import permutations, combinations
import math

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def find_all_odd_cycles(adj, n):
    """Find all directed odd cycles in tournament T."""
    cycles = []
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                if all(adj[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    # Canonical form: start with smallest vertex
                    min_idx = perm.index(min(perm))
                    canonical = perm[min_idx:] + perm[:min_idx]
                    if canonical not in cycles:
                        cycles.append(canonical)
                    break  # only need one orientation per vertex set
    return cycles

def find_all_odd_cycles_by_vset(adj, n):
    """Find odd directed cycles, keyed by vertex set."""
    cycles = {}  # frozenset(vertices) -> True if cycle exists
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                if all(adj[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    cycles[frozenset(combo)] = perm
                    break
    return cycles

def check_matroid_exchange(independent_sets):
    """Check the matroid exchange axiom on a collection of independent sets."""
    # For every pair I, J with |I| < |J|, check exists j in J\I with I+{j} independent.
    all_sets = [frozenset(s) for s in independent_sets]
    failures = 0
    total_checks = 0

    for I in all_sets:
        for J in all_sets:
            if len(I) >= len(J):
                continue
            total_checks += 1
            # Find j in J\I such that I + {j} is independent
            found = False
            for j in J - I:
                test = I | {j}
                if test in all_sets:
                    found = True
                    break
            if not found:
                failures += 1
                if failures <= 3:
                    print(f"    EXCHANGE FAIL: I={sorted(I)}, J={sorted(J)}")
                    print(f"      J\\I = {sorted(J-I)}")
                    for j in J - I:
                        test = I | {j}
                        print(f"      I+{j} = {sorted(test)}: {'indep' if test in all_sets else 'NOT indep'}")

    return failures, total_checks

# ============================================================
# TEST: Matroid property of vertex-disjoint odd cycles
# ============================================================
print("=" * 60)
print("MATROID TEST: vertex-disjoint odd cycles")
print("=" * 60)

for n in [5, 6, 7]:
    m = n*(n-1)//2
    num_tournaments = min(1 << m, 100 if n >= 7 else 1 << m)

    import random
    random.seed(42)

    matroid_pass = 0
    matroid_fail = 0

    for trial in range(num_tournaments):
        if n < 7:
            bits = trial
        else:
            bits = random.getrandbits(m)

        adj = tournament_from_bits(n, bits)
        cycles = find_all_odd_cycles_by_vset(adj, n)

        if not cycles:
            matroid_pass += 1
            continue

        # Find all independent sets (maximal collections of vertex-disjoint cycles)
        cycle_list = list(cycles.keys())
        num_cycles = len(cycle_list)

        # Find all collections of pairwise vertex-disjoint cycles
        indep_sets = set()
        indep_sets.add(frozenset())  # empty set

        def find_independent(current, start_idx, used_vertices):
            indep_sets.add(frozenset(current))
            for i in range(start_idx, num_cycles):
                c = cycle_list[i]
                if c & used_vertices:
                    continue
                current.append(i)
                find_independent(current, i+1, used_vertices | c)
                current.pop()

        find_independent([], 0, frozenset())

        # Convert to matroid ground set: elements are CYCLE INDICES
        # Independent sets: subsets of cycle indices that are vertex-disjoint
        failures, total = check_matroid_exchange(indep_sets)

        if failures == 0:
            matroid_pass += 1
        else:
            matroid_fail += 1
            if matroid_fail <= 2:
                print(f"\n  n={n} bits={bits}: MATROID FAILS ({failures}/{total})")
                print(f"    {num_cycles} cycles, {len(indep_sets)} independent sets")

    print(f"\n  n={n}: matroid pass={matroid_pass}, fail={matroid_fail}")

# ============================================================
# ALTERNATIVE: Matroid on VERTICES using cycle coverings
# ============================================================
print("\n" + "=" * 60)
print("PARTITION MATROID TEST")
print("=" * 60)

# A simpler question: for each vertex v, the cycles through v form a "partition" class.
# The independent sets are TRANSVERSALS — pick at most one cycle per "something".
# This would be a PARTITION matroid.
# But cycles can share multiple vertices, so it's not a clean partition.

# Let's check: is the maximum independent set size equal across all maximal independent sets?
# (necessary for matroid: all bases have the same size)

for n in [5, 6, 7]:
    m = n*(n-1)//2
    num_tournaments = min(1 << m, 200 if n >= 7 else 1 << m)

    import random
    random.seed(42)

    same_max = 0
    diff_max = 0

    for trial in range(num_tournaments):
        if n < 7:
            bits = trial
        else:
            bits = random.getrandbits(m)

        adj = tournament_from_bits(n, bits)
        cycles = find_all_odd_cycles_by_vset(adj, n)

        if not cycles:
            same_max += 1
            continue

        cycle_list = list(cycles.keys())
        num_cycles = len(cycle_list)

        # Find all maximal independent sets
        indep_sets = set()
        indep_sets.add(frozenset())

        def find_independent2(current, start_idx, used_vertices):
            indep_sets.add(frozenset(current))
            for i in range(start_idx, num_cycles):
                c = cycle_list[i]
                if c & used_vertices:
                    continue
                current.append(i)
                find_independent2(current, i+1, used_vertices | c)
                current.pop()

        find_independent2([], 0, frozenset())

        max_size = max(len(s) for s in indep_sets)
        # Are all maximal independent sets the same size?
        maximal = [s for s in indep_sets
                    if not any(s < t for t in indep_sets)]
        sizes = set(len(s) for s in maximal)

        if len(sizes) == 1:
            same_max += 1
        else:
            diff_max += 1
            if diff_max <= 3:
                print(f"  n={n} bits={bits}: maximal sizes={sorted(sizes)}, #maximal={len(maximal)}")

    print(f"  n={n}: all maximal same size={same_max}, different sizes={diff_max}")

# ============================================================
# BONUS: I(Omega, x) — characteristic polynomial structure
# ============================================================
print("\n" + "=" * 60)
print("INDEPENDENCE POLYNOMIAL STRUCTURE")
print("=" * 60)

# For small n, compute I(Omega(T), x) = sum_k a_k * x^k
# where a_k = # collections of k pairwise vertex-disjoint odd cycles

for n in [5, 6, 7]:
    m = n*(n-1)//2

    if n >= 7:
        import random
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(50)]
    else:
        iterator = range(1 << m)

    seen = set()
    I_polys = []

    for bits in iterator:
        adj = tournament_from_bits(n, bits)

        # Quick F computation for H
        F = [0]*n
        for P in permutations(range(n)):
            fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
            F[fwd] += 1
        H = sum(F)

        cycles = find_all_odd_cycles_by_vset(adj, n)
        cycle_list = list(cycles.keys())
        num_cycles = len(cycle_list)

        # Count by size of independent set
        a = [0] * (num_cycles + 1)

        def count_indep(start_idx, used_vertices, k):
            a[k] += 1
            for i in range(start_idx, num_cycles):
                c = cycle_list[i]
                if c & used_vertices:
                    continue
                count_indep(i+1, used_vertices | c, k+1)

        count_indep(0, frozenset(), 0)

        I_poly = tuple(ak for ak in a if ak > 0 or True)
        key = I_poly

        if key not in seen:
            seen.add(key)
            I_at_2 = sum(a[k] * 2**k for k in range(len(a)))
            I_list = [a[k] for k in range(len(a)) if k <= num_cycles][:10]  # truncate

            if n <= 6 or len(I_polys) < 10:
                print(f"  n={n}: alpha_1={num_cycles}, I(x)={I_list[:6]}, I(2)={I_at_2}, H={H}")
                if I_at_2 != H:
                    print(f"    *** OCF MISMATCH: I(2)={I_at_2} != H={H}")

            I_polys.append((I_list, H, num_cycles))

    print(f"  n={n}: {len(I_polys)} distinct I polynomials")
