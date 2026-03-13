#!/usr/bin/env python3
"""
cycle_disjointness_constraint.py -- kind-pasteur-2026-03-13-S61

Test whether the cycle-disjointness constraint from THM-165 holds more broadly.

At regular n=7: c5_dir + 2*disj_33 = 56 (constant)
This makes the c5 coefficient zero in the H bridge equation.

Questions:
1. Does a similar constraint hold at n=6 within score classes?
2. Does it hold at n=5?
3. What is the general form?

The constraint says: H = 1 + 2*alpha_1 + 4*alpha_2 where
alpha_1 = total directed odd cycles and alpha_2 = disjoint cycle pairs.
Within a score class, c3_dir is fixed, so:
  H = 1 + 2*c3_dir + 2*c5_dir + 2*c7_dir + 4*alpha_2
  H - 1 - 2*c3_dir = 2*(c5_dir + c7_dir) + 4*alpha_2

If H depends only on c_{n}_dir within the class, then:
  2*(c5_dir + c7_dir) + 4*alpha_2 must simplify.

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_subset(A, verts):
    k = len(verts)
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


# ========================================================================
# ANALYSIS 1: n=5 — FULL DECOMPOSITION
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: n=5 FULL DIRECTED CYCLE + DISJOINTNESS")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m

by_score = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    # Directed 3-cycles (from first vertex of each 3-subset)
    c3_dir = 0
    c3_vsets = []
    for a, b, c in combinations(range(n), 3):
        cnt = count_directed_ham_cycles_subset(A, [a, b, c])
        c3_dir += cnt
        if cnt > 0:
            c3_vsets.append(frozenset([a, b, c]))

    # Directed 5-cycles (Hamiltonian)
    c5_dir = count_directed_ham_cycles_subset(A, list(range(n)))

    # Total directed cycles
    alpha_1 = c3_dir + c5_dir

    # Disjoint pairs (vertex-disjoint directed cycles)
    # At n=5: only 3-3 pairs can be disjoint (3+3=6 > 5, so NO disjoint pairs!)
    # Wait: 3+3 = 6 > 5, so no two 3-cycles can be vertex-disjoint at n=5!
    # And a 5-cycle uses all vertices, so no disjoint pairs at all.
    alpha_2 = 0  # Always 0 at n=5

    # Verify H = 1 + 2*alpha_1 + 4*alpha_2
    H_check = 1 + 2 * alpha_1 + 4 * alpha_2

    by_score[scores].append({
        'H': H, 'c3_dir': c3_dir, 'c5_dir': c5_dir,
        'alpha_1': alpha_1, 'alpha_2': alpha_2,
        'H_check': H_check
    })

print(f"\nn={n}: alpha_2 = 0 always (no room for disjoint cycles)")
print(f"So H = 1 + 2*alpha_1 = 1 + 2*(c3_dir + c5_dir)")
print()

for sc in sorted(by_score.keys()):
    group = by_score[sc]
    Hs = set(d['H'] for d in group)
    c3s = set(d['c3_dir'] for d in group)
    c5s = set(d['c5_dir'] for d in group)
    a1s = set(d['alpha_1'] for d in group)
    checks = set(d['H_check'] for d in group)

    match = all(d['H'] == d['H_check'] for d in group)
    print(f"  Score {sc}: H={sorted(Hs)}, c3={sorted(c3s)}, c5={sorted(c5s)}, "
          f"alpha_1={sorted(a1s)}, H=1+2*a1? {match}")


# ========================================================================
# ANALYSIS 2: n=6 — DISJOINT PAIRS APPEAR
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: n=6 FULL DIRECTED CYCLE + DISJOINTNESS")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m

print(f"\nn={n}: 3+3=6=n, so disjoint 3-3 pairs become possible!")
print(f"But 3+5=8 > 6 and 5+5=10 > 6, so only 3-3 disjoint pairs exist.")
print()

by_score = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    # Directed 3-cycles
    c3_dir = 0
    c3_vsets = []
    for a, b, c in combinations(range(n), 3):
        cnt = count_directed_ham_cycles_subset(A, [a, b, c])
        c3_dir += cnt
        if cnt > 0:
            c3_vsets.append(frozenset([a, b, c]))

    # Directed 5-cycles (on each 5-subset)
    c5_dir = 0
    c5_vsets = []
    for subset in combinations(range(n), 5):
        cnt = count_directed_ham_cycles_subset(A, list(subset))
        c5_dir += cnt
        if cnt > 0:
            c5_vsets.append(frozenset(subset))

    # No directed 7-cycles (n=6 < 7)

    alpha_1 = c3_dir + c5_dir

    # Disjoint pairs: only 3-3 pairs (vertex-disjoint directed 3-cycles)
    # Each c3 vertex set has exactly 1 directed 3-cycle (from our counting method)
    # Two c3 vertex sets can be disjoint if they cover 6 vertices total
    # But the directed cycles themselves: if vset_i cap vset_j = empty,
    # then each directed cycle on vset_i is disjoint from each on vset_j.
    # Since each vset has exactly 1 directed cycle (from our counting), this gives 1 pair per disjoint vset pair.

    # Actually: from our enumerate, each 3-vertex set with a cycle gives 1 directed cycle.
    # But there are actually 2 directed 3-cycles per vertex set in a tournament!
    # Our count_directed_ham_cycles_subset counts from vertex 0 of the subset,
    # which gives 1 (since only one direction goes 0->next->next->0).
    # The other direction would be counted if we started from vertex 0 going backwards.

    # For the independence polynomial: we need to count ALL directed 3-cycles.
    # A 3-vertex tournament with a cycle has 2 directed cycles.
    # But wait, from the OCF verification, my code ALREADY gives the right answer.
    # Let me recheck: at n=5, alpha_1 = {5, 6, 7} and H = 1+2*alpha_1 = {11, 13, 15}.
    # From the output: H=11 has alpha_1=5, 1+2*5=11 YES.
    # H=15 has alpha_1=7 (for regular). nc=7 means 7 directed cycles total.
    # Regular T_5 has 5 three-cycles vertex sets (all C(5,3)=10 minus 5 transitive).
    # Wait no. C(5,3)=10, regular T_5 has c3 = 5*4*6/24 = 5 three-cycle vertex sets.
    # Each gives 1 directed cycle from our count, so c3_dir = 5.
    # Plus c5_dir = 2 (Hamiltonian cycles from vertex 0 in regular T_5).
    # alpha_1 = 5 + 2 = 7. H = 1 + 14 = 15. YES.

    # So for disjoint pair counting at n=6:
    # alpha_2 = number of pairs (i,j) where cycles[i] and cycles[j] are vertex-disjoint
    # Each cycle is represented by its vertex set and direction
    # Two cycles are disjoint iff their vertex sets are disjoint

    # Since each 3-vset gives 1 directed cycle, and 5-vsets give varying counts:
    # A disjoint pair of 3-3 type: 1*1 = 1 pair per pair of disjoint 3-vsets
    all_cycles = []  # (vertex_set, length)
    for vs in c3_vsets:
        all_cycles.append(vs)
    # For 5-cycles, each 5-vset contributes c5_count directed cycles
    for subset in combinations(range(n), 5):
        cnt = count_directed_ham_cycles_subset(A, list(subset))
        for _ in range(cnt):
            all_cycles.append(frozenset(subset))

    alpha_2 = 0
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if not (all_cycles[i] & all_cycles[j]):
                alpha_2 += 1

    # Higher alpha
    alpha_3 = 0
    # At n=6, max 3 disjoint 3-cycles would need 9 > 6 vertices, impossible
    # So alpha_3 = 0 always at n=6

    H_check = 1 + 2 * alpha_1 + 4 * alpha_2

    by_score[scores].append({
        'H': H, 'c3_dir': c3_dir, 'c5_dir': c5_dir,
        'alpha_1': alpha_1, 'alpha_2': alpha_2,
        'H_check': H_check
    })

# Check OCF
all_match = all(d['H'] == d['H_check'] for g in by_score.values() for d in g)
print(f"H = 1 + 2*alpha_1 + 4*alpha_2 for all n=6 tournaments? {all_match}")

if not all_match:
    # Find mismatches
    mismatch = 0
    for g in by_score.values():
        for d in g:
            if d['H'] != d['H_check']:
                mismatch += 1
                if mismatch <= 5:
                    print(f"  MISMATCH: H={d['H']}, check={d['H_check']}, "
                          f"a1={d['alpha_1']}, a2={d['alpha_2']}")
    print(f"Total mismatches: {mismatch}/{total}")

# Show constraint analysis within each score class
print(f"\nConstraint analysis: c5_dir + k*alpha_2 = const?")
print(f"{'Score':>25s} | {'c3':>4s} | {'c5 range':>12s} | {'a2 range':>12s} | {'c5+2a2':>10s} | {'H range':>10s}")
print(f"{'':->25s}-+-{'':->4s}-+-{'':->12s}-+-{'':->12s}-+-{'':->10s}-+-{'':->10s}")

for sc in sorted(by_score.keys()):
    group = by_score[sc]
    c3s = set(d['c3_dir'] for d in group)
    c5s = set(d['c5_dir'] for d in group)
    a2s = set(d['alpha_2'] for d in group)
    Hs = set(d['H'] for d in group)

    # Check c5 + 2*a2 = const
    combo = set(d['c5_dir'] + 2*d['alpha_2'] for d in group)

    c3_str = str(sorted(c3s)[0]) if len(c3s) == 1 else str(sorted(c3s))
    c5_str = f"{min(c5s)}-{max(c5s)}" if len(c5s) > 1 else str(sorted(c5s)[0])
    a2_str = f"{min(a2s)}-{max(a2s)}" if len(a2s) > 1 else str(sorted(a2s)[0])
    combo_str = str(sorted(combo)[0]) if len(combo) == 1 else str(sorted(combo))
    H_str = f"{min(Hs)}-{max(Hs)}" if len(Hs) > 1 else str(sorted(Hs)[0])

    const = len(combo) == 1
    marker = " *CONST*" if const and len(c5s) > 1 else ""
    print(f"  {str(sc):>23s} | {c3_str:>4s} | {c5_str:>12s} | {a2_str:>12s} | {combo_str:>10s} | {H_str:>10s}{marker}")

# Check: does H = f(c3_dir, c_n_dir) within each score class?
# n=6 is even, so c_n_dir doesn't exist (no n-cycles for odd n required by OCF)
# Actually, 6 is even so there are no directed 6-cycles in the OCF (only odd cycles).
# The longest odd cycle at n=6 is c5 (5-cycles).
print(f"\n\nDoes H = const + 2*c5_dir within score classes?")
for sc in sorted(by_score.keys()):
    group = by_score[sc]
    c5s = set(d['c5_dir'] for d in group)
    if len(c5s) <= 1:
        continue

    # Check H - 2*c5_dir = const
    residuals = set(d['H'] - 2*d['c5_dir'] for d in group)
    if len(residuals) == 1:
        print(f"  Score {sc}: YES! H - 2*c5_dir = {sorted(residuals)[0]}")
    else:
        print(f"  Score {sc}: NO. H - 2*c5_dir = {sorted(residuals)}")


# ========================================================================
# ANALYSIS 3: n=6 — DEEPER: THE FULL OCF RELATIONSHIP
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: n=6 FULL OCF CHECK")
print("=" * 70)

# At n=6: longest odd cycle is 5. These use 5 of 6 vertices.
# A 5-cycle is NOT Hamiltonian (n=6), so it doesn't conflict with everything.
# Specifically, a 5-cycle on {a,b,c,d,e} does NOT conflict with a 3-cycle on
# {f, ?, ?} if f is the excluded vertex and the other two are also not in {a,b,c,d,e}.
# But that's impossible: f is the only vertex outside, so any 3-cycle involving f
# must use two vertices from {a,b,c,d,e} and hence conflicts.
# So actually at n=6, a 5-cycle conflicts with ALL other cycles! Because any other
# cycle must use some vertex from the 5 used by the 5-cycle (only 1 vertex excluded,
# and a 3-cycle needs 3 vertices, at most 1 excluded, so at least 2 from the 5-cycle).

# Therefore at n=6, each 5-cycle also adds exactly 2 to H (same logic as Ham cycles).

# Check this:
for sc in sorted(by_score.keys()):
    group = by_score[sc]
    if len(set(d['c5_dir'] for d in group)) <= 1:
        continue

    # Partition by (c3_dir, alpha_2)
    subgroups = defaultdict(list)
    for d in group:
        subgroups[(d['c3_dir'], d['alpha_2'])].append(d)

    for key, sg in sorted(subgroups.items()):
        c5s = set(d['c5_dir'] for d in sg)
        Hs = set(d['H'] for d in sg)
        if len(c5s) <= 1:
            continue

        # Check if H - 2*c5_dir = const within this subgroup
        residuals = set(d['H'] - 2*d['c5_dir'] for d in sg)
        if len(residuals) == 1:
            print(f"  Score {sc}, c3={key[0]}, a2={key[1]}: H-2*c5 = {sorted(residuals)[0]} (const, {len(sg)} tours)")
        else:
            print(f"  Score {sc}, c3={key[0]}, a2={key[1]}: H-2*c5 = {sorted(residuals)} (NOT const, {len(sg)} tours)")


print("\n" + "=" * 70)
print("DONE.")
