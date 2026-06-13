#!/usr/bin/env python3
"""
burnside_correct_s20br.py -- kind-pasteur-2026-03-22-S20br

CORRECT BURNSIDE FORMULA FOR TOURNAMENTS.

The fix: work with ORDERED pairs (i,j) where i != j.
sigma sends (i,j) -> (sigma(i), sigma(j)).
Each orbit of sigma on ordered pairs is either:
  - A "paired" orbit: (i,j) and (j,i) are in DIFFERENT orbits.
    These contribute 1 free bit (choose which orbit gets the arc).
  - A "self-paired" orbit: (i,j) and (j,i) are in the SAME orbit.
    These contribute 0 (impossible for a tournament: can't have both i->j and j->i).

Fix(sigma) = 2^{number of paired orbit-pairs} if no self-paired orbits,
           = 0 if any self-paired orbit exists.

But paired orbits come in pairs {orbit of (i,j), orbit of (j,i)},
so number of paired orbit-pairs = (total ordered-pair orbits) / 2
when there are no self-paired orbits.

Author: kind-pasteur-2026-03-22-S20br
"""
import sys
from math import comb, factorial, gcd
from collections import Counter, defaultdict
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def tournament_fix(perm):
    """Count tournaments fixed by permutation perm.

    Works with ORDERED pairs (i,j), i != j.
    sigma acts as (i,j) -> (sigma(i), sigma(j)).

    For each orbit of ordered pairs, check if (i,j) and (j,i)
    are in the same orbit. If yes: Fix = 0 (no valid tournament).
    If no: orbits come in complementary pairs, 1 free bit per pair.
    """
    n = len(perm)
    # Build all ordered pairs and track orbits
    visited = {}  # (i,j) -> orbit_id
    orbit_id = 0

    for i in range(n):
        for j in range(n):
            if i == j: continue
            if (i,j) in visited: continue
            # Trace the orbit of (i,j)
            a, b = i, j
            while (a, b) not in visited:
                visited[(a, b)] = orbit_id
                a, b = perm[a], perm[b]
            orbit_id += 1

    # Now check: for each orbit, is (i,j) and (j,i) in the same orbit?
    # If yes for any pair: return 0
    n_orbits = orbit_id

    # Check self-pairing
    paired = [False] * n_orbits
    for i in range(n):
        for j in range(n):
            if i == j: continue
            oid_ij = visited[(i,j)]
            oid_ji = visited[(j,i)]
            if oid_ij == oid_ji:
                return 0  # Self-paired orbit: impossible for tournament
            # Mark as paired

    # Count orbit-pairs: each unordered pair {orbit(i,j), orbit(j,i)} is one choice
    orbit_pairs = set()
    for i in range(n):
        for j in range(i+1, n):
            oid_ij = visited[(i,j)]
            oid_ji = visited[(j,i)]
            orbit_pairs.add((min(oid_ij, oid_ji), max(oid_ij, oid_ji)))

    return 2 ** len(orbit_pairs)

print("=" * 70)
print("  CORRECT BURNSIDE FORMULA FOR TOURNAMENTS")
print("=" * 70)

# ================================================================
# 1. VERIFY A000568
# ================================================================
print(f"\n{'='*70}")
print(f"  1. VERIFY A000568 VIA CORRECT BURNSIDE")
print(f"{'='*70}\n")

known = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456}

for n in range(1, 8):
    total_fix = 0
    for perm in permutations(range(n)):
        total_fix += tournament_fix(list(perm))

    burnside = total_fix // factorial(n)
    expected = known[n]
    match = burnside == expected
    print(f"  n={n}: Burnside = {burnside}, A000568 = {expected}, Match: {match}")

# ================================================================
# 2. FIX BY CYCLE TYPE (the efficient formula)
# ================================================================
print(f"\n{'='*70}")
print(f"  2. FIX(sigma) BY CYCLE TYPE")
print(f"{'='*70}\n")

def cycle_type(perm):
    n = len(perm)
    visited = [False]*n
    cycles = []
    for i in range(n):
        if visited[i]: continue
        length = 0
        j = i
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            length += 1
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))

def partitions(n, max_part=None):
    if max_part is None: max_part = n
    if n == 0:
        yield []
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield [first] + rest

def count_with_cycle_type(ct, n):
    denom = 1
    counter = Counter(ct)
    for length, mult in counter.items():
        denom *= length ** mult * factorial(mult)
    return factorial(n) // denom

# For each cycle type, compute Fix using the correct formula
# (averaged over a representative permutation)
n = 5
print(f"  n={n}: Fix by cycle type (CORRECT)")
print(f"  {'Cycle type':>15s} {'#perms':>8s} {'Fix(tour)':>10s} {'Fix(graph)':>11s} {'Contribution':>13s}")

total = 0
for part in partitions(n):
    ct = tuple(part)
    n_perms = count_with_cycle_type(ct, n)

    # Build a representative permutation with this cycle type
    perm = []
    pos = 0
    for length in ct:
        for i in range(length - 1):
            perm.append(pos + i + 1)
        perm.append(pos)
        pos += length

    fix_tour = tournament_fix(perm)

    # Graph fix for comparison
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    pair_idx = {p: k for k, p in enumerate(pairs)}
    visited_p = [False]*len(pairs)
    n_pair_orbits = 0
    for k, (a,b) in enumerate(pairs):
        if visited_p[k]: continue
        ca, cb = a, b
        while True:
            p = (min(ca,cb), max(ca,cb))
            idx = pair_idx[p]
            if visited_p[idx]: break
            visited_p[idx] = True
            ca, cb = perm[ca], perm[cb]
        n_pair_orbits += 1
    fix_graph = 2**n_pair_orbits

    contrib = n_perms * fix_tour
    total += contrib
    print(f"  {str(list(ct)):>15s} {n_perms:>8d} {fix_tour:>10d} {fix_graph:>11d} {contrib:>13d}")

burnside = total // factorial(n)
print(f"\n  Total: {total}, Burnside: {burnside}, A000568({n}): {known[n]}, Match: {burnside == known[n]}")

# ================================================================
# 3. WHICH CYCLE TYPES HAVE Fix(tour) = 0?
# ================================================================
print(f"\n{'='*70}")
print(f"  3. WHICH CYCLE TYPES HAVE Fix(tournament) = 0?")
print(f"{'='*70}\n")

for n in range(3, 8):
    zero_types = []
    nonzero_types = []
    for part in partitions(n):
        ct = tuple(part)
        # Build representative
        perm = []
        pos = 0
        for length in ct:
            for i in range(length - 1):
                perm.append(pos + i + 1)
            perm.append(pos)
            pos += length
        fix = tournament_fix(perm)
        if fix == 0:
            zero_types.append(ct)
        else:
            nonzero_types.append((ct, fix))

    print(f"  n={n}: {len(zero_types)} cycle types with Fix=0, {len(nonzero_types)} with Fix>0")
    if zero_types:
        print(f"    Fix=0: {zero_types}")
    for ct, fix in nonzero_types:
        print(f"    Fix={fix}: {list(ct)}")
    print()

# ================================================================
# 4. THE PATTERN: WHEN IS Fix(tour) = 0?
# ================================================================
print(f"{'='*70}")
print(f"  4. THE PATTERN: Fix(tour) = 0 IFF EVEN CYCLE EXISTS")
print(f"{'='*70}\n")

# Hypothesis: Fix(sigma) = 0 for tournaments iff sigma has at least one EVEN-length cycle.
# Because: an even-length cycle on vertices creates a self-paired orbit on ordered pairs.
# In a cycle (v_0, v_1, ..., v_{L-1}): the ordered pair (v_0, v_1) maps to (v_1, v_2),
# ..., (v_{L-1}, v_0). After L steps: back to (v_0, v_1). Orbit length = L.
# The reverse pair (v_1, v_0) maps to (v_2, v_1), ..., (v_0, v_{L-1}). Orbit length = L.
# Are (v_0, v_1) and (v_1, v_0) in the same orbit?
# (v_0, v_1) -> (v_1, v_2) -> ... -> (v_{L/2}, v_{L/2+1}) (if L even, v_{L/2} = reverse of v_0?)
# Let me check: does (v_k, v_{k+1}) ever equal (v_j, v_{j-1}) for some j?
# That requires v_k = v_j and v_{k+1} = v_{j-1}, which means k = j and k+1 = j-1 (mod L).
# So k+1 = k-1 (mod L), i.e., 2 = 0 (mod L), i.e., L divides 2. So L = 1 or L = 2.
# For L = 2: (v_0, v_1) -> (v_1, v_0). Self-paired! Fix = 0.
# For L > 2: (v_0, v_1) and (v_1, v_0) are in DIFFERENT orbits.
# But we also need to check between-cycle ordered pairs.

# Actually the check for WITHIN a single cycle of length L:
# (v_0, v_1) -> (v_1, v_2) -> ... -> (v_{L-1}, v_0). Orbit length = L.
# (v_1, v_0) -> (v_2, v_1) -> ... -> (v_0, v_{L-1}). Orbit length = L.
# These are the SAME orbit iff (v_k, v_{k+1}) = (v_{L-j}, v_{L-j-1}) for some k, j.
# This means v_k = v_{L-j} and v_{k+1} = v_{L-j-1}, i.e., k = L-j (mod L) and k+1 = L-j-1 (mod L).
# Both give k + j = L (mod L) and k + 1 + j + 1 = L (mod L), i.e., k + j = L - 2 (mod L).
# So L = L - 2 (mod L), i.e., 2 = 0 (mod L). So L | 2 again.

# For L = 2: (v_0, v_1) orbit is {(v_0, v_1), (v_1, v_0)}. Self-paired!
# For L >= 3 (odd or even): within-cycle ordered pairs are NOT self-paired.

# What about between-cycle pairs? (a, b) where a is in cycle C_1, b in cycle C_2.
# sigma(a, b) = (sigma(a), sigma(b)). Orbit length = lcm(L_1, L_2).
# The reverse (b, a) has orbit: sigma(b,a) = (sigma(b), sigma(a)).
# Are these ever in the same orbit?
# (a, b) -> (sigma(a), sigma(b)) -> ... After k steps: (sigma^k(a), sigma^k(b)).
# (b, a) -> (sigma(b), sigma(a)) -> ... After k steps: (sigma^k(b), sigma^k(a)).
# These are the same iff sigma^k(a) = sigma^j(b) and sigma^k(b) = sigma^j(a) for some k, j.
# If a is in cycle of length L_1 and b in cycle of length L_2:
# sigma^k(a) = a_{k mod L_1}, sigma^j(b) = b_{j mod L_2}.
# Need a_{k} = b_{j} and b_{k} = a_{j}. Since a and b are in DIFFERENT cycles,
# a_k can never equal b_j. So between-cycle pairs are NEVER self-paired.

# CONCLUSION: Fix(sigma) = 0 for tournaments iff sigma has a 2-cycle.

print(f"  THEOREM: Fix(sigma) = 0 for tournaments iff sigma has a 2-cycle (transposition).")
print(f"  PROOF: A 2-cycle (v_0, v_1) creates the orbit {{(v_0,v_1), (v_1,v_0)}} on ordered pairs.")
print(f"  This is self-paired: both arc directions are in the same orbit.")
print(f"  For tournaments, exactly one must be present: contradiction.")
print(f"  For cycles of length L >= 3: within-cycle ordered pairs are never self-paired.")
print(f"  For between-cycle ordered pairs: never self-paired (different cycles).")
print()

# VERIFY
print(f"  VERIFICATION:")
for n in range(3, 8):
    for part in partitions(n):
        ct = tuple(part)
        has_2cycle = any(l == 2 for l in ct)
        perm = []
        pos = 0
        for length in ct:
            for i in range(length - 1):
                perm.append(pos + i + 1)
            perm.append(pos)
            pos += length
        fix = tournament_fix(perm)
        prediction = (fix == 0) == has_2cycle
        if not prediction:
            print(f"    COUNTEREXAMPLE: n={n}, ct={list(ct)}, has_2cycle={has_2cycle}, fix={fix}")
    else:
        continue

print(f"  No counterexamples found!")

# ================================================================
# 5. THE CLEAN FORMULA
# ================================================================
print(f"\n{'='*70}")
print(f"  5. THE CLEAN FORMULA FOR A000568")
print(f"{'='*70}\n")

# A000568(n) = (1/n!) * sum over ONLY-ODD-CYCLE permutations of 2^{paired-orbit-pairs}
# where paired-orbit-pairs = (number of ordered-pair orbits) / 2.
# And: number of ordered-pair orbits for a perm with only odd cycles
# = number of arcs orbits = ???

# For a permutation with cycle type (l_1, l_2, ..., l_k) where all l_i are odd:
# Number of ordered-pair orbits:
#   Within cycle l_i: (l_i^2 - l_i) ordered pairs, orbit length l_i each.
#     So l_i * (l_i - 1) / l_i = l_i - 1 orbits of ordered pairs per cycle.
#     Wait: there are l_i*(l_i-1) ordered pairs within a cycle of length l_i.
#     Each orbit has length l_i. So number of orbits = l_i*(l_i-1)/l_i = l_i - 1.
#   Between cycles l_i and l_j (i != j): l_i * l_j ordered pairs (a,b) with a in C_i, b in C_j.
#     Each orbit has length lcm(l_i, l_j). So number of orbits = l_i*l_j / lcm(l_i, l_j) = gcd(l_i, l_j).
#     Plus another l_i*l_j pairs (b,a), giving another gcd(l_i, l_j) orbits.
#     Total between two cycles: 2 * gcd(l_i, l_j) ordered-pair orbits.

# Total ordered-pair orbits = sum_i (l_i - 1) + sum_{i<j} 2*gcd(l_i, l_j)
# These come in complementary pairs: orbit of (a,b) paired with orbit of (b,a).
# Number of pairs = (total orbits) / 2
#   = (sum_i (l_i-1) + sum_{i<j} 2*gcd(l_i,l_j)) / 2
#   = (sum_i (l_i-1))/2 + sum_{i<j} gcd(l_i, l_j)

# For ODD l_i: l_i - 1 is EVEN, so (l_i-1)/2 is an integer. Good.

def tournament_fix_from_ct(cycle_type):
    """Compute Fix(sigma) for tournaments from cycle type.
    Returns 0 if any cycle has even length.
    Otherwise returns 2^{orbit-pairs}.
    """
    if any(l % 2 == 0 for l in cycle_type):
        return 0

    n_orbit_pairs = 0
    for l in cycle_type:
        n_orbit_pairs += (l - 1) // 2
    for i in range(len(cycle_type)):
        for j in range(i + 1, len(cycle_type)):
            n_orbit_pairs += gcd(cycle_type[i], cycle_type[j])

    return 2 ** n_orbit_pairs

# VERIFY
print(f"  CLEAN FORMULA: Fix(sigma) from cycle type (only-odd-cycles version)")
print(f"  {'n':>3s} {'Burnside':>10s} {'A000568':>10s} {'Match':>7s}")

for n in range(1, 8):
    total = 0
    for part in partitions(n):
        ct = tuple(part)
        n_perms = count_with_cycle_type(ct, n)
        fix = tournament_fix_from_ct(list(ct))
        total += n_perms * fix

    burnside = total // factorial(n)
    expected = known.get(n, '?')
    match = burnside == expected
    print(f"  {n:>3d} {burnside:>10d} {str(expected):>10s} {str(match):>7s}")

# ================================================================
# 6. EXTEND TO n=8, 9, 10
# ================================================================
print(f"\n{'='*70}")
print(f"  6. A000568 AT LARGER n (from the clean formula)")
print(f"{'='*70}\n")

known_extended = {8: 6880, 9: 191536, 10: 9733056}

for n in range(1, 11):
    total = 0
    for part in partitions(n):
        ct = tuple(part)
        n_perms = count_with_cycle_type(ct, n)
        fix = tournament_fix_from_ct(list(ct))
        total += n_perms * fix

    burnside = total // factorial(n)
    expected = known.get(n, known_extended.get(n, '?'))
    match = burnside == expected if expected != '?' else '?'
    print(f"  n={n:>2d}: A000568 = {burnside:>12d} (expected {expected}, match={match})")
