#!/usr/bin/env python3
"""
root_profile_5cycle_s18.py — kind-pasteur-2026-03-21-S18

At n=6, the 3-cycle root profile has 3 failures (H gap = 4 in each).
The gap of 4 = the contribution of one alpha_2 in OCF.

Hypothesis: adding 5-cycle information resolves the ambiguity.
Test: (3-cycle profile, alpha_2) determines H at n=6.
Also: (3-cycle profile, c5) determines H at n=6.
Also: does FULL odd-cycle root profile (all lengths) determine H?
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def directed_3cycle_profile(A, n):
    """Sorted profile: for each arc, count 3-cycles through it."""
    profile = []
    for i in range(n):
        for j in range(i+1, n):
            c = 0
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                if A[j][i] and A[i][k] and A[k][j]:
                    c += 1
            profile.append(c)
    return tuple(sorted(profile))

def count_directed_odd_cycles(A, n):
    """Count directed odd cycles by vertex set, return dict: frozenset -> count of directed cycles."""
    cycle_sets = defaultdict(int)
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            sub = list(subset)
            count = 0
            for perm in permutations(sub[1:]):
                ordering = [sub[0]] + list(perm)
                is_cycle = True
                for idx in range(length):
                    if not A[ordering[idx]][ordering[(idx+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    count += 1
            if count > 0:
                cycle_sets[frozenset(subset)] = count
    return cycle_sets

def compute_alpha(A, n):
    """Compute alpha_k (independence polynomial coefficients of conflict graph)."""
    cycles = count_directed_odd_cycles(A, n)
    cycle_list = list(cycles.keys())
    nc = len(cycle_list)

    # Build conflict adjacency
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycle_list[i] & cycle_list[j]:
                adj[i][j] = adj[j][i] = True

    # Count independent sets by size
    alpha = defaultdict(int)
    alpha[0] = 1
    for mask in range(1, 1 << nc):
        verts = [v for v in range(nc) if (mask >> v) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            alpha[len(verts)] += 1

    return dict(alpha)

print("=" * 72)
print("  ROOT PROFILE + 5-CYCLE RESOLUTION TEST (n=6)")
print("  kind-pasteur-2026-03-21-S18")
print("=" * 72)

n = 6
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Store data for the 3 failure profiles
failure_profiles = set()

# First pass: identify failure profiles
profile_to_H = defaultdict(set)
for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    H = count_hp(A, n)
    prof = directed_3cycle_profile(A, n)
    profile_to_H[prof].add(H)

for prof, Hs in profile_to_H.items():
    if len(Hs) > 1:
        failure_profiles.add(prof)

print(f"\n  Failure profiles: {len(failure_profiles)}")
for fp in sorted(failure_profiles):
    print(f"    {list(fp)} -> H in {sorted(profile_to_H[fp])}")

# Second pass: for failure profiles, compute alpha and c5
print(f"\n  Detailed analysis of failure cases:")
print()

# For each failure profile, analyze OCF structure
for fp in sorted(failure_profiles):
    print(f"  Profile {list(fp)}, H in {sorted(profile_to_H[fp])}:")

    # Find representative tournaments for each H value
    for target_H in sorted(profile_to_H[fp]):
        found = False
        for bits in range(2**m):
            A = np.zeros((n, n), dtype=int)
            for k, (i, j) in enumerate(pairs):
                if (bits >> k) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
            H = count_hp(A, n)
            if H != target_H:
                continue
            prof = directed_3cycle_profile(A, n)
            if prof != fp:
                continue

            # Found representative
            alpha = compute_alpha(A, n)

            # Count cycles by length
            cycles = count_directed_odd_cycles(A, n)
            c3 = sum(1 for s in cycles if len(s) == 3)
            c5 = sum(1 for s in cycles if len(s) == 5)

            # Verify OCF
            H_ocf = sum(alpha.get(k, 0) * (2**k) for k in range(10))

            print(f"    H={H}: alpha={dict(sorted(alpha.items()))}, c3={c3}, c5={c5}, OCF check: {H_ocf == H}")

            found = True
            break

        if not found:
            print(f"    H={target_H}: (no representative found)")
    print()

# Third pass: test (profile, c5) and (profile, alpha_2) resolution
print("=" * 72)
print("  TEST: (profile, c5) DETERMINES H?")
print("=" * 72)

profile_c5_to_H = defaultdict(set)
profile_a2_to_H = defaultdict(set)

processed = 0
for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1

    prof = directed_3cycle_profile(A, n)

    # Only compute expensive stuff for failure profiles (save time)
    # Actually, compute for ALL to check full resolution
    H = count_hp(A, n)

    # Count 5-cycles
    c5 = 0
    for subset in combinations(range(n), 5):
        sub = list(subset)
        has_cycle = False
        for perm in permutations(sub[1:]):
            ordering = [sub[0]] + list(perm)
            is_cycle = True
            for idx in range(5):
                if not A[ordering[idx]][ordering[(idx+1) % 5]]:
                    is_cycle = False
                    break
            if is_cycle:
                has_cycle = True
                break
        if has_cycle:
            c5 += 1

    profile_c5_to_H[(prof, c5)].add(H)

    processed += 1
    if processed % 100000 == 0:
        fails = sum(1 for v in profile_c5_to_H.values() if len(v) > 1)
        print(f"  ... {processed}/{2**m}: {len(profile_c5_to_H)} classes, {fails} failures")

determines_pc5 = all(len(v) == 1 for v in profile_c5_to_H.values())
print(f"\n  (profile, c5) classes: {len(profile_c5_to_H)}")
print(f"  Determines H: {determines_pc5}")

if not determines_pc5:
    fails = [(k, v) for k, v in profile_c5_to_H.items() if len(v) > 1]
    print(f"  Remaining failures: {len(fails)}")
    for (prof, c5_val), Hs in fails[:10]:
        print(f"    profile_sum={sum(prof)}, c5={c5_val}: H in {sorted(Hs)}")
else:
    print(f"  SUCCESS: 3-cycle root profile + c5 determines H at n=6!")

# Check: does alpha_1 alone (= c3) determine the profile?
print(f"\n  Bonus: does score sequence determine 3-cycle root profile?")
score_to_profiles = defaultdict(set)
for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    scores = tuple(sorted(A.sum(axis=1), reverse=True))
    prof = directed_3cycle_profile(A, n)
    score_to_profiles[scores].add(prof)

multi_score = sum(1 for v in score_to_profiles.values() if len(v) > 1)
print(f"  Score sequences with multiple profiles: {multi_score} / {len(score_to_profiles)}")

print("\n" + "=" * 72)
print("  SYNTHESIS")
print("=" * 72)
print("""
  At n=5: 3-cycle root profile determines H. (THEOREM)
  At n=6: 3-cycle root profile has 3 failures (H gap = 4 each).
          The gap is EXACTLY one alpha_2 contribution.
          Adding c5 (5-cycle count) should resolve.
  At n=7: 3-cycle root profile fails heavily (51/149).
          Need full odd-cycle structure.

  ROOT SYSTEM INTERPRETATION:
  The 3-cycle root profile encodes how the signed roots of A_{n-1}
  participate in triangular cycles. This captures alpha_1 perfectly
  but misses alpha_2 (vertex-disjoint cycle pairs) which first require
  5-cycles or 6+ vertices.

  The full root profile (all odd cycle lengths) should determine H
  via OCF, since H = I(Omega, 2) depends ONLY on the independence
  polynomial of the conflict graph, which is determined by the
  full cycle structure.
""")
