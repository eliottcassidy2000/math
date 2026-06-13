"""
forbidden7_exhaustive.py — opus-2026-04-04-S9

Exhaustively check ALL tournaments at n=7 with H divisible by 7.
For each: compute |Aut| and check if H/|Aut| = 7 ever.

Also: check the complete picture of when H ≡ 0 mod 7 first appears.
"""
import numpy as np
from itertools import permutations
from collections import Counter
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, make_tiles, tiling_to_tournament, count_hamiltonian_paths

def compute_aut_size(T, n):
    """Compute |Aut(T)| by checking all permutations."""
    count = 0
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[j]] == T[i][j] for i in range(n) for j in range(n)):
            count += 1
    return count

n = 7
tiles = make_tiles(n)
m = len(tiles)
N = 1 << m
print(f"n={n}, m={m}, N={N}")
print(f"Computing H for all tilings...")
H, _, _ = compute_all_H(n)

# Find all tilings with H divisible by 7
div7_tilings = [t for t in range(N) if int(H[t]) % 7 == 0]
print(f"Tilings with H divisible by 7: {len(div7_tilings)}")

# Group by H value
from collections import defaultdict
by_H = defaultdict(list)
for t in div7_tilings:
    by_H[int(H[t])].append(t)

print(f"\nH values divisible by 7: {sorted(by_H.keys())}")

# For each H value divisible by 7, check ALL tilings for |Aut|
tc7_found = False
tc21_found = False

for h in sorted(by_H.keys()):
    tilings = by_H[h]
    aut_dist = Counter()

    for t in tilings:
        T = tiling_to_tournament(n, tiles, t)
        aut = compute_aut_size(T, n)
        aut_dist[aut] += 1

        # Check if tiling count = 7 or 21
        tc = h // aut
        if tc == 7:
            tc7_found = True
            print(f"  *** TILING COUNT = 7 FOUND: H={h}, |Aut|={aut}, t={t}")
        if tc == 21:
            tc21_found = True
            print(f"  *** TILING COUNT = 21 FOUND: H={h}, |Aut|={aut}, t={t}")

    print(f"  H={h}: {len(tilings)} tilings, |Aut| dist = {dict(sorted(aut_dist.items()))}")

    # What tiling counts are achieved?
    achieved_tc = set()
    for aut, count in aut_dist.items():
        achieved_tc.add(h // aut)
    print(f"    Achieved tiling counts: {sorted(achieved_tc)}")

print(f"\n{'='*60}")
print(f"VERDICT: tiling count = 7 at n=7: {'FOUND' if tc7_found else 'NOT FOUND'}")
print(f"VERDICT: tiling count = 21 at n=7: {'FOUND' if tc21_found else 'NOT FOUND'}")

# Also: enumerate ALL tiling counts at n=7 via grouping by iso class
# Use a faster method: group by (H, score_sequence) then check Aut within groups
print(f"\n\nComputing ALL iso class tiling counts at n=7...")
print(f"(Using H + score_seq as iso class proxy, with Aut refinement)")

# Group by (H, sorted_scores)
groups = defaultdict(list)
for t in range(N):
    T = tiling_to_tournament(n, tiles, t)
    h = int(H[t])
    scores = tuple(sorted(int(sum(T[v])) for v in range(n)))
    groups[(h, scores)].append(t)

print(f"Groups (H, scores): {len(groups)}")

# Within each group, further split by canonical form
# For groups with size > 1, we need to check if tilings are truly isomorphic
# For now, use H/|Aut| for one representative per group
all_tiling_counts = []
for (h, scores), tilings in groups.items():
    # Compute |Aut| for ONE representative
    T = tiling_to_tournament(n, tiles, tilings[0])
    aut = compute_aut_size(T, n)
    tc = h // aut

    # The group might contain multiple iso classes with the same (H, scores)
    # For now, assume each group is one iso class (mostly true for distinct scores)
    # Number of tilings in this group should equal tc * (multiplicity of iso class copies)
    # If group size != tc: multiple iso classes in this group

    actual_count = len(tilings)
    expected_if_one_class = tc

    if actual_count == expected_if_one_class:
        all_tiling_counts.append(tc)
    else:
        # Multiple iso classes - check a second representative
        for t2 in tilings[1:2]:
            T2 = tiling_to_tournament(n, tiles, t2)
            aut2 = compute_aut_size(T2, n)
            tc2 = h // aut2
            if tc2 != tc:
                all_tiling_counts.append(tc)
                all_tiling_counts.append(tc2)
                break
        else:
            # Assume multiple copies of same class
            num_classes = actual_count // tc if tc > 0 else 1
            for _ in range(num_classes):
                all_tiling_counts.append(tc)

tc_set = sorted(set(all_tiling_counts))
print(f"\nApproximate tiling count set at n=7: {tc_set[:30]}...")
print(f"Total distinct values: {len(tc_set)}")
print(f"\nIs 7 in the set? {7 in tc_set}")
print(f"Is 21 in the set? {21 in tc_set}")

# Check all odd values up to max
max_tc = max(tc_set) if tc_set else 0
missing_odd = [k for k in range(1, min(max_tc+1, 200), 2) if k not in tc_set]
print(f"\nMissing odd tiling counts up to {min(max_tc, 200)}: {missing_odd[:30]}")
