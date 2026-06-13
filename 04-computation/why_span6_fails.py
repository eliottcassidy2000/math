#!/usr/bin/env python3
"""
Why can't 3 cyclic triples span 6 vertices with c3=3?

For each of the 480 3-triple configurations on 6 vertices,
check: what is the MINIMUM c3 achievable?

kind-pasteur-2026-03-06
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations

def get_cyclic_triples(T, n):
    result = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if T[a][b] and T[b][c] and T[c][a]:
            result.append(frozenset(combo))
        elif T[a][c] and T[c][b] and T[b][a]:
            result.append(frozenset(combo))
    return result

# Classify 3-triple configs by their intersection pattern
all_triples_6 = list(combinations(range(6), 3))

# Group by intersection pattern
from collections import Counter

pattern_min_c3 = {}

for idx_combo in combinations(range(len(all_triples_6)), 3):
    trips = [frozenset(all_triples_6[i]) for i in idx_combo]
    union = trips[0] | trips[1] | trips[2]
    if len(union) != 6:
        continue

    # Intersection pattern
    pw = []
    for i in range(3):
        for j in range(i+1, 3):
            pw.append(len(trips[i] & trips[j]))
    pattern = tuple(sorted(pw))

    # Find minimum c3 among tournaments containing all 3 as cyclic
    # Try all orientations of the 3 triples
    min_c3 = 100
    for orient in range(8):
        edges = {}
        conflict = False

        for tidx, trip in enumerate(trips):
            trip_sorted = sorted(trip)
            a, b, c = trip_sorted
            if orient & (1 << tidx):
                req = [(a,b),(b,c),(c,a)]
            else:
                req = [(a,c),(c,b),(b,a)]
            for (i,j) in req:
                if (i,j) in edges and edges[(i,j)] != 1:
                    conflict = True
                    break
                if (j,i) in edges and edges[(j,i)] != 0:
                    conflict = True
                    break
                edges[(i,j)] = 1
                edges[(j,i)] = 0
            if conflict:
                break
        if conflict:
            continue

        # Free edges
        all_pairs = [(i,j) for i in range(6) for j in range(i+1,6)]
        free = []
        for (i,j) in all_pairs:
            if (i,j) not in edges and (j,i) not in edges:
                free.append((i,j))

        for fmask in range(1 << len(free)):
            T = [[0]*6 for _ in range(6)]
            for (i,j), v in edges.items():
                T[i][j] = v
            for fidx, (i,j) in enumerate(free):
                if fmask & (1 << fidx):
                    T[i][j] = 1
                else:
                    T[j][i] = 1

            ok = True
            for i in range(6):
                for j in range(i+1,6):
                    if T[i][j] + T[j][i] != 1:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue

            triples_found = get_cyclic_triples(T, 6)
            c3 = len(triples_found)

            # Check that our 3 target triples are all present
            target_present = all(t in triples_found for t in trips)
            if target_present and c3 < min_c3:
                min_c3 = c3

    if pattern not in pattern_min_c3 or min_c3 < pattern_min_c3[pattern]:
        pattern_min_c3[pattern] = min_c3

print("Intersection pattern -> minimum c3 when all 3 triples are present:")
for pattern, mc3 in sorted(pattern_min_c3.items()):
    print(f"  {pattern}: min c3 = {mc3}")
    if mc3 <= 3:
        print(f"    ^^ This could give c3=3!")

print()
if all(mc3 >= 4 for mc3 in pattern_min_c3.values()):
    print("ALL patterns force c3 >= 4 when spanning 6 vertices!")
    print("This proves: at n=6, any 3 cyclic triples spanning 6 vertices")
    print("inevitably create at least 1 additional cyclic triple.")
else:
    print("Some patterns allow c3=3 -- need further investigation.")

print("\nDone.")
