#!/usr/bin/env python3
"""
Prove alpha_1=3 is impossible for all tournaments.

alpha_1 = total number of directed odd cycles in T.
At n<=4, max alpha_1 = 2. alpha_1=3 would need n>=5.
At n=5: alpha_1 = c3 + c5 where c3 = cyclic triples, c5 = directed 5-cycles.

For alpha_1=3: need c3 + c5 = 3.
Options: (c3=3,c5=0), (c3=2,c5=1), (c3=1,c5=2), (c3=0,c5=3)

Key constraint: c5 counts directed Hamiltonian cycles on 5-vertex subsets.
For n=5: c5 = directed Ham cycles of the whole tournament.
Score (2,2,2,2,2) regular: always c5=2.
Score (1,2,2,2,3) or similar: c5 can be 0 or other values.

For n>=6: alpha_1 = sum over all odd k of c_k. The gap at 3 might be explained
by some parity or modular constraint.

kind-pasteur-2026-03-06-S21
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import combinations

def count_cycles_detailed(T, n):
    """Return dict: length -> count of directed cycles of that length."""
    result = {}
    for k in range(3, n+1, 2):
        total = 0
        for combo in combinations(range(n), k):
            verts = list(combo)
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for vi in range(k):
                    if not (mask & (1 << vi)):
                        continue
                    c = dp.get((mask, vi), 0)
                    if c == 0:
                        continue
                    for ui in range(k):
                        if mask & (1 << ui):
                            continue
                        if T[verts[vi]][verts[ui]]:
                            key = (mask | (1 << ui), ui)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << k) - 1
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and T[verts[vi]][verts[0]]:
                    total += c
        result[k] = total
    return result

# At n=5: exhaustive search for alpha_1=3
n = 5
m = n*(n-1)//2
print(f"n=5: Searching for alpha_1=3...")
closest = []
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles = count_cycles_detailed(T, n)
    c3, c5 = cycles[3], cycles[5]
    a1 = c3 + c5
    if a1 == 3:
        print(f"  FOUND! bits={bits}, c3={c3}, c5={c5}")
    # Track tournaments with c3+c5 near 3
    if 2 <= a1 <= 4:
        closest.append((a1, c3, c5, bits))

print(f"\nTournaments with alpha_1 near 3:")
for a1, c3, c5, bits in sorted(closest):
    print(f"  alpha_1={a1}: c3={c3}, c5={c5}, bits={bits}")

# The constraint: at n=5, c5 is always 0 or 2 (for the single 5-vertex set)
# Actually: c5 = number of directed Hamiltonian cycles on all 5 vertices
# For a regular tournament: c5=2 always
# For non-regular: can be 0 or other
print(f"\nAll (c3, c5) pairs at n=5:")
pairs = set()
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles = count_cycles_detailed(T, n)
    pairs.add((cycles[3], cycles[5]))
for c3, c5 in sorted(pairs):
    print(f"  c3={c3}, c5={c5}, alpha_1={c3+c5}")

# At n=6: check (c3, c5) pairs where alpha_1=3
print(f"\nn=6: Searching for alpha_1=3...")
n = 6
m = n*(n-1)//2
found = False
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles = count_cycles_detailed(T, n)
    a1 = sum(cycles.values())
    if a1 == 3:
        found = True
        print(f"  FOUND! bits={bits}, cycles={cycles}")
        break
if not found:
    print(f"  NOT FOUND at n=6")

# Theoretical analysis: why can't alpha_1=3?
# At n=5: c5 is always even (0 or 2). So alpha_1 = c3 + c5 = c3 + even.
# c3 can be 0,1,2,3,4,5 at n=5.
# c5=0: alpha_1 = c3, achievable: 0,1,2,3,4,5
# c5=2: alpha_1 = c3 + 2, achievable: 2,3,4,5,6,7
# Wait - c5=2 with c3=1 would give alpha_1=3!
# But does (c3=1, c5=2) actually occur?

print(f"\nn=5: Does (c3=1, c5=2) occur?")
n = 5
m = n*(n-1)//2
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles = count_cycles_detailed(T, n)
    if cycles[3] == 1 and cycles[5] == 2:
        print(f"  YES! bits={bits}")
        break
else:
    print(f"  NO - (c3=1, c5=2) never occurs!")

print(f"\n  Does (c3=3, c5=0) occur?")
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles = count_cycles_detailed(T, n)
    if cycles[3] == 3 and cycles[5] == 0:
        print(f"  YES! bits={bits}")
        break
else:
    print(f"  NO - (c3=3, c5=0) never occurs!")

# So the gap is structural: certain (c3, c5) combinations never occur
# Let's check which (c3, c5) pairs are achievable
print(f"\nn=5: Which (c3, c5) pairs occur and how many tournaments for each?")
pair_counts = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles = count_cycles_detailed(T, n)
    pair = (cycles[3], cycles[5])
    pair_counts[pair] = pair_counts.get(pair, 0) + 1
for pair in sorted(pair_counts.keys()):
    c3, c5 = pair
    print(f"  (c3={c3}, c5={c5}): {pair_counts[pair]} tournaments, alpha_1={c3+c5}")

print("\nDone.")
