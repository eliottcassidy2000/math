#!/usr/bin/env python3
"""
Verify the (a1=12, a2=6, a3=1) case at n=9.
Does this triple actually occur? If so, real-rootedness FAILS at n=9.

Also: is (12, 6, 1) achievable? Turán says a2 <= a1^2/4 = 36, so a2=6 is fine.
But a1^2 - 4*a2 = 144 - 24 = 120 > 0, so C > 0 by Turán.
disc = 18*12*6*1 - 4*12^3*1 + 12^2*6^2 - 4*6^3 - 27*1^2
     = 1296 - 6912 + 5184 - 864 - 27 = -1323. Indeed negative!

But we need to check: does a tournament with exactly a1=12, a2=6, a3=1 exist?

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament
import numpy as np

def find_3cycles(T):
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i,j,k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i,j,k))
    return cycles

def count_a2_a3(cycles):
    m = len(cycles)
    cycle_sets = [frozenset(c) for c in cycles]
    disjoint_after = [[] for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if not (cycle_sets[i] & cycle_sets[j]):
                disjoint_after[i].append(j)
    a2 = sum(len(d) for d in disjoint_after)
    a3 = 0
    for i in range(m):
        ui = cycle_sets[i]
        for j in disjoint_after[i]:
            uij = ui | cycle_sets[j]
            for k in disjoint_after[j]:
                if k > j and not (cycle_sets[k] & uij):
                    a3 += 1
    return a2, a3

n = 9

# Search for tournaments with a1=12, a2=6, a3=1
print("Searching for tournament with (a1=12, a2=6, a3=1)...")
found = None
count = 0
for trial in range(500000):
    T = random_tournament(n)
    c3_list = find_3cycles(T)
    a1 = len(c3_list)
    if a1 == 12:
        a2, a3 = count_a2_a3(c3_list)
        if a2 == 6 and a3 == 1:
            found = T
            print(f"  FOUND at trial {trial}!")
            # Print adjacency matrix
            print("  Adjacency matrix:")
            for row in T:
                print("  ", row)
            # Print score sequence
            scores = sorted([sum(row) for row in T])
            print(f"  Score sequence: {scores}")
            # Print the 12 3-cycles
            print(f"  3-cycles ({len(c3_list)}):")
            for c in c3_list:
                print(f"    {c}")
            # Print disjoint pairs
            cycle_sets = [frozenset(c) for c in c3_list]
            pairs = []
            for i in range(len(c3_list)):
                for j in range(i+1, len(c3_list)):
                    if not (cycle_sets[i] & cycle_sets[j]):
                        pairs.append((i, j))
            print(f"  Disjoint pairs ({len(pairs)}):")
            for i, j in pairs:
                print(f"    {c3_list[i]} -- {c3_list[j]}")
            # Verify with numpy
            coeffs = [a3, a2, a1, 1]
            roots = np.roots(coeffs)
            print(f"  Polynomial: 1 + {a1}x + {a2}x^2 + {a3}x^3")
            print(f"  Roots: {roots}")
            print(f"  All real: {all(abs(r.imag) < 1e-10 for r in roots)}")
            print(f"  Discriminant: {18*a1*a2*a3 - 4*a1**3*a3 + a1**2*a2**2 - 4*a2**3 - 27*a3**2}")

            # CRITICAL: verify this is actually the independence polynomial
            # of Omega_3(T). The coefficients of I(G, x) are:
            # a_k = number of independent sets of size k in G.
            # Omega_3(T) has vertices = 3-cycles, edges = pairs sharing a vertex.
            # So a1 = #vertices = #3-cycles
            # a2 = #edges in complement = #pairs of DISJOINT 3-cycles...
            # WAIT. a2 should be the number of INDEPENDENT SETS of size 2.
            # An independent set in Omega_3 means pairwise NON-adjacent,
            # i.e., pairwise vertex-disjoint 3-cycles.
            # So a2 = number of pairs of disjoint 3-cycles. YES, this is correct.
            break
    count += 1
    if count % 100000 == 0:
        print(f"  ... tried {count}")

if not found:
    print(f"  NOT FOUND in {500000} trials")

    # Let's check what (12, *, *) triples actually exist
    print("\n  Searching for all (12, a2, a3) triples...")
    triples_12 = set()
    for trial in range(200000):
        T = random_tournament(n)
        c3_list = find_3cycles(T)
        if len(c3_list) == 12:
            a2, a3 = count_a2_a3(c3_list)
            triples_12.add((a2, a3))

    print(f"  Found {len(triples_12)} distinct (a2, a3) for a1=12:")
    for a2, a3 in sorted(triples_12):
        d = 18*12*a2*a3 - 4*12**3*a3 + 12**2*a2**2 - 4*a2**3 - 27*a3**2
        print(f"    a2={a2:3d}, a3={a3:2d}: disc={d:8d}")
