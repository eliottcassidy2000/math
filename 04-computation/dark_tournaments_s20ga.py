#!/usr/bin/env python3
"""
dark_tournaments_s20ga.py — Searching for "dark tournaments"
kind-pasteur-2026-03-25-S20ga

SETUP:
  |tournaments| = |even graphs| = A000568(n)
  |dark tournaments| = |odd graphs| = A000171(n) (to be verified)
  |even| + |odd| = |all graphs| = A000088(n)

So: |dark| = A000088(n) - A000568(n)

What combinatorial object has exactly |odd graphs| isomorphism classes?

CANDIDATES:
  1. Tournaments with one "special" arc (labeled arc) — but that gives C(n,2) × A000568
  2. Oriented graphs (not necessarily complete) — too many
  3. Graphs with a specific coloring property
  4. Tournaments on n vertices with one vertex "marked" — gives n × A000568 / avg_aut
  5. Semi-tournaments (some arcs undirected) — need specific count
  6. Anti-tournaments (complement of tournament = also a tournament)
  7. Signed tournaments (each arc has a sign ±1)
  8. Directed graphs that are "almost complete" (missing one arc)

Let me compute the actual counts and search for matches.
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DARK TOURNAMENTS: What object counts as A000088 - A000568?")
print("  kind-pasteur-2026-03-25-S20ga")
print("=" * 80)

# Known sequences
# A000568: tournaments  1, 1, 1, 2, 4, 12, 56, 456, 6880
# A000088: all graphs   1, 1, 2, 4, 11, 34, 156, 1044, 12346
# A000171: odd graphs (NOT in standard OEIS naming — let me compute)

# Even graphs = A000568 (Royle et al. 2023)
# Odd graphs = A000088 - A000568

A000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880}
A000088 = {1:1, 2:2, 3:4, 4:11, 5:34, 6:156, 7:1044, 8:12346}

print("\n  COUNTS:")
print(f"  {'n':>3} {'Graphs':>8} {'Even=Tour':>10} {'Odd=Dark':>10} {'Dark/Even':>10}")
for n in range(1, 9):
    total = A000088.get(n, '?')
    even = A000568.get(n, '?')
    if isinstance(total, int) and isinstance(even, int):
        odd = total - even
        ratio = odd / even if even > 0 else 0
        print(f"  {n:3d} {total:8d} {even:10d} {odd:10d} {ratio:10.4f}")

# The "dark tournament" count sequence
dark = {n: A000088[n] - A000568[n] for n in range(1, 9)}
print(f"\n  Dark tournament sequence: {[dark[n] for n in range(1, 9)]}")
print(f"  = 0, 1, 2, 7, 22, 100, 588, 5466")

# ================================================================
# SEARCH: What combinatorial objects have this count?
# ================================================================
print(f"\n  SEARCHING FOR MATCHES...")

# Candidate 1: Tournaments with one "forbidden" arc removed
# (complete digraphs minus one arc)
# = "almost-tournaments" = digraphs on n vertices with C(n,2)-1 arcs
# where every pair except one has exactly one direction.
# This is a tournament with one UNDIRECTED edge.

# Let's compute: how many iso classes of "tournaments with one undirected edge"?
# At n=3: take a tournament (2 classes), pick one of 3 arcs to "un-direct"
# Total labeled = 2^3 * 3... no, labeled tournament = 2^3 = 8, each with 3 arcs.
# Un-direct one arc: 8 * 3 = 24 labeled "semi-tournaments with one undirected edge"
# But undirecting arc (i,j) gives the same result as undirecting (j,i)... no, there's
# only one way to undirect an arc: remove the direction, making {i,j} undirected.

# Actually: a "tournament with one undirected edge" has:
# - C(n,2)-1 directed arcs (the tournament minus one)
# - 1 undirected edge
# The labeled count: C(n,2) * 2^{C(n,2)-1} (choose which edge, orient the rest)
# The iso class count: ?

# Let me compute directly at small n
for n in [3, 4, 5]:
    all_perms = list(permutations(range(n)))
    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(arcs)

    # Enumerate: (tournament with one undirected edge) = (orientation of m-1 arcs, 1 unoriented)
    canon_set = set()
    for undirected_idx in range(m):
        # The remaining m-1 arcs are oriented
        for mask in range(1 << (m-1)):
            # Build adjacency: orient all arcs except undirected_idx
            A = [[0]*n for _ in range(n)]
            bit_idx = 0
            for k in range(m):
                i, j = arcs[k]
                if k == undirected_idx:
                    A[i][j] = 2  # mark as undirected
                    A[j][i] = 2
                else:
                    if mask & (1 << bit_idx):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    bit_idx += 1

            # Canonicalize
            best = None
            for p in all_perms:
                s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
                if best is None or s < best:
                    best = s
            canon_set.add(best)

    print(f"  n={n}: 'tournament + 1 undirected edge' classes = {len(canon_set)}, dark target = {dark[n]}")

# Candidate 2: Graphs with odd number of edges where all degrees have same parity
# (this is related to the definition of odd graphs)
# Actually: an "even graph" is one where every vertex has even degree.
# An "odd graph" is one that is NOT an even graph (at least one odd-degree vertex).

# Candidate 3: Mixed graphs (some directed, some undirected)
# with a specific property count

# Candidate 4: Tournaments on n vertices with a MARKED VERTEX
# Count: sum over classes C of |Aut(C)|-orbits on vertices = sum of vertex orbits
# = total "perspectives" P(n)
# From the repo: P(3) = 4, P(4) = 12, P(5) = ?

# Let me compute rooted tournament count
for n in [3, 4, 5, 6]:
    if n > 5: continue  # too slow
    all_perms = list(permutations(range(n)))

    # Build tournament classes
    canon_map = {}
    for mask in range(1 << comb(n,2)):
        arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(arcs):
            if mask & (1 << k): A[i][j] = 1
            else: A[j][i] = 1

        best = None
        for p in all_perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        canon_map[mask] = best

    classes = set(canon_map.values())

    # For each class: count vertex orbits
    total_perspectives = 0
    for cn in classes:
        # Find a representative
        for mask, c in canon_map.items():
            if c == cn:
                arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
                A = [[0]*n for _ in range(n)]
                for k, (i,j) in enumerate(arcs):
                    if mask & (1 << k): A[i][j] = 1
                    else: A[j][i] = 1

                # Count vertex orbits under Aut
                aut = [p for p in all_perms
                       if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n))]
                vertex_orbits = len(set(
                    min(tuple(p[v] for p in aut) for p_unused in [0])
                    if False else
                    min(p[v] for p in aut)
                    for v in range(n)
                ))
                total_perspectives += vertex_orbits
                break

    print(f"  n={n}: perspectives (rooted tournaments) = {total_perspectives}, dark target = {dark[n]}")

# Candidate 5: Labeled forests / trees on n vertices?
# Cayley: labeled trees = n^{n-2}. Not matching.

# Candidate 6: Self-complementary graphs
# A000171: self-complementary graphs
# n=1:1, 2:0, 3:0, 4:1, 5:2, 6:1, 7:4, 8:10
# This is NOT what we want (too few)

# Let me just check OEIS for 0, 1, 2, 7, 22, 100, 588, 5466
print(f"\n  DARK SEQUENCE: 0, 1, 2, 7, 22, 100, 588, 5466")
print(f"  OEIS search needed for this sequence.")
print(f"  Note: starts from n=1. If from n=3: 2, 7, 22, 100, 588, 5466")

print("\nDONE.")
print("=" * 80)
