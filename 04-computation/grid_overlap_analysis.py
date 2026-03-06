#!/usr/bin/env python3
"""
Triangular Grid Overlap Structure
===================================
The tiling grid for n-vertex tournaments is a triangular grid with C(n,2) cells.
Key insight: it contains n copies of the (n-1)-grid (one per deleted vertex),
and C(n,2) copies of the (n-2)-grid (one per deleted pair).

Investigate how these overlapping subgrids relate to tournament structure.

kind-pasteur-2026-03-06-S18e
"""

import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles

# Visualize the triangular grid
for n in [4, 5, 6]:
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)  # = C(n,2)
    tiling_bits = m - (n - 1)  # non-path arcs

    print(f"\n{'='*60}")
    print(f"n={n}: {m} pairs, {n-1} path arcs, {tiling_bits} tiling bits")
    print(f"{'='*60}")

    # Show the grid as triangle
    print("Grid layout:")
    for i in range(n - 1):
        prefix = "  " * i
        row = " ".join(f"({i},{j})" for j in range(i + 1, n))
        print(f"  {prefix}{row}")

    # n copies of (n-1)-grid
    print(f"\n{n} copies of (n-1)={n-1} grid:")
    for v in range(n):
        sub = [(i, j) for i, j in pairs if i != v and j != v]
        print(f"  Delete v={v}: {len(sub)} cells = {sub}")

    # Overlap between copies
    print(f"\nOverlap between copy pairs:")
    for v1 in range(n):
        for v2 in range(v1 + 1, n):
            sub1 = set((i, j) for i, j in pairs if i != v1 and j != v1)
            sub2 = set((i, j) for i, j in pairs if i != v2 and j != v2)
            overlap = sub1 & sub2
            print(f"  Delete {v1} AND {v2}: {len(overlap)} overlap = (n-2)-grid on {n-2} verts")
            break
        break

    # Each pair (i,j) belongs to exactly n-2 copies (all except i and j)
    print(f"\nPair membership counts:")
    for p in pairs[:3]:  # show first 3
        copies = [v for v in range(n) if v not in p]
        print(f"  {p}: in {len(copies)} copies (v={copies})")
    print(f"  ... (each pair in exactly n-2 = {n-2} copies)")

# Now the key structural question from the user:
# The n=4 grid has 3 tiling bits.
# When we overlay 4 copies of the 3-tiling on the n=5 grid,
# which tiles overlap?

print(f"\n{'='*60}")
print("n=4 TILINGS ON n=5 GRID")
print(f"{'='*60}")

# n=4 Hamiltonian path: 3 path arcs out of 6 pairs
# Tiling = binary assignment of 3 non-path arcs
# For n=4, a specific path (say 0->1->2->3) uses arcs (0,1),(1,2),(2,3)
# Non-path arcs: (0,2),(0,3),(1,3) - these are the 3 tiling bits

# In the n=5 grid, delete v=4 gives the n=4 grid on {0,1,2,3}
# Delete v=0 gives the n=4 grid on {1,2,3,4}
# etc.

n = 5
pairs5 = [(i, j) for i in range(n) for j in range(i + 1, n)]
print(f"n=5 grid: {pairs5}")
print(f"n=5 non-path arcs for path 0->1->2->3->4:")
path_arcs = [(0, 1), (1, 2), (2, 3), (3, 4)]
nonpath = [p for p in pairs5 if p not in path_arcs]
print(f"  Non-path: {nonpath} ({len(nonpath)} tiles)")

# For each deleted vertex, what are the non-path arcs of the induced subpath?
print(f"\nDeleting each vertex from path 0->1->2->3->4:")
for v in range(5):
    # Subpath: remove v, reconnect
    subverts = [i for i in range(5) if i != v]
    subpairs = [(i, j) for i, j in pairs5 if i != v and j != v]

    # The induced path depends on where v was in the original path
    # If v was internal, deleting it breaks the path
    # But in the (n-1)-tournament T-v, there exists a different Hamiltonian path

    # For visualization, let's just look at the pair sets
    print(f"  v={v}: subgrid pairs = {subpairs}")

    # Which of the n=5 non-path arcs are in this subgrid?
    nonpath_in_sub = [p for p in nonpath if p in subpairs]
    print(f"          overlapping non-path arcs: {nonpath_in_sub}")

# The key insight: tiling bits of the n=5 grid that correspond to
# non-path arcs of T restricted to any 4-vertex subset are
# related to the tiling of the induced tournament

# Let me count: which n=5 tiling bits are in how many n=4 subgrids?
print(f"\nFor each n=5 tiling bit, membership in n=4 subgrids:")
for p in nonpath:
    copies = [v for v in range(5) if v not in p]
    print(f"  {p}: in {len(copies)} copies (delete {copies})")

# This is always n-2 = 3 copies
# So each of the 6 tiling bits appears in exactly 3 of the 5 copies

# STRUCTURAL OBSERVATION:
# The n=5 grid has 6 tiling bits. Each appears in 3/5 copies.
# Total = 6*3 = 18. Each copy has C(4,2)-(4-1) = 3 tiling bits.
# So 5*3 = 15... wait, that counts ALL bits in copies, not just non-path ones.
# Actually, the "tiling bits" of the n=4 subgrid depend on which path is used
# in the subtournament. The non-path arcs of T-v may differ from those of T.

print(f"\n{'='*60}")
print("PHYSICAL OVERLAY: 4 copies of n=3 grid on n=4 grid")
print(f"{'='*60}")
print("n=3 grid (3 cells):")
print("  (0,1) (0,2)")
print("        (1,2)")
print()
print("n=4 grid (6 cells):")
print("  (0,1) (0,2) (0,3)")
print("        (1,2) (1,3)")
print("              (2,3)")
print()
print("4 copies of n=3 (delete one vertex):")
print("  v=0: (1,2)(1,3)(2,3) = bottom-right triangle")
print("  v=1: (0,2)(0,3)(2,3) = skip column 1 / row 1")
print("  v=2: (0,1)(0,3)(1,3) = skip column 2 / row 2")
print("  v=3: (0,1)(0,2)(1,2) = top-left triangle")
print()
print("Each cell of n=4 grid is covered by exactly 2 copies:")
print(f"  (0,1): v=2,3")
print(f"  (0,2): v=1,3")
print(f"  (0,3): v=1,2")
print(f"  (1,2): v=0,3")
print(f"  (1,3): v=0,2")
print(f"  (2,3): v=0,1")

# Now for n+2: n=4 inside n=6
print(f"\n{'='*60}")
print("n=4 INSIDE n=6: delete 2 vertices")
print(f"{'='*60}")
n = 6
pairs6 = [(i, j) for i in range(n) for j in range(i + 1, n)]
print(f"n=6 grid: {len(pairs6)} cells")
print(f"C(6,2)=15 copies of n=4 grid (delete 2 vertices)")
for v1 in range(n):
    for v2 in range(v1 + 1, n):
        sub = [(i, j) for i, j in pairs6 if i not in (v1, v2) and j not in (v1, v2)]
        if v1 == 0:
            print(f"  Delete {{{v1},{v2}}}: {sub}")

# Each cell (i,j) is in C(n-2, 0) = ...
# Specifically, (i,j) is in the copy for (v1,v2) iff {v1,v2} and {i,j} are disjoint
# So (i,j) is in C(n-2, 2) copies... no.
# Delete {v1,v2}: keep cells not involving v1 or v2.
# Cell (i,j) is kept iff i,j not in {v1,v2}.
# Number of (v1,v2) pairs where i,j not in {v1,v2}: C(n-2, 2).
# For n=6: each cell is in C(4,2)=6 copies.

print(f"\nEach cell in n=6 grid is in C(4,2)=6 of the 15 copies")

# Now the even/odd parity connection:
# n=4 -> n=6 skips n=5 (even to even)
# The user says: "n+2 case is just as worthy due to even-odd stuff"
#
# For even n: sigma (anti-aut) is fpf, 3 orbits {a,sigma(a)} at n=6
# The 3-cycle disjoint pairs correspond to selecting one from each orbit
# This is a "transversal" of the orbit partition
#
# For the n=4 subgrids inside n=6: deleting 2 vertices from different orbits
# gives a 4-vertex tournament with specific symmetry properties

print(f"\n{'='*60}")
print("SIGMA-ORBIT CONNECTION")
print(f"{'='*60}")
# At n=6 with sigma = (1,0,5,4,3,2), orbits are {0,1}, {2,5}, {3,4}
sigma = (1, 0, 5, 4, 3, 2)
orbits = [(0, 1), (2, 5), (3, 4)]
print(f"Orbits: {orbits}")

# Deleting one from each orbit gives a 3-vertex tournament
# Deleting the sigma-complement gives the complementary 3-vertex tournament
print(f"\nDelete one from each orbit (3 vertex tournaments):")
from itertools import product
for choice in product(*orbits):
    complement = tuple(sigma[v] for v in choice)
    print(f"  Keep {sorted(choice)}, delete {sorted(complement)}")

# Deleting 2 vertices from SAME orbit
print(f"\nDelete entire orbit (4-vertex tournaments):")
for a, b in orbits:
    remaining = [v for v in range(6) if v not in (a, b)]
    print(f"  Delete orbit {{{a},{b}}}: keep {remaining}")

# These 3 copies of n=4 are special: they are T restricted to
# vertices from 2 orbits, which inherits sigma-symmetry!

print("\nDone.")
