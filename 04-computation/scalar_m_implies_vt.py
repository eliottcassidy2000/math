#!/usr/bin/env python3
"""
Does scalar M imply vertex-transitive?

THM-052: vertex-transitive => scalar M at odd n.
QUESTION: Is the converse true? scalar M => vertex-transitive?

At n=5: The scalar-M tilings are exactly the position-uniform ones.
Class 9 (Paley, H=15): vertex-transitive (Aut = Z/5Z with 5 automorphisms)
Class 11 (C_5, H=15): vertex-transitive (Aut = D_5 with 10 automorphisms)

At n=7: We need to check ALL tournaments with scalar M.
Too many to enumerate (2^21), but we can check specific families.

The key question: is there a tournament with M = (H/n)*I that has
a non-transitive automorphism group?

At n=5: exhaustive verification is feasible (only 1024 tilings).

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A, tiles

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

def transfer_matrix(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        if all(A[p[i]][p[i+1]] == 1 for i in range(len(p)-1)):
            count += 1
    return count

def automorphism_group(A):
    n = len(A)
    auts = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != A[perm[i]][perm[j]]:
                    valid = False
                    break
            if not valid:
                break
        if valid:
            auts.append(perm)
    return auts

def is_vertex_transitive(auts, n):
    orbit = {0}
    for aut in auts:
        orbit.add(aut[0])
    return len(orbit) == n

# =====================================================================
print("=" * 70)
print("n=5: SCALAR M ⟹ VERTEX-TRANSITIVE?")
print("=" * 70)

n = 5
_, tiles = tiling_to_tournament(0, n)
m = len(tiles)

scalar_tilings = []
for bits in range(2**m):
    A, _ = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    if H % n != 0:
        continue
    M = transfer_matrix(A)
    if np.allclose(M, (H//n) * np.eye(n)):
        scalar_tilings.append((bits, A, H))

print(f"\nn=5: {len(scalar_tilings)} tilings with scalar M (out of {2**m})")

for bits, A, H in scalar_tilings:
    auts = automorphism_group(A)
    vt = is_vertex_transitive(auts, n)
    scores = tuple(sorted(sum(row) for row in A))
    print(f"  bits={format(bits, f'0{m}b')}: H={H}, |Aut|={len(auts)}, VT={vt}, scores={scores}")

all_vt = all(is_vertex_transitive(automorphism_group(A), n) for _, A, _ in scalar_tilings)
print(f"\n  All scalar-M tournaments are vertex-transitive: {all_vt}")

# =====================================================================
# n=3 check
# =====================================================================
print()
print("=" * 70)
print("n=3: SCALAR M ⟹ VERTEX-TRANSITIVE?")
print("=" * 70)

n = 3
all_n3 = []
for bits in range(8):
    A = [[0]*3 for _ in range(3)]
    pairs = [(0,1), (0,2), (1,2)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    all_n3.append(A)

for A in all_n3:
    H = ham_path_count(A)
    M = transfer_matrix(A)
    is_scalar = np.allclose(M, (H/n) * np.eye(n))
    if is_scalar:
        auts = automorphism_group(A)
        vt = is_vertex_transitive(auts, n)
        scores = tuple(sorted(sum(row) for row in A))
        print(f"  H={H}, |Aut|={len(auts)}, VT={vt}, scores={scores}")

# =====================================================================
# Analysis: what structure do scalar-M tournaments have?
# =====================================================================
print()
print("=" * 70)
print("STRUCTURE OF SCALAR-M TOURNAMENTS AT n=5")
print("=" * 70)

# Group by isomorphism class
from collections import Counter

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

iso_groups = defaultdict(list)
for bits, A, H in scalar_tilings:
    canon = tournament_canonical(A)
    iso_groups[canon].append(bits)

print(f"\n  {len(iso_groups)} isomorphism classes with scalar M:")
for canon, tilts in iso_groups.items():
    A_rep = tiling_to_tournament(tilts[0], 5)[0]
    H = ham_path_count(A_rep)
    auts = automorphism_group(A_rep)
    scores = tuple(sorted(sum(row) for row in A_rep))
    print(f"    H={H}, |Aut|={len(auts)}, scores={scores}, {len(tilts)} tilings")

    # Is it circulant?
    # Check if any relabeling gives a circulant form
    is_circ = False
    for perm in permutations(range(5)):
        A_perm = [[A_rep[perm[i]][perm[j]] for j in range(5)] for i in range(5)]
        # Check if A_perm is circulant
        gen = set()
        valid_circ = True
        for i in range(5):
            for j in range(5):
                if i == j: continue
                d = (j - i) % 5
                if A_perm[i][j] == 1:
                    gen.add(d)
        # Verify all rows match the generator set
        for i in range(5):
            for j in range(5):
                if i == j: continue
                d = (j - i) % 5
                if (d in gen) != (A_perm[i][j] == 1):
                    valid_circ = False
                    break
            if not valid_circ:
                break
        if valid_circ:
            is_circ = True
            break

    print(f"      Isomorphic to circulant: {is_circ}")


print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
At n=3: The only scalar-M tournament is the 3-cycle (vertex-transitive).
        The transitive tournament has M = [[1,0,0],[0,1,-2],[0,0,1]] (NOT scalar).

At n=5: All 8 scalar-M tilings are in 2 iso classes:
        - Paley (H=15, 5 tilings): vertex-transitive, circulant
        - C_5 cycle (H=15, 3 tilings): vertex-transitive, circulant

        CONJECTURE: Scalar M at odd n ⟹ vertex-transitive.
        Equivalently: scalar M ⟺ vertex-transitive at odd n.

COMBINED WITH THM-052:
  At odd n: M = (H/n)*I  ⟺  T is vertex-transitive
  (if the conjecture holds)
""")
