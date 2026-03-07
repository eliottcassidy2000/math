#!/usr/bin/env python3
"""
CONVERSE OF THM-052 IS FALSE: Scalar M does NOT imply vertex-transitive.

At n=5: There are tournaments with M = (H/n)*I = 3*I that have |Aut|=3
and are NOT vertex-transitive. They have score sequence (1,2,2,2,3).

This script analyzes these counterexamples in detail:
1. What is the automorphism group?
2. What is the tournament structure?
3. What property DO they have that gives scalar M?

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

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

# =====================================================================
print("=" * 70)
print("ANALYSIS OF SCALAR-M NON-VERTEX-TRANSITIVE TOURNAMENTS AT n=5")
print("=" * 70)

n = 5
_, tiles = tiling_to_tournament(0, n)
m = len(tiles)

# Find all scalar-M tournaments
scalar_tilings = []
for bits in range(2**m):
    A, _ = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    if H % n != 0:
        continue
    M = transfer_matrix(A)
    if np.allclose(M, (H//n) * np.eye(n)):
        scalar_tilings.append((bits, A, H))

# Separate VT and non-VT
vt_list = []
non_vt_list = []
for bits, A, H in scalar_tilings:
    auts = automorphism_group(A)
    if is_vertex_transitive(auts, n):
        vt_list.append((bits, A, H, auts))
    else:
        non_vt_list.append((bits, A, H, auts))

print(f"\nScalar-M tournaments: {len(scalar_tilings)}")
print(f"  Vertex-transitive: {len(vt_list)}")
print(f"  NOT vertex-transitive: {len(non_vt_list)}")

# =====================================================================
print("\n" + "=" * 70)
print("NON-VT SCALAR-M TOURNAMENTS — DETAILED STRUCTURE")
print("=" * 70)

for bits, A, H, auts in non_vt_list:
    print(f"\n  bits={format(bits, f'0{m}b')}: H={H}")
    print(f"  Adjacency matrix:")
    for row in A:
        print(f"    {row}")

    scores = [sum(row) for row in A]
    print(f"  Out-degrees (scores): {scores}")

    print(f"  Automorphisms ({len(auts)}):")
    for aut in auts:
        cycle = []
        visited = set()
        for start in range(n):
            if start in visited:
                continue
            c = []
            x = start
            while x not in visited:
                visited.add(x)
                c.append(x)
                x = aut[x]
            if len(c) > 1:
                cycle.append(tuple(c))
            elif len(c) == 1:
                cycle.append((c[0],))
        print(f"    {aut} -> cycles: {cycle}")

    # Vertex orbits
    orbits = []
    assigned = set()
    for v in range(n):
        if v in assigned:
            continue
        orb = {v}
        for aut in auts:
            orb.add(aut[v])
        orbits.append(sorted(orb))
        assigned |= orb
    print(f"  Vertex orbits: {orbits}")

    # Check: is there an anti-automorphism?
    print(f"  Anti-automorphisms (T -> T^op):")
    anti_auts = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != A[perm[j]][perm[i]]:  # note: swapped!
                    valid = False
                    break
            if not valid:
                break
        if valid:
            anti_auts.append(perm)
    print(f"    Found {len(anti_auts)} anti-automorphisms")
    for aa in anti_auts:
        print(f"      {aa}")

    # Check: does reversal composed with any permutation preserve?
    # i.e., is there tau such that A[tau(i)][tau(j)] = A[j][i] for all i,j?
    # This is what the palindromic proof needs.

    # N(d,j) analysis
    print(f"\n  N(v, j) matrix (rows=vertices, cols=positions):")
    ham_paths = [p for p in permutations(range(n)) if all(A[p[i]][p[i+1]] == 1 for i in range(n-1))]
    N = np.zeros((n, n), dtype=int)
    for p in ham_paths:
        for j, v in enumerate(p):
            N[v][j] += 1
    for v in range(n):
        print(f"    vertex {v}: {list(N[v])}")

    # Check if N[v,j] is palindromic for each v
    print(f"  Palindromic N[v,*]?")
    for v in range(n):
        row = list(N[v])
        is_pal = all(row[j] == row[n-1-j] for j in range(n))
        print(f"    vertex {v}: {row} -> palindromic={is_pal}")

    # The alternating-sign vector for each vertex
    print(f"  Alternating position sums (should all be H/n = {H//n}):")
    for v in range(n):
        alt_sum = sum((-1)**j * N[v][j] for j in range(n))
        print(f"    vertex {v}: sum_j (-1)^j N[{v},j] = {alt_sum}")

    break  # Just analyze the first one in detail

# =====================================================================
print("\n" + "=" * 70)
print("ISOMORPHISM CLASS ANALYSIS")
print("=" * 70)

iso_groups = defaultdict(list)
for bits, A, H, auts in non_vt_list:
    canon = tournament_canonical(A)
    iso_groups[canon].append(bits)

print(f"\n  {len(iso_groups)} non-VT iso class(es) with scalar M")
for canon, tilts in iso_groups.items():
    print(f"  {len(tilts)} tilings in this class")

    # What tournament is this? Check known classes
    A_rep, _ = tiling_to_tournament(tilts[0], n)
    scores = tuple(sorted(sum(row) for row in A_rep))
    print(f"  Score sequence: {scores}")

    # Is this the "almost-regular" tournament?
    # At n=5, score (1,2,2,2,3) means one source-like vertex (out-degree 3)
    # and one sink-like vertex (out-degree 1)

    # Check: is it a "near-circulant" — circulant with one edge reversed?
    # Try all circulant tournaments and see if any differ in exactly 1 edge
    half = list(range(1, (n+1)//2))
    gen_sets = set()
    for mask in range(1 << len(half)):
        gs = set()
        for k, d in enumerate(half):
            if mask & (1 << k):
                gs.add(d)
            else:
                gs.add(n - d)
        gen_sets.add(frozenset(gs))

    min_diff = n*n
    for gs in gen_sets:
        A_circ = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j and (j-i)%n in gs:
                    A_circ[i][j] = 1
        # Check all relabelings of A_rep vs A_circ
        for perm in permutations(range(n)):
            A_perm = [[A_rep[perm[i]][perm[j]] for j in range(n)] for i in range(n)]
            diff = sum(1 for i in range(n) for j in range(n) if A_perm[i][j] != A_circ[i][j])
            min_diff = min(min_diff, diff)

    print(f"  Minimum edge difference from any circulant: {min_diff//2} edges")


# =====================================================================
# Now check: what do ALL scalar-M tournaments have in common?
# =====================================================================
print("\n" + "=" * 70)
print("COMMON PROPERTY: WHAT MAKES M SCALAR?")
print("=" * 70)

# All 12 iso classes at n=5
all_tilings_by_class = defaultdict(list)
for bits in range(2**m):
    A, _ = tiling_to_tournament(bits, n)
    canon = tournament_canonical(A)
    if not all_tilings_by_class[canon]:
        all_tilings_by_class[canon].append((bits, A))

print(f"\n  Total iso classes at n=5: {len(all_tilings_by_class)}")

# For each iso class, check: scalar M? VT? regular scores? doubly regular?
print("\n  Class | H | Scalar M? | VT? | Scores | |Aut| | Doubly regular?")
print("  " + "-" * 70)

for canon in sorted(all_tilings_by_class.keys()):
    bits, A = all_tilings_by_class[canon][0]
    H = ham_path_count(A)
    M = transfer_matrix(A)
    is_scalar = np.allclose(M, (H/n) * np.eye(n))
    auts = automorphism_group(A)
    vt = is_vertex_transitive(auts, n)
    scores = tuple(sorted(sum(row) for row in A))

    # Doubly regular: every pair of vertices has same number of common out-neighbors
    doubly_reg = True
    common_counts = set()
    for i in range(n):
        for j in range(i+1, n):
            common = sum(1 for k in range(n) if k != i and k != j and A[i][k] and A[j][k])
            common_counts.add(common)
    doubly_reg = len(common_counts) <= 1 if scores == (2,2,2,2,2) else False

    print(f"  {scores} | H={H:2d} | scalar={str(is_scalar):5s} | VT={str(vt):5s} | |Aut|={len(auts):2d} | DR={doubly_reg}")

# =====================================================================
# Key question: what weaker-than-VT property do the non-VT ones have?
# =====================================================================
print("\n" + "=" * 70)
print("POSITION-UNIFORM PROPERTY")
print("=" * 70)

# A tournament is "position-uniform" if every vertex appears the same number
# of times in each position across all Hamiltonian paths.
# Formally: N[v,j] = H/n for all v,j.

print("\n  Checking position-uniformity for scalar-M tournaments:")
for bits, A, H, auts in non_vt_list[:1]:
    ham_paths = [p for p in permutations(range(n)) if all(A[p[i]][p[i+1]] == 1 for i in range(n-1))]
    N = np.zeros((n, n), dtype=int)
    for p in ham_paths:
        for j, v in enumerate(p):
            N[v][j] += 1

    is_pos_uniform = np.all(N == H // n)
    print(f"  Non-VT example: N matrix =")
    for row in N:
        print(f"    {list(row)}")
    print(f"  Position-uniform (N[v,j] = {H//n} for all v,j)? {is_pos_uniform}")

for bits, A, H, auts in vt_list[:1]:
    ham_paths = [p for p in permutations(range(n)) if all(A[p[i]][p[i+1]] == 1 for i in range(n-1))]
    N = np.zeros((n, n), dtype=int)
    for p in ham_paths:
        for j, v in enumerate(p):
            N[v][j] += 1

    is_pos_uniform = np.all(N == H // n)
    print(f"\n  VT example: N matrix =")
    for row in N:
        print(f"    {list(row)}")
    print(f"  Position-uniform? {is_pos_uniform}")

# =====================================================================
# Even weaker: alternating-position-uniform
# M[a,a] = sum_j (-1)^j N(a,j) = H/n for all a => diagonal is scalar
# M[a,b] = sum_j (-1)^j N_{a->b}(j) = 0 for a != b => off-diagonal is zero
# =====================================================================
print("\n" + "=" * 70)
print("ALTERNATING-POSITION ANALYSIS")
print("=" * 70)

print("\n  For scalar M, we need:")
print("  (1) sum_j (-1)^j N(v,j) = H/n for all v")
print("  (2) Off-diagonal M[a,b] = 0 for a != b")
print()

# Check (1) for all tournaments
for bits in range(2**m):
    A, _ = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    ham_paths = [p for p in permutations(range(n)) if all(A[p[i]][p[i+1]] == 1 for i in range(n-1))]
    N = np.zeros((n, n), dtype=int)
    for p in ham_paths:
        for j, v in enumerate(p):
            N[v][j] += 1

    alt_sums = [sum((-1)**j * N[v][j] for j in range(n)) for v in range(n)]
    if len(set(alt_sums)) == 1:
        # All alternating sums equal
        M = transfer_matrix(A)
        is_scalar = np.allclose(M, (H/n) * np.eye(n))
        if not is_scalar:
            scores = tuple(sorted(sum(row) for row in A))
            print(f"  INTERESTING: Equal alt-sums but NOT scalar M! bits={format(bits, f'0{m}b')}, scores={scores}")
            print(f"    alt_sums={alt_sums}, M diagonal={[M[i][i] for i in range(n)]}")

print("  (Check complete: any tournament with equal alt-sums but non-scalar M shown above)")

# =====================================================================
print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
THM-052 CONVERSE IS FALSE:
  At n=5, there exist tournaments with M = 3*I that are NOT vertex-transitive.
  They have |Aut| = 3, score sequence (1,2,2,2,3).

The correct characterization of scalar M is:
  vertex-transitive => scalar M (THM-052, proved)
  scalar M =/=> vertex-transitive (disproved at n=5)

QUESTION: What IS the exact characterization of scalar M?
  - Not "regular" (the non-VT ones have scores (1,2,2,2,3))
  - Not "vertex-transitive"
  - Perhaps: "the alternating position distribution is uniform"?
  - Or: "the signed IE decomposition happens to cancel"?
""")
