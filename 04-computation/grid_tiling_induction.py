#!/usr/bin/env python3
"""
Grid Tiling Induction: The user's insight about overlapping subgrids.

The tiling grid for n-vertex tournaments is a triangular grid with C(n,2) cells.
Each cell (i,j) with i<j represents the arc direction between vertices i and j.

For a Hamiltonian path P = v_0 -> v_1 -> ... -> v_{n-1}:
- "Path arcs": (v_k, v_{k+1}) for k=0,...,n-2. These are FIXED.
- "Tiling bits": the remaining C(n,2) - (n-1) = (n-1)(n-2)/2 arcs.

The key insight: when we delete a vertex v to get the (n-1)-subgrid,
the Hamiltonian paths of T-v correspond to tilings of the (n-1)-grid.
The tiling bits of T that fall in the (n-1)-subgrid (excluding v) are a
SUBSET of the T tiling bits plus some path arcs that become tiling bits
in the restricted tournament.

This analysis connects H(T) and H(T-v) via the tiling structure.

QUESTION: For the H-maximizer, how do the n overlapping (n-1)-tilings
relate? Does the maximizer achieve a "balanced" overlap pattern?

kind-pasteur-2026-03-06-S18f
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import permutations

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def find_all_ham_paths(T):
    """Return list of all Hamiltonian paths as tuples of vertex indices."""
    n = len(T)
    paths = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            paths.append(perm)
    return paths

def path_arcs(path):
    """Return set of directed arcs used by the path."""
    return {(path[i], path[i+1]) for i in range(len(path)-1)}

def tiling_bits(T, path):
    """Return the tiling (non-path arc directions) for a given path."""
    n = len(T)
    p_arcs = path_arcs(path)
    tiling = {}
    for i in range(n):
        for j in range(i+1, n):
            if (i, j) not in p_arcs and (j, i) not in p_arcs:
                tiling[(i, j)] = T[i][j]
    return tiling

# ============================================================
# Analyze the tiling structure for small n
# ============================================================
print("=" * 70)
print("TILING STRUCTURE ANALYSIS")
print("=" * 70)

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    n_tiling = m - (n - 1)

    print(f"\n{'='*60}")
    print(f"n={n}: {m} arcs, {n-1} path arcs, {n_tiling} tiling bits")
    print(f"{'='*60}")

    # Find the maximizer
    max_h = 0
    max_bits = None
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h > max_h:
            max_h = h
            max_bits = bits

    T = tournament_from_bits(n, max_bits)
    paths = find_all_ham_paths(T)

    print(f"Maximizer: bits={max_bits}, H={max_h}")
    print(f"Score: {score_seq(T)}")

    # For each path, compute the tiling Hamming weight (number of 1-bits in tiling)
    hw_counts = {}
    for p in paths:
        tb = tiling_bits(T, p)
        hw = sum(tb.values())
        hw_counts[hw] = hw_counts.get(hw, 0) + 1

    print(f"Tiling Hamming weight distribution:")
    for hw in sorted(hw_counts.keys()):
        print(f"  HW={hw}: {hw_counts[hw]} paths")

    # For the maximizer, how does the tiling of a path in T
    # relate to the tilings in T-v?
    if n <= 5:
        print(f"\nPath -> deletion mapping:")
        for p in paths[:5]:  # first 5 paths
            p_arcs_set = path_arcs(p)
            print(f"  Path: {'->'.join(map(str, p))}")

            for v_pos in range(n):
                v = p[v_pos]
                # Remove v from path
                sub_path = tuple(x for x in p if x != v)
                # Check if sub_path is valid in T-v
                valid = True
                if v_pos > 0 and v_pos < n-1:
                    # Internal deletion: need arc p[v_pos-1] -> p[v_pos+1]
                    if not T[p[v_pos-1]][p[v_pos+1]]:
                        valid = False

                # Reindex: vertices in T-v are {0,...,n-1}\{v}
                verts = [i for i in range(n) if i != v]
                reindex = {old: new for new, old in enumerate(verts)}

                if valid:
                    sub_path_reindexed = tuple(reindex[x] for x in sub_path)
                    # This IS a Ham path in T-v
                    # What are the tiling bits?
                    T_v = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
                    tb_sub = tiling_bits(T_v, sub_path_reindexed)
                    hw_sub = sum(tb_sub.values())
                    status = f"valid, HW(T-v)={hw_sub}"
                else:
                    status = "INVALID (arc goes wrong way)"

                if v_pos <= 1 or v_pos >= n-2:  # show endpoints and near-endpoints
                    print(f"    del v={v} (pos {v_pos}): {status}")

# ============================================================
# The overlap structure: for each tiling bit of the n-grid,
# which (n-1)-subgrids contain it, and what role does it play?
# ============================================================
print(f"\n{'='*70}")
print("TILING BIT OVERLAP STRUCTURE")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    n_tiling = m - (n - 1)

    # For path 0->1->...->n-1:
    path = tuple(range(n))
    p_arcs = path_arcs(path)

    # Non-path arcs (tiling bits)
    tiling_cells = [(i, j) for i in range(n) for j in range(i+1, n)
                     if (i, j) not in p_arcs and (j, i) not in p_arcs]

    print(f"\nn={n}: path 0->1->...->n-1")
    print(f"Path arcs: {sorted(p_arcs)}")
    print(f"Tiling bits: {tiling_cells}")
    print(f"Number of tiling bits: {len(tiling_cells)}")

    # For each tiling bit, which vertex deletions keep it?
    print(f"\nTiling bit membership in (n-1)-subgrids:")
    for (i, j) in tiling_cells:
        containing = [v for v in range(n) if v != i and v != j]
        # In the (n-1)-subgrid from deleting v, arc (i,j) exists.
        # Is (i,j) still a tiling bit in T-v, or does it become a path arc?
        # That depends on which Hamiltonian path we use in T-v.
        # But from the GRID perspective, (i,j) is always in the subgrid for v not in {i,j}.
        print(f"  ({i},{j}): in {len(containing)} copies (del {containing})")

    # KEY: For each pair of (n-1)-subgrids (del v1, del v2),
    # what tiling bits do they share?
    print(f"\nShared tiling bits between pairs of (n-1)-subgrids:")
    for v1 in range(n):
        for v2 in range(v1+1, n):
            shared = [(i,j) for (i,j) in tiling_cells
                      if v1 not in (i,j) and v2 not in (i,j)]
            if v1 == 0:
                print(f"  del {v1}, del {v2}: {len(shared)} shared tiling bits")

# ============================================================
# The crucial question: for the MAXIMIZER, what is the
# "tiling constraint" that makes it optimal?
# ============================================================
print(f"\n{'='*70}")
print("MAXIMIZER TILING CONSTRAINT")
print("=" * 70)

# At n=5, maximizer has H=15 and 3 tiling bits.
# 2^3 = 8 possible tilings, but only 15 paths exist (not 2^3 * something).
# The tiling Hamming weight distribution should be symmetric (bell curve).
n = 5
T = tournament_from_bits(n, 76)  # bits=76 is regular Paley T_5
h = hamiltonian_path_count(T)
paths = find_all_ham_paths(T)
print(f"\nn=5 Paley (bits=76): H={h}")

# For EACH path, what are the 3 tiling bits?
tiling_data = []
for p in paths:
    tb = tiling_bits(T, p)
    hw = sum(tb.values())
    tiling_data.append((p, tb, hw))

hw_dist = {}
for _, _, hw in tiling_data:
    hw_dist[hw] = hw_dist.get(hw, 0) + 1
print(f"HW distribution: {dict(sorted(hw_dist.items()))}")
print(f"Expected symmetric: HW k and HW (3-k) should match")

# Now: the grid-overlap question.
# The n=5 grid has 10 arcs. Path uses 4, leaving 6 tiling bits.
# Wait, n=5 has C(5,2)=10 arcs, path uses 4, so 6 non-path arcs.
# But "tiling bits" depends on WHICH path we're looking at.
# The path 0->1->2->3->4 has path arcs (0,1),(1,2),(2,3),(3,4).
# Non-path: (0,2),(0,3),(0,4),(1,3),(1,4),(2,4) = 6 tiling bits.

# Hmm, earlier I said 3. Let me recheck.
# C(5,2) - 4 = 10 - 4 = 6. That's 6 tiling bits, not 3.
# The "(n-1)(n-2)/2" formula gives 4*3/2 = 6. Correct.
print(f"\nCorrection: n=5 has {10-4}=6 tiling bits (not 3)")

# At n=4: C(4,2)=6 arcs, 3 path arcs, 3 tiling bits. Correct.
# At n=5: C(5,2)=10, 4 path arcs, 6 tiling bits. 2^6 = 64.
# H=15 means 15/64 = 23.4% of tilings correspond to paths.

# The key structural question: among the 2^6 = 64 possible
# "tiling assignments" (each non-path arc can go either way),
# which ones correspond to actual Hamiltonian paths?

# For a fixed starting path order, the tiling bits are determined
# by the tournament. But DIFFERENT paths have different "path arcs"
# and thus different tiling bits.

# Let me think about this differently.
# The TOURNAMENT T determines ALL arc directions.
# Different Hamiltonian paths P of T select different sets of n-1 arcs as "path arcs."
# The remaining C(n,2)-(n-1) arcs are "tiling bits" for that particular path.
# These are not free to vary -- they're determined by T.

# So the "tiling" is NOT a binary choice per cell.
# The tournament T IS the tiling -- it's the fixed assignment of all arc directions.
# Different Hamiltonian paths just SELECT different (n-1)-subsets as their path arcs.

# The user's grid insight is about the PHYSICAL structure of the grid itself,
# not about varying tiling bits. The n-grid has n copies of (n-1)-grid,
# and the arc directions in the overlap region are SHARED.

# This means: if we know the tournament T on n vertices, the subtournament T-v
# on n-1 vertices is completely determined. The "tiling" of T-v is the restriction
# of T's arc directions to the (n-1)-subgrid.

# So the inductive question is: given that T-v has certain structural properties
# (e.g., max H in its class), what does that force on T?

# For the MAXIMIZER T at odd n: ALL vertex deletions give the (n-1)-maximizer.
# This means the entire tournament is constrained by n overlapping optimality conditions.
# Each (n-1)-subgrid must contain the maximizer tiling.

# How many degrees of freedom remain?
# At n=7 (regular, 240 maximizers): the n=6 subgrids have H=45.
# Each n=6 subgrid has C(6,2)=15 arcs, and the n=7 grid has C(7,2)=21 arcs.
# Each arc is in n-2=5 of the 7 subgrids.
# So 7 * 15 = 105 constraints (with n-2=5 sharing) = 21 independent arcs.
# This means the 21 arc directions are (in principle) fully determined
# by the 7 subgrid optimality conditions!

print(f"\n{'='*70}")
print("DEGREES OF FREEDOM ANALYSIS")
print("=" * 70)

for n in range(4, 9):
    total_arcs = n * (n-1) // 2
    sub_arcs = (n-1) * (n-2) // 2
    total_constraints = n * sub_arcs  # n subgrids, each with sub_arcs constraints
    sharing = n - 2  # each arc is in n-2 subgrids
    independent_constraints = total_constraints // sharing  # = total_arcs
    print(f"n={n}: {total_arcs} arcs, {n} subgrids of {sub_arcs} arcs each, "
          f"sharing factor {sharing}, "
          f"ratio = {total_constraints}/{total_arcs} = {total_constraints/total_arcs:.1f}")

# The ratio is n*(n-1)(n-2)/2 / (n(n-1)/2) = n-2 = sharing factor.
# So the system is exactly determined (each constraint counted n-2 times).
# The n overlapping optimality conditions + consistency = unique tournament?

# At n=7: there are 240 maximizers, not 1. So the conditions don't
# uniquely determine T. But 240/2^21 ≈ 0.01% is very small.

for n in range(4, 9):
    total = 1 << (n*(n-1)//2)
    max_count = {4: 24, 5: 64, 6: 480, 7: 240, 8: '?'}
    mc = max_count.get(n, '?')
    if isinstance(mc, int):
        frac = mc / total * 100
        print(f"n={n}: {mc}/{total} = {frac:.4f}% are maximizers")

print("\nDone.")
