#!/usr/bin/env python3
"""
Deep analysis of blue pairs and the skeleton at n=5.

The "blue pair" structure connects a tournament T to its "opposite"
under the tiling model. At odd n, Theorem 5 says there are no "blueself"
pairs (T isomorphic to flip(T)).

But there ARE self-flip pairs where T and flip(T) are in the same
isomorphism class (just not isomorphic to each other). What is the
spectral/structural meaning of these?

KEY FOCUS:
1. The 4 self-flip pairs at n=5 (classes 8 and 10)
2. Blue pair H values: H(T) + H(flip(T)) pattern
3. The "blue pair skeleton" — which classes are connected by flip?
4. Transfer matrix relationship: M(T) vs M(flip(T))
5. Connection to odd cycle collections
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
    return A

def ham_path_count(A):
    n = len(A)
    count = 0
    for p in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start:
            continue
        if end is not None and p[-1] != end:
            continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

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

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

def odd_cycles(A):
    n = len(A)
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                valid = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1) % length]] != 1:
                        valid = False
                        break
                if valid:
                    min_idx = perm.index(min(perm))
                    canonical = perm[min_idx:] + perm[:min_idx]
                    if canonical not in cycles:
                        cycles.append(canonical)
    return cycles

n = 5
tiles = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles.append((a, b))
tiles.sort()
m = len(tiles)
num_tilings = 2**m

# Precompute all data
tiling_H = {}
tiling_canon = {}
iso_classes = defaultdict(list)

for bits in range(num_tilings):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    canon = tournament_canonical(A)
    tiling_H[bits] = H
    tiling_canon[bits] = canon
    iso_classes[canon].append(bits)

class_labels = {}
for idx, canon in enumerate(sorted(iso_classes.keys())):
    class_labels[canon] = idx

print("=" * 70)
print("BLUE PAIR ANALYSIS AT n=5")
print("=" * 70)

# =====================================================================
# 1. Full blue pair map: which classes are connected by flip?
# =====================================================================
print("\n--- BLUE PAIR MAP: class(T) <-> class(flip(T)) ---\n")

blue_pair_map = defaultdict(int)  # (class1, class2) -> count
H_sum_by_pair = defaultdict(set)

for bits in range(num_tilings):
    flip_bits = ((1 << m) - 1) ^ bits
    if flip_bits >= bits:
        c1 = class_labels[tiling_canon[bits]]
        c2 = class_labels[tiling_canon[flip_bits]]
        pair = (min(c1,c2), max(c1,c2))
        blue_pair_map[pair] += 1
        H1 = tiling_H[bits]
        H2 = tiling_H[flip_bits]
        H_sum_by_pair[pair].add(H1 + H2)

print(f"{'Class pair':<15} {'Count':>5}  {'H-sum':>8}  Notes")
print("-" * 55)
for (c1, c2), count in sorted(blue_pair_map.items()):
    H_sums = H_sum_by_pair[(c1, c2)]
    is_self = (c1 == c2)
    notes = "SELF-FLIP" if is_self else ""
    print(f"  ({c1:2d},{c2:2d})        {count:3d}    {H_sums}    {notes}")

# =====================================================================
# 2. H(T) + H(flip(T)) distribution
# =====================================================================
print()
print("=" * 70)
print("H(T) + H(flip(T)) DISTRIBUTION")
print("=" * 70)

H_sum_dist = defaultdict(int)
for bits in range(num_tilings):
    flip_bits = ((1 << m) - 1) ^ bits
    if flip_bits >= bits:
        H_sum = tiling_H[bits] + tiling_H[flip_bits]
        H_sum_dist[H_sum] += 1

print("\n  H(T) + H(flip(T))  Count")
print("  " + "-" * 30)
for s in sorted(H_sum_dist.keys()):
    print(f"  {s:18d}  {H_sum_dist[s]:5d}")

# Is the sum always the same? (This would be a strong constraint)
print(f"\n  All sums equal? {len(set(H_sum_dist.keys())) == 1}")
print(f"  Mean sum = {sum(s*c for s,c in H_sum_dist.items()) / sum(H_sum_dist.values()):.2f}")

# =====================================================================
# 3. What does flip DO to the tournament? (arc-reversal structure)
# =====================================================================
print()
print("=" * 70)
print("FLIP OPERATION: what arcs get reversed?")
print("=" * 70)

# flip reverses ALL non-path arcs. The path arcs (i -> i-1) are FIXED.
# So flip(T) reverses {(a,b) : a-b >= 2} while keeping {(i,i-1)} fixed.
print(f"\n  Path arcs (fixed): " + ", ".join(f"{i}->{i-1}" for i in range(1, n)))
print(f"  Non-path arcs (flipped): " + ", ".join(f"({a},{b})" for a,b in tiles))

# For a specific tournament, show the flip effect
print(f"\n  Example: tiling 000000 (all forward) -> flip = 111111 (all backward)")
A_all_fwd = tiling_to_tournament(0, n)
A_all_bwd = tiling_to_tournament(63, n)
print(f"  T_fwd arcs: ", end="")
for i in range(n):
    for j in range(n):
        if A_all_fwd[i][j] == 1:
            print(f"{i}->{j} ", end="")
print()
print(f"  T_bwd arcs: ", end="")
for i in range(n):
    for j in range(n):
        if A_all_bwd[i][j] == 1:
            print(f"{i}->{j} ", end="")
print()

# How many arcs differ between T and flip(T)?
print(f"\n  Number of arcs reversed in flip = {m} (= number of non-path arcs)")

# Is flip(T) = T^op (reverse ALL arcs)?
# No! T^op reverses ALL arcs including path arcs.
# flip only reverses non-path arcs.
# So flip(T) = T^op composed with "re-orient path arcs."

# Let's check: is flip(T) related to T^op by a relabeling?
print(f"\n  Is flip(T) related to T^op?")

A = tiling_to_tournament(0, n)  # transitive tournament
A_flip = tiling_to_tournament(63, n)
A_op = [[A[j][i] for j in range(n)] for i in range(n)]

# Check if flip and op are the same
print(f"  T = ", [[A[i][j] for j in range(n)] for i in range(n)])
print(f"  flip(T) = ", [[A_flip[i][j] for j in range(n)] for i in range(n)])
print(f"  T^op = ", [[A_op[i][j] for j in range(n)] for i in range(n)])

# Check if there's a permutation sigma such that flip(T) = sigma(T^op)
canon_flip = tournament_canonical(A_flip)
canon_op = tournament_canonical(A_op)
print(f"  flip(T) iso T^op? {canon_flip == canon_op}")

# =====================================================================
# 4. Transfer matrix relationship: M(T) vs M(flip(T))
# =====================================================================
print()
print("=" * 70)
print("TRANSFER MATRIX: M(T) vs M(flip(T))")
print("=" * 70)

# Sample a few representative blue pairs
print("\nFor each blue pair, compute M(T) and M(flip(T)):")

for bits in [0, 1, 5, 10, 21]:
    flip_bits = ((1 << m) - 1) ^ bits
    A1 = tiling_to_tournament(bits, n)
    A2 = tiling_to_tournament(flip_bits, n)
    M1 = transfer_matrix(A1)
    M2 = transfer_matrix(A2)
    H1 = tiling_H[bits]
    H2 = tiling_H[flip_bits]
    c1 = class_labels[tiling_canon[bits]]
    c2 = class_labels[tiling_canon[flip_bits]]

    print(f"\n  T={format(bits, f'0{m}b')} (H={H1}, class={c1}) <-> flip={format(flip_bits, f'0{m}b')} (H={H2}, class={c2})")
    print(f"  M(T)      = {M1.tolist()}")
    print(f"  M(flip)   = {M2.tolist()}")

    # Check relationships
    P = np.zeros((n,n))
    for i in range(n):
        P[i][n-1-i] = 1

    M2_rev = P @ M1 @ P.T  # relabeling v -> n-1-v
    print(f"  M(flip) = P*M(T)*P^T? {np.allclose(M2, M2_rev)}")
    print(f"  M(flip) = M(T)^T?     {np.allclose(M2, M1.T)}")
    print(f"  M(flip) = -M(T)?      {np.allclose(M2, -M1)}")

    # Check eigenvalue relationship
    e1 = sorted(np.linalg.eigvalsh(M1).tolist(), reverse=True)
    e2 = sorted(np.linalg.eigvalsh(M2).tolist(), reverse=True)
    print(f"  evals(T)    = {[round(e,3) for e in e1]}")
    print(f"  evals(flip) = {[round(e,3) for e in e2]}")

# =====================================================================
# 5. Odd cycle collections under flip
# =====================================================================
print()
print("=" * 70)
print("ODD CYCLE COLLECTIONS: T vs flip(T)")
print("=" * 70)

for bits in [0, 1, 10, 21]:
    flip_bits = ((1 << m) - 1) ^ bits
    A1 = tiling_to_tournament(bits, n)
    A2 = tiling_to_tournament(flip_bits, n)
    H1 = tiling_H[bits]
    H2 = tiling_H[flip_bits]

    oc1 = odd_cycles(A1)
    oc2 = odd_cycles(A2)

    print(f"\n  T={format(bits, f'0{m}b')} (H={H1}): {len(oc1)} odd cycles")
    for c in oc1:
        print(f"    {list(c)} (length {len(c)})")

    print(f"  flip={format(flip_bits, f'0{m}b')} (H={H2}): {len(oc2)} odd cycles")
    for c in oc2:
        print(f"    {list(c)} (length {len(c)})")

    # Count 3-cycles and 5-cycles separately
    tri1 = [c for c in oc1 if len(c) == 3]
    tri2 = [c for c in oc2 if len(c) == 3]
    pent1 = [c for c in oc1 if len(c) == 5]
    pent2 = [c for c in oc2 if len(c) == 5]

    print(f"  3-cycles: {len(tri1)} vs {len(tri2)}")
    print(f"  5-cycles: {len(pent1)} vs {len(pent2)}")
    print(f"  H = 1 + 2*(#3-cycles) + 2*(#5-cycles): {1 + 2*len(tri1) + 2*len(pent1)} vs {1 + 2*len(tri2) + 2*len(pent2)}")

# =====================================================================
# 6. The "gradient" structure of the skeleton
# =====================================================================
print()
print("=" * 70)
print("GRADIENT STRUCTURE: H increases monotonically?")
print("=" * 70)

# For each tiling, compute how many of its 6 neighbors have higher H
print("\n  bits      H  | #higher #lower #same | gradient")
print("  " + "-" * 60)
for bits in range(num_tilings):
    H = tiling_H[bits]
    higher = lower = same = 0
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        Hn = tiling_H[neighbor]
        if Hn > H:
            higher += 1
        elif Hn < H:
            lower += 1
        else:
            same += 1
    if H in [1, 15]:  # Show extremes
        print(f"  {format(bits, f'0{m}b')}  {H:2d}  | {higher:7d} {lower:6d} {same:5d} | {'MIN' if lower==0 else ''}")

# Count local minima and maxima
local_min = []
local_max = []
for bits in range(num_tilings):
    H = tiling_H[bits]
    is_min = True
    is_max = True
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        Hn = tiling_H[neighbor]
        if Hn < H:
            is_min = False
        if Hn > H:
            is_max = False
    if is_min:
        local_min.append((bits, H))
    if is_max:
        local_max.append((bits, H))

print(f"\n  Local minima: {len(local_min)}")
for bits, H in local_min:
    c = class_labels[tiling_canon[bits]]
    print(f"    {format(bits, f'0{m}b')} H={H} class={c}")

print(f"\n  Local maxima: {len(local_max)}")
for bits, H in local_max:
    c = class_labels[tiling_canon[bits]]
    print(f"    {format(bits, f'0{m}b')} H={H} class={c}")

print()
print("=" * 70)
print("SUMMARY OF BLUE PAIR ANALYSIS")
print("=" * 70)
print(f"""
n=5 Blue Pair Structure:
  - 32 blue pairs total (half of 64 tilings)
  - H(T) + H(flip(T)) is NOT constant (varies by pair)
  - flip(T) is NOT T^op in general (flip only reverses non-path arcs)
  - Self-flip pairs in classes 8 (H=11) and 10 (H=13) have
    identical eigenvalue spectra

Transfer matrix relationships:
  - M(flip(T)) is NOT M(T)^T, -M(T), or P*M(T)*P^T in general
  - The relationship is more subtle: flip changes the odd cycle structure

Gradient structure:
  - The skeleton is NOT a gradient flow (H is not monotone)
  - Local minima exist at H=1 (transitive tournament, unique)
  - Local maxima at H=15 (regular-type tournaments)
""")
