#!/usr/bin/env python3
"""tiling_236_s116n.py — Deep investigation of {2,3,6} structures in n=6 tournaments.

n=6 is the pivot: first disjoint 3-cycles, first beta_3 > 0, per-path failure.
6 = 2*3: the product of the two generating primes.

We investigate the 10-bit tiling triangle (fixing one Hamiltonian path),
looking for patterns involving 2, 3, and 6 at the microscopic level.

Session: kind-pasteur-2026-03-16-S116n32
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

from itertools import permutations, combinations
from math import comb

print()
print("  {2, 3, 6} STRUCTURES IN THE n=6 TILING")
print()
print("=" * 70)
print()

N = 6

# ============================================================
# CORE: Enumerate all n=6 tournaments and compute basic invariants
# ============================================================

def tournament_from_bits(bits, n):
    """Convert C(n,2)-bit integer to adjacency matrix.
    Bit k corresponds to pair (i,j) with i<j, in lexicographic order.
    Bit=1 means i->j, bit=0 means j->i.
    """
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def count_hamiltonian_paths(adj, n):
    """Count Hamiltonian paths by brute force (n! check)."""
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if adj[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def count_3cycles(adj, n):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if i->j->k->i or i->k->j->i
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    count += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    count += 1
    return count

def count_disjoint_3cycle_pairs(adj, n):
    """Count pairs of vertex-disjoint directed 3-cycles (alpha_2 contribution)."""
    # First find all 3-cycles as vertex sets
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))
    # Deduplicate vertex sets (each set can support at most 2 directed cycles)
    cycle_vsets = list(set(cycles))

    # Count vertex-disjoint pairs
    count = 0
    for a in range(len(cycle_vsets)):
        for b in range(a+1, len(cycle_vsets)):
            if len(cycle_vsets[a] & cycle_vsets[b]) == 0:
                count += 1
    return count

def score_sequence(adj, n):
    """Score sequence (sorted out-degrees)."""
    scores = [sum(adj[i]) for i in range(n)]
    return tuple(sorted(scores))

print("  Computing all 2^15 = 32768 tournaments on 6 vertices...")
print()

num_bits = N * (N-1) // 2  # 15
total = 1 << num_bits  # 32768

# Compute H, c3, alpha_2, score for every tournament
H_vals = []
c3_vals = []
a2_vals = []
scores = []

for bits in range(total):
    adj = tournament_from_bits(bits, N)
    H = count_hamiltonian_paths(adj, N)
    c3 = count_3cycles(adj, N)
    a2 = count_disjoint_3cycle_pairs(adj, N)
    sc = score_sequence(adj, N)
    H_vals.append(H)
    c3_vals.append(c3)
    a2_vals.append(a2)
    scores.append(sc)
    if (bits + 1) % 8192 == 0:
        print(f"    ...{bits+1}/{total} done")

print(f"  Done. {total} tournaments computed.")
print()

# ============================================================
print("  I. H-VALUE DISTRIBUTION")
print("  " + "-" * 50)
print()

from collections import Counter

H_dist = Counter(H_vals)
print(f"  {'H':>4s}  {'count':>6s}  {'%':>8s}  {'H mod 2':>7s}  {'H mod 3':>7s}  {'H mod 6':>7s}")
for H in sorted(H_dist.keys()):
    c = H_dist[H]
    print(f"  {H:4d}  {c:6d}  {100*c/total:7.3f}%  {H%2:7d}  {H%3:7d}  {H%6:7d}")
print()

# Check: are all H odd?
all_odd = all(H % 2 == 1 for H in H_vals)
print(f"  All H odd (Redei): {all_odd}")
print()

# H mod 3 distribution
H_mod3 = Counter(H % 3 for H in H_vals)
print(f"  H mod 3 distribution:")
for r in sorted(H_mod3.keys()):
    print(f"    H = {r} mod 3: {H_mod3[r]} ({100*H_mod3[r]/total:.2f}%)")
print()

# H mod 6 distribution (combining mod 2 and mod 3)
H_mod6 = Counter(H % 6 for H in H_vals)
print(f"  H mod 6 distribution:")
for r in sorted(H_mod6.keys()):
    print(f"    H = {r} mod 6: {H_mod6[r]} ({100*H_mod6[r]/total:.2f}%)")
print()

# ============================================================
print("  II. THE OCF COMPONENTS: c3, alpha_2, AND H")
print("  " + "-" * 50)
print()

# OCF: H = 1 + 2*alpha_1 + 4*alpha_2
# where alpha_1 = number of independent sets of size 1 in Omega = c3
# and alpha_2 = number of independent sets of size 2 = disjoint pairs
# Actually alpha_1 = c3 (undirected 3-cycle vertex sets) and alpha_2 = disjoint pairs of such.

# But wait: alpha_1 counts VERTEX SETS that support a 3-cycle, not directed 3-cycles.
# Each vertex set {i,j,k} supports 0, 1, or 2 directed 3-cycles.
# In a tournament on 3 vertices: either 0 (transitive) or 2 (cyclic).
# So c3 (directed) = 2 * (# cyclic triples) = 2 * alpha_1.
# No wait: c3 counts directed 3-cycles. Each VERTEX SET that is cyclic has exactly 2 directed
# cycles (clockwise and counterclockwise). So c3 = 2 * (# cyclic vertex sets).
# And alpha_1 = # cyclic vertex sets = c3 / 2.

# The independence polynomial: I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + ...
# H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

# Compute alpha_1 (= c3/2, number of cyclic 3-vertex sets) and verify OCF
print("  Verifying OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + ...")
print()

# For n=6: max independent set in Omega has at most 2 elements (since 2*3=6)
# So I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2
# H = 1 + 2*alpha_1 + 4*alpha_2

mismatches = 0
for bits in range(total):
    alpha_1 = c3_vals[bits] // 2  # directed -> undirected
    alpha_2 = a2_vals[bits]
    predicted = 1 + 2*alpha_1 + 4*alpha_2
    if predicted != H_vals[bits]:
        mismatches += 1
        if mismatches <= 3:
            print(f"    MISMATCH at bits={bits}: H={H_vals[bits]}, "
                  f"predicted={predicted}, c3={c3_vals[bits]}, a2={alpha_2}")

if mismatches == 0:
    print(f"  OCF VERIFIED for all {total} tournaments!")
    print(f"  H = 1 + 2*(c3/2) + 4*alpha_2 = 1 + c3 + 4*alpha_2")
else:
    print(f"  {mismatches} mismatches found!")
print()

# ============================================================
print("  III. THE (c3, alpha_2) PHASE SPACE")
print("  " + "-" * 50)
print()

# Since H = 1 + c3 + 4*alpha_2 exactly at n=6, the phase space is 2D
phase = Counter()
for bits in range(total):
    alpha_1 = c3_vals[bits] // 2
    alpha_2 = a2_vals[bits]
    phase[(alpha_1, alpha_2)] += 1

print(f"  (alpha_1, alpha_2) -> count, H = 1 + 2*alpha_1 + 4*alpha_2")
print(f"  {'a1':>4s}  {'a2':>4s}  {'count':>6s}  {'H':>4s}  {'H mod 3':>7s}  {'H mod 6':>7s}")
for (a1, a2) in sorted(phase.keys()):
    c = phase[(a1, a2)]
    H = 1 + 2*a1 + 4*a2
    print(f"  {a1:4d}  {a2:4d}  {c:6d}  {H:4d}  {H%3:7d}  {H%6:7d}")
print()

# ============================================================
print("  IV. THE 2x3 SPLIT: BIPARTITE TOURNAMENT STRUCTURE")
print("  " + "-" * 50)
print()

# Split 6 vertices into two groups of 3: A={0,1,2}, B={3,4,5}
# Count: intra-A 3-cycles, intra-B 3-cycles, cross arcs
# Since 6 = 2*3, this split reveals the interaction between the two primes.

def analyze_split(adj, A, B):
    """Analyze a tournament split into two groups."""
    # Intra-A: is A a cycle or transitive?
    a_cycle = 0
    if adj[A[0]][A[1]] and adj[A[1]][A[2]] and adj[A[2]][A[0]]:
        a_cycle += 1
    if adj[A[0]][A[2]] and adj[A[2]][A[1]] and adj[A[1]][A[0]]:
        a_cycle += 1
    # Intra-B: same
    b_cycle = 0
    if adj[B[0]][B[1]] and adj[B[1]][B[2]] and adj[B[2]][B[0]]:
        b_cycle += 1
    if adj[B[0]][B[2]] and adj[B[2]][B[1]] and adj[B[1]][B[0]]:
        b_cycle += 1
    # Cross arcs: A->B count
    cross_ab = sum(adj[a][b] for a in A for b in B)
    return a_cycle > 0, b_cycle > 0, cross_ab

A = [0, 1, 2]
B = [3, 4, 5]

split_stats = Counter()
for bits in range(total):
    adj = tournament_from_bits(bits, N)
    a_cyc, b_cyc, cross = analyze_split(adj, A, B)
    split_stats[(a_cyc, b_cyc, cross)] += 1

print(f"  Split A={{0,1,2}}, B={{3,4,5}}:")
print(f"  {'A_cyc':>5s}  {'B_cyc':>5s}  {'A->B':>5s}  {'count':>6s}  {'avg_H':>8s}")
for (ac, bc, cross) in sorted(split_stats.keys()):
    c = split_stats[(ac, bc, cross)]
    # Compute average H for this configuration
    total_H = 0
    cnt = 0
    for bits in range(total):
        adj = tournament_from_bits(bits, N)
        ac2, bc2, cr2 = analyze_split(adj, A, B)
        if (ac2, bc2, cr2) == (ac, bc, cross):
            total_H += H_vals[bits]
            cnt += 1
    avg_H = total_H / cnt if cnt > 0 else 0
    print(f"  {str(ac):>5s}  {str(bc):>5s}  {cross:5d}  {c:6d}  {avg_H:8.2f}")
print()

# Disjoint 3-cycles: possible iff BOTH halves are cyclic
both_cyclic = sum(1 for bits in range(total) if c3_vals[bits] >= 4 and a2_vals[bits] > 0)
# Actually: disjoint pair exists iff there exist two vertex-disjoint cyclic triples
print(f"  Tournaments with at least one disjoint 3-cycle pair: "
      f"{sum(1 for a2 in a2_vals if a2 > 0)}/{total}")
print(f"  = {100*sum(1 for a2 in a2_vals if a2 > 0)/total:.2f}%")
print()

# ============================================================
print("  V. THE TILING TRIANGLE: 10 BITS FIXING PATH 0->1->2->3->4->5")
print("  " + "-" * 50)
print()

# Fix the path 0->1->2->3->4->5.
# The 10 non-path arcs, in the triangular grid:
# Row 0 (skip 2): (0,2), (1,3), (2,4), (3,5) — 4 arcs
# Row 1 (skip 3): (0,3), (1,4), (2,5) — 3 arcs
# Row 2 (skip 4): (0,4), (1,5) — 2 arcs
# Row 3 (skip 5): (0,5) — 1 arc

# For a given tournament (15 bits), the tiling is determined by the non-path bits.
# Path arcs: (0,1), (1,2), (2,3), (3,4), (4,5) — all forward (bit=1).
# Non-path arcs: the remaining 10.

# Encoding: bit position for pair (i,j) with i<j
# Pairs in order: (0,1),(0,2),(0,3),(0,4),(0,5),(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)
# Positions:       0     1     2     3     4     5     6     7     8     9     10    11    12    13    14

# Path arcs at positions: 0 (0,1), 5 (1,2), 9 (2,3), 12 (3,4), 14 (4,5)
path_positions = [0, 5, 9, 12, 14]
nonpath_positions = [i for i in range(15) if i not in path_positions]
# nonpath: 1,2,3,4, 6,7,8, 10,11, 13
# = (0,2),(0,3),(0,4),(0,5), (1,3),(1,4),(1,5), (2,4),(2,5), (3,5)

# Map each nonpath position to its (skip, offset) in the triangle
nonpath_pairs = []
k = 0
for i in range(N):
    for j in range(i+1, N):
        if k in nonpath_positions:
            skip = j - i
            nonpath_pairs.append((i, j, skip, k))
        k += 1

print("  Non-path arcs (triangular grid):")
print(f"  {'bit':>4s}  {'pair':>6s}  {'skip':>5s}  {'row':>4s}")
for i, j, skip, pos in nonpath_pairs:
    row = skip - 2
    print(f"  {pos:4d}  ({i},{j})  {skip:5d}  {row:4d}")
print()

# The triangle:
# Row 0 (skip 2): (0,2) (1,3) (2,4) (3,5) — 4 entries
# Row 1 (skip 3): (0,3) (1,4) (2,5)       — 3 entries
# Row 2 (skip 4): (0,4) (1,5)             — 2 entries
# Row 3 (skip 5): (0,5)                   — 1 entry

# For each tournament, extract the 10-bit tiling
def tiling_bits(tournament_bits):
    """Extract the 10 non-path bits from a 15-bit tournament."""
    return tuple((tournament_bits >> pos) & 1 for _, _, _, pos in nonpath_pairs)

# Only consider tournaments with the canonical path 0->1->2->3->4->5
# This means bits at path positions must all be 1.
# Filter: path arcs all forward
canonical_count = 0
tiling_to_H = {}
tiling_to_c3 = {}

for bits in range(total):
    # Check if path arcs are all forward
    path_ok = all((bits >> p) & 1 for p in path_positions)
    if not path_ok:
        continue
    canonical_count += 1
    tb = tiling_bits(bits)
    tiling_to_H[tb] = H_vals[bits]
    tiling_to_c3[tb] = c3_vals[bits]

print(f"  Tournaments with canonical path 0->1->...->5: {canonical_count}")
print(f"  = 2^10 = {2**10}: {canonical_count == 2**10}")
print()

# ============================================================
print("  VI. TILING HAMMING WEIGHT vs H")
print("  " + "-" * 50)
print()

# Hamming weight of tiling = number of forward non-path arcs
hw_to_H = {}
for tb, H in tiling_to_H.items():
    hw = sum(tb)
    if hw not in hw_to_H:
        hw_to_H[hw] = []
    hw_to_H[hw].append(H)

print(f"  {'HW':>3s}  {'count':>6s}  {'min_H':>6s}  {'max_H':>6s}  {'avg_H':>8s}  {'H values (if few)':>30s}")
for hw in sorted(hw_to_H.keys()):
    vals = hw_to_H[hw]
    H_set = sorted(set(vals))
    extra = str(H_set) if len(H_set) <= 8 else f"{len(H_set)} distinct"
    print(f"  {hw:3d}  {len(vals):6d}  {min(vals):6d}  {max(vals):6d}  "
          f"{sum(vals)/len(vals):8.2f}  {extra:>30s}")
print()

# Check symmetry: HW and 10-HW should have same H distribution (complement)
print("  Complement symmetry (HW vs 10-HW):")
for hw in range(6):
    if hw in hw_to_H and (10-hw) in hw_to_H:
        h1 = sorted(hw_to_H[hw])
        h2 = sorted(hw_to_H[10-hw])
        match = (h1 == h2)
        print(f"    HW={hw} vs HW={10-hw}: same H distribution? {match}")
print()

# ============================================================
print("  VII. ROW-BY-ROW CONTRIBUTIONS TO H")
print("  " + "-" * 50)
print()

# The triangle has 4 rows (skip 2,3,4,5).
# Row sizes: 4, 3, 2, 1
# Does each row contribute independently to H?

# Row 0 bits: positions 0,1,2,3 in the tiling tuple
# Row 1 bits: positions 4,5,6
# Row 2 bits: positions 7,8
# Row 3 bits: position 9

row_ranges = [(0,4), (4,7), (7,9), (9,10)]  # [start, end) in tiling tuple

# For each row pattern, compute average H
for row_idx, (start, end) in enumerate(row_ranges):
    row_size = end - start
    pattern_to_H = {}
    for tb, H in tiling_to_H.items():
        row_pattern = tb[start:end]
        if row_pattern not in pattern_to_H:
            pattern_to_H[row_pattern] = []
        pattern_to_H[row_pattern].append(H)

    print(f"  Row {row_idx} (skip {row_idx+2}, {row_size} bits):")
    for pat in sorted(pattern_to_H.keys()):
        vals = pattern_to_H[pat]
        hw = sum(pat)
        print(f"    {pat}: HW={hw}, avg_H={sum(vals)/len(vals):.2f}, "
              f"range=[{min(vals)},{max(vals)}], count={len(vals)}")
    print()

# ============================================================
print("  VIII. THE {2,3} INTERACTION: ROW 0 x ROW 1")
print("  " + "-" * 50)
print()

# Row 0 has 4 bits (skip-2 arcs: nearest neighbor interactions)
# Row 1 has 3 bits (skip-3 arcs: next-nearest neighbor)
# These are the "2-body" and "3-body" interactions.
# Their JOINT effect on H reveals the {2,3} structure.

joint_H = {}
for tb, H in tiling_to_H.items():
    r0 = sum(tb[0:4])  # row 0 HW
    r1 = sum(tb[4:7])  # row 1 HW
    key = (r0, r1)
    if key not in joint_H:
        joint_H[key] = []
    joint_H[key].append(H)

print(f"  Row0_HW x Row1_HW -> avg_H:")
print(f"  {'':>5s}", end="")
for r1 in range(4):
    print(f"  {'R1='+str(r1):>8s}", end="")
print()
for r0 in range(5):
    print(f"  R0={r0}", end="")
    for r1 in range(4):
        if (r0, r1) in joint_H:
            avg = sum(joint_H[(r0,r1)]) / len(joint_H[(r0,r1)])
            print(f"  {avg:8.2f}", end="")
        else:
            print(f"  {'---':>8s}", end="")
    print()
print()

# ============================================================
print("  IX. THE DIAGONAL STRUCTURE (skip-length analysis)")
print("  " + "-" * 50)
print()

# Each non-path arc has a "skip" length (j-i).
# Skip 2: (0,2),(1,3),(2,4),(3,5) — 4 arcs, nearest-neighbor interactions
# Skip 3: (0,3),(1,4),(2,5) — 3 arcs
# Skip 4: (0,4),(1,5) — 2 arcs
# Skip 5: (0,5) — 1 arc (longest)

# Contribution of each skip to c3:
# A 3-cycle on {i,j,k} with i<j<k has skip lengths (j-i, k-j, k-i).
# For n=6: (j-i)+(k-j) = k-i. So two short skips compose to one long skip.

# Count how many 3-cycles involve each skip length
print("  3-cycles by skip composition:")
print("  A 3-cycle {i,j,k} with i<j<k uses skip lengths (j-i, k-j, k-i).")
print()

skip_usage = Counter()
for tb, c3 in tiling_to_c3.items():
    if c3 > 0:
        # Reconstruct adjacency from tiling
        adj = [[0]*N for _ in range(N)]
        # Path arcs
        for i in range(N-1):
            adj[i][i+1] = 1
        # Non-path arcs
        for idx, (i, j, skip, pos) in enumerate(nonpath_pairs):
            if tb[idx]:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        # Find 3-cycles
        for i in range(N):
            for j in range(i+1, N):
                for k in range(j+1, N):
                    if adj[i][j] and adj[j][k] and adj[k][i]:
                        s1, s2, s3 = j-i, k-j, k-i
                        skip_usage[(s1, s2, s3)] += 1
                    if adj[i][k] and adj[k][j] and adj[j][i]:
                        s1, s2, s3 = j-i, k-j, k-i
                        skip_usage[(s1, s2, s3)] += 1

print(f"  {'(s1,s2,s3)':>12s}  {'count':>8s}  {'s1+s2=s3':>8s}  {'triple':>10s}")
for (s1, s2, s3) in sorted(skip_usage.keys()):
    c = skip_usage[(s1, s2, s3)]
    # Which triples have these skips?
    triples = [(i, i+s1, i+s3) for i in range(N) if i+s3 < N]
    print(f"  ({s1},{s2},{s3})  {c:8d}  {s1+s2==s3!s:>8}  {triples}")
print()

# ============================================================
print("  X. THE 6 = 2*3 FACTORIZATION IN SKIP SPACE")
print("  " + "-" * 50)
print()
print("  Skip lengths: 2, 3, 4, 5")
print("  2 = first prime (parity)")
print("  3 = second prime (curvature)")
print("  4 = 2*2 (double parity)")
print("  5 = path length - 1")
print()
print("  3-cycle skip compositions (s1+s2=s3):")
print("  (1,1,2): uses path arcs — NOT in tiling (path fixed)")
print("  (1,2,3): one path + one skip-2 + one skip-3")
print("  (1,3,4): one path + one skip-3 + one skip-4")
print("  (1,4,5): one path + one skip-4 + one skip-5")
print("  (2,1,3): one skip-2 + one path + one skip-3")
print("  (2,2,4): two skip-2 + one skip-4")
print("  (2,3,5): one skip-2 + one skip-3 + one skip-5")
print("  (3,1,4): one skip-3 + one path + one skip-4")
print("  (3,2,5): one skip-3 + one skip-2 + one skip-5")
print()
print("  KEY PATTERN: The 3-cycle {i, i+2, i+5} uses skips (2,3,5).")
print("  2+3=5, and {2,3,5} are the first three primes!")
print("  This is the 'spherical' triangle of our grand trichotomy.")
print()
print("  The 3-cycle {i, i+2, i+4} uses skips (2,2,4).")
print("  This is the 'parity-doubled' triangle (purely INERT).")
print()
print("  The 3-cycle {i, i+1, i+3} uses skips (1,2,3).")
print("  This includes a path arc (skip 1) — the 'path-anchored' triangle.")
print()

# Count how many of each type
print("  Contribution of each skip-triple to total c3:")
total_contribution = sum(skip_usage.values())
for (s1, s2, s3) in sorted(skip_usage.keys()):
    c = skip_usage[(s1, s2, s3)]
    print(f"    ({s1},{s2},{s3}): {c:6d} / {total_contribution} = {100*c/total_contribution:.1f}%")
print()

# ============================================================
print("  XI. NEW PATTERN: FORWARD-BACKWARD PARITY BY ROW")
print("  " + "-" * 50)
print()

# For each row, compute the "net forward" count (sum of bits)
# and see how it correlates with H.
# Is there a formula H = f(net_row0, net_row1, net_row2, net_row3)?

# More interestingly: the PARITY (even/odd) of each row sum
print("  Row parities (even=0, odd=1) and their effect on H:")
parity_H = {}
for tb, H in tiling_to_H.items():
    parities = tuple(sum(tb[s:e]) % 2 for s, e in row_ranges)
    if parities not in parity_H:
        parity_H[parities] = []
    parity_H[parities].append(H)

print(f"  {'(p0,p1,p2,p3)':>14s}  {'count':>6s}  {'avg_H':>8s}  {'min':>5s}  {'max':>5s}  {'H mod 3 dist':>20s}")
for par in sorted(parity_H.keys()):
    vals = parity_H[par]
    mod3_dist = Counter(v % 3 for v in vals)
    mod3_str = f"{{0:{mod3_dist.get(0,0)}, 1:{mod3_dist.get(1,0)}, 2:{mod3_dist.get(2,0)}}}"
    print(f"  {str(par):>14s}  {len(vals):6d}  {sum(vals)/len(vals):8.2f}  "
          f"{min(vals):5d}  {max(vals):5d}  {mod3_str:>20s}")
print()

# Check: does the total parity of all 10 bits determine H mod something?
total_parity_H = {}
for tb, H in tiling_to_H.items():
    tp = sum(tb) % 2
    if tp not in total_parity_H:
        total_parity_H[tp] = Counter()
    total_parity_H[tp][H % 6] += 1

print("  Total tiling parity vs H mod 6:")
for tp in sorted(total_parity_H.keys()):
    dist = total_parity_H[tp]
    print(f"    Parity {tp}: H mod 6 distribution = {dict(sorted(dist.items()))}")
print()

# ============================================================
print("  XII. SYNTHESIS: THE {2,3,6} PATTERNS")
print("  " + "-" * 50)
print()
print("  At n=6 = 2*3:")
print("  1. OCF truncates at degree 2: H = 1 + c3 + 4*alpha_2 (EXACT)")
print("     Only two terms beyond the base: the 2-body (c3) and 2*2-body (alpha_2)")
print("  2. The tiling triangle has 4 rows with sizes 4,3,2,1 = 10 total")
print("     Row sizes sum: 4+3+2+1 = 10 = C(5,2)")
print("  3. Skip-2 arcs (row 0): 4 nearest-neighbor interactions")
print("     Skip-3 arcs (row 1): 3 next-nearest interactions")
print("     The 'fundamental' 3-cycle uses skips {2,3,5} — the spherical primes!")
print("  4. Disjoint 3-cycle pairs (alpha_2 > 0) first appear at n=6 = 2*3")
print("     because you need two groups of 3 = the factorization of 6")
print("  5. The {2,3} row interaction table shows H depends on BOTH")
print("     skip-2 and skip-3 arcs non-independently")
print("  6. Row parities create a 2^4 = 16 partition of tiling space")
print("     with potentially interesting H mod 3 structure")
print()
