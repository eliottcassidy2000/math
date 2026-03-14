#!/usr/bin/env python3
"""
block_transitive_family.py — opus-2026-03-14-S71g

The block-transitive tournament family: exact formulas for H.

DEFINITION: A block-transitive tournament T(b₁, b₂, ..., b_k) has:
  - k blocks of sizes b₁, ..., b_k (Σ bᵢ = n)
  - Within each block: some specific tournament (typically a directed cycle)
  - Between blocks: transitive ordering (block i beats block j for i < j)

THEOREM (from our n=9 computation):
  If each block is a directed 3-cycle (bᵢ = 3 for all i):
  T(3, 3, ..., 3) has H = 3^k.

KEY QUESTION: What happens with other block types?
  - b = 1 (single vertex): trivially transitive within, contributes factor 1
  - b = 3 (directed 3-cycle): contributes factor 3
  - b = 5 (tournament on 5 vertices): depends on the specific tournament
  - b = 2 (two vertices): one dominates, contributes factor 1

The "brick product" interpretation:
  H(T(b₁,...,b_k)) = ∏ H_internal(block_i)
  where H_internal is the contribution from each block.

Is this true? Let's verify.
"""

from itertools import permutations, combinations
from math import factorial

def count_hp(A, n):
    """Held-Karp DP."""
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def build_block_tournament(blocks, internal_tournaments=None):
    """Build block-transitive tournament.

    blocks: list of block sizes
    internal_tournaments: list of adjacency matrices for each block
                         (if None, use directed cycle)
    """
    n = sum(blocks)
    A = [[0]*n for _ in range(n)]

    # Build block positions
    positions = []
    pos = 0
    for b in blocks:
        positions.append(list(range(pos, pos + b)))
        pos += b

    # Within-block tournaments
    for idx, (block_pos, b) in enumerate(zip(positions, blocks)):
        if internal_tournaments and internal_tournaments[idx] is not None:
            int_A = internal_tournaments[idx]
            for i in range(b):
                for j in range(b):
                    if i != j:
                        A[block_pos[i]][block_pos[j]] = int_A[i][j]
        else:
            # Default: directed cycle
            for i in range(b):
                A[block_pos[i]][block_pos[(i+1) % b]] = 1

    # Between blocks: transitive (earlier beats later)
    for i in range(len(blocks)):
        for j in range(i+1, len(blocks)):
            for a in positions[i]:
                for b_v in positions[j]:
                    A[a][b_v] = 1

    return A, n

# ============================================================
# Part 1: Pure 3-cycle blocks
# ============================================================

print("=" * 70)
print("BLOCK-TRANSITIVE FAMILY WITH 3-CYCLE BLOCKS")
print("=" * 70)

for k in range(1, 7):
    blocks = [3] * k
    A, n = build_block_tournament(blocks)
    if n <= 18:
        H = count_hp(A, n)
        print(f"  k={k}, n={3*k}: H={H}, 3^k={3**k}, match={H==3**k}")
    else:
        print(f"  k={k}, n={3*k}: too large for DP")

# ============================================================
# Part 2: Mixed block sizes — what's H_internal?
# ============================================================

print(f"\n{'='*70}")
print("INTERNAL H CONTRIBUTION PER BLOCK SIZE")
print(f"{'='*70}")

# H_internal(b) = H of the b-vertex tournament (when embedded in block-transitive)
# For a directed cycle on b vertices:

# Build single-block tournament
for b in range(1, 8):
    # Directed cycle on b vertices
    A_int = [[0]*b for _ in range(b)]
    for i in range(b):
        A_int[i][(i+1) % b] = 1

    if b >= 2:
        H_int = count_hp(A_int, b)
    else:
        H_int = 1

    print(f"  b={b}: H_internal(directed {b}-cycle) = {H_int}")

# Now verify: H(T(b₁,...,b_k)) = ∏ H_internal(bᵢ)?
print(f"\n  Verification: H(block-transitive) = ∏ H_internal?")

test_cases = [
    [1, 3],
    [3, 1],
    [1, 1, 3],
    [3, 3],
    [1, 3, 3],
    [3, 5],
    [5, 3],
    [1, 1, 1, 3],
    [5, 5],
    [3, 3, 3],
    [1, 5],
    [7],
    [3, 7],
    [1, 3, 5],
]

# First get H_internal for each block size
h_internal = {}
for b in range(1, 8):
    A_int = [[0]*b for _ in range(b)]
    for i in range(b):
        A_int[i][(i+1) % b] = 1
    h_internal[b] = count_hp(A_int, b) if b >= 2 else 1

print(f"\n  H_internal: {h_internal}")

for blocks in test_cases:
    n = sum(blocks)
    if n > 15:
        continue
    A, n = build_block_tournament(blocks)
    H = count_hp(A, n)
    predicted = 1
    for b in blocks:
        predicted *= h_internal[b]
    match = H == predicted
    blocks_str = "+".join(str(b) for b in blocks)
    print(f"  T({blocks_str}), n={n}: H={H}, ∏H_int={predicted}, match={match}")

# ============================================================
# Part 3: The product formula — WHY does it work?
# ============================================================

print(f"\n{'='*70}")
print("PRODUCT FORMULA FOR BLOCK-TRANSITIVE TOURNAMENTS")
print(f"{'='*70}")

print("""
THEOREM (Block-Transitive Product Formula):

  For T = T(b₁, ..., b_k) with transitive inter-block ordering:

  H(T) = ∏ᵢ H(Bᵢ) × multinomial adjustment

  Wait — is it really just the product, or does the transitive
  inter-block structure add a multinomial factor?

  Key observation: In a Hamiltonian path of T, all vertices of block i
  must appear BEFORE all vertices of block j if i < j... NO! That's wrong.
  The HP can interleave vertices from different blocks.

  Example: T(1,3) with block 1 = {v} and block 2 = {a,b,c} cycle.
  v beats a, b, c. A path could be: v→a→b→c or v→b→c→a or v→c→a→b.
  Or... can it start from a vertex in block 2? b→c→?→v: c doesn't beat v.
  Only v can be the starting point (it beats everyone in block 2).

  Actually: in block-transitive with transitive inter-block:
  Every vertex in block i beats every vertex in block j (j>i).
  So the HP must visit all vertices of block k (the last) last...
  NO: it can interleave! v1→a1→v2→a2→... if there are arcs.

  But v (block 1) beats a,b,c (block 2). The arcs WITHIN block 2 are
  the cycle a→b→c→a. So from v we can go to a, b, or c.
  v→a→b→c: valid HP. H should be related to internal HPs.

  For block 2 (cycle a→b→c→a): HPs of length 3 = 3 (a→b→c, b→c→a, c→a→b).
  For block 1 (single vertex v): HP of length 1 = 1.
  H(T) = 1 × 3 = 3? Let's check.
""")

A, n = build_block_tournament([1, 3])
H = count_hp(A, n)
print(f"  T(1,3): H = {H}")
print(f"  Product: H_internal(1) × H_internal(3) = {h_internal[1]} × {h_internal[3]} = {h_internal[1]*h_internal[3]}")

# Actually, the paths in T(1,3) must start at v (it beats everyone).
# Then v→a, v→b, or v→c. From each, we traverse the internal cycle.
# v→a→b→c: ✓ (a→b✓, b→c✓)
# v→b→c→a: ✓
# v→c→a→b: ✓
# That's 3. And h_internal(3) = 3. So the product works!

# What about T(3,1)? Block 1 = {a,b,c} cycle, block 2 = {v} (beaten by all).
A, n = build_block_tournament([3, 1])
H = count_hp(A, n)
print(f"  T(3,1): H = {H}")
print(f"  Product: {h_internal[3]} × {h_internal[1]} = {h_internal[3]*h_internal[1]}")

# The paths in T(3,1): must end at v (v is beaten by everyone in block 1).
# a→b→c→v: ✓ (a→b✓, b→c✓, c→v✓)
# b→c→a→v: ✓
# c→a→b→v: ✓
# That's 3. Product works!

# T(3,3): block 1 = {a,b,c} cycle, block 2 = {d,e,f} cycle.
# All of block 1 beats all of block 2.
A, n = build_block_tournament([3, 3])
H = count_hp(A, n)
print(f"  T(3,3): H = {H}")
print(f"  Product: {h_internal[3]}² = {h_internal[3]**2}")

# PROOF of product formula:
# In a block-transitive tournament with blocks B₁ > B₂ > ... > B_k:
# A Hamiltonian path must "eventually" leave each block.
# But actually, can it interleave? In T(3,3):
# Can we go a→d→b→e→c→f? Check: a→d (block1→block2 ✓).
# But d→b? d is in block 2, b is in block 1. Block 1 beats block 2!
# So b→d, not d→b. INVALID.
#
# AH HA! Because block 1 beats block 2, we can go FROM block 1 TO block 2
# but NOT back. So once we leave block 1, we can NEVER return.
# Similarly for any pair of blocks.
#
# This means: the HP must visit ALL of block 1 first, then ALL of block 2, etc.
# The internal ordering within each block is an independent HP.
# So H(T) = ∏ H(Bᵢ). QED!

print(f"""
PROOF OF PRODUCT FORMULA:

In a block-transitive tournament T(B₁, ..., B_k) where block i
beats block j for i < j:

  Every arc from block j to block i goes FROM i TO j (not j to i).
  So once the HP leaves block i, it can NEVER return to block i.

  Therefore: the HP visits ALL vertices of some block first,
  then ALL of another, etc. The inter-block ordering must be
  1, 2, ..., k (since block 1 beats all others, block 2 beats
  blocks 3,...,k, etc.).

  Within each block, the vertices can be traversed in any valid
  HP order for that block's internal tournament.

  The internal HPs are independent.

  Therefore: H(T) = ∏ᵢ H(Bᵢ).  QED.

This explains:
  - H(T(3,...,3)) = 3^k (each 3-cycle has H=3)
  - H(T(b₁,...,b_k)) = ∏ H(Bᵢ)
  - The "brick product" in OCF terms: Ω decomposes into independent blocks
""")

# ============================================================
# Part 4: Which H values are achievable by block-transitive?
# ============================================================

print(f"{'='*70}")
print("ACHIEVABLE H VALUES FROM BLOCK-TRANSITIVE TOURNAMENTS")
print(f"{'='*70}")

# The achievable H values are products of H(B) for various block tournaments.
# H(B) for b-vertex directed cycle:
# b=1: H=1
# b=2: H=1 (one arc, one path)
# b=3: H=3
# b=4: H=? (directed 4-cycle, not a tournament — need to fill in)

# Actually, for b=4: directed 4-cycle 0→1→2→3→0.
# Need to fill in 0-2, 1-3 arcs. Two choices: 0→2 or 2→0, 1→3 or 3→1.
# This gives 4 different tournaments, each potentially different H.

# Let me just enumerate all tournaments on small blocks:
print(f"\n  All H values for b-vertex tournaments:")
for b in range(1, 7):
    h_vals = set()
    if b <= 1:
        h_vals.add(1)
    else:
        total_edges = b * (b-1) // 2
        for bits in range(2**total_edges):
            A = [[0]*b for _ in range(b)]
            idx = 0
            for i in range(b):
                for j in range(i+1, b):
                    if bits & (1 << idx):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    idx += 1
            H = count_hp(A, b)
            h_vals.add(H)
    print(f"    b={b}: H ∈ {sorted(h_vals)}")

# Products of these values: which positive integers are NOT achievable?
# From b=1: {1}
# From b=2: {1}
# From b=3: {1, 3}
# From b=4: {1, 4, 8}... let me check

all_h_by_size = {}
for b in range(1, 7):
    h_vals = set()
    if b <= 1:
        h_vals.add(1)
    else:
        total_edges = b * (b-1) // 2
        for bits in range(2**total_edges):
            A = [[0]*b for _ in range(b)]
            idx = 0
            for i in range(b):
                for j in range(i+1, b):
                    if bits & (1 << idx):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    idx += 1
            H = count_hp(A, b)
            h_vals.add(H)
    all_h_by_size[b] = sorted(h_vals)

# Generate all products up to some bound
from functools import reduce
from operator import mul

# Available factors: union of all H values for blocks of any size
available_factors = set()
for b, vals in all_h_by_size.items():
    available_factors.update(vals)

print(f"\n  Available factors (H values for single blocks ≤ 6): {sorted(available_factors)}")

# Generate products up to 200
achievable = {1}
MAX = 200
for factor in sorted(available_factors):
    if factor <= 1:
        continue
    new_achievable = set()
    for a in achievable:
        prod = a * factor
        while prod <= MAX:
            new_achievable.add(prod)
            prod *= factor
    achievable.update(new_achievable)

# Keep extending with products
changed = True
while changed:
    changed = False
    new_set = set()
    for a in achievable:
        for f in sorted(available_factors):
            if f <= 1:
                continue
            p = a * f
            if p <= MAX and p not in achievable:
                new_set.add(p)
                changed = True
    achievable.update(new_set)

# Which odd numbers ≤ 200 are not achievable?
missing_odd = [v for v in range(1, MAX+1, 2) if v not in achievable]
print(f"\n  Odd H values NOT achievable by block-transitive ≤ {MAX}:")
print(f"    {missing_odd[:30]}")
print(f"    Total: {len(missing_odd)} out of {MAX//2 + 1} odd values")

# Which of these are also NOT achievable by ANY tournament?
# (From our H-spectrum data)
print(f"\n  Note: H=7 and H=21 are known to be impossible for ALL tournaments")
print(f"  Many of the above are achievable by non-block-transitive tournaments")

print(f"\n{'='*70}")
print("ANALYSIS COMPLETE")
print(f"{'='*70}")
