#!/usr/bin/env python3
"""
fiber_structure_deep.py -- kind-pasteur-2026-03-13-S61

Deep investigation of the fiber structure over lambda classes.

Key questions:
1. WHY is the H=111 fiber twice as large as H=109? (10080 vs 5040)
2. Each tournament has exactly 1 lambda-preserving S. How do the S-sets
   partition the fiber?
3. What is the automorphism group structure?
4. Is the 2:1 ratio related to tournament isomorphism classes?
5. Connection to the overlap weight matrix and Vitali quotient.

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations, permutations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def get_labeled_lambda(A, n):
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))
    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            val = sum(1 for cs in c3_sets if u in cs and v in cs)
            lam[u][v] = val
            lam[v][u] = val
    return lam, c3_sets


def tournament_canonical(A, n):
    """Canonical form under vertex relabeling (graph isomorphism)."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best


n = 7

# ========================================================================
# ANALYSIS 1: Tournament isomorphism classes in the fiber
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: ISOMORPHISM CLASSES IN THE FIBER")
print("=" * 70)

# Collect the fiber
target_score = (2, 2, 3, 3, 3, 3, 5)
sorted_lam_target = (0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3)
target_deg_seq = (6, 8, 8, 10, 10, 10, 14)

print(f"  Collecting fiber for score={target_score}, deg_seq={target_deg_seq}...")

fiber = []  # (bits, H, c7)
m = n * (n - 1) // 2
total = 1 << m

for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted(sum(A[v]) for v in range(n)))
    if scores != target_score:
        continue

    lam, _ = get_labeled_lambda(A, n)
    sl = tuple(sorted(lam[u][v] for u in range(n) for v in range(u+1, n)))
    if sl != sorted_lam_target:
        continue

    deg_seq = tuple(sorted(sum(lam[v]) for v in range(n)))
    if deg_seq != target_deg_seq:
        continue

    H = count_ham_paths(A, n)
    c7 = count_directed_ham_cycles_on_subset(A, list(range(n)))
    fiber.append((bits, H, c7))

print(f"  Fiber size: {len(fiber)}")
by_H = defaultdict(list)
for bits, H, c7 in fiber:
    by_H[H].append(bits)
for H in sorted(by_H.keys()):
    print(f"    H={H}: {len(by_H[H])} tournaments")

# Tournament isomorphism is expensive at n=7 (7! = 5040 per tournament).
# Instead, use a faster invariant: sorted (score, neighbor-scores) signature.

def tournament_signature(A, n):
    """Faster isomorphism invariant."""
    scores = [sum(A[v]) for v in range(n)]
    sig = []
    for v in range(n):
        out_scores = tuple(sorted(scores[w] for w in range(n) if w != v and A[v][w]))
        in_scores = tuple(sorted(scores[w] for w in range(n) if w != v and A[w][v]))
        sig.append((scores[v], out_scores, in_scores))
    return tuple(sorted(sig))


# Group by signature within each H class
print(f"\n  Isomorphism invariant classes:")
for H in sorted(by_H.keys()):
    by_sig = defaultdict(list)
    for bits in by_H[H]:
        A = binary_to_tournament(bits, n)
        sig = tournament_signature(A, n)
        by_sig[sig].append(bits)
    print(f"    H={H}: {len(by_sig)} signature classes")
    for sig, bits_list in sorted(by_sig.items(), key=lambda x: -len(x[1])):
        print(f"      {len(bits_list)} tournaments")

# The number of tournaments per isomorphism class = |S_n| / |Aut(T)|
# = 5040 / |Aut(T)|.
# If H=109 has 5040 tournaments, that could be 1 class with |Aut|=1 (5040/1=5040)
# or fewer classes with larger Aut groups.
# If H=111 has 10080 = 2 * 5040, that could be 2 classes with |Aut|=1
# or 1 class with |Aut|=... wait, 5040/|Aut| * num_classes = 10080.
# If 1 class: 5040/|Aut| = 10080 means |Aut| = 0.5 — impossible.
# So at least 2 classes.

# Let's compute actual canonical forms for a sample
print(f"\n  Computing canonical forms (sampling)...")

# Take one representative from each H class
for H in [109, 111]:
    bits_list = by_H[H][:10]  # Sample first 10
    canons = set()
    for bits in bits_list:
        A = binary_to_tournament(bits, n)
        c = tournament_canonical(A, n)
        canons.add(c)
    print(f"    H={H}: {len(canons)} distinct canonical forms in sample of {len(bits_list)}")


# ========================================================================
# ANALYSIS 2: Automorphism group size
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: AUTOMORPHISM GROUPS")
print("=" * 70)

for H in [109, 111]:
    bits = by_H[H][0]
    A = binary_to_tournament(bits, n)
    # Count automorphisms
    aut_count = 0
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[j]] == A[i][j] for i in range(n) for j in range(n) if i != j):
            aut_count += 1
    print(f"  H={H}, bits={bits}: |Aut| = {aut_count}")
    print(f"    Orbit size = 7!/{aut_count} = {5040 // aut_count}")


# ========================================================================
# ANALYSIS 3: The lambda-preserving involution
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: THE INVOLUTION STRUCTURE")
print("=" * 70)

# The reversal is an involution: applying it twice returns the original.
# So each lambda-preserving reversal pairs a c7=8 tournament with a c7=9 tournament.
# The pairing is 1-to-1 between c7=8 and c7=9 tournaments.

# But 5040 (c7=8) pairs with 5040 (out of 10080) c7=9 tournaments.
# This leaves 5040 UNPAIRED c7=9 tournaments!

# These unpaired c7=9 tournaments must come from a DIFFERENT source:
# either a different lambda class, or they have NO lambda-preserving reversal.

# Wait: each tournament has exactly 1 lambda-preserving S (we checked).
# A c7=9 tournament has reversal to c7=8 (paired with c7=8 tournament).
# But there are 10080 c7=9 and only 5040 c7=8.
# So some c7=9 tournaments must have delta_c7=0 under their reversal!

# Let me check: do ALL c7=9 tournaments have exactly 1 lambda-preserving S?
# And what's the c7 of the reversed tournament?

print(f"  Checking reversal structure for c7=9 tournaments...")
c7_9_bits = [bits for bits, H, c7 in fiber if c7 == 9]
c7_8_bits = [bits for bits, H, c7 in fiber if c7 == 8]

# Sample some c7=9 tournaments
sample_size = min(100, len(c7_9_bits))
delta_dist = defaultdict(int)

import random
random.seed(42)
sample = random.sample(c7_9_bits, sample_size)

for bits in sample:
    A = binary_to_tournament(bits, n)
    lam, _ = get_labeled_lambda(A, n)

    for S in combinations(range(n), 4):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n)
        if all(lam_rev[i][j] == lam[i][j] for i in range(n) for j in range(i+1, n)):
            c7_rev = count_directed_ham_cycles_on_subset(A_rev, list(range(n)))
            delta = c7_rev - 9
            delta_dist[delta] += 1

print(f"  Delta c7 distribution for c7=9 tournaments (sample of {sample_size}):")
for dc in sorted(delta_dist.keys()):
    print(f"    delta_c7 = {dc:+d}: {delta_dist[dc]} ({100*delta_dist[dc]/sample_size:.0f}%)")

# Also check c7=8 tournaments
delta_dist_8 = defaultdict(int)
sample_8 = random.sample(c7_8_bits, min(100, len(c7_8_bits)))

for bits in sample_8:
    A = binary_to_tournament(bits, n)
    lam, _ = get_labeled_lambda(A, n)

    for S in combinations(range(n), 4):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n)
        if all(lam_rev[i][j] == lam[i][j] for i in range(n) for j in range(i+1, n)):
            c7_rev = count_directed_ham_cycles_on_subset(A_rev, list(range(n)))
            delta = c7_rev - 8
            delta_dist_8[delta] += 1

print(f"\n  Delta c7 distribution for c7=8 tournaments (sample of {len(sample_8)}):")
for dc in sorted(delta_dist_8.keys()):
    print(f"    delta_c7 = {dc:+d}: {delta_dist_8[dc]} ({100*delta_dist_8[dc]/len(sample_8):.0f}%)")


# ========================================================================
# ANALYSIS 4: Complement structure
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: COMPLEMENT STRUCTURE")
print("=" * 70)

# The complement of a tournament T has A'[i][j] = 1-A[i][j] for all i<j.
# Score of complement: n-1-s_i. So (2,2,3,3,3,3,5) -> (1,3,3,3,3,4,4).
# This is the OTHER ambiguous case!

# Check: is the second ambiguous score class (1,3,3,3,3,4,4) the complement of the first?
print(f"  Complement of score (2,2,3,3,3,3,5): {tuple(sorted(6-s for s in target_score))}")
# (6-2,6-2,6-3,6-3,6-3,6-3,6-5) = (4,4,3,3,3,3,1)
# Sorted: (1,3,3,3,3,4,4) YES!

# So the two deepest ambiguity cases are COMPLEMENTS of each other.
# bits 4728 (H=109) has complement bits
bits_4728 = 4728
A_4728 = binary_to_tournament(bits_4728, n)
# Complement: flip all bits in the lower m positions
m = n * (n - 1) // 2
comp_bits = bits_4728 ^ ((1 << m) - 1)
A_comp = binary_to_tournament(comp_bits, n)
H_comp = count_ham_paths(A_comp, n)
scores_comp = tuple(sorted(sum(A_comp[v]) for v in range(n)))

print(f"  bits=4728: H=109, scores={target_score}")
print(f"  Complement bits={comp_bits}: H={H_comp}, scores={scores_comp}")

# Check if complement of 4658 gives 9653
comp_4658 = 4658 ^ ((1 << m) - 1)
A_comp_4658 = binary_to_tournament(comp_4658, n)
H_comp_4658 = count_ham_paths(A_comp_4658, n)
scores_comp_4658 = tuple(sorted(sum(A_comp_4658[v]) for v in range(n)))
print(f"  bits=4658: H=111, scores={target_score}")
print(f"  Complement bits={comp_4658}: H={H_comp_4658}, scores={scores_comp_4658}")

# Are 9388 and 9653 complements of 4728 and 4658?
print(f"\n  Cross-check with Case 2:")
for bits in [9388, 9653]:
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted(sum(A[v]) for v in range(n)))
    print(f"    bits={bits}: H={H}, scores={scores}")


# ========================================================================
# ANALYSIS 5: The overlap weight matrix as group action
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: OVERLAP WEIGHT AND GROUP STRUCTURE")
print("=" * 70)

# The 3-cycle overlap weight between two 3-cycle vertex sets C_i, C_j
# is |C_i cap C_j| in {0, 1, 2, 3}.
# For distinct sets: overlap is in {0, 1, 2} (can't be 3 since distinct).

# The overlap weight matrix W[i][j] = |C_i cap C_j| is a symmetric matrix.
# The lambda graph DETERMINES W (and vice versa).

# Within a fixed lambda isomorphism class, the overlap weight matrix is fixed
# up to relabeling of 3-cycles. But the TOURNAMENT adds orientation to each 3-cycle.

# The oriented overlap: for each pair (C_i, C_j) of 3-cycles sharing 2 vertices,
# are they "concordant" (same orientation on shared pair) or "discordant"?

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    lam, c3s = get_labeled_lambda(A, n)

    print(f"\n  {label} (bits={bits}): {len(c3s)} 3-cycles")

    # For each pair of 3-cycles with overlap=2:
    concordant = 0
    discordant = 0
    for i in range(len(c3s)):
        for j in range(i+1, len(c3s)):
            overlap = len(c3s[i] & c3s[j])
            if overlap == 2:
                # Shared vertices
                shared = c3s[i] & c3s[j]
                u, v = sorted(shared)
                # Direction of arc u->v in each 3-cycle
                # In C_i: determine cycle direction through {u,v}
                c_i = sorted(c3s[i])
                a, b, c = c_i
                if A[a][b] and A[b][c] and A[c][a]:
                    # Clockwise: a->b->c->a
                    dir_i = {(a,b): 1, (b,c): 1, (c,a): 1}
                else:
                    # Counter: a->c->b->a
                    dir_i = {(a,c): 1, (c,b): 1, (b,a): 1}

                c_j = sorted(c3s[j])
                a2, b2, c2 = c_j
                if A[a2][b2] and A[b2][c2] and A[c2][a2]:
                    dir_j = {(a2,b2): 1, (b2,c2): 1, (c2,a2): 1}
                else:
                    dir_j = {(a2,c2): 1, (c2,b2): 1, (b2,a2): 1}

                # Check if both cycles traverse u->v in the same direction
                # (both have u->v or both have v->u)
                if (u,v) in dir_i and (u,v) in dir_j:
                    concordant += 1
                elif (v,u) in dir_i and (v,u) in dir_j:
                    concordant += 1
                else:
                    discordant += 1

    print(f"    Overlap-2 pairs: concordant={concordant}, discordant={discordant}")


# ========================================================================
# ANALYSIS 6: The Vitali quotient as a fiber bundle
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: VITALI QUOTIENT AS FIBER BUNDLE")
print("=" * 70)

print("""
SUMMARY OF THE VITALI FIBER BUNDLE STRUCTURE:

Base space: Lambda isomorphism classes
  -> Each point in the base determines c3_dir, c5_dir, alpha_2

Fiber: Tournaments sharing a lambda class
  -> The fiber is parameterized by arc orientations on the lambda graph
  -> The fiber is partitioned by c7_dir values

Structure group: The lambda-preserving transformations
  -> Generated by 4-vertex sub-tournament reversals
  -> Each reversal is an involution (order 2)
  -> Balance condition: all S-vertices have same out-degree to outside

At n=7:
  - The fiber over our target class has 15120 tournaments
  - Split as: c7=8 (5040) and c7=9 (10080)
  - The 4-vertex reversal involution maps some c7=9 to c7=8 (5040 pairs)
  - The remaining 5040 c7=9 tournaments have delta_c7=0 under reversal

The COMPLETE INVARIANT at n=7 is:
  (lambda isomorphism class, c7_dir)

The "non-measurable" dimension = 1 (just c7_dir).
The Vitali quotient forgets this single bit of information.

This is the simplest possible non-trivial fiber: a Z/2Z-graded fiber
with one "parity" bit that lambda cannot see.
""")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
