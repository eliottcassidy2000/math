#!/usr/bin/env python3
"""
Test self-converse for Cayley tournaments on the Frobenius group F_21 = Z/7 x| Z/3.

This is the SMALLEST non-abelian group of odd order (order 21).
For abelian groups, inversion always gives self-converse.
For F_21, inversion works ONLY if S is a union of conjugacy classes.

We test whether ALL Cayley tournaments on F_21 are self-converse,
including those with non-normal connection sets.

Group structure: F_21 = <a, b | a^7 = b^3 = 1, bab^{-1} = a^2>
Elements: (i, j) for i in Z/7, j in Z/3
Multiplication: (i1, j1)(i2, j2) = (i1 + 2^{j1} * i2, j1 + j2)
Inverse: (i, j)^{-1} = (-2^{-j} * i, -j) = (-2^{3-j} * i, 3-j) all mod respective

kind-pasteur-2026-03-06-S25e
"""

import sys
from itertools import combinations

# Group operations for F_21 = Z/7 x| Z/3
# Elements: (i, j) for i in Z_7, j in Z_3
# Multiply: (i1,j1)(i2,j2) = (i1 + 2^j1 * i2 mod 7, j1+j2 mod 3)
# Note: 2^0=1, 2^1=2, 2^2=4 (mod 7). 2^3=8=1 mod 7.

POW2 = [1, 2, 4]  # 2^j mod 7 for j=0,1,2

def mul(a, b):
    """Multiply two elements of F_21."""
    return ((a[0] + POW2[a[1]] * b[0]) % 7, (a[1] + b[1]) % 3)

def inv(a):
    """Inverse of element in F_21."""
    j_inv = (-a[1]) % 3
    i_inv = (-POW2[j_inv] * a[0]) % 7
    return (i_inv, j_inv)

# Enumerate all 21 elements
elements = [(i, j) for j in range(3) for i in range(7)]
idx = {e: k for k, e in enumerate(elements)}
e_id = (0, 0)

# Verify group axioms
for a in elements:
    assert mul(a, e_id) == a, f"Right identity fails at {a}"
    assert mul(e_id, a) == a, f"Left identity fails at {a}"
    assert mul(a, inv(a)) == e_id, f"Right inverse fails at {a}"
    assert mul(inv(a), a) == e_id, f"Left inverse fails at {a}"

# Compute conjugacy classes
def conjugate(g, h):
    """h g h^{-1}"""
    return mul(mul(h, g), inv(h))

conj_classes = []
classified = set()
for g in elements:
    if g in classified:
        continue
    cls = set()
    for h in elements:
        cls.add(conjugate(g, h))
    conj_classes.append(frozenset(cls))
    classified.update(cls)

print("Conjugacy classes of F_21:")
for i, cls in enumerate(conj_classes):
    print(f"  C{i}: {sorted(cls)} (size {len(cls)})")

# Pair up non-identity elements: {g, g^{-1}}
nonid = [e for e in elements if e != e_id]
pairs = []
seen = set()
for g in nonid:
    if g not in seen:
        gi = inv(g)
        pairs.append((g, gi))
        seen.add(g)
        seen.add(gi)

print(f"\nInverse pairs ({len(pairs)} total):")
for g, gi in pairs:
    print(f"  {g} <-> {gi}")

# Check which pairs have both elements in the same conjugacy class
# (If they do, choosing one from the pair breaks the conjugacy class)
for g, gi in pairs:
    same_class = any(g in cls and gi in cls for cls in conj_classes)
    if same_class:
        print(f"  WARNING: {g} and {gi} are in the same conjugacy class!")

# Total connection sets: 2^10 = 1024
# (10 pairs, choose one from each)
n = 21
total_sets = 2 ** len(pairs)
print(f"\nTotal connection sets: {total_sets}")

# For each connection set, check if it's normal (union of conj classes)
# and if the tournament is self-converse.
# For n=21, brute-force isomorphism (21! permutations) is impossible.
# Instead, we use a smarter approach:
# 1. Build the adjacency matrix
# 2. Use a canonical form based on sorted neighbor lists
# 3. Compare T and T^op canonical forms

def adj_matrix(S_set):
    """Build adjacency matrix for right Cayley tournament with connection set S."""
    A = [[0]*n for _ in range(n)]
    for a in elements:
        for b in elements:
            if a == b:
                continue
            diff = mul(inv(a), b)  # a^{-1} b
            if diff in S_set:
                A[idx[a]][idx[b]] = 1
    return A

def is_normal_set(S_set):
    """Check if S is a union of conjugacy classes."""
    for cls in conj_classes:
        # Either all or none of cls should be in S
        intersection = cls & S_set
        if len(intersection) > 0 and len(intersection) < len(cls):
            return False
    return True

def canonical_hash(A):
    """Compute a hash that is invariant under isomorphism.
    Uses sorted out-neighbor multiset (degree sequence, then higher-order invariants).
    NOTE: This is NOT a complete isomorphism invariant, but a necessary condition.
    """
    n = len(A)
    # Level 1: out-degrees (all same for regular tournaments)
    out_deg = tuple(sorted(sum(row) for row in A))

    # Level 2: for each vertex, count 3-cycles
    three_cycles = []
    for i in range(n):
        cnt = 0
        for j in range(n):
            if A[i][j] == 0:
                continue
            for k in range(n):
                if A[j][k] == 1 and A[k][i] == 1:
                    cnt += 1
        three_cycles.append(cnt)
    tc = tuple(sorted(three_cycles))

    # Level 3: for each vertex, in-neighbor set structure
    in_out = []
    for i in range(n):
        out_set = frozenset(j for j in range(n) if A[i][j] == 1)
        in_set = frozenset(j for j in range(n) if A[j][i] == 1)
        # How many edges within out-set?
        edges_in_out = sum(1 for j in out_set for k in out_set if A[j][k] == 1)
        in_out.append(edges_in_out)
    io = tuple(sorted(in_out))

    return (out_deg, tc, io)

def transpose_matrix(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]

# Check a few specific connection sets
print("\n" + "=" * 70)
print("SELF-CONVERSE CHECK FOR F_21 CAYLEY TOURNAMENTS")
print("=" * 70)

# Normal sets: those that are unions of conjugacy classes
# Non-identity conjugacy classes: C1, C2, ..., C_k
# For each class C with C^{-1} = C', choose C or C'.
# First, pair up conjugacy classes under inversion
conj_class_pairs = []
conj_seen = set()
for i, cls in enumerate(conj_classes):
    if i in conj_seen or e_id in cls:
        continue
    inv_cls = frozenset(inv(g) for g in cls)
    j = next(k for k, c in enumerate(conj_classes) if c == inv_cls)
    conj_class_pairs.append((i, j))
    conj_seen.add(i)
    conj_seen.add(j)

print(f"\nConjugacy class inverse pairs:")
for i, j in conj_class_pairs:
    print(f"  C{i} (size {len(conj_classes[i])}) <-> C{j} (size {len(conj_classes[j])})")

# Generate normal connection sets
normal_count = 0
normal_hashes = set()
for bits in range(2**len(conj_class_pairs)):
    S = set()
    for k, (i, j) in enumerate(conj_class_pairs):
        if (bits >> k) & 1:
            S.update(conj_classes[i])
        else:
            S.update(conj_classes[j])
    if len(S) == 10:  # (21-1)/2 = 10
        normal_count += 1
        A = adj_matrix(S)
        h = canonical_hash(A)
        normal_hashes.add(h)

print(f"\nNormal connection sets with |S|=10: {normal_count}")

# Now test: for NON-normal S, is T still self-converse?
# We can't check all 2^10 = 1024 sets exhaustively with isomorphism,
# but we can check the canonical hash as a NECESSARY condition.

# If T and T^op have DIFFERENT canonical hashes, they are NOT isomorphic.
# If they have the SAME hash, they MIGHT be isomorphic (need further check).

not_self_converse_candidates = 0
definitely_self_converse = 0
possibly_self_converse = 0
non_normal_count = 0

for bits in range(total_sets):
    S = set()
    for k, (g, gi) in enumerate(pairs):
        if (bits >> k) & 1:
            S.add(g)
        else:
            S.add(gi)

    normal = is_normal_set(S)
    if not normal:
        non_normal_count += 1

    A = adj_matrix(S)
    AT = transpose_matrix(A)

    h_A = canonical_hash(A)
    h_AT = canonical_hash(AT)

    if h_A != h_AT:
        not_self_converse_candidates += 1
        if not normal:
            print(f"  NOT SELF-CONVERSE (by hash): S={sorted(S)}, normal={normal}")
            print(f"    hash(T)  = {h_A}")
            print(f"    hash(T') = {h_AT}")
            if not_self_converse_candidates >= 5:
                print("  ... (stopping after 5 examples)")
                break

print(f"\n  Summary:")
print(f"    Total connection sets: {total_sets}")
print(f"    Non-normal sets: {non_normal_count}")
print(f"    Hash mismatch (T vs T^op): {not_self_converse_candidates}")
if not_self_converse_candidates == 0:
    print(f"    => ALL hashes match! Consistent with all being self-converse.")
    print(f"    (But hash match is necessary, not sufficient for isomorphism)")
else:
    print(f"    => {not_self_converse_candidates} have DIFFERENT hashes => NOT self-converse!")
    print(f"    THIS WOULD DISPROVE 'all VT tournaments are self-converse'!")

# For a stronger test: use the automorphism group structure
# T is self-converse iff there exists sigma in Sym(21) with
# A[sigma(i)][sigma(j)] = A[j][i] for all i,j.
# Since G acts regularly by left multiplication, we can restrict search
# to: sigma = L_g composed with group automorphism alpha, composed with inversion.
# But actually sigma can be ANY permutation.

# A smarter approach: check if T and T^op have the same "orbit partition"
# under the automorphism group.

print("\n" + "=" * 70)
print("DEEPER STRUCTURE TEST")
print("=" * 70)

# Check if ALL F_21 Cayley tournaments have the same 3-cycle count
# (which would indicate they might all be isomorphic to each other)
three_cycle_counts = set()
for bits in range(min(total_sets, 100)):  # Check first 100
    S = set()
    for k, (g, gi) in enumerate(pairs):
        if (bits >> k) & 1:
            S.add(g)
        else:
            S.add(gi)
    A = adj_matrix(S)
    tc = sum(1 for i in range(n) for j in range(n) for k in range(n)
             if A[i][j] == 1 and A[j][k] == 1 and A[k][i] == 1)
    three_cycle_counts.add(tc)

print(f"  Distinct 3-cycle counts (first 100 sets): {sorted(three_cycle_counts)}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
