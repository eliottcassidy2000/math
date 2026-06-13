#!/usr/bin/env python3
"""
h21_correct_structural.py — Correct structural analysis for H=21.

Using DIRECTED odd cycles as Ω vertices (not vertex sets).
H=21 requires α₁=10, α₂=0 (all pairs intersect).
Key: 5-element sets can contribute MULTIPLE directed cycles.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def fast_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if not c or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def count_directed_odd_cycles(A, n):
    """Count directed odd cycles, grouped by vertex set."""
    groups = defaultdict(int)
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    groups[frozenset(verts)] += 1
    return groups

print("=" * 70)
print("CORRECT H=21 STRUCTURAL ANALYSIS (DIRECTED CYCLES)")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: α₁ distribution at n=5,6 (CORRECT: counting directed cycles)
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: Correct α₁ (directed cycle count) vs H ---")

for n in [5, 6]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    ne = len(edges)

    alpha1_to_H = defaultdict(set)
    H_to_alpha1 = defaultdict(set)

    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        H = fast_hp(A, n)
        groups = count_directed_odd_cycles(A, n)
        alpha1 = sum(groups.values())  # TOTAL directed cycles

        alpha1_to_H[alpha1].add(H)
        H_to_alpha1[H].add(alpha1)

    print(f"\n  n={n}:")
    # H near 21
    for h in [15, 17, 19, 21, 23, 25]:
        if h in H_to_alpha1:
            print(f"    H={h} ← α₁ ∈ {sorted(H_to_alpha1[h])}")
        else:
            print(f"    H={h}: NOT ACHIEVED")

    # Does α₁=10 exist?
    if 10 in alpha1_to_H:
        print(f"    α₁=10 → H ∈ {sorted(alpha1_to_H[10])}")
    else:
        print(f"    α₁=10: NOT ACHIEVED")

# ═══════════════════════════════════════════════════════════════════
# Part 2: At n=5, detailed I.P. structure
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: n=5 detailed --- ")
print("  Key: 5-element vertex set can have up to 12 directed Hamiltonian cycles")
print("  (actually at most 3 for a tournament on 5 vertices)")

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

# How many directed Ham cycles does each tournament on 5 vertices have?
hc_counts = Counter()
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Count directed 5-cycles
    v0 = 0
    count_5 = 0
    for perm in permutations(range(1, n)):
        cycle = (0,) + perm
        ok = True
        for i in range(n):
            if A[cycle[i]][cycle[(i+1) % n]] != 1:
                ok = False
                break
        if ok:
            count_5 += 1

    hc_counts[count_5] += 1

print(f"  n=5 tournaments by #directed 5-cycles: {dict(sorted(hc_counts.items()))}")

# ═══════════════════════════════════════════════════════════════════
# Part 3: At n=5, does H=7 become possible with correct counting?
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: H=7 with correct directed-cycle α₁ ---")
print("  H=7 requires α₁=3, all pairwise intersecting.")
print("  With directed cycles: 3 could be:")
print("    (a) 3 directed 3-cycles (same as vertex-set analysis)")
print("    (b) 2 directed 3-cycles + 1 directed 5-cycle")
print("    (c) 1 directed 3-cycle + 2 directed 5-cycles on same vertex set")
print("    (d) 3 directed 5-cycles on same vertex set")
print()
print("  Case (c) and (d): at n=5, a subtournament has ≤3 directed 5-cycles.")
print("  3 cycles on same 5-element set are all pairwise adjacent (same vertices).")
print("  This gives α₁=3, α₂=0, I=7. BUT:")
print("  If the tournament has ANY 3-cycle, α₁ increases!")
print()
print("  So H=7 requires: 3 directed cycles with NO additional cycles.")
print("  With directed counting: a 5-vertex set with 3 Ham cycles + 0 three-cycles")
print("  would need the tournament to have no cyclic triples.")

# Check: at n=5, does any tournament have exactly 3 directed 5-cycles and 0 three-cycles?
found = 0
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Count 3-cycles
    dc3 = 0
    for v0, v1, v2 in combinations(range(n), 3):
        if A[v0][v1] and A[v1][v2] and A[v2][v0]:
            dc3 += 1
        elif A[v0][v2] and A[v2][v1] and A[v1][v0]:
            dc3 += 1

    if dc3 > 0:
        continue

    # Count 5-cycles
    dc5 = 0
    for perm in permutations(range(1, n)):
        cycle = (0,) + perm
        ok = True
        for i in range(n):
            if A[cycle[i]][cycle[(i+1) % n]] != 1:
                ok = False
                break
        if ok:
            dc5 += 1

    if dc5 > 0:
        found += 1
        H = fast_hp(A, n)
        print(f"  bits={bits}: dc3=0, dc5={dc5}, H={H}")

if found == 0:
    print("  NO tournament at n=5 has 0 three-cycles with any 5-cycles!")
    print("  (A tournament without 3-cycles is transitive → no odd cycles at all)")

# ═══════════════════════════════════════════════════════════════════
# Part 4: H=21 requires α₁=10 with correct counting
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 4: What (α₁_correct, α₂) gives H=21? ---")
print("  Still need 10 = α₁ + 2α₂ + 4α₃ + ...")
print("  Where α₁ now counts DIRECTED cycles.")
print()
print("  INSIGHT: Two directed cycles on the SAME vertex set are always")
print("  adjacent (share all vertices). So they can never be in the same")
print("  independent set. This means the multiplicity d(S) of a vertex set")
print("  S contributes a FACTOR of d(S) to the count of independent sets")
print("  that include exactly one cycle from S.")
print()
print("  In other words: I(Ω, x) = sum over collections of pairwise-disjoint")
print("  vertex sets, product of d(S_i) × x^k where k = number of sets.")
print()
print("  So H = I(Ω,2) = sum_{k≥0} 2^k × Σ_{pairwise-disjoint k-tuples} Π d(S_i)")
print("  = sum_{k≥0} 2^k × β_k")
print("  where β_k = Σ Π d(S_i) over all k-tuples of pairwise-disjoint vertex sets.")
print()
print("  For H=21: 1 + 2β₁ + 4β₂ + ... = 21")
print("  20 = 2β₁ + 4β₂ + ...")
print("  10 = β₁ + 2β₂ + ...")
print()
print("  β₁ = Σ d(S) over all vertex sets S with an odd cycle")
print("  β₂ = Σ_{disjoint pairs (S,T)} d(S)·d(T)")

# ═══════════════════════════════════════════════════════════════════
# Part 5: Compute β₁ distribution at n=5,6
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 5: β₁ distribution ---")

for n in [5, 6]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    ne = len(edges)

    beta1_to_H = defaultdict(set)
    H_to_beta1 = defaultdict(set)

    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        H = fast_hp(A, n)
        groups = count_directed_odd_cycles(A, n)
        beta1 = sum(groups.values())  # = α₁ in directed cycle sense

        beta1_to_H[beta1].add(H)
        H_to_beta1[H].add(beta1)

    print(f"\n  n={n}:")
    # What β₁ gives H near 21?
    for h in [19, 21, 23, 25]:
        if h in H_to_beta1:
            print(f"    H={h} ← β₁ ∈ {sorted(H_to_beta1[h])}")

    # Does β₁=10 exist?
    if 10 in beta1_to_H:
        print(f"    β₁=10 → H ∈ {sorted(beta1_to_H[10])}")
    else:
        print(f"    β₁=10: NOT ACHIEVED")

    # What β₁ values exist near 10?
    for b in [8, 9, 10, 11, 12]:
        if b in beta1_to_H:
            print(f"    β₁={b} → H ∈ {sorted(beta1_to_H[b])}")

# ═══════════════════════════════════════════════════════════════════
# Part 6: Full (β₁,β₂) at n=5,6 near H=21
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 6: (β₁,β₂) → H near 21 ---")

for n in [5, 6]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    ne = len(edges)

    b1b2_to_H = defaultdict(set)

    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        H = fast_hp(A, n)
        if H < 15 or H > 27:
            continue

        groups = count_directed_odd_cycles(A, n)
        vs_list = list(groups.items())  # (vertex_set, multiplicity)
        n_vs = len(vs_list)

        beta1 = sum(groups.values())

        # β₂ = sum over disjoint pairs of vertex sets of d(S)×d(T)
        beta2 = 0
        for i in range(n_vs):
            for j in range(i+1, n_vs):
                if not (vs_list[i][0] & vs_list[j][0]):  # disjoint
                    beta2 += vs_list[i][1] * vs_list[j][1]

        b1b2_to_H[(beta1, beta2)].add(H)

    print(f"\n  n={n}: (β₁,β₂) → H for 15 ≤ H ≤ 27:")
    for key in sorted(b1b2_to_H.keys()):
        h_set = sorted(b1b2_to_H[key])
        # Check if H=21 possible
        expected = 1 + 2*key[0] + 4*key[1]
        match = "✓" if expected in h_set else "✗"
        marker = " ← TARGET" if expected == 21 else ""
        print(f"    ({key[0]:2d},{key[1]:2d}): H ∈ {h_set}, expected 1+2·{key[0]}+4·{key[1]} = {expected} {match}{marker}")
