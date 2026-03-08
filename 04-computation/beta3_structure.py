#!/usr/bin/env python3
"""
beta3_structure.py - Understand beta3=1 tournament structure at n=6

The 320 beta3=1 tournaments split into:
- 80 with scores (1,1,1,4,4,4): c3=2
- 240 with scores (2,2,2,3,3,3): c3=8

What is special about these tournaments?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os
import numpy as np
from collections import Counter, defaultdict
from itertools import permutations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    c3 += 1
    return c3


# ============================================================
# Find all beta3=1 tournaments at n=6
# ============================================================
print("=" * 70)
print("BETA3=1 TOURNAMENT STRUCTURE: n=6")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)

b3_tours = []  # (bits, scores, c3, A)
for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    if betti[3] > 0:
        scores = tuple(sorted([sum(row) for row in A]))
        c3 = count_c3(A, n)
        b3_tours.append((bits, scores, c3, [row[:] for row in A]))

print(f"Total beta3>0: {len(b3_tours)}")

# ============================================================
# Score (1,1,1,4,4,4) analysis
# ============================================================
print(f"\n--- Score (1,1,1,4,4,4): 80 tournaments ---")

tours_114 = [(b, sc, c3, A) for b, sc, c3, A in b3_tours if sc == (1,1,1,4,4,4)]
print(f"Count: {len(tours_114)}")

# All have c3=2
c3_vals = set(c3 for _, _, c3, _ in tours_114)
print(f"c3 values: {c3_vals}")

# Show a representative
bits0, sc0, c30, A0 = tours_114[0]
print(f"\nExample: bits={bits0}")
for i in range(n):
    print(f"  {A0[i]}  (out-deg={sum(A0[i])})")

# Identify the structure: 3 "sinks" and 3 "sources"
low_deg = [i for i in range(n) if sum(A0[i]) <= 1]
high_deg = [i for i in range(n) if sum(A0[i]) >= 4]
print(f"\nLow-outdeg vertices: {low_deg}, High-outdeg: {high_deg}")

# What's the pattern between low and high?
print(f"Edges from high to low:")
for h in high_deg:
    for l in low_deg:
        print(f"  {h}->{l}: {A0[h][l]}")

print(f"Edges within low:")
for i in low_deg:
    for j in low_deg:
        if i != j:
            print(f"  {i}->{j}: {A0[i][j]}")

print(f"Edges within high:")
for i in high_deg:
    for j in high_deg:
        if i != j:
            print(f"  {i}->{j}: {A0[i][j]}")

# Check: is this the "blow-up" of a 3-cycle?
# i.e., partition {0,...,5} into 3 pairs, with edges between pairs following a 3-cycle pattern
# Low-deg vertices: each beats exactly 1 other vertex
# High-deg: each beats 4

# Actually with scores (1,1,1,4,4,4):
# Each low-deg vertex beats exactly 1 vertex
# Each high-deg vertex beats 4 vertices (all except one)

# Check complement: flip all edges
A0_comp = [[1 - A0[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
comp_scores = tuple(sorted([sum(row) for row in A0_comp]))
print(f"\nComplement scores: {comp_scores}")
# Should be (1,1,1,4,4,4) reversed = (1,1,1,4,4,4)
# Actually complement of (s_i) is (n-1-s_i) = (4,4,4,1,1,1) which sorts to (1,1,1,4,4,4)
# So these are SELF-COMPLEMENTARY scores!


# ============================================================
# Score (2,2,2,3,3,3) analysis
# ============================================================
print(f"\n\n--- Score (2,2,2,3,3,3): 240 tournaments ---")

tours_223 = [(b, sc, c3, A) for b, sc, c3, A in b3_tours if sc == (2,2,2,3,3,3)]
print(f"Count: {len(tours_223)}")

c3_vals = set(c3 for _, _, c3, _ in tours_223)
print(f"c3 values: {c3_vals}")

bits1, sc1, c31, A1 = tours_223[0]
print(f"\nExample: bits={bits1}")
for i in range(n):
    print(f"  {A1[i]}  (out-deg={sum(A1[i])})")

# These have c3=8 (maximum is C(6,3)/2 = 10 for regular, but (2,2,2,3,3,3) is almost-regular)
# Total (2,2,2,3,3,3) tournaments: 2640 (from our earlier data)
# Of these, 240 have beta3=1 (9.1%)

# What makes the beta3=1 ones special?
# Let's check the 2640 total and find what distinguishes the 240

all_223_bits = []
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    if scores == (2,2,2,3,3,3):
        all_223_bits.append(bits)

print(f"\nTotal (2,2,2,3,3,3) tournaments: {len(all_223_bits)}")
b3_bits = set(b for b, _, _, _ in tours_223)

# Compute distinguishing properties
print("\nDistinguishing properties:")
# Check c3, c5
c3_b3 = Counter()
c3_nob3 = Counter()
for bits in all_223_bits:
    A = build_adj(n, bits)
    c3 = count_c3(A, n)
    if bits in b3_bits:
        c3_b3[c3] += 1
    else:
        c3_nob3[c3] += 1

print(f"  c3 for beta3=1: {dict(sorted(c3_b3.items()))}")
print(f"  c3 for beta3=0: {dict(sorted(c3_nob3.items()))}")


# ============================================================
# Check isomorphism classes
# ============================================================
print(f"\n{'='*70}")
print("ISOMORPHISM CLASSES")
print("=" * 70)

def canonical_form(A, n):
    """Simple canonical form by trying all permutations."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if best is None or form < best:
            best = form
    return best

# For (1,1,1,4,4,4):
print("\n(1,1,1,4,4,4) isomorphism classes:")
iso_114 = defaultdict(int)
for bits, _, _, A in tours_114:
    cf = canonical_form(A, n)
    iso_114[cf] += 1
print(f"  Number of iso classes: {len(iso_114)}")
for cf, cnt in iso_114.items():
    print(f"  Size {cnt}")

# For (2,2,2,3,3,3) with beta3=1:
print("\n(2,2,2,3,3,3) beta3=1 iso classes:")
iso_223_b3 = defaultdict(int)
for bits, _, _, A in tours_223[:30]:  # Sample to avoid N! computation
    cf = canonical_form(A, n)
    iso_223_b3[cf] += 1
print(f"  Iso classes from first 30: {len(iso_223_b3)}")
for cf, cnt in sorted(iso_223_b3.items(), key=lambda x: -x[1]):
    print(f"  Size {cnt}")


# ============================================================
# Is (1,1,1,4,4,4) a "blow-up" of C3?
# ============================================================
print(f"\n{'='*70}")
print("STRUCTURAL ANALYSIS")
print("=" * 70)

# A tournament on {0,1,2,3,4,5} with scores (1,1,1,4,4,4)
# The 3 low-deg vertices form a specific pattern

# Check: is the induced subtournament on high-deg vertices a 3-cycle?
print("\n(1,1,1,4,4,4) induced subtournaments:")
for bits, _, _, A in tours_114[:5]:
    low = [i for i in range(n) if sum(A[i]) <= 1]
    high = [i for i in range(n) if sum(A[i]) >= 4]

    # Induced on high
    high_adj = [[A[high[i]][high[j]] for j in range(3)] for i in range(3)]
    high_c3 = 1 if (high_adj[0][1] and high_adj[1][2] and high_adj[2][0]) or \
                    (high_adj[1][0] and high_adj[2][1] and high_adj[0][2]) else 0

    # Induced on low
    low_adj = [[A[low[i]][low[j]] for j in range(3)] for i in range(3)]
    low_c3 = 1 if (low_adj[0][1] and low_adj[1][2] and low_adj[2][0]) or \
                   (low_adj[1][0] and low_adj[2][1] and low_adj[0][2]) else 0

    # All high->low edges
    h2l = sum(A[h][l] for h in high for l in low)
    l2h = sum(A[l][h] for l in low for h in high)

    print(f"  bits={bits}: high_c3={high_c3}, low_c3={low_c3}, h->l={h2l}, l->h={l2h}")


# ============================================================
# Connection to complement / reversal
# ============================================================
print(f"\n{'='*70}")
print("COMPLEMENT AND REVERSAL")
print("=" * 70)

# Check if these tournaments are self-complementary or self-converse
for bits, sc, c3, A in tours_114[:3]:
    # Complement (swap 0<->1 in adjacency, not diagonal)
    A_comp = [[1 - A[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    comp_scores = tuple(sorted([sum(row) for row in A_comp]))

    # Converse (transpose)
    A_conv = [[A[j][i] for j in range(n)] for i in range(n)]
    conv_scores = tuple(sorted([sum(row) for row in A_conv]))

    print(f"  bits={bits}: comp_scores={comp_scores}, conv_scores={conv_scores}")

    # Check if complement or converse is in our beta3=1 set
    comp_bits = None
    conv_bits = None
    for b2, _, _, A2 in b3_tours:
        if A2 == A_comp:
            comp_bits = b2
        if A2 == A_conv:
            conv_bits = b2

    print(f"    complement in beta3=1: {comp_bits is not None}")
    print(f"    converse in beta3=1: {conv_bits is not None}")


print("\n\nDone.")
