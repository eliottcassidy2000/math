#!/usr/bin/env python3
"""
Design Theory + Hereditary Maximizer Connection

The cyclic triples of Paley T_p form a 2-(p, 3, (p+1)/4) BIBD.
At p=7: two Fano planes, 7 disjoint pairs, alpha_2=7.
At p=11: 2-(11,3,3), 55 triples, 495 disjoint pairs.

Question: When we delete a vertex v from T_p, the surviving triples
form a "residual design." Is this residual design structure related
to why T_p - v maximizes H at n=p-1?

Also: What is the independence polynomial of the Paley maximizer?
How does it connect to the design parameters?

kind-pasteur-2026-03-06-S18f
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph

def indep_poly_coeffs(adj):
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def build_paley(p):
    """Build Paley tournament T_p for prime p = 3 mod 4."""
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

# ============================================================
# Paley T_7: Full cycle structure and I.P.
# ============================================================
print("=" * 70)
print("PALEY T_7: CYCLE STRUCTURE AND INDEPENDENCE POLYNOMIAL")
print("=" * 70)

T7 = build_paley(7)
h7 = hamiltonian_path_count(T7)
cycles7 = find_odd_cycles(T7)
cg7 = conflict_graph(cycles7)
ip7 = indep_poly_coeffs(cg7)

c3 = [c for c in cycles7 if len(c) == 3]
c5 = [c for c in cycles7 if len(c) == 5]
c7 = [c for c in cycles7 if len(c) == 7]

print(f"H(T_7) = {h7}")
print(f"Cycles: c3={len(c3)}, c5={len(c5)}, c7={len(c7)}, total={len(cycles7)}")
print(f"I.P. = {ip7}")
print(f"Check: H = {sum(c * (2**k) for k, c in enumerate(ip7))}")

# Identify the disjoint pairs (alpha_2)
print(f"\nalpha_2 = {ip7[2] if len(ip7) > 2 else 0}")
print(f"Disjoint 3-cycle pairs:")
disjoint_3pairs = []
for i in range(len(c3)):
    for j in range(i+1, len(c3)):
        if not cg7[cycles7.index(c3[i])][cycles7.index(c3[j])]:
            v1 = set(c3[i])
            v2 = set(c3[j])
            if v1.isdisjoint(v2):
                disjoint_3pairs.append((c3[i], c3[j]))
                print(f"  {c3[i]} + {c3[j]} (vertices: {sorted(v1)} + {sorted(v2)})")

print(f"Total disjoint 3-cycle pairs: {len(disjoint_3pairs)}")

# ============================================================
# Deletion analysis: T_7 - v
# ============================================================
print(f"\n{'='*70}")
print("T_7 - v: CYCLE STRUCTURE OF VERTEX DELETIONS")
print("=" * 70)

# By vertex-transitivity, just check v=0
v = 0
verts = [i for i in range(7) if i != v]
T6 = [[T7[verts[i]][verts[j]] for j in range(6)] for i in range(6)]
h6 = hamiltonian_path_count(T6)
cycles6 = find_odd_cycles(T6)
cg6 = conflict_graph(cycles6)
ip6 = indep_poly_coeffs(cg6)

c3_6 = [c for c in cycles6 if len(c) == 3]
c5_6 = [c for c in cycles6 if len(c) == 5]

print(f"T_7 - v=0: H={h6}")
print(f"Cycles: c3={len(c3_6)}, c5={len(c5_6)}, total={len(cycles6)}")
print(f"I.P. = {ip6}")
print(f"Check: H = {sum(c * (2**k) for k, c in enumerate(ip6))}")

# How do the surviving cycles relate to the original?
# A cycle in T_7 that doesn't include v survives in T_7-v.
surviving = [c for c in cycles7 if v not in c]
print(f"\nSurviving cycles from T_7 (not through v=0): {len(surviving)}")
for c in surviving:
    # Re-index to T_7-v indices
    reindexed = tuple(verts.index(x) for x in c)
    in_t6 = reindexed in cycles6 or tuple(reversed(reindexed)) in cycles6
    # Actually cycles are canonicalized differently, let me check by vertex sets
    found = False
    for c6 in cycles6:
        if set(c6) == set(reindexed):
            found = True
            break
    kind = f"{len(c)}-cycle"
    print(f"  {c} -> {reindexed} ({kind}, in T_7-v: {found})")

# New cycles in T_7-v that weren't in T_7 (after reindexing)
print(f"\nNew cycles in T_7-v (not from T_7):")
surviving_sets = [frozenset(verts.index(x) for x in c) for c in surviving]
for c6 in cycles6:
    if frozenset(c6) not in surviving_sets:
        print(f"  {c6} ({len(c6)}-cycle)")

# ============================================================
# The Residual Design
# ============================================================
print(f"\n{'='*70}")
print("RESIDUAL DESIGN: Triples of T_7 - v")
print("=" * 70)

# In T_7, the 14 cyclic triples form a 2-(7,3,2) design (two Fano planes).
# Deleting v=0: the 8 triples not containing 0 form the "residual design."
# But we only have 14 total, so 14 - 6 = 8 triples survive? No...
# Each triple contains 3 vertices. Through v=0: triples with 0 as a vertex.
# v=0 is in exactly lambda*(p-1)/(3-1) = 2*6/2 = 6 triples.
# So 14 - 6 = 8 triples survive.

triples_through_0 = [c for c in c3 if 0 in c]
triples_not_through_0 = [c for c in c3 if 0 not in c]
print(f"3-cycles through 0: {len(triples_through_0)}")
print(f"3-cycles not through 0: {len(triples_not_through_0)}")

for c in triples_not_through_0:
    reindexed = tuple(verts.index(x) for x in c)
    print(f"  {c} -> {reindexed} in T_7-v")

# These 8 triples on 6 points: is this a known design?
# 2-(7,3,2) with one point deleted: residual design is 2-(6,3,2)
# But wait: a residual design of a BIBD is obtained by deleting a point
# and keeping only blocks that contained it, removing the deleted point.
# That's different from what we want (keeping blocks NOT containing it).

# What we have: 8 triples on {1,2,3,4,5,6}, each triple = 3-cycle in T_7-v.
# T_7-v has these as its 3-cycles (c3=8 confirmed above).

# Structure: are these the 3-cycles of the n=6 maximizer?
# The n=6 maximizer (T_7-v) has H=45 = max.
# Its I.P. is either [1,14,4] (Type A) or [1,20,1] (Type B).
# From above: ip6 = [1,14,4] (for this specific deletion).
# So T_7-v is Type A (fewer cycles, more disjoint pairs).

print(f"\nT_7 - v=0: I.P. = {ip6}")
print(f"Type A (I.P.=[1,14,4], del=11) or Type B (I.P.=[1,20,1], del=13)?")
if ip6 == [1, 14, 4]:
    print("  TYPE A: fewer cycles, more disjoint pairs")
elif ip6 == [1, 20, 1]:
    print("  TYPE B: more cycles, fewer disjoint pairs")

# ============================================================
# Check the disjoint pair structure in T_7-v
# ============================================================
print(f"\n{'='*70}")
print("DISJOINT PAIRS IN T_7 - v")
print("=" * 70)

# Count disjoint pairs among ALL cycles of T_7-v
disjoint_all = 0
disjoint_33 = 0
disjoint_35 = 0
disjoint_55 = 0
for i in range(len(cycles6)):
    for j in range(i+1, len(cycles6)):
        if not cg6[i][j]:
            li, lj = len(cycles6[i]), len(cycles6[j])
            disjoint_all += 1
            if li == 3 and lj == 3:
                disjoint_33 += 1
            elif li == 5 and lj == 5:
                disjoint_55 += 1
            else:
                disjoint_35 += 1

print(f"Disjoint pairs by type:")
print(f"  3-3 pairs: {disjoint_33}")
print(f"  3-5 pairs: {disjoint_35}")
print(f"  5-5 pairs: {disjoint_55}")
print(f"  Total: {disjoint_all}")
print(f"  alpha_2 = {ip6[2] if len(ip6) > 2 else 0}")

# Which 3-cycle pairs are disjoint?
for i in range(len(c3_6)):
    for j in range(i+1, len(c3_6)):
        ci, cj = c3_6[i], c3_6[j]
        idx_i = cycles6.index(ci)
        idx_j = cycles6.index(cj)
        if not cg6[idx_i][idx_j]:
            print(f"  Disjoint 3-3: {ci} + {cj}")

# ============================================================
# The KEY question: Why is T_7-v the maximizer at n=6?
# Hypothesis: The Fano plane structure of T_7 forces T_7-v
# to have the maximum number of disjoint pairs among all
# n=6 tournaments in its score class.
# ============================================================
print(f"\n{'='*70}")
print("FANO STRUCTURE AND MAXIMALITY")
print("=" * 70)

# In T_7, the 7 disjoint 3-cycle pairs correspond to deleting one vertex.
# T_7 - v has 8 surviving 3-cycles. These 8 should contain 4 disjoint pairs
# (from the Fano structure: 7 pairs total, 6 through v, leaving... wait,
# the disjoint pairs aren't "through v" — they're pairs of vertex-disjoint cycles.

# Let me list ALL 7 disjoint pairs in T_7 and check which survive in T_7-v.
print(f"All disjoint 3-cycle pairs in T_7:")
for c1, c2 in disjoint_3pairs:
    survives = (0 not in c1) and (0 not in c2)
    print(f"  {c1} + {c2}: both survive in T_7-v? {survives}")

surviving_pairs = [(c1, c2) for c1, c2 in disjoint_3pairs
                   if (0 not in c1) and (0 not in c2)]
print(f"\nSurviving disjoint 3-cycle pairs: {len(surviving_pairs)}")

# But alpha_2 = 4 for T_7-v. So there are 4 independent pairs.
# 4 = alpha_2 counts independent SETS of size 2, which for cycles means
# vertex-disjoint pairs from ANY cycle type, not just 3-3.
# The 3-3 disjoint pairs that survive might be fewer than 4.

# Actually, alpha_2 from the I.P. [1,14,4] means there are 4 independent
# sets of size 2 in Omega. These could be any combination of cycle types.

print(f"\nalpha_2 = {ip6[2]}, breakdown: 3-3={disjoint_33}, 3-5={disjoint_35}, 5-5={disjoint_55}")

# ============================================================
# Compare: What would a RANDOM tournament with score (2,2,2,3,3,3)
# have for alpha_2?
# ============================================================
print(f"\n{'='*70}")
print("COMPARISON: alpha_2 DISTRIBUTION AT n=6 SCORE (2,2,2,3,3,3)")
print("=" * 70)

import random
random.seed(42)

n = 6
target = (2, 2, 2, 3, 3, 3)
alpha2_counts = {}
h_counts = {}
samples = 0
for bits in range(1 << 15):
    T = tournament_from_bits(n, bits)
    s = tuple(sorted(sum(T[i]) for i in range(n)))
    if s != target:
        continue
    samples += 1
    cycles = find_odd_cycles(T)
    if cycles:
        cg = conflict_graph(cycles)
        ip = indep_poly_coeffs(cg)
        a2 = ip[2] if len(ip) > 2 else 0
        h = sum(c * (2**k) for k, c in enumerate(ip))
    else:
        a2 = 0
        h = 1
    alpha2_counts[a2] = alpha2_counts.get(a2, 0) + 1
    h_counts[h] = h_counts.get(h, 0) + 1

print(f"Tournaments with score {target}: {samples}")
print(f"\nalpha_2 distribution:")
for a2 in sorted(alpha2_counts.keys()):
    print(f"  alpha_2={a2}: {alpha2_counts[a2]} ({100*alpha2_counts[a2]/samples:.1f}%)")
print(f"\nH distribution (top 5):")
for h in sorted(h_counts.keys(), reverse=True)[:5]:
    print(f"  H={h}: {h_counts[h]}")

# KEY: alpha_2=4 is the MAXIMUM alpha_2 in this score class.
# And it's achieved by T_7-v (the maximizer).
# alpha_2=1 is achieved by the Type B maximizer (I.P.=[1,20,1]).

print(f"\nMax alpha_2 = {max(alpha2_counts.keys())}")
print(f"Max H = {max(h_counts.keys())}")

print("\nDone.")
