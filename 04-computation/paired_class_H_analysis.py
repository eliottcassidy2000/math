#!/usr/bin/env python3
"""
Paired vs Self-Converse Isomorphism Classes: H-value Analysis
=============================================================
Hypothesis: Non-self-converse (NSC) tournaments come in "paired" iso classes
{C, C^op} where C != C^op. These paired classes should have LOWER H values
than self-converse (SC) classes.

Stronger hypothesis: Paired classes MINIMIZE H(T).

Also investigate:
- Distribution of H values across SC vs paired classes
- Whether blueself classes (at even n) achieve the MAXIMUM H
- Relationship between class pairing and independence polynomial

kind-pasteur-2026-03-06-S18
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, hamiltonian_path_count,
                             opposite_tournament, find_odd_cycles,
                             conflict_graph)
from itertools import permutations
from collections import defaultdict

def canonical_form(T):
    """Return a canonical form for tournament T (smallest adj matrix under vertex permutations)."""
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def is_self_converse(T):
    """Check if T is isomorphic to T^op."""
    Top = opposite_tournament(T)
    return canonical_form(T) == canonical_form(Top)

def score_sequence(T):
    """Sorted score sequence (descending)."""
    n = len(T)
    scores = sorted([sum(T[i]) for i in range(n)], reverse=True)
    return tuple(scores)

def independence_polynomial_coefficients(adj):
    """Compute independence polynomial coefficients from adjacency matrix."""
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

# ══════════════════════════════════════════════════════════════
# n=4: Exhaustive
# ══════════════════════════════════════════════════════════════
print("=" * 70)
print("PAIRED vs SELF-CONVERSE CLASS ANALYSIS")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    # Group tournaments by isomorphism class
    class_map = {}  # canonical_form -> list of (bits, T, H)
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        cf = canonical_form(T)
        if cf not in class_map:
            class_map[cf] = []
        class_map[cf].append((bits, h))

    # Classify each class as SC or paired
    sc_classes = []
    paired_classes = []  # store as pairs
    seen_paired = set()

    for cf in class_map:
        if cf in seen_paired:
            continue
        # Reconstruct tournament from first member
        bits0 = class_map[cf][0][0]
        T0 = tournament_from_bits(n, bits0)
        Top = opposite_tournament(T0)
        cf_op = canonical_form(Top)

        h_val = class_map[cf][0][1]  # H is same for all class members
        ss = score_sequence(T0)

        if cf == cf_op:
            sc_classes.append((cf, h_val, ss, len(class_map[cf])))
        else:
            paired_classes.append((cf, cf_op, h_val, ss, len(class_map[cf])))
            seen_paired.add(cf)
            seen_paired.add(cf_op)

    print(f"\nn={n}: {len(class_map)} iso classes, {len(sc_classes)} SC, {len(paired_classes)} paired")
    print(f"  SC classes (H, score_seq, class_size):")
    for cf, h, ss, sz in sorted(sc_classes, key=lambda x: x[1]):
        print(f"    H={h}, scores={ss}, size={sz}")

    print(f"  Paired classes (H, score_seq, class_size):")
    for cf, cf_op, h, ss, sz in sorted(paired_classes, key=lambda x: x[2]):
        print(f"    H={h}, scores={ss}, size={sz}")

    sc_h = [h for _, h, _, _ in sc_classes]
    paired_h = [h for _, _, h, _, _ in paired_classes]

    if sc_h:
        print(f"  SC H values: {sorted(sc_h)}, mean={sum(sc_h)/len(sc_h):.2f}")
    if paired_h:
        print(f"  Paired H values: {sorted(paired_h)}, mean={sum(paired_h)/len(paired_h):.2f}")

    # Check: does paired always have lower H?
    if sc_h and paired_h:
        min_sc = min(sc_h)
        max_paired = max(paired_h)
        print(f"  Min SC H = {min_sc}, Max Paired H = {max_paired}")
        print(f"  Paired ALWAYS < SC? {max_paired < min_sc}")
        print(f"  Paired MINIMIZE H? {max_paired <= min_sc}")

# ══════════════════════════════════════════════════════════════
# n=6: Exhaustive (2^15 = 32768 tournaments, 56 classes)
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("n=6 EXHAUSTIVE ANALYSIS")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
print(f"n={n}, m={m}, enumerating {1 << m} tournaments...")

class_map6 = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    cf = canonical_form(T)
    if cf not in class_map6:
        class_map6[cf] = {'h': h, 'count': 0, 'bits0': bits}
    class_map6[cf]['count'] += 1

print(f"Found {len(class_map6)} iso classes")

sc6 = []
paired6 = []
seen6 = set()

for cf, info in class_map6.items():
    if cf in seen6:
        continue
    T0 = tournament_from_bits(n, info['bits0'])
    Top = opposite_tournament(T0)
    cf_op = canonical_form(Top)
    ss = score_sequence(T0)

    # Compute independence polynomial
    cycles = find_odd_cycles(T0)
    if cycles:
        cg = conflict_graph(cycles)
        ip = independence_polynomial_coefficients(cg)
    else:
        ip = [1]

    if cf == cf_op:
        sc6.append((info['h'], ss, info['count'], ip))
    else:
        paired6.append((info['h'], ss, info['count'], ip))
        seen6.add(cf)
        seen6.add(cf_op)

sc6.sort(key=lambda x: x[0])
paired6.sort(key=lambda x: x[0])

print(f"\n{len(sc6)} SC classes, {len(paired6)} paired classes")

print("\nSC classes:")
for h, ss, cnt, ip in sc6:
    print(f"  H={h:3d}, scores={ss}, size={cnt:4d}, I.P.={ip}")

print("\nPaired classes:")
for h, ss, cnt, ip in paired6:
    print(f"  H={h:3d}, scores={ss}, size={cnt:4d}, I.P.={ip}")

sc_h6 = [h for h, _, _, _ in sc6]
paired_h6 = [h for h, _, _, _ in paired6]

print(f"\nSC H values: {sc_h6}")
print(f"Paired H values: {paired_h6}")
print(f"SC mean H = {sum(sc_h6)/len(sc_h6):.2f}")
print(f"Paired mean H = {sum(paired_h6)/len(paired_h6):.2f}")

min_sc6 = min(sc_h6)
max_paired6 = max(paired_h6)
print(f"Min SC H = {min_sc6}, Max Paired H = {max_paired6}")
print(f"Paired ALWAYS < min SC? {max_paired6 < min_sc6}")

# Check overlap
overlap = set(sc_h6) & set(paired_h6)
print(f"H values in BOTH SC and Paired: {sorted(overlap) if overlap else 'NONE'}")

# Who achieves the global max and min?
all_h = sc_h6 + paired_h6
print(f"\nGlobal max H = {max(all_h)}, achieved by SC? {max(all_h) in sc_h6}")
print(f"Global min H = {min(all_h)}, achieved by Paired? {min(all_h) in paired_h6}")

# ══════════════════════════════════════════════════════════════
# n=7: Can we do this? 456 classes. Let's try sampling.
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("n=7 ANALYSIS (random sampling for iso class identification)")
print("=" * 70)

import random
random.seed(42)

n = 7
# We can't enumerate all 2^21 = 2M tournaments for canonical form.
# Instead, use score sequence + H as a proxy for class identification,
# and just check SC vs NSC for each tournament.

sc_h7 = []
nsc_h7 = []
sc_scores7 = defaultdict(list)
nsc_scores7 = defaultdict(list)

# Sample 5000 random tournaments
for _ in range(5000):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    h = hamiltonian_path_count(T)
    sc = is_self_converse(T)
    ss = score_sequence(T)
    if sc:
        sc_h7.append(h)
        sc_scores7[ss].append(h)
    else:
        nsc_h7.append(h)
        nsc_scores7[ss].append(h)

print(f"Sampled 5000: {len(sc_h7)} SC, {len(nsc_h7)} NSC")
if sc_h7:
    print(f"SC: mean H={sum(sc_h7)/len(sc_h7):.2f}, min={min(sc_h7)}, max={max(sc_h7)}")
if nsc_h7:
    print(f"NSC: mean H={sum(nsc_h7)/len(nsc_h7):.2f}, min={min(nsc_h7)}, max={max(nsc_h7)}")

# H value ranges
print(f"\nSC H range: [{min(sc_h7)}, {max(sc_h7)}]")
print(f"NSC H range: [{min(nsc_h7)}, {max(nsc_h7)}]")

# Unique H values
sc_unique = sorted(set(sc_h7))
nsc_unique = sorted(set(nsc_h7))
print(f"SC unique H values: {sc_unique}")
print(f"NSC unique H values: {nsc_unique}")
both = set(sc_unique) & set(nsc_unique)
sc_only = set(sc_unique) - set(nsc_unique)
nsc_only = set(nsc_unique) - set(sc_unique)
print(f"H values in BOTH: {sorted(both)}")
print(f"H values SC-only: {sorted(sc_only)}")
print(f"H values NSC-only: {sorted(nsc_only)}")

# ══════════════════════════════════════════════════════════════
# Key test: At each n, is the global H-maximizer ALWAYS SC?
# And is the global H-minimizer ALWAYS paired/NSC?
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("H-MAXIMIZER / H-MINIMIZER ANALYSIS")
print("=" * 70)

for n in [3, 4, 5]:
    m_tiles = n * (n - 1) // 2
    max_h = 0
    min_h = float('inf')
    max_sc = None
    min_sc = None

    for bits in range(1 << m_tiles):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        sc = is_self_converse(T)
        if h > max_h:
            max_h = h
            max_sc = sc
        if h < min_h:
            min_h = h
            min_sc = sc

    print(f"n={n}: max H={max_h} (SC={max_sc}), min H={min_h} (SC={min_sc})")

# n=6
max_h6 = max(sc_h6 + paired_h6)
min_h6 = min(sc_h6 + paired_h6)
max_is_sc = max_h6 in sc_h6
min_is_paired = min_h6 in paired_h6
print(f"n=6: max H={max_h6} (SC={max_is_sc}), min H={min_h6} (Paired={min_is_paired}, SC={min_h6 in sc_h6})")

# ══════════════════════════════════════════════════════════════
# H distribution by score sequence regularity
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("SCORE SEQUENCE REGULARITY vs H (n=6)")
print("=" * 70)

def score_regularity(ss):
    """Variance of score sequence — lower = more regular."""
    n = len(ss)
    mean = sum(ss) / n
    return sum((s - mean)**2 for s in ss) / n

for h, ss, cnt, ip in sc6:
    reg = score_regularity(ss)
    print(f"  SC  H={h:3d}, scores={ss}, regularity={reg:.3f}, I.P.={ip}")

for h, ss, cnt, ip in paired6:
    reg = score_regularity(ss)
    print(f"  PAI H={h:3d}, scores={ss}, regularity={reg:.3f}, I.P.={ip}")

# ══════════════════════════════════════════════════════════════
# n=8 sampling: Does the pattern hold?
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("n=8 SAMPLING (1000 random)")
print("=" * 70)

n = 8
sc_h8 = []
nsc_h8 = []

for _ in range(1000):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    h = hamiltonian_path_count(T)
    sc = is_self_converse(T)
    if sc:
        sc_h8.append(h)
    else:
        nsc_h8.append(h)

print(f"Sampled 1000: {len(sc_h8)} SC, {len(nsc_h8)} NSC")
if sc_h8:
    print(f"SC: mean H={sum(sc_h8)/len(sc_h8):.2f}, min={min(sc_h8)}, max={max(sc_h8)}")
if nsc_h8:
    print(f"NSC: mean H={sum(nsc_h8)/len(nsc_h8):.2f}, min={min(nsc_h8)}, max={max(nsc_h8)}")

# Histograms
from collections import Counter
print(f"\nSC H histogram (top 10):")
for h, c in Counter(sc_h8).most_common(10):
    print(f"  H={h}: {c}")
print(f"NSC H histogram (top 10):")
for h, c in Counter(nsc_h8).most_common(10):
    print(f"  H={h}: {c}")

print("\nDone.")
