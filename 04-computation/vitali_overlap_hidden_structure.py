#!/usr/bin/env python3
"""
vitali_overlap_hidden_structure.py -- kind-pasteur-2026-03-13-S61

Deep exploration: overlap weights, Vitali sets, and hidden higher-dimensional
structure in small tournaments.

The key insight: the {2,1,0} overlap weight between cycle pairs carries
TERNARY information that refines the binary Omega adjacency. This creates
a hierarchy:
  - Level 0: Score sequence (measurable, determines c3)
  - Level 1: Binary Omega adjacency (determines H via OCF)
  - Level 2: Ternary overlap weights (refines Omega, extra info)
  - Level 3: Full tournament structure

The "Vitali set" analogy: just as Vitali's construction shows sets that
are not Lebesgue measurable (requiring the axiom of choice), the cycle
structure beyond the score sequence requires "non-constructive" information
about the tournament that cannot be deduced from the score alone.

This script:
1. Exhaustive n=5,6: overlap weight matrix eigenstructure
2. Overlap weight as a metric on cycle space
3. n=7 regular: ternary code determines H-class?
4. Hidden dimension counting: how many DOF are in the overlap structure?
5. Vitali quotient: what happens when we mod out by score class?
6. Connection to Fourier degree: overlap weight vs Fourier coefficient

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import combinations
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


def enumerate_directed_cycles(A, n, max_k=None):
    if max_k is None:
        max_k = n
    cycles = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            nc = count_directed_ham_cycles_on_subset(A, verts)
            for _ in range(nc):
                cycles.append((frozenset(subset), k))
    return cycles


def overlap_weight_vector(cycles):
    """Return the full overlap weight distribution as a sorted tuple.
    This is the ternary fingerprint of the tournament's cycle structure."""
    n = len(cycles)
    weights = []
    for i in range(n):
        for j in range(i+1, n):
            ov = len(cycles[i][0] & cycles[j][0])
            weights.append(ov)
    return tuple(sorted(weights))


def overlap_weight_histogram(cycles):
    """Return histogram of overlap weights."""
    n = len(cycles)
    hist = defaultdict(int)
    for i in range(n):
        for j in range(i+1, n):
            ov = len(cycles[i][0] & cycles[j][0])
            hist[ov] += 1
    return dict(sorted(hist.items()))


def overlap_weight_matrix_eigenvalues(cycles):
    """Compute eigenvalues of the overlap weight matrix W[i,j] = |V(C_i) cap V(C_j)|.
    Returns sorted eigenvalues (approximate, using power iteration for small n)."""
    nc = len(cycles)
    if nc == 0:
        return []

    # Build W matrix
    W = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        W[i][i] = len(cycles[i][0])  # self-overlap = size of cycle
        for j in range(i+1, nc):
            ov = len(cycles[i][0] & cycles[j][0])
            W[i][j] = ov
            W[j][i] = ov

    # For small nc, compute characteristic polynomial... too expensive.
    # Instead, compute trace and Frobenius norm as spectral invariants.
    trace = sum(W[i][i] for i in range(nc))
    frob_sq = sum(W[i][j]**2 for i in range(nc) for j in range(nc))

    return trace, frob_sq, nc


# ========================================================================
# ANALYSIS 1: EXHAUSTIVE n=5 — OVERLAP WEIGHT AS TOURNAMENT INVARIANT
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: OVERLAP WEIGHT STRUCTURE AT n=5 (EXHAUSTIVE)")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m

# Group tournaments by overlap weight vector
by_ov_vec = defaultdict(list)
by_score_and_H = defaultdict(list)

for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    cycles = enumerate_directed_cycles(A, n)
    nc = len(cycles)

    ov_hist = overlap_weight_histogram(cycles)
    # Use histogram as fingerprint (cheaper than full sorted vector)
    ov_key = tuple(sorted(ov_hist.items()))

    by_ov_vec[ov_key].append({'bits': bits, 'H': H, 'scores': scores, 'nc': nc})
    by_score_and_H[(scores, H)].append({'bits': bits, 'ov_key': ov_key, 'nc': nc})

print(f"\nTotal tournaments: {total}")
print(f"Distinct overlap weight histograms: {len(by_ov_vec)}")
print(f"Distinct (score, H) classes: {len(by_score_and_H)}")

# Key question: does the overlap histogram determine more than just (score, H)?
print(f"\nOverlap histogram classes:")
for ov_key in sorted(by_ov_vec.keys()):
    group = by_ov_vec[ov_key]
    H_set = set(d['H'] for d in group)
    sc_set = set(d['scores'] for d in group)
    print(f"  {dict(ov_key)}: count={len(group)}, H={sorted(H_set)}, scores={sorted(sc_set)}")

# Within each (score, H) class, how many distinct overlap histograms?
print(f"\nRefinement: within each (score, H) class, how many overlap histograms?")
for (scores, H) in sorted(by_score_and_H.keys()):
    group = by_score_and_H[(scores, H)]
    ov_keys = set(d['ov_key'] for d in group)
    if len(ov_keys) > 1:
        print(f"  Score {scores}, H={H}: {len(ov_keys)} distinct overlap histograms!")
        for ok in sorted(ov_keys):
            cnt = sum(1 for d in group if d['ov_key'] == ok)
            print(f"    {dict(ok)}: {cnt} tournaments")
    else:
        print(f"  Score {scores}, H={H}: 1 overlap histogram ({len(group)} tournaments)")


# ========================================================================
# ANALYSIS 2: THE {2,1,0} TERNARY CODE AT n=5
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: THE {2,1,0} TERNARY CODE AT n=5")
print("=" * 70)

# At n=5, 3-cycles have 3 vertices each, so overlap is 0, 1, or 2
# There's no overlap=3 (3-cycles on same vertex set counted separately)
# For directed cycles: at n=5 we also have 5-cycles (using all 5 vertices)
# 5-cycles overlap with everything by at least 3 vertices

# Focus on 3-cycle pairs only: the ternary structure is {0,1,2}
print("\n3-cycle pair overlap at n=5:")
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    H = count_ham_paths(A, n)
    cycles = enumerate_directed_cycles(A, n)

    # 3-cycle only
    c3 = [(fs, k) for fs, k in cycles if k == 3]
    c5 = [(fs, k) for fs, k in cycles if k == 5]

    if len(c3) < 2:
        continue

    # Ternary code: overlap weights between 3-cycle pairs
    nc3 = len(c3)
    ternary = []
    for i in range(nc3):
        for j in range(i+1, nc3):
            ov = len(c3[i][0] & c3[j][0])
            ternary.append(ov)

    ternary_code = tuple(sorted(ternary))

    if bits < 200 or H in [9, 11, 13, 15]:
        print(f"  bits={bits}: H={H}, c3={nc3}, c5={len(c5)}, "
              f"3-cycle ternary code={ternary_code}")
        break  # Just show first few

# Comprehensive: group by ternary code within score class (1,2,2,2,3)
print("\nGrouping score (1,2,2,2,3) by 3-cycle ternary code:")
ternary_groups = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    H = count_ham_paths(A, n)
    cycles = enumerate_directed_cycles(A, n)
    c3 = [(fs, k) for fs, k in cycles if k == 3]
    c5 = [(fs, k) for fs, k in cycles if k == 5]
    nc3 = len(c3)

    ternary = []
    for i in range(nc3):
        for j in range(i+1, nc3):
            ov = len(c3[i][0] & c3[j][0])
            ternary.append(ov)

    ternary_code = tuple(sorted(ternary))
    ternary_groups[ternary_code].append({'H': H, 'c3': nc3, 'c5': len(c5), 'bits': bits})

for tc in sorted(ternary_groups.keys()):
    group = ternary_groups[tc]
    H_set = set(d['H'] for d in group)
    c5_set = set(d['c5'] for d in group)
    print(f"  Ternary code {tc}: count={len(group)}, H={sorted(H_set)}, c5={sorted(c5_set)}")


# ========================================================================
# ANALYSIS 3: SPECTRAL INVARIANTS OF OVERLAP MATRIX
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: SPECTRAL INVARIANTS OF OVERLAP MATRIX AT n=5")
print("=" * 70)

# For each tournament, compute trace(W) and ||W||_F^2 of the overlap weight matrix
spectral_data = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    cycles = enumerate_directed_cycles(A, n)

    if len(cycles) == 0:
        spectral_data[(scores, H)].append((0, 0, 0))
        continue

    trace, frob_sq, nc = overlap_weight_matrix_eigenvalues(cycles)
    spectral_data[(scores, H)].append((trace, frob_sq, nc))

print(f"\nSpectral invariants by (score, H):")
for (scores, H) in sorted(spectral_data.keys()):
    group = spectral_data[(scores, H)]
    traces = set(t for t, f, nc in group)
    frobs = set(f for t, f, nc in group)
    ncs = set(nc for t, f, nc in group)
    print(f"  Score {scores}, H={H}: nc={sorted(ncs)}, "
          f"trace(W)={sorted(traces)}, ||W||_F^2={sorted(frobs)}")


# ========================================================================
# ANALYSIS 4: n=6 — OVERLAP WEIGHT AND ALPHA_2
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: OVERLAP WEIGHT AND ALPHA_2 AT n=6")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m

# For n=6, alpha_2 matters. Does overlap weight histogram determine alpha_2?
alpha2_by_ov = defaultdict(list)
count = 0
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    cycles = enumerate_directed_cycles(A, n)

    nc = len(cycles)
    c3 = [(fs, k) for fs, k in cycles if k == 3]
    c5 = [(fs, k) for fs, k in cycles if k == 5]

    # Compute alpha_2: number of vertex-disjoint cycle pairs
    disj_pairs = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if not (cycles[i][0] & cycles[j][0]):
                disj_pairs += 1

    # Overlap histogram (3-cycle pairs only for speed)
    nc3 = len(c3)
    ov_hist_3 = defaultdict(int)
    for i in range(nc3):
        for j in range(i+1, nc3):
            ov = len(c3[i][0] & c3[j][0])
            ov_hist_3[ov] += 1

    ov_key = (len(c3), len(c5), tuple(sorted(ov_hist_3.items())))
    alpha2_by_ov[ov_key].append({
        'H': H, 'scores': scores, 'alpha2': disj_pairs,
        'c3': len(c3), 'c5': len(c5), 'nc': nc
    })
    count += 1

print(f"Processed {count} tournaments")
print(f"Distinct overlap profiles: {len(alpha2_by_ov)}")

# Check: does overlap profile determine alpha_2?
ambig_count = 0
for ov_key in sorted(alpha2_by_ov.keys()):
    group = alpha2_by_ov[ov_key]
    a2_set = set(d['alpha2'] for d in group)
    if len(a2_set) > 1:
        ambig_count += 1
        H_set = set(d['H'] for d in group)
        print(f"  AMBIGUOUS: c3={ov_key[0]}, c5={ov_key[1]}, "
              f"ov_hist={dict(ov_key[2])}, alpha2={sorted(a2_set)}, "
              f"H={sorted(H_set)}, count={len(group)}")

if ambig_count == 0:
    print("  3-cycle overlap histogram DETERMINES alpha_2!")
else:
    print(f"  {ambig_count} ambiguous overlap profiles for alpha_2")


# ========================================================================
# ANALYSIS 5: HIDDEN DIMENSION — DOF IN CYCLE STRUCTURE BEYOND SCORES
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: HIDDEN DIMENSION COUNT")
print("=" * 70)

# At n=5: score has 3 distinct classes. H has 7 values.
# How many bits of information are in:
# (a) Score sequence
# (b) H value
# (c) Overlap weight histogram
# (d) Full tournament (up to isomorphism)

# Information = log2(number of equivalence classes)
n = 5
m = n * (n - 1) // 2
total = 1 << m

scores_n5 = set()
H_n5 = set()
ov_hists_n5 = set()
full_orbits_n5 = set()

for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    cycles = enumerate_directed_cycles(A, n)
    ov_hist = tuple(sorted(overlap_weight_histogram(cycles).items()))

    scores_n5.add(scores)
    H_n5.add(H)
    ov_hists_n5.add(ov_hist)

print(f"n=5:")
print(f"  Total tournaments: {total}")
print(f"  Score classes: {len(scores_n5)} ({math.log2(len(scores_n5)):.2f} bits)")
print(f"  H classes: {len(H_n5)} ({math.log2(len(H_n5)):.2f} bits)")
print(f"  Overlap histogram classes: {len(ov_hists_n5)} ({math.log2(len(ov_hists_n5)):.2f} bits)")
print(f"  Full tournament space: {total} ({m:.2f} bits)")

# Also n=6
n = 6
m = n * (n - 1) // 2
total = 1 << m

scores_n6 = set()
H_n6 = set()
ov_hists_n6 = set()

for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    cycles = enumerate_directed_cycles(A, n)
    ov_hist = tuple(sorted(overlap_weight_histogram(cycles).items()))

    scores_n6.add(scores)
    H_n6.add(H)
    ov_hists_n6.add(ov_hist)

print(f"\nn=6:")
print(f"  Total tournaments: {total}")
print(f"  Score classes: {len(scores_n6)} ({math.log2(len(scores_n6)):.2f} bits)")
print(f"  H classes: {len(H_n6)} ({math.log2(len(H_n6)):.2f} bits)")
print(f"  Overlap histogram classes: {len(ov_hists_n6)} ({math.log2(len(ov_hists_n6)):.2f} bits)")
print(f"  Full tournament space: {total} ({m:.2f} bits)")


# ========================================================================
# ANALYSIS 6: VITALI QUOTIENT — NON-MEASURABLE CYCLE STRUCTURE
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: VITALI QUOTIENT — STRUCTURE BEYOND SCORES")
print("=" * 70)

# The "Vitali set" analogy: within a score class, different tournaments have
# different cycle structures. The score class is like a measurable set.
# The finer cycle structure is like a Vitali partition — it exists but
# requires additional information (like the axiom of choice for tournaments).

# At n=5: focus on score class (1,2,2,2,3)
print("\nScore class (1,2,2,2,3) at n=5:")
n = 5
m = n * (n - 1) // 2
total = 1 << m

class_data = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    H = count_ham_paths(A, n)
    cycles = enumerate_directed_cycles(A, n)
    c3 = sum(1 for _, k in cycles if k == 3)
    c5 = sum(1 for _, k in cycles if k == 5)

    # Compute full cycle vertex set incidence
    c3_sets = [fs for fs, k in cycles if k == 3]
    c5_sets = [fs for fs, k in cycles if k == 5]

    # The "hidden structure": given c3 (fixed by score), c5 varies
    # This variation is the "non-measurable" part
    class_data.append({
        'bits': bits, 'H': H, 'c3': c3, 'c5': c5,
        'c3_sets': c3_sets, 'c5_sets': c5_sets
    })

# Group by c5 (the non-score-determined quantity)
by_c5 = defaultdict(list)
for d in class_data:
    by_c5[d['c5']].append(d)

print(f"  Total in class: {len(class_data)}")
print(f"  c3 is CONSTANT = {class_data[0]['c3']} (score-determined)")
print(f"  c5 varies:")
for c5 in sorted(by_c5.keys()):
    group = by_c5[c5]
    H_set = set(d['H'] for d in group)
    print(f"    c5={c5}: count={len(group)}, H={sorted(H_set)}")

# Key question: what property of the tournament determines c5?
# Check if any simple structural feature separates c5 values
print(f"\n  What determines c5? Let's check vertex participation:")
for c5 in sorted(by_c5.keys()):
    group = by_c5[c5]
    d = group[0]
    A = binary_to_tournament(d['bits'], n)

    # Which vertex has out-degree 1 and which has out-degree 3?
    out_degrees = [sum(A[v]) for v in range(n)]
    sink = out_degrees.index(min(out_degrees))
    source = out_degrees.index(max(out_degrees))

    # Check if 5-cycles pass through both sink and source
    c5_through_sink = sum(1 for fs in d['c5_sets'] if sink in fs)
    c5_through_source = sum(1 for fs in d['c5_sets'] if source in fs)

    # 3-cycle sets: which contain the source/sink
    c3_through_sink = sum(1 for fs in d['c3_sets'] if sink in fs)
    c3_through_source = sum(1 for fs in d['c3_sets'] if source in fs)

    print(f"    c5={c5}: out_degrees={out_degrees}, "
          f"c3_through_sink={c3_through_sink}, c3_through_source={c3_through_source}, "
          f"c5_through_sink={c5_through_sink}, c5_through_source={c5_through_source}")


# ========================================================================
# ANALYSIS 7: OVERLAP WEIGHT AS A METRIC ON CYCLE SPACE
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: OVERLAP WEIGHT METRIC PROPERTIES")
print("=" * 70)

# The overlap weight W(C_i, C_j) = |V(C_i) cap V(C_j)| is NOT a metric
# (it increases with overlap, not decreases). But d(C_i, C_j) = k - W(C_i, C_j)
# where k = max(|C_i|, |C_j|) gives a pseudo-metric.

# For 3-cycle pairs: d = 3 - W in {1, 2, 3} (3-overlap=same cycle doesn't happen
# because we count directed cycles, not vertex sets)
# Actually W in {0, 1, 2} so d in {1, 2, 3}

# Check if the overlap-distance satisfies triangle inequality for 3-cycles
n = 5
m = n * (n - 1) // 2
total = 1 << m

triangle_violations = 0
triangle_total = 0

for bits in range(min(200, total)):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    cycles = enumerate_directed_cycles(A, n)
    c3 = [(fs, k) for fs, k in cycles if k == 3]
    nc3 = len(c3)

    if nc3 < 3:
        continue

    # Check triangle inequality for d = 3 - W
    for i in range(nc3):
        for j in range(i+1, nc3):
            for k_idx in range(j+1, nc3):
                wij = len(c3[i][0] & c3[j][0])
                wik = len(c3[i][0] & c3[k_idx][0])
                wjk = len(c3[j][0] & c3[k_idx][0])
                dij = 3 - wij
                dik = 3 - wik
                djk = 3 - wjk
                triangle_total += 1
                if dij > dik + djk or dik > dij + djk or djk > dij + dik:
                    triangle_violations += 1
    break  # Just check first tournament in class

print(f"Triangle inequality check (d = 3 - overlap weight) on 3-cycles:")
print(f"  Total triples: {triangle_total}")
print(f"  Violations: {triangle_violations}")
if triangle_violations == 0:
    print(f"  d = 3 - W is a VALID metric on 3-cycles!")
else:
    print(f"  d = 3 - W violates triangle inequality")


# ========================================================================
# ANALYSIS 8: n=7 REGULAR — TERNARY STRUCTURE DETERMINES H-CLASS
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 8: n=7 REGULAR — TERNARY OVERLAP STRUCTURE")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m

print(f"Scanning all {total} tournaments for regular ones...")

regular_data = []
checked = 0
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = [sum(A[v]) for v in range(n)]
    if any(s != 3 for s in scores):
        continue

    H = count_ham_paths(A, n)

    # Get 3-cycle structure
    c3 = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3.append(frozenset([a, b, c_v]))
    c3 = list(set(c3))  # vertex sets

    # 3-cycle overlap histogram
    nc3 = len(c3)
    ov_hist = defaultdict(int)
    for i in range(nc3):
        for j in range(i+1, nc3):
            ov = len(c3[i] & c3[j])
            ov_hist[ov] += 1

    ov_key = tuple(sorted(ov_hist.items()))

    # Count disjoint pairs
    disj = ov_hist.get(0, 0)

    regular_data.append({
        'bits': bits, 'H': H, 'c3': nc3, 'ov_key': ov_key,
        'disj33': disj, 'ov_hist': dict(ov_hist)
    })
    checked += 1

print(f"Found {checked} regular tournaments")

# Group by H
by_H = defaultdict(list)
for d in regular_data:
    by_H[d['H']].append(d)

print(f"\nOverlap structure by H-class:")
for H in sorted(by_H.keys()):
    group = by_H[H]
    ov_keys = set(d['ov_key'] for d in group)
    disj_set = set(d['disj33'] for d in group)
    print(f"  H={H}: count={len(group)}, disj_33={sorted(disj_set)}")
    for ok in sorted(ov_keys):
        cnt = sum(1 for d in group if d['ov_key'] == ok)
        print(f"    overlap hist {dict(ok)}: {cnt} tournaments")

# Key: does the overlap histogram UNIQUELY determine the H-class?
print(f"\nDoes 3-cycle overlap histogram determine H-class?")
by_ov_key = defaultdict(list)
for d in regular_data:
    by_ov_key[d['ov_key']].append(d)

for ok in sorted(by_ov_key.keys()):
    group = by_ov_key[ok]
    H_set = set(d['H'] for d in group)
    if len(H_set) > 1:
        print(f"  AMBIGUOUS: {dict(ok)} -> H = {sorted(H_set)}")
    else:
        print(f"  {dict(ok)} -> H = {sorted(H_set)[0]} ({len(group)} tournaments)")


# ========================================================================
# ANALYSIS 9: THE 5-CYCLE STRUCTURE IN REGULAR n=7
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 9: 5-CYCLE OVERLAP STRUCTURE IN REGULAR n=7")
print("=" * 70)

# We know c5_dir + 2*disj_33 = 56 for regular n=7.
# But what about the overlap structure of 5-cycles themselves?

# Sample one tournament from each H-class
for H_target in sorted(by_H.keys()):
    d = by_H[H_target][0]
    bits = d['bits']
    A = binary_to_tournament(bits, n)

    # Get 5-cycles (directed)
    c5 = []
    for subset in combinations(range(n), 5):
        verts = list(subset)
        nc5 = count_directed_ham_cycles_on_subset(A, verts)
        for _ in range(nc5):
            c5.append(frozenset(subset))

    print(f"\n  H={H_target}: {len(c5)} directed 5-cycles")

    # 5-cycle vertex set distribution
    c5_vsets = list(set(c5))
    nc5_vsets = len(c5_vsets)
    directed_per_vset = defaultdict(int)
    for fs in c5:
        directed_per_vset[fs] += 1

    dpv_dist = defaultdict(int)
    for fs, cnt in directed_per_vset.items():
        dpv_dist[cnt] += 1
    print(f"    Vertex sets: {nc5_vsets}")
    print(f"    Directed cycles per vertex set: {dict(sorted(dpv_dist.items()))}")

    # Overlap between 5-cycle vertex sets
    c5_ov_hist = defaultdict(int)
    for i in range(nc5_vsets):
        for j in range(i+1, nc5_vsets):
            ov = len(c5_vsets[i] & c5_vsets[j])
            c5_ov_hist[ov] += 1
    print(f"    5-cycle vertex set overlap: {dict(sorted(c5_ov_hist.items()))}")

    # 3-5 cross-overlap: how do 3-cycles and 5-cycles interact?
    c3_vsets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_vsets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_vsets.append(frozenset([a, b, c_v]))
    c3_vsets = list(set(c3_vsets))

    # Each 3-cycle must be a subset of some 5-cycle vertex set (or not)
    containment = 0
    total_35 = 0
    for c3_vs in c3_vsets:
        for c5_vs in c5_vsets:
            total_35 += 1
            if c3_vs.issubset(c5_vs):
                containment += 1

    print(f"    3-in-5 containment: {containment}/{total_35} "
          f"({100*containment/total_35:.1f}%)")


# ========================================================================
# ANALYSIS 10: WEIGHTED INDEPENDENCE POLYNOMIAL
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 10: OVERLAP-WEIGHTED INDEPENDENCE POLYNOMIAL")
print("=" * 70)

# Instead of I(Omega, x), consider the WEIGHTED version where each
# independent set S of size k contributes x^k * prod_{distinct pairs (i,j) in S} f(W(i,j))
# where f(0) = 1 (disjoint pair), f(w) = 0 for w > 0 (overlapping pair excluded)
# This is just the standard independence polynomial.
#
# But we can define a SOFT version: f(w) = q^w for some parameter q in [0,1]
# At q=0: standard I.P. (only truly disjoint sets contribute)
# At q=1: all pairs allowed (I = (1+x)^n trivially)
# The transition from q=0 to q=1 reveals the "softness" of the conflict structure.

print("Soft independence: I_q(Omega, x) at n=5, score (1,2,2,2,3)")
n = 5
m = n * (n - 1) // 2
total = 1 << m

for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    H = count_ham_paths(A, n)
    cycles = enumerate_directed_cycles(A, n)
    nc = len(cycles)

    # For each q in {0, 0.25, 0.5, 0.75, 1}, compute I_q(Omega, 2)
    # I_q = sum over subsets S of cycles: prod_{pairs in S} q^{W(i,j)} * 2^|S|
    # For independent sets, W(i,j) = 0 so q^0 = 1.
    # For overlapping pairs, q^w with w > 0.

    # Precompute overlap weights
    W = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            ov = len(cycles[i][0] & cycles[j][0])
            W[i][j] = ov
            W[j][i] = ov

    qs = [0.0, 0.25, 0.5, 0.75, 1.0]
    results = []
    for q in qs:
        # Enumerate all 2^nc subsets
        I_q = 0.0
        for mask in range(1 << nc):
            # Get elements in this subset
            elems = [i for i in range(nc) if mask & (1 << i)]
            k = len(elems)

            # Product of q^{W(i,j)} over all pairs
            prod_val = 1.0
            for ii in range(len(elems)):
                for jj in range(ii+1, len(elems)):
                    w = W[elems[ii]][elems[jj]]
                    if q == 0.0 and w > 0:
                        prod_val = 0.0
                        break
                    elif q > 0.0:
                        prod_val *= q ** w
                if prod_val == 0.0:
                    break

            I_q += prod_val * (2 ** k)
        results.append(I_q)

    print(f"  bits={bits}, H={H}: I_q(Omega,2) = "
          + ", ".join(f"q={q}:{r:.1f}" for q, r in zip(qs, results)))
    if H == 15:  # Show a few examples
        break

# Show for each H value
print(f"\n  Per H-value (one representative each):")
seen_H = set()
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    H = count_ham_paths(A, n)
    if H in seen_H:
        continue
    seen_H.add(H)

    cycles = enumerate_directed_cycles(A, n)
    nc = len(cycles)

    W = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            ov = len(cycles[i][0] & cycles[j][0])
            W[i][j] = ov
            W[j][i] = ov

    # Compute for q = 0, 0.5, 1
    for q in [0.0, 0.5, 1.0]:
        I_q = 0.0
        for mask in range(1 << nc):
            elems = [i for i in range(nc) if mask & (1 << i)]
            k = len(elems)
            prod_val = 1.0
            for ii in range(len(elems)):
                for jj in range(ii+1, len(elems)):
                    w = W[elems[ii]][elems[jj]]
                    if q == 0.0 and w > 0:
                        prod_val = 0.0
                        break
                    elif q > 0:
                        prod_val *= q ** w
                if prod_val == 0:
                    break
            I_q += prod_val * (2 ** k)

        if q == 0.0:
            print(f"    H={H}: I_0 = {I_q:.0f} (= H, standard OCF)", end="")
        elif q == 0.5:
            print(f", I_0.5 = {I_q:.1f}", end="")
        else:
            print(f", I_1 = {I_q:.0f} (= 2^nc={2**nc})")


print(f"\n\n{'='*70}")
print("DONE.")
print("=" * 70)
