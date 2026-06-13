#!/usr/bin/env python3
"""
beta1_discriminant_S72.py — opus-2026-03-15-S72

THE MISSING INVARIANT: What distinguishes β₁=0 from β₁=1?

We know: (score_seq, t3) does NOT determine β₁ at n≥5.
The ambiguity is always β₁ ∈ {0, 1}.

CANDIDATES for the discriminant:
  1. Score sequence of 3-cycle participation (c3v = per-vertex 3-cycle count)
  2. Number of strongly-connected 4-subtournaments (SC4)
  3. Number of 4-cycles
  4. Second eigenvalue of adjacency matrix
  5. Total transitive triple count (= C(n,3) - t3)
  6. Per-vertex dominance pattern (in-out-degree sequences)

If we find it → prove it → INSTANT β₁ computation → skip the biggest rank computation

ALSO: Exploit β₂=0 (PROVED) for speedup:
  β₁ = C(n,2) - (n-1) - dim(Ω₂) + rank(∂₃)  [seesaw formula]
  So computing rank(∂₃) gives β₁ for FREE once we know dim(Ω₂).
  But we still need dim(Ω₂), which requires the constraint matrix.

  HOWEVER: if dim(Ω₂) has a closed form... we'd skip EVERYTHING.
"""

import numpy as np
from math import comb
from collections import Counter, defaultdict
import sys
sys.stdout.reconfigure(line_buffering=True)

PRIME = 2**31 - 1

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def score_seq(A, n):
    return tuple(sorted(int(np.sum(A[i])) for i in range(n)))

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i,j] and A[j,k] and A[k,i]) or \
                   (A[i,k] and A[k,j] and A[j,i]):
                    c3 += 1
    return c3

def per_vertex_3cycles(A, n):
    """Count 3-cycles through each vertex. Returns sorted tuple."""
    c3v = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i,j] and A[j,k] and A[k,i]) or \
                   (A[i,k] and A[k,j] and A[j,i]):
                    c3v[i] += 1
                    c3v[j] += 1
                    c3v[k] += 1
    return tuple(sorted(c3v))

def count_4cycles(A, n):
    """Count directed 4-cycles."""
    c4 = 0
    for i in range(n):
        for j in range(n):
            if j == i or A[i,j] == 0: continue
            for k in range(n):
                if k in (i,j) or A[j,k] == 0: continue
                for l in range(n):
                    if l in (i,j,k) or A[k,l] == 0 or A[l,i] == 0: continue
                    c4 += 1
    return c4 // 4

def count_sc4(A, n):
    """Count 4-vertex strongly-connected subtournaments."""
    sc4 = 0
    for verts in __import__('itertools').combinations(range(n), 4):
        # A 4-tournament is SC iff it has at least one 3-cycle
        # (equivalently: NOT transitive)
        sub = np.zeros((4, 4), dtype=np.int8)
        for i, vi in enumerate(verts):
            for j, vj in enumerate(verts):
                sub[i, j] = A[vi, vj]
        # Check: is sub transitive? (score seq = 0,1,2,3)
        scores = sorted(int(np.sum(sub[i])) for i in range(4))
        if scores != [0, 1, 2, 3]:
            sc4 += 1
    return sc4

def dominance_graph_edges(A, n):
    """Count edges in the 'dominance graph': i dominates j if i beats
    all of j's out-neighbors that j doesn't beat."""
    # Simpler: just use eigenvalues
    eigs = np.sort(np.real(np.linalg.eigvals(A.astype(float))))[::-1]
    return round(eigs[1], 6)

def compute_beta1(A, n):
    """Fast β₁ computation."""
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                edges.append((i, j))
                edge_idx[(i, j)] = len(edges) - 1

    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue
            for k in range(n):
                if k in (i, j) or A[j][k] == 0: continue
                if A[i][k] == 1:
                    triples.append((i, j, k))

    if not triples:
        return comb(n, 2) - (n - 1)

    bd = np.zeros((len(edges), len(triples)), dtype=np.int64)
    for t, (i, j, k) in enumerate(triples):
        if (j, k) in edge_idx:
            bd[edge_idx[(j, k)], t] = (bd[edge_idx[(j, k)], t] + 1) % PRIME
        if (i, k) in edge_idx:
            bd[edge_idx[(i, k)], t] = (bd[edge_idx[(i, k)], t] - 1) % PRIME
        if (i, j) in edge_idx:
            bd[edge_idx[(i, j)], t] = (bd[edge_idx[(i, j)], t] + 1) % PRIME

    M = bd % PRIME
    nrows, ncols = M.shape
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0: continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), PRIME - 2, PRIME)
        M[rank] = M[rank] * inv % PRIME
        factors = M[:, col].copy()
        factors[rank] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[rank])) % PRIME
        rank += 1

    return max(0, comb(n, 2) - (n - 1) - rank)

print("=" * 72)
print("β₁ DISCRIMINANT SEARCH")
print("opus-2026-03-15-S72")
print("=" * 72)

# ============================================================
# STEP 1: Collect all tournament data at n=5
# ============================================================

n = 5
m = n * (n - 1) // 2
total = 2**m

print(f"\nCollecting invariants for all {total} tournaments at n={n}...")

data = []
for bits in range(total):
    A = bits_to_adj(bits, n)
    b1 = compute_beta1(A, n)
    ss = score_seq(A, n)
    t3 = count_3cycles(A, n)
    c3v = per_vertex_3cycles(A, n)
    c4 = count_4cycles(A, n)
    sc4 = count_sc4(A, n)
    eig2 = dominance_graph_edges(A, n)
    data.append({
        'bits': bits, 'A': A, 'b1': b1, 'ss': ss, 't3': t3,
        'c3v': c3v, 'c4': c4, 'sc4': sc4, 'eig2': eig2
    })

print(f"Done. β₁ distribution: {Counter(d['b1'] for d in data)}")

# ============================================================
# STEP 2: Find which invariant resolves the ambiguity
# ============================================================

print(f"\n" + "=" * 72)
print("STEP 2: Which invariant resolves β₁ ∈ {{0,1}} ambiguity?")
print("=" * 72)

# Identify ambiguous groups: same (ss, t3) but different β₁
ambiguous_groups = defaultdict(list)
for d in data:
    ambiguous_groups[(d['ss'], d['t3'])].append(d)

amb_keys = [k for k, v in ambiguous_groups.items()
            if len(set(d['b1'] for d in v)) > 1]

print(f"\n  {len(amb_keys)} ambiguous (score_seq, t3) groups:")

for key in amb_keys:
    ss, t3 = key
    group = ambiguous_groups[key]
    b0 = [d for d in group if d['b1'] == 0]
    b1 = [d for d in group if d['b1'] == 1]
    print(f"\n  score_seq={ss}, t3={t3}: {len(b0)} with β₁=0, {len(b1)} with β₁=1")

    # Test each candidate invariant
    candidates = {
        'c3v (per-vertex 3-cycles)': lambda d: d['c3v'],
        'c4 (4-cycle count)': lambda d: d['c4'],
        'sc4 (SC 4-subtournaments)': lambda d: d['sc4'],
        'eig2 (2nd eigenvalue)': lambda d: d['eig2'],
    }

    for cand_name, cand_fn in candidates.items():
        vals_b0 = Counter(cand_fn(d) for d in b0)
        vals_b1 = Counter(cand_fn(d) for d in b1)
        overlap = set(vals_b0.keys()) & set(vals_b1.keys())

        if not overlap:
            print(f"    ★ {cand_name} RESOLVES the ambiguity!")
            print(f"      β₁=0: {dict(vals_b0)}")
            print(f"      β₁=1: {dict(vals_b1)}")
        else:
            print(f"    ✗ {cand_name}: overlap at {overlap}")
            print(f"      β₁=0: {dict(vals_b0)}")
            print(f"      β₁=1: {dict(vals_b1)}")

# ============================================================
# STEP 3: Test the winning discriminant at n=6
# ============================================================

print(f"\n" + "=" * 72)
print("STEP 3: Verify winning discriminant at n=6")
print("=" * 72)

n6 = 6
m6 = n6 * (n6 - 1) // 2
total6 = 2**m6

print(f"\nCollecting data for all {total6} tournaments at n={n6}...")
import time
t0 = time.time()

# Focus on just the invariants that worked
data6 = []
count = 0
for bits in range(total6):
    A = bits_to_adj(bits, n6)
    b1 = compute_beta1(A, n6)
    ss = score_seq(A, n6)
    t3 = count_3cycles(A, n6)
    c3v = per_vertex_3cycles(A, n6)
    sc4 = count_sc4(A, n6)
    data6.append({
        'b1': b1, 'ss': ss, 't3': t3, 'c3v': c3v, 'sc4': sc4
    })
    count += 1
    if count % 5000 == 0:
        print(f"  {count}/{total6} ({time.time()-t0:.1f}s)")

print(f"Done in {time.time()-t0:.1f}s. β₁ distribution: {Counter(d['b1'] for d in data6)}")

# Check if (ss, t3, winning_invariant) determines β₁
# First check c3v
groups_c3v = defaultdict(set)
for d in data6:
    groups_c3v[(d['ss'], d['t3'], d['c3v'])].add(d['b1'])
violations_c3v = sum(1 for v in groups_c3v.values() if len(v) > 1)
print(f"\n  (ss, t3, c3v) determines β₁ at n=6: {'YES' if violations_c3v==0 else f'NO ({violations_c3v} violations)'}")

# Check sc4
groups_sc4 = defaultdict(set)
for d in data6:
    groups_sc4[(d['ss'], d['t3'], d['sc4'])].add(d['b1'])
violations_sc4 = sum(1 for v in groups_sc4.values() if len(v) > 1)
print(f"  (ss, t3, sc4) determines β₁ at n=6: {'YES' if violations_sc4==0 else f'NO ({violations_sc4} violations)'}")

# Check just c3v alone (without ss and t3)
groups_c3v_only = defaultdict(set)
for d in data6:
    groups_c3v_only[d['c3v']].add(d['b1'])
violations_c3v_only = sum(1 for v in groups_c3v_only.values() if len(v) > 1)
print(f"  c3v ALONE determines β₁ at n=6: {'YES' if violations_c3v_only==0 else f'NO ({violations_c3v_only} violations)'}")

# Check sc4 alone
groups_sc4_only = defaultdict(set)
for d in data6:
    groups_sc4_only[d['sc4']].add(d['b1'])
violations_sc4_only = sum(1 for v in groups_sc4_only.values() if len(v) > 1)
print(f"  sc4 ALONE determines β₁ at n=6: {'YES' if violations_sc4_only==0 else f'NO ({violations_sc4_only} violations)'}")

# What about just (c3v, sc4)?
groups_combo = defaultdict(set)
for d in data6:
    groups_combo[(d['c3v'], d['sc4'])].add(d['b1'])
violations_combo = sum(1 for v in groups_combo.values() if len(v) > 1)
print(f"  (c3v, sc4) determines β₁ at n=6: {'YES' if violations_combo==0 else f'NO ({violations_combo} violations)'}")

# ============================================================
# STEP 4: Synthesize the formula
# ============================================================

print(f"\n" + "=" * 72)
print("STEP 4: Synthesize β₁ formula")
print("=" * 72)

# At n=5, which invariant worked best?
# Let me just check: is β₁ always 0 when transitive ratio is high enough?
# Or: β₁ = 1 iff ???

# Simple test: at n=5, what's the relationship between sc4 and β₁?
sc4_vs_b1 = defaultdict(lambda: Counter())
for d in data:
    sc4_vs_b1[d['sc4']][d['b1']] += 1

print(f"\n  n=5: SC4 count vs β₁:")
for sc4, dist in sorted(sc4_vs_b1.items()):
    total_here = sum(dist.values())
    print(f"    sc4={sc4}: β₁=0: {dist.get(0,0)}, β₁=1: {dist.get(1,0)} "
          f"(P(β₁=1) = {dist.get(1,0)/total_here:.3f})")

# At n=5, is β₁ fully determined by c3v?
c3v_vs_b1_5 = defaultdict(set)
for d in data:
    c3v_vs_b1_5[d['c3v']].add(d['b1'])
violations_c3v_5 = sum(1 for v in c3v_vs_b1_5.values() if len(v) > 1)
print(f"\n  n=5: c3v determines β₁: {'YES' if violations_c3v_5==0 else f'NO ({violations_c3v_5} violations)'}")

if violations_c3v_5 == 0:
    print(f"  ★ β₁ IS determined by per-vertex 3-cycle distribution at n=5!")
    for c3v, b1_set in sorted(c3v_vs_b1_5.items()):
        print(f"    c3v={c3v} → β₁={list(b1_set)[0]}")

print("\nDONE — beta1_discriminant_S72.py complete")
