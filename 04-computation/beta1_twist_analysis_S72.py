#!/usr/bin/env python3
"""
beta1_twist_analysis_S72.py — opus-2026-03-15-S72

THE MÖBIUS TWIST: β₁=1 with full edge coverage
=================================================

At n=5: 144 tournaments have β₁=1 despite ALL edges being covered
by transitive triples. These are the "twisted" cases.

Questions:
1. What distinguishes twisted from untwisted β₁=1 cases?
2. What is the actual obstruction preventing cycle-filling?
3. Can we characterize the twist algebraically?
"""

import numpy as np
from itertools import combinations
import sys
import time
sys.stdout.reconfigure(line_buffering=True)

def adj_from_bits(n, bits):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

PRIME = 32749

def rank_mod(M, p):
    M = M.copy() % p
    nrows, ncols = M.shape
    r = 0
    for col in range(ncols):
        nz = np.where(M[r:, col] != 0)[0]
        if len(nz) == 0:
            continue
        pivot = nz[0] + r
        if pivot != r:
            M[[r, pivot]] = M[[pivot, r]]
        inv = pow(int(M[r, col]), p - 2, p)
        M[r] = M[r] * inv % p
        factors = M[:, col].copy()
        factors[r] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[r])) % p
        r += 1
    return r

def nullspace_mod(M, p):
    """Return basis for null space of M mod p."""
    M = M.copy() % p
    nrows, ncols = M.shape
    pivot_cols = []
    r = 0
    for col in range(ncols):
        nz = np.where(M[r:, col] != 0)[0]
        if len(nz) == 0:
            continue
        pivot = nz[0] + r
        if pivot != r:
            M[[r, pivot]] = M[[pivot, r]]
        inv = pow(int(M[r, col]), p - 2, p)
        M[r] = M[r] * inv % p
        factors = M[:, col].copy()
        factors[r] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[r])) % p
        pivot_cols.append(col)
        r += 1

    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    if not free_cols:
        return np.zeros((0, ncols), dtype=np.int64)
    basis = np.zeros((len(free_cols), ncols), dtype=np.int64)
    for i, fc in enumerate(free_cols):
        basis[i, fc] = 1
        for rr, pc in enumerate(pivot_cols):
            if rr < M.shape[0]:
                basis[i, pc] = (-M[rr, fc]) % p
    return basis

def full_analysis(A, n):
    """Compute β₁ and extract cycle representatives, coverage info."""
    edges = []
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j] == 1:
                edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    ne = len(edges)

    # Transitive triples
    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[j][k] == 1 and A[i][k] == 1:
                    triples.append((i, j, k))
    nt = len(triples)

    # ∂₁: edges → vertices
    d1 = np.zeros((n, ne), dtype=np.int64)
    for idx, (i, j) in enumerate(edges):
        d1[j, idx] += 1
        d1[i, idx] -= 1

    # ∂₂: triples → edges
    d2 = np.zeros((ne, nt), dtype=np.int64)
    for t_idx, (i, j, k) in enumerate(triples):
        if (j, k) in edge_idx:
            d2[edge_idx[(j, k)], t_idx] += 1
        if (i, k) in edge_idx:
            d2[edge_idx[(i, k)], t_idx] -= 1
        if (i, j) in edge_idx:
            d2[edge_idx[(i, j)], t_idx] += 1

    rank_d1 = rank_mod(d1 % PRIME, PRIME)
    rank_d2 = rank_mod(d2 % PRIME, PRIME)
    ker_d1 = ne - rank_d1
    beta1 = ker_d1 - rank_d2

    # Edge coverage
    covered = set()
    for (i, j, k) in triples:
        covered.add((j, k))
        covered.add((i, k))
        covered.add((i, j))
    uncovered = set(edges) - covered

    # Find cycle representatives if β₁ > 0
    cycles = []
    if beta1 > 0:
        # ker(d1) basis
        ker_basis = nullspace_mod(d1 % PRIME, PRIME)
        # For each ker vector, check if it's in im(d2)
        for vec in ker_basis:
            augmented = np.hstack([d2, vec.reshape(-1, 1)])
            if rank_mod(augmented % PRIME, PRIME) > rank_d2:
                # This vector is NOT in im(d2) — it's a cycle representative
                cycle_edges = []
                for idx in range(ne):
                    coeff = int(vec[idx]) % PRIME
                    if coeff != 0:
                        if coeff > PRIME // 2:
                            coeff -= PRIME
                        cycle_edges.append((edges[idx], coeff))
                cycles.append(cycle_edges)

    return {
        'beta1': beta1,
        'edges': edges,
        'triples': triples,
        'uncovered': uncovered,
        'rank_d2': rank_d2,
        'nt': nt,
        'cycles': cycles,
    }

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def t3_count(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                fwd = A[i][j] * A[j][k] * A[k][i]
                bwd = A[j][i] * A[i][k] * A[k][j]
                count += fwd + bwd
    return count

print("=" * 72)
print("THE MÖBIUS TWIST — β₁=1 with full coverage")
print("opus-2026-03-15-S72")
print("=" * 72)

n = 5
num_edges = n * (n - 1) // 2
total = 1 << num_edges

# Classify all tournaments
twisted = []     # β₁=1, full coverage
untwisted = []   # β₁=1, has uncovered edges
clean = []       # β₁=0

t0 = time.time()
for bits in range(total):
    A = adj_from_bits(n, bits)
    info = full_analysis(A, n)

    if info['beta1'] == 0:
        clean.append((bits, info))
    elif len(info['uncovered']) == 0:
        twisted.append((bits, info))
    else:
        untwisted.append((bits, info))

elapsed = time.time() - t0
print(f"\nn=5: {len(clean)} clean, {len(twisted)} twisted, {len(untwisted)} untwisted ({elapsed:.1f}s)")

# Analyze twisted cases
print(f"\n{'='*60}")
print(f"TWISTED CASES (β₁=1, all edges covered)")
print(f"{'='*60}")

tw_scores = {}
tw_t3 = {}
tw_nt = {}
for bits, info in twisted:
    A = adj_from_bits(n, bits)
    s = score_seq(A, n)
    t = t3_count(A, n)
    tw_scores[s] = tw_scores.get(s, 0) + 1
    tw_t3[t] = tw_t3.get(t, 0) + 1
    tw_nt[info['nt']] = tw_nt.get(info['nt'], 0) + 1

print(f"\n  Score sequences: {tw_scores}")
print(f"  t3 distribution: {tw_t3}")
print(f"  #triples distribution: {tw_nt}")

# Compare with untwisted
print(f"\n{'='*60}")
print(f"UNTWISTED CASES (β₁=1, has uncovered edges)")
print(f"{'='*60}")

utw_scores = {}
utw_t3 = {}
utw_nt = {}
for bits, info in untwisted:
    A = adj_from_bits(n, bits)
    s = score_seq(A, n)
    t = t3_count(A, n)
    utw_scores[s] = utw_scores.get(s, 0) + 1
    utw_t3[t] = utw_t3.get(t, 0) + 1
    utw_nt[info['nt']] = utw_nt.get(info['nt'], 0) + 1

print(f"\n  Score sequences: {utw_scores}")
print(f"  t3 distribution: {utw_t3}")
print(f"  #triples distribution: {utw_nt}")

# Deep dive into a twisted example
print(f"\n{'='*60}")
print(f"DEEP DIVE: One twisted tournament")
print(f"{'='*60}")

bits, info = twisted[0]
A = adj_from_bits(n, bits)
print(f"\n  bits={bits}, score_seq={score_seq(A, n)}, t3={t3_count(A, n)}")
print(f"  Adjacency:")
for i in range(n):
    print(f"    {list(A[i])}")
print(f"  #edges={len(info['edges'])}, #triples={info['nt']}, rank(∂₂)={info['rank_d2']}")
print(f"  Uncovered edges: {info['uncovered']}")

if info['cycles']:
    print(f"\n  Unbounding cycle(s):")
    for ci, cycle in enumerate(info['cycles']):
        print(f"    Cycle {ci}: {cycle}")
        # Which edges in the cycle?
        cycle_edge_set = set(e for e, c in cycle)
        print(f"    Edges in cycle: {cycle_edge_set}")

        # For each edge in cycle, which triples cover it?
        for e, coeff in cycle:
            covering = []
            for t in info['triples']:
                i, j, k = t
                if e == (j, k) or e == (i, k) or e == (i, j):
                    covering.append(t)
            print(f"      Edge {e} (coeff={coeff}): covered by {len(covering)} triples: {covering}")

# Now the key question: WHY can't these triples fill the cycle?
print(f"\n{'='*60}")
print(f"WHY CAN'T TRIPLES FILL THE CYCLE?")
print(f"{'='*60}")

print(f"""
  The cycle is a vector c ∈ ker(∂₁) with c ∉ im(∂₂).
  Even though every edge in c is covered by some triple,
  no LINEAR COMBINATION of triple boundaries equals c.

  This is like: every side of a Möbius strip is locally
  a piece of paper, but you can't tile it consistently.

  Let's see the actual linear algebra obstruction:
""")

bits, info = twisted[0]
A = adj_from_bits(n, bits)
edges = info['edges']
triples = info['triples']
edge_idx = {e: i for i, e in enumerate(edges)}
ne = len(edges)
nt = len(triples)

d2 = np.zeros((ne, nt), dtype=np.int64)
for t_idx, (i, j, k) in enumerate(triples):
    if (j, k) in edge_idx:
        d2[edge_idx[(j, k)], t_idx] += 1
    if (i, k) in edge_idx:
        d2[edge_idx[(i, k)], t_idx] -= 1
    if (i, j) in edge_idx:
        d2[edge_idx[(i, j)], t_idx] += 1

d1 = np.zeros((n, ne), dtype=np.int64)
for idx, (i, j) in enumerate(edges):
    d1[j, idx] += 1
    d1[i, idx] -= 1

ker_basis = nullspace_mod(d1 % PRIME, PRIME)
print(f"  dim(ker ∂₁) = {len(ker_basis)}")
print(f"  rank(∂₂) = {info['rank_d2']}")
print(f"  β₁ = {len(ker_basis)} - {info['rank_d2']} = {len(ker_basis) - info['rank_d2']}")

# Show the cycle and what im(∂₂) looks like
print(f"\n  Boundary matrix ∂₂ ({ne}×{nt}):")
print(f"  Edges: {edges}")
print(f"  Triples: {triples}")

# Show the cycle
for vec in ker_basis:
    augmented = np.hstack([d2, vec.reshape(-1, 1)])
    if rank_mod(augmented % PRIME, PRIME) > info['rank_d2']:
        print(f"\n  Unbounding cycle vector:")
        for idx in range(ne):
            v = int(vec[idx]) % PRIME
            if v > PRIME // 2:
                v -= PRIME
            if v != 0:
                print(f"    {v:+d} × {edges[idx]}")

        # Try to express it as sum of triple boundaries
        print(f"\n  Triple boundaries:")
        for t_idx, (i, j, k) in enumerate(triples):
            boundary = []
            for e_idx in range(ne):
                v = int(d2[e_idx, t_idx])
                if v != 0:
                    boundary.append(f"{v:+d}×{edges[e_idx]}")
            print(f"    ∂({i},{j},{k}) = {' '.join(boundary)}")

        break

# The crucial structure: look at 3-cycles vs transitive triples
print(f"\n{'='*60}")
print(f"3-CYCLES vs TRANSITIVE TRIPLES")
print(f"{'='*60}")

for label, cases in [("twisted", twisted[:3]), ("untwisted", untwisted[:3])]:
    print(f"\n  --- {label} examples ---")
    for bits, info in cases:
        A = adj_from_bits(n, bits)

        # Count 3-cycles
        three_cycles = []
        for i in range(n):
            for j in range(n):
                if i == j or A[i][j] == 0:
                    continue
                for k in range(n):
                    if k == i or k == j:
                        continue
                    if A[j][k] == 1 and A[k][i] == 1:
                        cyc = tuple(sorted([i, j, k]))
                        if cyc not in [tuple(sorted(c)) for c in three_cycles]:
                            three_cycles.append((i, j, k))

        tt = info['nt']
        tc = len(three_cycles)
        print(f"    bits={bits}: t3={t3_count(A,n)}, #trans_triples={tt}, 3-cycles={tc}, ratio={tt/max(tc,1):.2f}")

# Comprehensive: is there a simple invariant separating twisted from clean+untwisted?
print(f"\n{'='*60}")
print(f"SEARCHING FOR TWIST INVARIANT")
print(f"{'='*60}")

# For each tournament, compute rank(∂₂), dim(ker ∂₁), β₁, and some new candidates
print(f"\n  Comparing rank(∂₂) distributions:")

clean_ranks = {}
twisted_ranks = {}
untw_ranks = {}

for bits, info in clean:
    r = info['rank_d2']
    clean_ranks[r] = clean_ranks.get(r, 0) + 1
for bits, info in twisted:
    r = info['rank_d2']
    twisted_ranks[r] = twisted_ranks.get(r, 0) + 1
for bits, info in untwisted:
    r = info['rank_d2']
    untw_ranks[r] = untw_ranks.get(r, 0) + 1

print(f"  Clean (β₁=0):    rank(∂₂) dist = {dict(sorted(clean_ranks.items()))}")
print(f"  Twisted (β₁=1):  rank(∂₂) dist = {dict(sorted(twisted_ranks.items()))}")
print(f"  Untwisted (β₁=1): rank(∂₂) dist = {dict(sorted(untw_ranks.items()))}")

# What about the PARITY of #triples?
print(f"\n  Comparing #triples distributions:")
for label, cases in [("Clean", clean), ("Twisted", twisted), ("Untwisted", untwisted)]:
    nt_dist = {}
    for bits, info in cases:
        nt = info['nt']
        nt_dist[nt] = nt_dist.get(nt, 0) + 1
    print(f"  {label:12s}: #triples dist = {dict(sorted(nt_dist.items()))}")

# The definitive test: what about CYCLE-TRANSITIVITY RATIO?
# For each directed 3-cycle, how many of its edges participate in transitive triples?
print(f"\n{'='*60}")
print(f"CYCLE-EDGE PARTICIPATION IN TRANSITIVE TRIPLES")
print(f"{'='*60}")

for label, cases in [("Twisted", twisted[:5]), ("Untwisted", untwisted[:5])]:
    print(f"\n  --- {label} ---")
    for bits, info in cases:
        A = adj_from_bits(n, bits)
        # Find 3-cycles
        three_cyc_edges = set()
        for i in range(n):
            for j in range(n):
                if i == j or A[i][j] == 0: continue
                for k in range(n):
                    if k == i or k == j: continue
                    if A[j][k] == 1 and A[k][i] == 1:
                        three_cyc_edges.add((i,j))
                        three_cyc_edges.add((j,k))
                        three_cyc_edges.add((k,i))

        # How many 3-cycle edges are ALSO in transitive triples?
        trans_edges = set()
        for (i, j, k) in info['triples']:
            trans_edges.add((i,j))
            trans_edges.add((j,k))
            trans_edges.add((i,k))

        both = three_cyc_edges & trans_edges
        only_cyc = three_cyc_edges - trans_edges

        print(f"    bits={bits}: 3-cyc edges={len(three_cyc_edges)}, trans edges={len(trans_edges)}, "
              f"both={len(both)}, only-cycle={len(only_cyc)}")

# The real question: what LINEAR DEPENDENCE in ∂₂ causes the twist?
print(f"\n{'='*60}")
print(f"LINEAR DEPENDENCE STRUCTURE")
print(f"{'='*60}")

print(f"\n  For a twisted tournament, rank(∂₂) = 5 (one less than the 6 it needs)")
print(f"  This means there's a linear relation among the triple boundaries.")
print(f"  What IS this relation?")

bits, info = twisted[0]
A = adj_from_bits(n, bits)
edges = info['edges']
triples = info['triples']
edge_idx = {e: i for i, e in enumerate(edges)}
ne = len(edges)
nt = len(triples)

d2 = np.zeros((ne, nt), dtype=np.int64)
for t_idx, (i, j, k) in enumerate(triples):
    if (j, k) in edge_idx:
        d2[edge_idx[(j, k)], t_idx] += 1
    if (i, k) in edge_idx:
        d2[edge_idx[(i, k)], t_idx] -= 1
    if (i, j) in edge_idx:
        d2[edge_idx[(i, j)], t_idx] += 1

# Find the null space of d2^T (left null space of d2 = cokernel)
# No wait — we want ker(d2^T): linear relations among rows of d2
# Actually: ker of d2 as a map R^nt → R^ne
# That tells us which combinations of triples have zero boundary
ns = nullspace_mod(d2.T % PRIME, PRIME)
print(f"\n  Left nullspace of ∂₂ (relations among rows):")
print(f"    dim = {len(ns)} (these are the relations)")

# Actually we want right nullspace: ker(d2) = relations among triple columns
ns_right = nullspace_mod(d2 % PRIME, PRIME)
print(f"\n  Right nullspace of ∂₂ (2-cycles = relations among triples):")
print(f"    dim = {len(ns_right)}")
for vi, vec in enumerate(ns_right):
    rel = []
    for t_idx in range(nt):
        v = int(vec[t_idx]) % PRIME
        if v > PRIME // 2:
            v -= PRIME
        if v != 0:
            rel.append(f"{v:+d}×∂({triples[t_idx]})")
    if rel:
        print(f"    Relation {vi}: {' '.join(rel)} = 0")
        # What does this mean? These triples form a closed surface (2-cycle)
        # Their boundaries cancel. This is a 2-cycle!

        # Which vertices are involved?
        verts = set()
        for t_idx in range(nt):
            v = int(vec[t_idx]) % PRIME
            if v > PRIME // 2:
                v -= PRIME
            if v != 0:
                for x in triples[t_idx]:
                    verts.add(x)
        print(f"           Vertices: {sorted(verts)}")

# Compare with clean tournament
bits_clean, info_clean = clean[0]
A_clean = adj_from_bits(n, bits_clean)
triples_clean = info_clean['triples']
edge_idx_clean = {e: i for i, e in enumerate(info_clean['edges'])}
ne_c = len(info_clean['edges'])
nt_c = len(triples_clean)

d2_c = np.zeros((ne_c, nt_c), dtype=np.int64)
for t_idx, (i, j, k) in enumerate(triples_clean):
    if (j, k) in edge_idx_clean:
        d2_c[edge_idx_clean[(j, k)], t_idx] += 1
    if (i, k) in edge_idx_clean:
        d2_c[edge_idx_clean[(i, k)], t_idx] -= 1
    if (i, j) in edge_idx_clean:
        d2_c[edge_idx_clean[(i, j)], t_idx] += 1

ns_c = nullspace_mod(d2_c % PRIME, PRIME)
print(f"\n  For comparison, clean tournament bits={bits_clean}:")
print(f"    #triples={nt_c}, rank(∂₂)={info_clean['rank_d2']}")
print(f"    2-cycles (right nullspace dim) = {len(ns_c)}")

print(f"\n{'='*60}")
print(f"FINAL SYNTHESIS")
print(f"{'='*60}")
print(f"""
  At n=5:
  - 720 tournaments have β₁=0 (all edges covered, all 1-cycles bound)
  - 160 tournaments have β₁=1 with uncovered edge (obvious obstruction)
  - 144 tournaments have β₁=1 with full coverage (THE TWIST)

  The 144 twisted cases have:
  - All edges covered by transitive triples
  - But the triple system has a "topological twist"
  - The boundaries of triples can cover every edge individually
    but cannot SIMULTANEOUSLY fill a 1-cycle
  - This is exactly like a Möbius strip: locally orientable,
    globally non-orientable

  The 2-cycles (right nullspace of ∂₂) reveal the structure:
  - Relations among triples show which surfaces are "closed"
  - In twisted tournaments, these closed surfaces don't
    span all of ker(∂₁)
""")

print("DONE — beta1_twist_analysis_S72.py complete")
