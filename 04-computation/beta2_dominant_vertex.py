#!/usr/bin/env python3
r"""beta2_dominant_vertex.py - Dominant vertex characterization of b1

KEY THEOREM (proved algebraically):
  If vertex d dominates all of {a,b,c} (d->a, d->b, d->c) or is
  dominated by all (a->d, b->d, c->d), then for 3-cycle C=(a,b,c):
    z(a,b) + z(b,c) + z(c,a) = 0
  for any TT-cocycle z.

PROOF: If d->a, d->b, d->c and C = a->b->c->a, then:
  TT(d,a,b): z(d,a)+z(a,b)=z(d,b)
  TT(d,b,c): z(d,b)+z(b,c)=z(d,c)
  TT(d,c,a): z(d,c)+z(c,a)=z(d,a)
  Sum: z(a,b)+z(b,c)+z(c,a) = 0.
  Similarly for d dominated by all.

COROLLARY: A 3-cycle has free cycle sum only if NO external vertex
dominates or is dominated by all 3 of its vertices.

Define: A 3-cycle C is "free" if no external vertex dominates/is-dominated.

CONJECTURE: b1 = #{cycle graph components consisting entirely of free cycles}

This script verifies the conjecture at n=4,5,6.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_3cycles(A, n):
    """Directed 3-cycles as (i,j,k) with i->j->k->i."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i, j, k))
                elif A[i][k] and A[k][j] and A[j][i]:
                    cycles.append((i, k, j))
    return cycles


def is_dominated_cycle(A, n, cycle_verts):
    """Check if any external vertex dominates or is dominated by all cycle verts."""
    a, b, c = cycle_verts
    for d in range(n):
        if d in (a, b, c):
            continue
        # d dominates all?
        if A[d][a] and A[d][b] and A[d][c]:
            return True
        # d dominated by all?
        if A[a][d] and A[b][d] and A[c][d]:
            return True
    return False


def cycle_graph_components(A, n, cycles):
    """Connected components of 3-cycle graph (edge-sharing adjacency)."""
    c3 = len(cycles)
    if c3 == 0:
        return []

    # Build edge sets
    cycle_edge_sets = []
    for (i, j, k) in cycles:
        es = set()
        for x, y in [(i,j),(j,k),(k,i),(j,i),(k,j),(i,k)]:
            if A[x][y]:
                es.add((x,y))
        cycle_edge_sets.append(es)

    # Adjacency
    adj = [[] for _ in range(c3)]
    for p in range(c3):
        for q in range(p+1, c3):
            if cycle_edge_sets[p] & cycle_edge_sets[q]:
                adj[p].append(q)
                adj[q].append(p)

    # Components
    visited = set()
    components = []
    for start in range(c3):
        if start in visited:
            continue
        comp = []
        stack = [start]
        visited.add(start)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in visited:
                    visited.add(v)
                    stack.append(v)
        components.append(comp)

    return components


def count_free_components(A, n, cycles, components):
    """Count components where ALL cycles are free (non-dominated)."""
    free_count = 0
    for comp in components:
        all_free = True
        for idx in comp:
            cycle_verts = cycles[idx]
            if is_dominated_cycle(A, n, cycle_verts):
                all_free = False
                break
        if all_free:
            free_count += 1
    return free_count


def compute_b1_tt(A, n):
    """Compute b1 via TT boundary rank in Z_1."""
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    num_edges = len(edges)

    D1 = np.zeros((n, num_edges))
    for idx, (a, b) in enumerate(edges):
        D1[b, idx] += 1
        D1[a, idx] -= 1

    U, S, Vt = np.linalg.svd(D1, full_matrices=True)
    rk_d1 = int(sum(s > 1e-8 for s in S))
    Z1_basis = Vt[rk_d1:]
    dim_Z1 = Z1_basis.shape[0]

    tt_triples = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b or not A[b][c]:
                    continue
                if A[a][c]:
                    tt_triples.append((a, b, c))

    if not tt_triples:
        return dim_Z1

    M_TT = np.zeros((num_edges, len(tt_triples)))
    for j, (a, b, c) in enumerate(tt_triples):
        M_TT[edge_idx[(b, c)], j] += 1
        M_TT[edge_idx[(a, c)], j] -= 1
        M_TT[edge_idx[(a, b)], j] += 1

    M_TT_Z1 = Z1_basis @ M_TT
    sv = np.linalg.svd(M_TT_Z1, compute_uv=False)
    rk = int(sum(s > 1e-8 for s in sv))

    return dim_Z1 - rk


# ============================================================
# Part 1: Verify b1 = #free_components at n=4,5
# ============================================================
print("=" * 70)
print("b1 vs #FREE COMPONENTS")
print("=" * 70)

for n_val in [4, 5]:
    total = 2 ** (n_val*(n_val-1)//2)
    match = 0
    mismatch = 0
    mismatch_examples = []

    for bits in range(total):
        A = [[0] * n_val for _ in range(n_val)]
        idx = 0
        for i in range(n_val):
            for j in range(i + 1, n_val):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        b1 = compute_b1_tt(A, n_val)
        cycles = get_3cycles(A, n_val)
        components = cycle_graph_components(A, n_val, cycles)
        fc = count_free_components(A, n_val, cycles, components)

        if b1 == fc:
            match += 1
        else:
            mismatch += 1
            if len(mismatch_examples) < 5:
                mismatch_examples.append((bits, b1, fc, len(cycles), len(components)))

    print(f"\nn={n_val}: {match} match, {mismatch} mismatch")
    if mismatch_examples:
        for bits, b1, fc, c3, nc in mismatch_examples:
            print(f"  bits={bits}: b1={b1}, #free_comp={fc}, c3={c3}, #comp={nc}")


# ============================================================
# Part 2: At n=6 (exhaustive)
# ============================================================
print(f"\n{'=' * 70}")
print("n=6 EXHAUSTIVE")
print("=" * 70)

n = 6
total = 2 ** (n*(n-1)//2)
match = 0
mismatch = 0
mismatch_detail = Counter()
t0 = time.time()

for bits in range(total):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    b1 = compute_b1_tt(A, n)
    cycles = get_3cycles(A, n)
    components = cycle_graph_components(A, n, cycles)
    fc = count_free_components(A, n, cycles, components)

    if b1 == fc:
        match += 1
    else:
        mismatch += 1
        mismatch_detail[(b1, fc)] += 1

elapsed = time.time() - t0
print(f"n=6 ({elapsed:.0f}s): {match} match, {mismatch} mismatch")
if mismatch_detail:
    for (b1, fc), cnt in sorted(mismatch_detail.items()):
        print(f"  b1={b1}, #free_comp={fc}: {cnt} cases")


# ============================================================
# Part 3: Sampled at n=7,8
# ============================================================
print(f"\n{'=' * 70}")
print("SAMPLED n=7,8")
print("=" * 70)

for n in [7, 8]:
    random.seed(42)
    trials = 500 if n <= 7 else 200
    match = 0
    mismatch = 0
    t0 = time.time()

    for _ in range(trials):
        A = random_tournament(n)
        b1 = compute_b1_tt(A, n)
        cycles = get_3cycles(A, n)
        components = cycle_graph_components(A, n, cycles)
        fc = count_free_components(A, n, cycles, components)

        if b1 == fc:
            match += 1
        else:
            mismatch += 1

    elapsed = time.time() - t0
    print(f"n={n} ({elapsed:.0f}s, {trials} trials): {match} match, {mismatch} mismatch")


# ============================================================
# Part 4: When free components exist, how many?
# ============================================================
print(f"\n{'=' * 70}")
print("FREE COMPONENT COUNT DISTRIBUTION")
print("=" * 70)

for n in [5, 6]:
    total = 2 ** (n*(n-1)//2)
    fc_dist = Counter()

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        cycles = get_3cycles(A, n)
        components = cycle_graph_components(A, n, cycles)
        fc = count_free_components(A, n, cycles, components)
        fc_dist[fc] += 1

    print(f"\nn={n}: #free_components distribution: {dict(sorted(fc_dist.items()))}")
    print(f"  Max free components: {max(fc_dist.keys())}")


# ============================================================
# Part 5: Can we prove max 1 free component?
# ============================================================
print(f"\n{'=' * 70}")
print("FREE COMPONENT ANALYSIS: WHY AT MOST 1?")
print("=" * 70)

# For a free 3-cycle C = (a,b,c), every external vertex d has
# score 1 or 2 in {a,b,c}. I.e., d beats exactly 1 or 2 of them.
#
# For two free cycles to coexist in different components:
# they must be in different components AND both free.
#
# Let's check: for free cycles, what is the vertex structure?

n = 5
free_cycle_info = Counter()
for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    cycles = get_3cycles(A, n)
    for cyc in cycles:
        a, b, c = cyc
        if not is_dominated_cycle(A, n, cyc):
            # Free cycle. What are the external vertex scores?
            ext_scores = []
            for d in range(n):
                if d in (a, b, c):
                    continue
                s = sum(A[d][v] for v in (a, b, c))
                ext_scores.append(s)
            free_cycle_info[tuple(sorted(ext_scores))] += 1

print(f"n=5: External vertex score patterns for free cycles:")
for pattern, cnt in sorted(free_cycle_info.items()):
    print(f"  External scores {pattern}: {cnt}")


# ============================================================
# Part 6: For n=6, check if disjoint free cycles can exist
# ============================================================
print(f"\n{'=' * 70}")
print("DISJOINT FREE CYCLES AT n=6")
print("=" * 70)

n = 6
total = 2 ** (n*(n-1)//2)
disjoint_free_count = 0
total_checked = 0

for bits in range(total):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    cycles = get_3cycles(A, n)
    free_cycles = [c for c in cycles if not is_dominated_cycle(A, n, c)]

    if len(free_cycles) < 2:
        continue

    total_checked += 1

    # Check for vertex-disjoint free cycles
    for i in range(len(free_cycles)):
        for j in range(i+1, len(free_cycles)):
            v1 = set(free_cycles[i])
            v2 = set(free_cycles[j])
            if not v1 & v2:
                disjoint_free_count += 1
                if disjoint_free_count <= 3:
                    print(f"  DISJOINT FREE: bits={bits}")
                    print(f"    C1={free_cycles[i]}, C2={free_cycles[j]}")
                    # Check components
                    components = cycle_graph_components(A, n, cycles)
                    # Find component of each
                    idx_i = cycles.index(free_cycles[i])
                    idx_j = cycles.index(free_cycles[j])
                    for k, comp in enumerate(components):
                        if idx_i in comp:
                            print(f"    C1 in component {k} (size {len(comp)})")
                        if idx_j in comp:
                            print(f"    C2 in component {k} (size {len(comp)})")
                    b1 = compute_b1_tt(A, n)
                    print(f"    b1={b1}")
                break
        else:
            continue
        break

print(f"\nn=6: {disjoint_free_count} tournaments with disjoint free cycles")
print(f"  (out of {total_checked} with >= 2 free cycles)")


print("\n\nDone.")
