#!/usr/bin/env python3
r"""beta2_tt_boundary_span.py - TT boundaries span im(d_2)

KEY DISCOVERY: TT (transitive triple) boundaries span ALL of im(d_2).
This means b_1 = dim(Z_1) - rk(TT boundary matrix).

A TT triple (a,b,c) has a->b->c AND a->c.
boundary(a,b,c) = (b,c) - (a,c) + (a,b) in R^{edges}.

The TT boundary matrix M_TT has:
- Rows indexed by edges (directed)
- Columns indexed by TT triples
- M_TT[(a,b), (a,b,c)] = +1
- M_TT[(b,c), (a,b,c)] = +1
- M_TT[(a,c), (a,b,c)] = -1

b_1 = dim(Z_1) - rk(M_TT restricted to Z_1)
    = (n-1)(n-2)/2 - rk(projection of M_TT columns into Z_1)

CLAIM: rk(M_TT in Z_1) >= (n-1)(n-2)/2 - 1 for all tournaments.

Can we prove this? The TT triples form a very structured set.
For a tournament on n vertices: the number of TT triples equals
the number of "transitive triples" = C(n,3) - c_3(T) where c_3 is directed 3-cycles.

At n=5: C(5,3) = 10. c_3 in {0,1,...,5}. So #TT in {5,...,10}.
dim(Z_1) = 6. So we have at least 5 TT triples, and need rk >= 5.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

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


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c3 += 1
    return c3


def is_sc(A, n):
    if n <= 1:
        return True
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[u][v] and v not in visited:
                visited.add(v)
                stack.append(v)
    if len(visited) < n:
        return False
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[v][u] and v not in visited:
                visited.add(v)
                stack.append(v)
    return len(visited) == n


def compute_tt_boundary_rank_in_Z1(A, n):
    """Compute rank of TT boundary matrix projected into Z_1."""
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    num_edges = len(edges)  # = n*(n-1)/2... wait, n*(n-1)/2 for tournament

    # Actually for a tournament, num_edges = C(n,2) = n(n-1)/2.

    # d_1 matrix: maps edges to vertices. d_1(a,b) = (b) - (a).
    D1 = np.zeros((n, num_edges))
    for idx, (a, b) in enumerate(edges):
        D1[b, idx] += 1
        D1[a, idx] -= 1

    # Z_1 = ker(D1) in R^{edges}
    U, S, Vt = np.linalg.svd(D1, full_matrices=True)
    rk_d1 = int(sum(s > 1e-8 for s in S))
    Z1_basis = Vt[rk_d1:]  # shape (dim_Z1, num_edges)
    dim_Z1 = Z1_basis.shape[0]

    # TT triples
    tt_triples = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b or not A[b][c]:
                    continue
                if A[a][c]:  # TT: a->c
                    tt_triples.append((a, b, c))

    if not tt_triples:
        return dim_Z1, 0, len(tt_triples)

    # TT boundary matrix in edge coordinates
    M_TT = np.zeros((num_edges, len(tt_triples)))
    for j, (a, b, c) in enumerate(tt_triples):
        M_TT[edge_idx[(b, c)], j] += 1
        M_TT[edge_idx[(a, c)], j] -= 1
        M_TT[edge_idx[(a, b)], j] += 1

    # Project into Z_1
    M_TT_Z1 = Z1_basis @ M_TT  # shape (dim_Z1, num_tt)
    sv = np.linalg.svd(M_TT_Z1, compute_uv=False)
    rk = int(sum(s > 1e-8 for s in sv))

    return dim_Z1, rk, len(tt_triples)


# ============================================================
# Part 1: Verify at small n
# ============================================================
print("=" * 70)
print("TT BOUNDARY RANK IN Z_1")
print("=" * 70)

for n in [3, 4, 5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    rank_deficit_dist = Counter()
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

        dim_Z1, rk, num_tt = compute_tt_boundary_rank_in_Z1(A, n)
        deficit = dim_Z1 - rk  # = b1
        c3 = count_c3(A, n)
        sc = is_sc(A, n)
        rank_deficit_dist[(deficit, sc, c3)] += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): dim_Z1={(n-1)*(n-2)//2}")
    print(f"  (deficit, SC, c3) => count")
    for (d, sc, c3), cnt in sorted(rank_deficit_dist.items()):
        print(f"  deficit={d}, SC={sc}, c3={c3}: {cnt}")


# ============================================================
# Part 2: What determines b1? Is it c3 and SC?
# ============================================================
print(f"\n{'=' * 70}")
print("b1 DETERMINANTS")
print("=" * 70)

# From the data:
# n=3: deficit=0 iff c3=0 (transitive), deficit=1 iff c3=1 (3-cycle)
# n=4: deficit=0 iff SC=False or SC=True with specific structure
# n=5: deficit=0 for all non-SC, and for SC with c3 in {3,4}
#       deficit=1 for SC with c3 in {3,4,5}
# So c3 alone doesn't determine b1.

# Let me check if #TT triples vs dim_Z1 matters.
print("\nAt n=5: num_TT vs deficit")
n = 5
tt_vs_deficit = Counter()
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

    dim_Z1, rk, num_tt = compute_tt_boundary_rank_in_Z1(A, n)
    deficit = dim_Z1 - rk
    tt_vs_deficit[(num_tt, deficit)] += 1

for (tt, d), cnt in sorted(tt_vs_deficit.items()):
    print(f"  #TT={tt}, deficit={d}: {cnt}")


# ============================================================
# Part 3: Sampled at larger n
# ============================================================
print(f"\n{'=' * 70}")
print("SAMPLED AT LARGER n")
print("=" * 70)

for n in [7, 8, 9, 10, 15]:
    random.seed(42)
    trials = max(10, min(100, 500 // n))
    deficit_dist = Counter()

    for trial in range(trials):
        A = random_tournament(n)
        dim_Z1, rk, num_tt = compute_tt_boundary_rank_in_Z1(A, n)
        deficit = dim_Z1 - rk
        deficit_dist[deficit] += 1

    print(f"\nn={n}: {trials} trials, dim_Z1={(n-1)*(n-2)//2}, C(n,3)={n*(n-1)*(n-2)//6}")
    print(f"  deficit distribution: {dict(sorted(deficit_dist.items()))}")
    print(f"  b1 <= 1: {all(d <= 1 for d in deficit_dist.keys())}")


# ============================================================
# Part 4: When #TT >= dim_Z1, is rk always = dim_Z1?
# ============================================================
print(f"\n{'=' * 70}")
print("WHEN #TT >= dim_Z1, IS rk = dim_Z1?")
print("=" * 70)

# #TT triples = 3 * (C(n,3) - c_3) for ordered TT, or (C(n,3) - c_3) for unordered.
# Wait: for each unordered triple {a,b,c} with a transitive arrangement (say a->b->c, a->c),
# there are 2 ordered TT triples: (a,b,c) (via a->b->c) is one TT path.
# Actually for a TT triple a->b->c with a->c, the allowed 2-path is just (a,b,c).
# But for the same 3-set, b->c->a is NOT an allowed 2-path (needs b->c AND c->a, but a->c).
# So each unordered transitive triple gives exactly 1 TT 2-path? No:
# If a->b, b->c, a->c: TT path is (a,b,c). That's the only one.
# But if also b->a... no, tournament: exactly one direction.

# Actually for ordered triples: for a triple (a,b,c), it's TT iff a->b, b->c, a->c.
# For 3-set {i,j,k} with one transitive arrangement (say i->j->k, i->k):
# The TT 2-paths are (i,j,k). Just one.
# So #TT = C(n,3) - c_3.

n = 5
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

    c3 = count_c3(A, n)
    tt_count = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b or not A[b][c] or not A[a][c]:
                    continue
                tt_count += 1

    expected = n*(n-1)*(n-2)//6 - c3
    # Wait, this gives unordered count. The ordered TT count should be different.
    # For each unordered transitive triple, there's exactly 1 ordered TT 2-path.
    # C(n,3) - c_3 = number of unordered transitive triples.
    # But I'm counting ORDERED TT 2-paths above...
    # For the triple a->b->c->a (3-cycle), NO TT path exists (none has both a->b->c and a->c).
    # For the triple a->b, b->c, a->c: one TT path (a,b,c).
    # For the triple a->b, b->c, c->a: 3-cycle (no TT).
    # So #(ordered TT 2-paths) = C(n,3) - c_3 = #(transitive triples).

    if tt_count != n*(n-1)*(n-2)//6 - c3:
        print(f"  MISMATCH: tt={tt_count}, expected={n*(n-1)*(n-2)//6 - c3}")
        break
else:
    print(f"n=5: #TT = C(n,3) - c_3 verified for all tournaments")


# ============================================================
# Part 5: For SC tournaments, #TT >= dim_Z1 iff c3 <= ?
# ============================================================
print(f"\n{'=' * 70}")
print("#TT vs dim_Z1")
print("=" * 70)

for n in [5, 6, 7]:
    dim_Z1 = (n-1)*(n-2)//2
    cn3 = n*(n-1)*(n-2)//6
    max_c3_for_enough_tt = cn3 - dim_Z1
    print(f"\nn={n}: dim_Z1={dim_Z1}, C(n,3)={cn3}")
    print(f"  #TT = C(n,3) - c_3 >= dim_Z1 iff c_3 <= {max_c3_for_enough_tt}")
    print(f"  Max c_3 = C(n,3)/4 = {cn3/4:.1f} (regular tournament)")


# ============================================================
# Part 6: Sum_v b1(T\v) analysis
# ============================================================
print(f"\n{'=' * 70}")
print("Sum_v b1(T\\v) ANALYSIS")
print("=" * 70)

for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    sum_b1_dist = Counter()
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

        s = 0
        for v in range(n):
            others = [x for x in range(n) if x != v]
            m = n - 1
            B = [[0] * m for _ in range(m)]
            for i in range(m):
                for j in range(m):
                    B[i][j] = A[others[i]][others[j]]

            dim_Z1, rk, _ = compute_tt_boundary_rank_in_Z1(B, m)
            b1v = dim_Z1 - rk
            s += b1v

        sum_b1_dist[s] += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): Sum_v b1(T\\v)")
    for s, cnt in sorted(sum_b1_dist.items()):
        print(f"  Sum={s}: {cnt}")
    print(f"  Max sum: {max(sum_b1_dist.keys())}")
    print(f"  Sum <= 3: {all(s <= 3 for s in sum_b1_dist.keys())}")


print("\n\nDone.")
