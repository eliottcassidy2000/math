"""
sum_b1_deletion_analysis.py — Prove HYP-282: Sum_v b1(T\\v) <= 3

Key observations from previous sessions:
- b1(T) in {0, 1} always (THM-103)
- b1=1 iff strongly connected (at n>=4, equivalent conditions)
- When b1(T)=0: most vertices are "good" (b1(T\\v) <= b1(T) = 0)
- Sum_v b1(T\\v) <= 3 verified n<=10

Can we characterize EXACTLY which vertices v have b1(T\\v)=1?
Recall b1=1 iff all 3-cycles are in a single connected component
of the "shared directed edge" graph AND that component is free
(not dominated by any external vertex).

At T\\v: vertex v is removed. Which 3-cycles of T survive to T\\v?
A 3-cycle (a,b,c) in T survives in T\\v iff v not in {a,b,c}.

So b1(T\\v) depends on the 3-cycle structure of T after removing v.

Key idea: v is "bad" (b1(T\\v) = 1) iff the 3-cycle connectivity
of T\\v allows a free component.

Author: kind-pasteur-S46 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
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

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_3cycles(A, n):
    """Find all directed 3-cycles."""
    cycles = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] == 1 and A[k][i] == 1:
                    # Directed cycle i->j->k->i
                    canon = min((i,j,k), (j,k,i), (k,i,j))
                    cycles.append(canon)
    return list(set(cycles))

def is_strongly_connected(A, n):
    """Check if digraph is strongly connected."""
    # BFS from 0
    visited = set([0])
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[v][u] == 1 and u not in visited:
                visited.add(u)
                queue.append(u)
    if len(visited) < n:
        return False
    # Reverse BFS from 0
    visited = set([0])
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[u][v] == 1 and u not in visited:
                visited.add(u)
                queue.append(u)
    return len(visited) == n

def compute_b1_fast(A, n):
    """Compute b1 using the characterization: b1 = 1 iff the tournament is
    strongly connected AND has a free 3-cycle component.

    Actually, for speed, use the rank formula:
    b1 = dim(Z_1) - rank(d_2) where dim(Z_1) = C(n,2) - (n-1)

    d_2 is the boundary map restricted to TT triples only (THM-103 says
    TT boundaries span all of im(d_2)).
    """
    # Build TT triples
    tts = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[i][k] == 1 and A[j][k] == 1:
                    # i->j, i->k, j->k: transitive triple (i,j,k)
                    tts.append((i, j, k))

    # All directed edges (= Omega_1 = A_1)
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    num_edges = len(edges)

    if not tts:
        # No TT => rank(d_2) = 0 => b1 = dim(Z_1) - 0 = C(n,2) - (n-1)
        # But this is wrong for n>=3; if no TT triples then Omega_2 = A_2 = all 2-paths
        # Actually we need the full computation. For tournaments, TT exists for n>=3.
        return 0

    # Build boundary matrix for TT triples
    bd = np.zeros((num_edges, len(tts)))
    for j, (a, b, c) in enumerate(tts):
        # d_2(a,b,c) = (b,c) - (a,c) + (a,b)
        if (b, c) in edge_idx:
            bd[edge_idx[(b, c)], j] += 1
        if (a, c) in edge_idx:
            bd[edge_idx[(a, c)], j] -= 1
        if (a, b) in edge_idx:
            bd[edge_idx[(a, b)], j] += 1

    sv = np.linalg.svd(bd, compute_uv=False)
    rank_d2 = int(sum(s > 1e-8 for s in sv))
    dim_Z1 = comb(n, 2) - (n - 1)
    b1 = dim_Z1 - rank_d2
    return max(b1, 0)

def compute_b1_deletion(A, n, v):
    """Compute b1(T\\v) by removing vertex v."""
    # Create (n-1) x (n-1) adjacency matrix
    vertices = [i for i in range(n) if i != v]
    n1 = len(vertices)
    A1 = np.zeros((n1, n1), dtype=int)
    for i, vi in enumerate(vertices):
        for j, vj in enumerate(vertices):
            A1[i][j] = A[vi][vj]
    return compute_b1_fast(A1, n1)


def main():
    print("=" * 70)
    print("SUM b1(T\\v) DELETION ANALYSIS — HYP-282")
    print("=" * 70)

    # Part 1: Exhaustive at n=5,6 — distribution of sum_v b1(T\\v)
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        sum_counter = Counter()
        bad_vertex_count = Counter()
        b1_counter = Counter()

        print(f"\n--- n={n} exhaustive ({N} tournaments) ---")

        for bits in range(N):
            A = bits_to_adj(bits, n)
            b1 = compute_b1_fast(A, n)
            b1_counter[b1] += 1

            sum_b1 = 0
            bad_verts = []
            for v in range(n):
                b1v = compute_b1_deletion(A, n, v)
                if b1v > b1:
                    bad_verts.append(v)
                    sum_b1 += b1v

            total_sum = sum(compute_b1_deletion(A, n, v) for v in range(n))
            sum_counter[total_sum] += 1
            bad_vertex_count[len(bad_verts)] += 1

            if total_sum > 3:
                scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
                print(f"  VIOLATION: bits={bits}, scores={scores}, sum={total_sum}, b1={b1}")

            if (bits + 1) % 5000 == 0:
                print(f"  n={n}: {bits+1}/{N} done", flush=True)

        print(f"\n  b1 distribution: {dict(sorted(b1_counter.items()))}")
        print(f"  Sum_v b1(T\\v) distribution:")
        for s in sorted(sum_counter.keys()):
            cnt = sum_counter[s]
            pct = 100*cnt/N
            print(f"    sum={s}: {cnt} ({pct:.1f}%)")
        print(f"  #bad vertices distribution:")
        for nb in sorted(bad_vertex_count.keys()):
            print(f"    #bad={nb}: {bad_vertex_count[nb]}")
        print(f"  Max sum: {max(sum_counter.keys())}")

    # Part 2: n=7 sampled
    print("\n--- Part 2: n=7 sampled ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 500
    sum_counter7 = Counter()
    b1_by_sum = defaultdict(Counter)

    for trial in range(N):
        A = random_tournament(n, rng)
        b1 = compute_b1_fast(A, n)
        total_sum = sum(compute_b1_deletion(A, n, v) for v in range(n))
        sum_counter7[total_sum] += 1
        b1_by_sum[b1][total_sum] += 1

        if total_sum > 3:
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            print(f"  VIOLATION at n=7: trial={trial}, scores={scores}, sum={total_sum}, b1={b1}")

        if (trial + 1) % 100 == 0:
            print(f"  n=7: {trial+1}/{N} done", flush=True)

    print(f"\n  Sum_v b1(T\\v) distribution (n=7):")
    for s in sorted(sum_counter7.keys()):
        cnt = sum_counter7[s]
        pct = 100*cnt/N
        print(f"    sum={s}: {cnt} ({pct:.1f}%)")
    for b1_val in sorted(b1_by_sum.keys()):
        print(f"  When b1={b1_val}: {dict(sorted(b1_by_sum[b1_val].items()))}")

    # Part 3: n=8 sampled
    print("\n--- Part 3: n=8 sampled ---")
    n = 8
    rng8 = np.random.RandomState(99)
    N8 = 200
    sum_counter8 = Counter()

    for trial in range(N8):
        A = random_tournament(n, rng8)
        b1 = compute_b1_fast(A, n)
        total_sum = sum(compute_b1_deletion(A, n, v) for v in range(n))
        sum_counter8[total_sum] += 1

        if total_sum > 3:
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            print(f"  VIOLATION at n=8: trial={trial}, scores={scores}, sum={total_sum}")

        if (trial + 1) % 50 == 0:
            print(f"  n=8: {trial+1}/{N8} done", flush=True)

    print(f"\n  Sum_v b1(T\\v) distribution (n=8):")
    for s in sorted(sum_counter8.keys()):
        cnt = sum_counter8[s]
        pct = 100*cnt/N8
        print(f"    sum={s}: {cnt} ({pct:.1f}%)")

    # Part 4: Characterize bad vertices
    print("\n--- Part 4: Bad vertex characterization at n=6 ---")
    n = 6
    N = 2**(n*(n-1)//2)
    bad_score_dist = Counter()
    bad_c3_dist = Counter()

    for bits in range(N):
        A = bits_to_adj(bits, n)
        b1 = compute_b1_fast(A, n)
        if b1 > 0:
            continue  # Skip b1=1 case (all vertices are good)

        scores = [int(sum(A[i])) for i in range(n)]
        cycles = find_3cycles(A, n)

        for v in range(n):
            b1v = compute_b1_deletion(A, n, v)
            if b1v == 1:  # v is bad
                # Score of bad vertex
                bad_score_dist[scores[v]] += 1
                # Number of 3-cycles through v
                c3v = sum(1 for c in cycles if v in c)
                bad_c3_dist[c3v] += 1

    print(f"  Bad vertex score distribution (when b1(T)=0):")
    for sc in sorted(bad_score_dist.keys()):
        print(f"    score={sc}: {bad_score_dist[sc]}")
    print(f"  Bad vertex c3(v) distribution:")
    for c3 in sorted(bad_c3_dist.keys()):
        print(f"    c3(v)={c3}: {bad_c3_dist[c3]}")

    # Part 5: Connection between bad vertices and strongly connected components
    print("\n--- Part 5: Bad vertices and SC structure at n=6 ---")
    # When b1(T)=0 and b1(T\\v)=1: is T\\v always strongly connected?
    n = 6
    count_sc = 0
    count_not_sc = 0

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        b1 = compute_b1_fast(A, n)
        if b1 > 0:
            continue

        for v in range(n):
            b1v = compute_b1_deletion(A, n, v)
            if b1v == 1:
                # Check if T\\v is strongly connected
                vertices = [i for i in range(n) if i != v]
                n1 = len(vertices)
                A1 = np.zeros((n1, n1), dtype=int)
                for i, vi in enumerate(vertices):
                    for j, vj in enumerate(vertices):
                        A1[i][j] = A[vi][vj]
                if is_strongly_connected(A1, n1):
                    count_sc += 1
                else:
                    count_not_sc += 1

    print(f"  When b1(T)=0, b1(T\\v)=1:")
    print(f"    T\\v strongly connected: {count_sc}")
    print(f"    T\\v NOT strongly connected: {count_not_sc}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
