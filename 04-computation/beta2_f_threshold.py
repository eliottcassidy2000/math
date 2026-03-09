"""
Verify: f(v) >= n-3 is necessary for v to be bad.
f(v) = #{cycles C not containing v with Dom(C) subset {v}}
     = #{free cycles not through v} + #{uniquely-v-dominated cycles}

If f(v) >= n-3 is necessary for bad, and we can show at most 3 vertices
have f >= n-3, then |BAD| <= 3.
"""
import numpy as np
from itertools import combinations
from collections import defaultdict

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

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

def find_3cycles(A, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def shared_directed_edge(c1, c2):
    edges1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
    edges2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
    return len(edges1 & edges2) > 0

def cycle_graph_components(cycles):
    if not cycles:
        return []
    nc = len(cycles)
    adj = [[] for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if shared_directed_edge(cycles[i], cycles[j]):
                adj[i].append(j)
                adj[j].append(i)
    visited = [False]*nc
    components = []
    for start in range(nc):
        if visited[start]:
            continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]:
                continue
            visited[v] = True
            comp.append(v)
            for u in adj[v]:
                if not visited[u]:
                    stack.append(u)
        components.append(comp)
    return components

def is_dominated(A, n, cyc):
    for d in range(n):
        if d in cyc: continue
        a, b, c = cyc
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            return True
    return False

def compute_b1(A, n):
    cycles = find_3cycles(A, n)
    if not cycles:
        return 0
    comps = cycle_graph_components(cycles)
    return sum(1 for comp in comps if all(not is_dominated(A, n, cycles[ci]) for ci in comp))

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def compute_f(A, n):
    """For each vertex v, compute f(v) = #{cycles freed by deleting v}."""
    cycles = find_3cycles(A, n)
    f = [0]*n
    for cyc in cycles:
        doms = []
        for d in range(n):
            if d in cyc: continue
            a, b, c = cyc
            if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                doms.append(d)
        if len(doms) == 0:
            for v in range(n):
                if v not in cyc:
                    f[v] += 1
        elif len(doms) == 1:
            f[doms[0]] += 1
    return f

def main():
    rng = np.random.RandomState(42)

    print("=" * 60)
    print("f(v) THRESHOLD TEST: f(bad) >= n-3 > f(good)?")
    print("=" * 60)

    for n in [5, 6, 7, 8]:
        threshold = n - 3

        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample_range = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 5000, 8: 2000}[n]
            sample_range = range(nsamp)
            method = f"sampled({nsamp})"

        n_b1_0 = 0
        bad_below_thresh = 0  # bad vertex with f < threshold (VIOLATION)
        good_at_or_above = 0  # good vertex with f >= threshold
        total_bad = 0
        total_good = 0
        min_bad_f = float('inf')
        max_good_f = 0

        for idx in sample_range:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue
            n_b1_0 += 1

            f = compute_f(A, n)
            for v in range(n):
                B = delete_vertex(A, n, v)
                is_bad = compute_b1(B, n-1) > 0
                if is_bad:
                    total_bad += 1
                    if f[v] < threshold:
                        bad_below_thresh += 1
                    min_bad_f = min(min_bad_f, f[v])
                else:
                    total_good += 1
                    if f[v] >= threshold:
                        good_at_or_above += 1
                    max_good_f = max(max_good_f, f[v])

        print(f"\nn={n} ({method}), threshold=n-3={threshold}")
        print(f"  b1=0 tournaments: {n_b1_0}")
        print(f"  Total bad vertices: {total_bad}, good: {total_good}")
        print(f"  Bad with f < {threshold}: {bad_below_thresh} (VIOLATIONS)")
        print(f"  Good with f >= {threshold}: {good_at_or_above}")
        print(f"  min f(bad) = {min_bad_f if min_bad_f < float('inf') else 'N/A'}")
        print(f"  max f(good) = {max_good_f}")
        if min_bad_f < float('inf'):
            print(f"  STRICT SEPARATION: {min_bad_f > max_good_f}")

    # DETAILED: At n=7, show f(v) distribution for bad vs good
    print(f"\n{'='*60}")
    print("f(v) DISTRIBUTION at n=7")
    print(f"{'='*60}")

    n = 7
    bad_f_dist = defaultdict(int)
    good_f_dist = defaultdict(int)

    for _ in range(5000):
        A = random_tournament(n, rng)
        if compute_b1(A, n) != 0:
            continue
        f = compute_f(A, n)
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad_f_dist[f[v]] += 1
            else:
                good_f_dist[f[v]] += 1

    print("  BAD vertex f distribution:")
    for fv in sorted(bad_f_dist.keys()):
        print(f"    f={fv}: {bad_f_dist[fv]}")
    print("  GOOD vertex f distribution:")
    for fv in sorted(good_f_dist.keys()):
        print(f"    f={fv}: {good_f_dist[fv]}")

    # Can we prove: |{v: f(v) >= n-3}| <= 3?
    print(f"\n{'='*60}")
    print("COUNT OF VERTICES WITH f >= n-3")
    print(f"{'='*60}")

    for n in [5, 6, 7, 8]:
        threshold = n - 3
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample_range = range(total)
        else:
            nsamp = {7: 5000, 8: 2000}[n]
            sample_range = range(nsamp)

        count_above_dist = defaultdict(int)
        n_b1_0 = 0

        for idx in sample_range:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue
            n_b1_0 += 1

            f = compute_f(A, n)
            count_above = sum(1 for v in range(n) if f[v] >= threshold)
            count_above_dist[count_above] += 1

        print(f"\n  n={n}, threshold={threshold}:")
        for k in sorted(count_above_dist.keys()):
            pct = 100*count_above_dist[k]/n_b1_0
            print(f"    {k} vertices with f >= {threshold}: {count_above_dist[k]} ({pct:.1f}%)")

if __name__ == '__main__':
    main()
