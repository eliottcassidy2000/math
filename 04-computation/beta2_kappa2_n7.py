"""
Quick check: kappa>=2 structure at n=7.
Also: attempt to prove |BAD| < n via counting.
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
        if d in cyc:
            continue
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

def is_sc(A, n):
    for start in range(n):
        visited = set()
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited: continue
            visited.add(v)
            for u in range(n):
                if A[v][u] and u not in visited:
                    stack.append(u)
        if len(visited) < n:
            return False
    return True

def main():
    rng = np.random.RandomState(123)

    # Test 1: Count all-bad at various n
    print("TEST 1: Can ALL vertices be bad?")
    print("="*50)

    for n in [5, 6, 7, 8]:
        nsamp = {5: 1024, 6: 32768, 7: 5000, 8: 2000}[n]
        all_bad = 0
        n_b1_0 = 0

        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample_range = range(total)
        else:
            sample_range = range(nsamp)

        for idx in sample_range:
            if n <= 6:
                A = np.zeros((n, n), dtype=int)
                bits = idx
                k = 0
                for i in range(n):
                    for j in range(i+1, n):
                        if bits & (1 << k):
                            A[i][j] = 1
                        else:
                            A[j][i] = 1
                        k += 1
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue
            n_b1_0 += 1

            all_are_bad = True
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) == 0:
                    all_are_bad = False
                    break
            if all_are_bad:
                all_bad += 1

        print(f"  n={n}: {all_bad}/{n_b1_0} have all vertices bad")

    # Test 2: For each free cycle C in T (when b1=0), how many external
    # vertices dominate it? (Distribution of dom(C))
    print(f"\nTEST 2: Domination count per cycle")
    print("="*50)

    for n in [5, 6]:
        total = 2**(n*(n-1)//2)
        dom_dist = defaultdict(int)
        total_cycles = 0

        for bits in range(total):
            A = np.zeros((n, n), dtype=int)
            k = 0
            for i in range(n):
                for j in range(i+1, n):
                    if bits & (1 << k):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    k += 1

            if compute_b1(A, n) != 0:
                continue

            cycles = find_3cycles(A, n)
            for cyc in cycles:
                dom_count = 0
                for d in range(n):
                    if d in cyc: continue
                    a, b, c = cyc
                    if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                        dom_count += 1
                dom_dist[dom_count] += 1
                total_cycles += 1

        print(f"\n  n={n} (total cycles in b1=0 tournaments: {total_cycles}):")
        for dc in sorted(dom_dist.keys()):
            pct = 100*dom_dist[dc]/total_cycles
            print(f"    dom={dc}: {dom_dist[dc]} ({pct:.1f}%)")

    # Test 3: Crucial counting argument
    # For each vertex v, define:
    # f(v) = #{cycles C not containing v with dom(C) = 0 or dom(C) = {v}}
    # = #{cycles that become free when v is deleted}
    # Sum_v f(v) = 3*#{free cycles} + #{uniquely-dominated cycles}
    # If sum_v f(v) < n * (some threshold), then not all can be bad
    print(f"\nTEST 3: Counting argument for |BAD| < n")
    print("="*50)

    n = 6
    total = 2**(n*(n-1)//2)

    for bits in range(total):
        A = np.zeros((n, n), dtype=int)
        k = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << k):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                k += 1

        if compute_b1(A, n) != 0:
            continue

        bad = set()
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.add(v)

        if len(bad) < 3:
            continue

        cycles = find_3cycles(A, n)

        # For each vertex: count cycles that become free upon deletion
        f = [0]*n
        for cyc in cycles:
            doms = []
            for d in range(n):
                if d in cyc: continue
                a, b, c = cyc
                if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                    doms.append(d)
            if len(doms) == 0:
                # Free cycle: becomes free for any deleted vertex not in it
                for v in range(n):
                    if v not in cyc:
                        f[v] += 1
            elif len(doms) == 1:
                f[doms[0]] += 1  # Only freed when deleting the unique dom

        # Print f values for bad vs good
        bad_f = [f[v] for v in bad]
        good_f = [f[v] for v in range(n) if v not in bad]

        # Just print first few
        if bits < 1000:
            print(f"  bits={bits}: bad_f={sorted(bad_f)}, good_f={sorted(good_f)}, sum_f={sum(f)}")

    # Test 4: Does v being bad require f(v) >= threshold?
    print(f"\nTEST 4: f(v) for bad vs good vertices")
    print("="*50)

    n = 6
    total = 2**(n*(n-1)//2)
    bad_f_dist = defaultdict(int)
    good_f_dist = defaultdict(int)

    for bits in range(total):
        A = np.zeros((n, n), dtype=int)
        k = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << k):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                k += 1

        if compute_b1(A, n) != 0:
            continue

        bad = set()
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.add(v)

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

        for v in range(n):
            if v in bad:
                bad_f_dist[f[v]] += 1
            else:
                good_f_dist[f[v]] += 1

    print(f"  f(v) for BAD vertices:")
    for fv in sorted(bad_f_dist.keys()):
        print(f"    f={fv}: {bad_f_dist[fv]}")
    print(f"  f(v) for GOOD vertices:")
    for fv in sorted(good_f_dist.keys()):
        print(f"    f={fv}: {good_f_dist[fv]}")

    # Is f(bad) > f(good) always? Or f(bad) >= some threshold?
    all_f_vals = sorted(set(list(bad_f_dist.keys()) + list(good_f_dist.keys())))
    min_bad_f = min(bad_f_dist.keys()) if bad_f_dist else None
    max_good_f = max(good_f_dist.keys()) if good_f_dist else None
    print(f"\n  min f(bad) = {min_bad_f}, max f(good) = {max_good_f}")
    if min_bad_f is not None and max_good_f is not None:
        print(f"  Strict separation? {min_bad_f > max_good_f}")

if __name__ == '__main__':
    main()
