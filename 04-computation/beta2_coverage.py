"""
Coverage analysis: for each vertex u, compute
coverage(u) = |{v != u : u is in some cycle C with Dom(C) <= {v}, v not in C}|

If coverage(u) < n-1 for some u, then not all vertices can be bad.
"""
import numpy as np
from itertools import combinations
from collections import defaultdict

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
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def get_dom(A, n, cyc):
    doms = []
    a, b, c = cyc
    for d in range(n):
        if d in {a,b,c}: continue
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            doms.append(d)
    return doms

def compute_coverage(A, n):
    """For each vertex u, compute coverage(u) =
    |{v != u : exists cycle C containing u with Dom(C) <= {v}, v not in C}|.

    A vertex v "covers" u if u is in some freed cycle in F_v.
    """
    cycles = find_3cycles(A, n)

    # For each cycle, compute dom
    cycle_doms = []
    for cyc in cycles:
        cycle_doms.append(set(get_dom(A, n, cyc)))

    # For each vertex u: which v's cover it?
    coverage = [set() for _ in range(n)]

    for i, cyc in enumerate(cycles):
        dom = cycle_doms[i]
        verts = set(cyc)
        if len(dom) == 0:
            # Free cycle. Covers u for all v not in cycle.
            for u in verts:
                for v in range(n):
                    if v != u and v not in verts:
                        coverage[u].add(v)
        elif len(dom) == 1:
            # Uniquely dominated by d = dom[0].
            d = list(dom)[0]
            # Covers u for v = d (if d not in cycle, which is true by definition)
            for u in verts:
                coverage[u].add(d)
        # dom >= 2: removing any single vertex doesn't free this cycle

    return [len(c) for c in coverage]

def main():
    rng = np.random.RandomState(42)

    print("=" * 60)
    print("COVERAGE ANALYSIS: Can all vertices have coverage = n-1?")
    print("=" * 60)

    for n in [5, 6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample_range = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 5000, 8: 2000}[n]
            sample_range = range(nsamp)
            method = f"sampled({nsamp})"

        n_b1_0 = 0
        max_coverage = defaultdict(int)  # distribution of max(coverage(u))
        all_full_coverage = 0  # tournaments where ALL vertices have coverage = n-1
        min_coverage_dist = defaultdict(int)  # min over vertices

        for idx in sample_range:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            cycles = find_3cycles(A, n)
            if not cycles:
                continue
            # Quick b1 check
            from itertools import combinations as comb
            def shared_de(c1, c2):
                e1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
                e2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
                return len(e1 & e2) > 0
            nc = len(cycles)
            adj = [[] for _ in range(nc)]
            for i in range(nc):
                for j in range(i+1, nc):
                    if shared_de(cycles[i], cycles[j]):
                        adj[i].append(j)
                        adj[j].append(i)
            visited = [False]*nc
            components = []
            for start in range(nc):
                if visited[start]: continue
                comp = []
                stack = [start]
                while stack:
                    vv = stack.pop()
                    if visited[vv]: continue
                    visited[vv] = True
                    comp.append(vv)
                    for uu in adj[vv]:
                        if not visited[uu]: stack.append(uu)
                components.append(comp)
            def is_dom(cyc):
                return len(get_dom(A, n, cyc)) > 0
            b1 = sum(1 for comp in components if all(not is_dom(cycles[ci]) for ci in comp))
            if b1 != 0:
                continue
            n_b1_0 += 1

            cov = compute_coverage(A, n)
            min_c = min(cov)
            max_c = max(cov)
            max_coverage[max_c] += 1
            min_coverage_dist[min_c] += 1
            if min_c == n - 1:
                all_full_coverage += 1

        print(f"\nn={n} ({method}): {n_b1_0} with b1=0")
        print(f"  Tournaments where ALL vertices have coverage = {n-1}: {all_full_coverage}/{n_b1_0}")
        print(f"  Min coverage distribution:")
        for k in sorted(min_coverage_dist.keys()):
            pct = 100*min_coverage_dist[k]/n_b1_0 if n_b1_0 > 0 else 0
            print(f"    min_coverage = {k}: {min_coverage_dist[k]} ({pct:.1f}%)")
        print(f"  Max coverage distribution:")
        for k in sorted(max_coverage.keys()):
            pct = 100*max_coverage[k]/n_b1_0 if n_b1_0 > 0 else 0
            print(f"    max_coverage = {k}: {max_coverage[k]} ({pct:.1f}%)")

    # At n=6: characterize the min-coverage vertex
    print(f"\n{'='*60}")
    print("n=6: ALWAYS some vertex with coverage < n-1?")
    print(f"{'='*60}")

    n = 6
    total = 2**(n*(n-1)//2)
    always_some_low = True

    for bits in range(total):
        A = bits_to_adj(bits, n)
        cycles = find_3cycles(A, n)
        if not cycles: continue
        nc = len(cycles)
        adj = [[] for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if shared_de(cycles[i], cycles[j]):
                    adj[i].append(j)
                    adj[j].append(i)
        visited = [False]*nc
        components = []
        for start in range(nc):
            if visited[start]: continue
            comp = []
            stack = [start]
            while stack:
                vv = stack.pop()
                if visited[vv]: continue
                visited[vv] = True
                comp.append(vv)
                for uu in adj[vv]:
                    if not visited[uu]: stack.append(uu)
            components.append(comp)
        def is_dom(cyc):
            return len(get_dom(A, n, cyc)) > 0
        b1 = sum(1 for comp in components if all(not is_dom(cycles[ci]) for ci in comp))
        if b1 != 0: continue

        cov = compute_coverage(A, n)
        if min(cov) == n - 1:
            always_some_low = False
            print(f"  ALL coverage = {n-1}! bits={bits}, scores={[sum(A[i]) for i in range(n)]}")
            break

    if always_some_low:
        print(f"  CONFIRMED: Always some vertex with coverage < {n-1}")
    else:
        print(f"  FAILED: Some tournament has all vertices with full coverage")

if __name__ == '__main__':
    main()
