"""
When exactly 3 bad vertices exist, do they form a 3-cycle?
Also: prove that non-SC tournaments always have b1=0 (all deletions good).
And characterize the kappa>=2 case.
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
    verts = set(cyc)
    for d in range(n):
        if d in verts:
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
    """Check if tournament is strongly connected."""
    for start in range(n):
        visited = set()
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            for u in range(n):
                if A[v][u] and u not in visited:
                    stack.append(u)
        if len(visited) < n:
            return False
    return True

def vertex_connectivity(A, n):
    """Compute vertex connectivity of tournament."""
    if not is_sc(A, n):
        return 0
    if n <= 2:
        return n - 1
    for v in range(n):
        B = delete_vertex(A, n, v)
        if not is_sc(B, n-1):
            return 1
    return 2  # at least 2, good enough for our purposes

def main():
    print("=" * 70)
    print("BAD TRIPLE STRUCTURE + NON-SC VERIFICATION")
    print("=" * 70)

    # PART 1: Verify non-SC => all deletions good
    print("\nPART 1: Non-SC tournaments always have b1=0 for all deletions")
    n = 6
    total = 2**(n*(n-1)//2)
    non_sc_count = 0
    non_sc_all_good = 0
    non_sc_b1_nonzero = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if is_sc(A, n):
            continue
        non_sc_count += 1
        b1 = compute_b1(A, n)
        if b1 != 0:
            non_sc_b1_nonzero += 1
            continue
        all_good = True
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                all_good = False
                break
        if all_good:
            non_sc_all_good += 1

    print(f"  n=6: {non_sc_count} non-SC tournaments")
    print(f"  b1 != 0 among non-SC: {non_sc_b1_nonzero}")
    print(f"  b1=0 and all deletions good: {non_sc_all_good}/{non_sc_count - non_sc_b1_nonzero}")

    # PART 2: When 3 bad vertices exist, do they form a 3-cycle?
    print(f"\n{'='*70}")
    print("PART 2: Structure of 3 bad vertices")
    print(f"{'='*70}")

    for n in [5, 6]:
        total = 2**(n*(n-1)//2)
        triple_is_cycle = 0
        triple_total = 0
        triple_is_trans = 0
        triple_score_patterns = defaultdict(int)
        triple_in_same_scc = 0

        for bits in range(total):
            A = bits_to_adj(bits, n)
            if compute_b1(A, n) != 0:
                continue

            bad = []
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad.append(v)

            if len(bad) != 3:
                continue
            triple_total += 1

            v1, v2, v3 = bad
            # Check if {v1,v2,v3} forms a 3-cycle
            if (A[v1][v2] and A[v2][v3] and A[v3][v1]) or \
               (A[v1][v3] and A[v3][v2] and A[v2][v1]):
                triple_is_cycle += 1

            # Check if transitive
            if (A[v1][v2] and A[v2][v3] and A[v1][v3]) or \
               (A[v1][v3] and A[v3][v2] and A[v1][v2]) or \
               (A[v2][v1] and A[v1][v3] and A[v2][v3]) or \
               (A[v2][v3] and A[v3][v1] and A[v2][v1]) or \
               (A[v3][v1] and A[v1][v2] and A[v3][v2]) or \
               (A[v3][v2] and A[v2][v1] and A[v3][v1]):
                triple_is_trans += 1

            scores = [sum(A[i]) for i in range(n)]
            bad_scores = tuple(sorted([scores[v] for v in bad]))
            triple_score_patterns[bad_scores] += 1

        print(f"\nn={n}: {triple_total} tournaments with exactly 3 bad vertices")
        print(f"  3 bad vertices form a 3-cycle: {triple_is_cycle}/{triple_total}")
        print(f"  3 bad vertices form transitive: {triple_is_trans}/{triple_total}")
        print(f"  Score patterns:")
        for pat, cnt in sorted(triple_score_patterns.items(), key=lambda x: -x[1]):
            print(f"    {pat}: {cnt}")

    # PART 3: For bad pairs, check domination structure
    print(f"\n{'='*70}")
    print("PART 3: Bad pairs — arc direction and domination")
    print(f"{'='*70}")

    n = 6
    pair_both_up = 0  # both dominate from above
    pair_both_down = 0
    pair_up_down = 0  # one up one down
    pair_total = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)

        if len(bad) < 2:
            continue

        for v, w in combinations(bad, 2):
            pair_total += 1
            # How does v dominate cycles? Check if v beats all of some uniquely-dominated cycle
            cycles = find_3cycles(A, n)
            v_up = False; v_down = False; w_up = False; w_down = False
            for cyc in cycles:
                if v in cyc or w in cyc:
                    continue
                a, b, c = cyc
                # v dominates from above?
                if A[v][a] and A[v][b] and A[v][c]:
                    # unique?
                    other = False
                    for d in range(n):
                        if d in {v, a, b, c}:
                            continue
                        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                            other = True; break
                    if not other:
                        v_up = True
                if A[a][v] and A[b][v] and A[c][v]:
                    other = False
                    for d in range(n):
                        if d in {v, a, b, c}:
                            continue
                        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                            other = True; break
                    if not other:
                        v_down = True
                if A[w][a] and A[w][b] and A[w][c]:
                    other = False
                    for d in range(n):
                        if d in {w, a, b, c}:
                            continue
                        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                            other = True; break
                    if not other:
                        w_up = True
                if A[a][w] and A[b][w] and A[c][w]:
                    other = False
                    for d in range(n):
                        if d in {w, a, b, c}:
                            continue
                        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                            other = True; break
                    if not other:
                        w_down = True

    print(f"  Total bad pairs: {pair_total}")

    # PART 4: Kappa analysis
    print(f"\n{'='*70}")
    print("PART 4: Connectivity vs bad count")
    print(f"{'='*70}")

    n = 6
    kappa_bad = defaultdict(lambda: defaultdict(int))
    sc_count = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue
        if not is_sc(A, n):
            continue  # already handled
        sc_count += 1

        kap = vertex_connectivity(A, n)

        bad_count = 0
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad_count += 1

        kappa_bad[kap][bad_count] += 1

    print(f"\nn=6 SC tournaments with b1=0: {sc_count}")
    for kap in sorted(kappa_bad.keys()):
        print(f"  kappa={kap}:")
        for bc in sorted(kappa_bad[kap].keys()):
            print(f"    {bc} bad: {kappa_bad[kap][bc]}")

    # PART 5: CRITICAL — Do bad vertices form EXACTLY the non-dominant vertices?
    print(f"\n{'='*70}")
    print("PART 5: Bad vertex = vertex that uniquely dominates a cycle?")
    print(f"{'='*70}")

    n = 6
    # For each vertex: does it uniquely dominate at least one cycle?
    perfect_match = 0
    mismatch = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        # Compute bad vertices
        bad = set()
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.add(v)

        # Vertices that uniquely dominate at least one cycle
        uniq_dom = set()
        for cyc in cycles:
            doms = []
            for d in range(n):
                if d in cyc:
                    continue
                a, b, c = cyc
                if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                    doms.append(d)
            if len(doms) == 1:
                uniq_dom.add(doms[0])

        if bad == uniq_dom:
            perfect_match += 1
        else:
            mismatch += 1
            if mismatch <= 3:
                print(f"  MISMATCH: bad={bad}, uniq_dom={uniq_dom}")
                scores = [sum(A[i]) for i in range(n)]
                print(f"    scores={scores}")

    total_b1_0_sc = sum(sum(kappa_bad[k].values()) for k in kappa_bad)
    non_sc_b1_0 = 27968 - total_b1_0_sc  # approximate

    print(f"\n  SC b1=0 tournaments: {sc_count}")
    print(f"  bad == uniq_dom: {perfect_match}")
    print(f"  mismatch: {mismatch}")

if __name__ == '__main__':
    main()
