"""
Verify: bad vertices always form a TRANSITIVE set.
Test at n=7,8 (sampled). Also check: can there be 4+ bad vertices at higher n?
Explore WHY transitive and WHY at most 3.
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

def is_transitive_set(A, verts):
    """Check if the subtournament on verts is transitive."""
    k = len(verts)
    if k <= 1:
        return True
    # Check: can we order them v1,...,vk with vi->vj for i<j?
    from itertools import permutations
    for perm in permutations(verts):
        ok = True
        for i in range(k):
            for j in range(i+1, k):
                if not A[perm[i]][perm[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False

def main():
    print("=" * 70)
    print("TRANSITIVE BAD SET: Verification at n=7,8")
    print("=" * 70)

    rng = np.random.RandomState(42)

    for n in [7, 8]:
        nsamp = {7: 5000, 8: 2000}[n]
        n_b1_0 = 0
        bad_count_dist = defaultdict(int)
        transitive_count = 0
        cycle_count = 0
        total_with_bad = 0

        # Track transitive ordering: does the "top" bad vertex have max score?
        top_max_score = 0
        bot_min_score = 0

        for _ in range(nsamp):
            A = random_tournament(n, rng)
            if compute_b1(A, n) != 0:
                continue
            n_b1_0 += 1

            bad = []
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad.append(v)

            bad_count_dist[len(bad)] += 1
            if not bad:
                continue
            total_with_bad += 1

            if len(bad) >= 2:
                trans = is_transitive_set(A, bad)
                if trans:
                    transitive_count += 1
                else:
                    cycle_count += 1

            if len(bad) == 3:
                # Check transitive ordering
                v1, v2, v3 = bad
                scores = [sum(A[i]) for i in range(n)]
                # Find the ordering
                for perm in [(v1,v2,v3),(v1,v3,v2),(v2,v1,v3),(v2,v3,v1),(v3,v1,v2),(v3,v2,v1)]:
                    if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                        top = perm[0]
                        bot = perm[2]
                        if scores[top] == max(scores[v] for v in bad):
                            top_max_score += 1
                        if scores[bot] == min(scores[v] for v in bad):
                            bot_min_score += 1
                        break

        n_with_multi_bad = sum(v for k, v in bad_count_dist.items() if k >= 2)
        print(f"\nn={n} (sampled {nsamp}): {n_b1_0} with b1=0")
        print(f"  Bad count distribution: {dict(bad_count_dist)}")
        print(f"  Tournaments with >=2 bad: {n_with_multi_bad}")
        print(f"  Bad set transitive: {transitive_count}/{n_with_multi_bad}")
        print(f"  Bad set has 3-cycle: {cycle_count}/{n_with_multi_bad}")
        n_triple = bad_count_dist.get(3, 0)
        if n_triple > 0:
            print(f"  Among 3-bad: top=max_score {top_max_score}/{n_triple}, bot=min_score {bot_min_score}/{n_triple}")

    # PART 2: WHY transitive? Proof sketch via domination directions
    print(f"\n{'='*70}")
    print("WHY TRANSITIVE: Domination direction analysis")
    print(f"{'='*70}")

    n = 6
    total = 2**(n*(n-1)//2)

    # For each bad pair (v,w) with v->w:
    # How does v dominate cycles? How does w dominate?
    pair_pattern = defaultdict(int)

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

        cycles = find_3cycles(A, n)

        for v in bad:
            for w in bad:
                if v == w or not A[v][w]:
                    continue
                # v->w. How does v dominate cycles not containing v?
                # Specifically cycles in the free comp of T\v that don't contain w
                v_dom_above = 0  # v->a,b,c
                v_dom_below = 0  # a,b,c->v
                w_dom_above = 0
                w_dom_below = 0

                for cyc in cycles:
                    if v in cyc:
                        continue
                    a, b, c = cyc
                    # Check if uniquely dominated by v
                    if A[v][a] and A[v][b] and A[v][c]:
                        other = False
                        for d in range(n):
                            if d in {v,a,b,c}: continue
                            if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                                other = True; break
                        if not other:
                            v_dom_above += 1
                    if A[a][v] and A[b][v] and A[c][v]:
                        other = False
                        for d in range(n):
                            if d in {v,a,b,c}: continue
                            if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                                other = True; break
                        if not other:
                            v_dom_below += 1

                for cyc in cycles:
                    if w in cyc:
                        continue
                    a, b, c = cyc
                    if A[w][a] and A[w][b] and A[w][c]:
                        other = False
                        for d in range(n):
                            if d in {w,a,b,c}: continue
                            if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                                other = True; break
                        if not other:
                            w_dom_above += 1
                    if A[a][w] and A[b][w] and A[c][w]:
                        other = False
                        for d in range(n):
                            if d in {w,a,b,c}: continue
                            if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                                other = True; break
                        if not other:
                            w_dom_below += 1

                # v->w. v dominates from? w dominates from?
                pattern = (f"v:above={v_dom_above},below={v_dom_below}",
                          f"w:above={w_dom_above},below={w_dom_below}")
                pair_pattern[pattern] += 1

    print("\nBad pair (v->w) domination patterns:")
    for pat, cnt in sorted(pair_pattern.items(), key=lambda x: -x[1])[:15]:
        print(f"  {pat}: {cnt}")

    # PART 3: The key structural theorem to prove
    print(f"\n{'='*70}")
    print("STRUCTURAL ANALYSIS: Why ≤ 3 bad vertices?")
    print(f"{'='*70}")

    # Hypothesis: Bad vertices form a chain v1->v2->...->vk (transitive).
    # The "top" vertex v1 dominates cycles from ABOVE (v1->all 3 vertices).
    # The "bottom" vertex vk dominates from BELOW (all 3 vertices->vk).
    # This means:
    #   - out(v1) contains cycles uniquely dominated by v1
    #   - in(vk) contains cycles uniquely dominated by vk
    #   - For the middle vertex v2 (if exists): it dominates from both?

    # Let's check: in the transitive ordering, does the top always dominate from above?
    n = 6
    top_above_only = 0
    bot_below_only = 0
    mid_pattern = defaultdict(int)
    n_triple = 0

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
        n_triple += 1

        # Find transitive order
        v1, v2, v3 = None, None, None
        for perm in [(bad[0],bad[1],bad[2]),(bad[0],bad[2],bad[1]),
                     (bad[1],bad[0],bad[2]),(bad[1],bad[2],bad[0]),
                     (bad[2],bad[0],bad[1]),(bad[2],bad[1],bad[0])]:
            if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                v1, v2, v3 = perm
                break

        cycles = find_3cycles(A, n)

        for label, v in [('top', v1), ('mid', v2), ('bot', v3)]:
            above = 0; below = 0
            for cyc in cycles:
                if v in cyc: continue
                a, b, c = cyc
                if A[v][a] and A[v][b] and A[v][c]:
                    other = False
                    for d in range(n):
                        if d in {v,a,b,c}: continue
                        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                            other = True; break
                    if not other: above += 1
                elif A[a][v] and A[b][v] and A[c][v]:
                    other = False
                    for d in range(n):
                        if d in {v,a,b,c}: continue
                        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                            other = True; break
                    if not other: below += 1

            if label == 'top':
                if above > 0 and below == 0:
                    top_above_only += 1
            elif label == 'bot':
                if below > 0 and above == 0:
                    bot_below_only += 1
            else:
                mid_pattern[(above, below)] += 1

    print(f"\nn=6: {n_triple} tournaments with exactly 3 bad vertices")
    print(f"  Top vertex dominates ONLY from above: {top_above_only}/{n_triple}")
    print(f"  Bottom vertex dominates ONLY from below: {bot_below_only}/{n_triple}")
    print(f"  Middle vertex pattern (above, below):")
    for pat, cnt in sorted(mid_pattern.items(), key=lambda x: -x[1]):
        print(f"    {pat}: {cnt}")

if __name__ == '__main__':
    main()
