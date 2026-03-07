#!/usr/bin/env python3
"""
At n=9, if t3=3 and every vertex in a 3-cycle,
then the 3 three-cycles partition [9] into 3 disjoint triples.
What is the minimum alpha_1?

With 3 disjoint triples: 9 inter-block arcs per pair of triples,
total 27 inter-block arcs. These create additional cycles.

Let's enumerate: fix blocks {0,1,2}, {3,4,5}, {6,7,8}.
Try all tournaments with exactly 1 directed 3-cycle per block
and check alpha_1.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations, product
from collections import Counter
from math import comb

def find_directed_cycles_dp(adj, n, k):
    result = []
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        full = (1 << k) - 1
        for S in range(1, full + 1):
            for i in range(k):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(k):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, k):
            if (full, j) in dp and adj[v[j]][v[0]]:
                count += dp[(full, j)]
        if count > 0:
            result.append((frozenset(verts), count))
    return result


def main():
    n = 9
    blocks = [{0,1,2}, {3,4,5}, {6,7,8}]
    block_verts = [[0,1,2], [3,4,5], [6,7,8]]

    # Each block has 3 arcs. A directed 3-cycle on {a,b,c} requires
    # exactly one of the two cyclic orientations: a->b->c->a or a->c->b->a.
    # For block {a,b,c}: 2^3 = 8 tournaments on 3 vertices,
    # of which 2 are 3-cycles and 6 are transitive.
    # We need t3=1 per block, so each block is a 3-cycle (2 orientations).
    # Total intra-block: 2^3 = 8 choices.

    # Inter-block arcs: for each pair of blocks (3*3=9 arcs per pair),
    # there are 3 pairs of blocks = 27 inter-block arcs.
    # Each arc has 2 orientations: 2^27 possibilities. Too many!

    # Let's sample instead.
    import random
    random.seed(42)

    print(f"=== n=9: 3 disjoint 3-cycles partitioning vertices ===")
    print(f"Blocks: {blocks}")

    alpha_vals = []
    comp_vals = []

    for trial in range(100000):
        adj = [[0]*n for _ in range(n)]

        # Set each block to a directed 3-cycle
        for blk in block_verts:
            if random.random() < 0.5:
                # a->b->c->a
                adj[blk[0]][blk[1]] = 1
                adj[blk[1]][blk[2]] = 1
                adj[blk[2]][blk[0]] = 1
            else:
                # a->c->b->a
                adj[blk[0]][blk[2]] = 1
                adj[blk[2]][blk[1]] = 1
                adj[blk[1]][blk[0]] = 1

        # Random inter-block arcs
        for bi in range(3):
            for bj in range(3):
                if bi == bj:
                    continue
                for u in block_verts[bi]:
                    for v in block_verts[bj]:
                        if u < v and adj[u][v] == 0 and adj[v][u] == 0:
                            if random.random() < 0.5:
                                adj[u][v] = 1
                            else:
                                adj[v][u] = 1

        # Count t3
        scores = [sum(adj[i]) for i in range(n)]
        t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)

        # Count all cycles
        alpha1 = 0
        comp = {}
        for k in range(3, n+1, 2):
            ck = 0
            for vs, d in find_directed_cycles_dp(adj, n, k):
                ck += d
            comp[k] = ck
            alpha1 += ck

        alpha_vals.append(alpha1)
        comp_vals.append(tuple(sorted(comp.items())))

    print(f"Trials: {len(alpha_vals)}")
    print(f"alpha_1 range: [{min(alpha_vals)}, {max(alpha_vals)}]")
    print(f"alpha_1 distribution (top 20):")
    for val, cnt in Counter(alpha_vals).most_common(20):
        print(f"  alpha_1={val}: {cnt}")

    if min(alpha_vals) > 10:
        print(f"\n*** 3 disjoint 3-cycles at n=9 => alpha_1 >= {min(alpha_vals)} > 10 ***")
        print("*** So t3=3 with every-vertex-in-3-cycle at n=9 is INCOMPATIBLE with H=21 ***")

    # Also check: can t3 STAY at 3?
    t3_dist = Counter()
    for trial in range(len(alpha_vals)):
        # Recompute - actually let me just check t3 in the above
        pass

    # Redo with t3 tracking
    print(f"\n=== Also checking t3 in these tournaments ===")
    t3_vals = []
    random.seed(42)
    for trial in range(10000):
        adj = [[0]*n for _ in range(n)]
        for blk in block_verts:
            if random.random() < 0.5:
                adj[blk[0]][blk[1]] = 1
                adj[blk[1]][blk[2]] = 1
                adj[blk[2]][blk[0]] = 1
            else:
                adj[blk[0]][blk[2]] = 1
                adj[blk[2]][blk[1]] = 1
                adj[blk[1]][blk[0]] = 1
        for bi in range(3):
            for bj in range(3):
                if bi == bj:
                    continue
                for u in block_verts[bi]:
                    for v in block_verts[bj]:
                        if u < v and adj[u][v] == 0 and adj[v][u] == 0:
                            if random.random() < 0.5:
                                adj[u][v] = 1
                            else:
                                adj[v][u] = 1
        scores = [sum(adj[i]) for i in range(n)]
        t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)
        t3_vals.append(t3)

    print(f"t3 range: [{min(t3_vals)}, {max(t3_vals)}]")
    print(f"t3 distribution (top 10):")
    for val, cnt in Counter(t3_vals).most_common(10):
        print(f"  t3={val}: {cnt}")

    if min(t3_vals) > 3:
        print(f"  NOTE: inter-block arcs force ADDITIONAL 3-cycles (t3 > 3)")


if __name__ == "__main__":
    main()
