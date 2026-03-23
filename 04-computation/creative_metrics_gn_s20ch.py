#!/usr/bin/env python3
"""
creative_metrics_gn_s20ch.py -- kind-pasteur-2026-03-22-S20ch

CREATIVE METRICS FOR THE META-GRAPH G_n.

New metrics that reveal hidden structure:
1. Distance-to-source and distance-to-sink for each class
2. Tiling entropy S(C) = log(H/|Aut|)
3. The bottleneck H (which H-level is widest?)
4. Sink-to-sink distance
5. H-convexity (is every H-interval connected?)
6. Degree-H correlation (assortativity by H)
7. Max flow from source to sinks

Author: kind-pasteur-2026-03-22-S20ch
"""
import sys
import numpy as np
from math import comb, factorial, log2
from collections import defaultdict
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

print("=" * 70)
print("  CREATIVE METRICS FOR G_n")
print("=" * 70)

for n in [5, 6]:
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}\n")

    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    # Build classes
    canon_map = defaultdict(list)
    H_map = {}
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        H = count_hp(A, n)
        H_map[bits] = H
        cf = canonical(A, n)
        canon_map[cf].append(bits)

    classes = []
    cf_to_id = {}
    for cf, members in sorted(canon_map.items(), key=lambda x: H_map[x[1][0]]):
        cid = len(classes)
        cf_to_id[cf] = cid
        classes.append({'id': cid, 'H': H_map[members[0]], 'size': len(members),
                       'members': set(members)})

    bits_to_class = {}
    for c in classes:
        for b in c['members']:
            bits_to_class[b] = c['id']

    N = len(classes)

    # Build adjacency
    adj = np.zeros((N, N), dtype=int)
    for c in classes:
        T = list(c['members'])[0]
        for k in range(m):
            nb = T ^ (1 << k)
            if bits_to_class[nb] != c['id']:
                adj[c['id']][bits_to_class[nb]] = 1
                adj[bits_to_class[nb]][c['id']] = 1

    # Floyd-Warshall distances
    dist = np.full((N, N), N+1)
    np.fill_diagonal(dist, 0)
    for i in range(N):
        for j in range(N):
            if adj[i][j]: dist[i][j] = 1
    for k in range(N):
        for i in range(N):
            for j in range(N):
                if dist[i][k] + dist[k][j] < dist[i][j]:
                    dist[i][j] = dist[i][k] + dist[k][j]

    # Identify source and sinks
    H_min = min(c['H'] for c in classes)
    H_max = max(c['H'] for c in classes)
    sources = [c['id'] for c in classes if c['H'] == H_min]
    sinks = [c['id'] for c in classes if c['H'] == H_max]

    # ---- METRIC 1: Distance to source and sink ----
    print(f"  METRIC 1: DISTANCE TO SOURCE/SINK")
    print(f"  {'Class':>6s} {'H':>5s} {'d_src':>6s} {'d_sink':>7s} {'|Aut|':>6s} {'deg':>4s}")
    for c in classes:
        d_src = min(dist[c['id']][s] for s in sources)
        d_sink = min(dist[c['id']][s] for s in sinks)
        aut = factorial(n) // c['size']
        deg = sum(adj[c['id']])
        if n <= 5 or c['id'] < 10 or c['id'] >= N - 5:
            print(f"  {c['id']:>6d} {c['H']:>5d} {d_src:>6d} {d_sink:>7d} {aut:>6d} {deg:>4d}")

    # ---- METRIC 2: Tiling entropy ----
    print(f"\n  METRIC 2: TILING ENTROPY S(C) = log2(H/|Aut|)")
    entropies = []
    for c in classes:
        aut = factorial(n) // c['size']
        tilings = c['H'] / aut
        if tilings > 0:
            ent = log2(tilings)
        else:
            ent = 0
        entropies.append(ent)
        if n <= 5:
            print(f"  Class {c['id']:>3d} (H={c['H']:>3d}, |Aut|={aut}): S = {ent:.3f}")

    print(f"  Entropy range: [{min(entropies):.3f}, {max(entropies):.3f}]")
    print(f"  Average entropy: {np.mean(entropies):.3f}")

    # ---- METRIC 3: Bottleneck H ----
    H_levels = defaultdict(list)
    for c in classes:
        H_levels[c['H']].append(c['id'])

    widths = {h: len(cids) for h, cids in H_levels.items()}
    bottleneck_H = max(widths.keys(), key=lambda h: widths[h])
    max_width = widths[bottleneck_H]
    print(f"\n  METRIC 3: BOTTLENECK H = {bottleneck_H} (width = {max_width})")
    print(f"  Width prediction C(n-2, (n-2)//2) = {comb(n-2, (n-2)//2)}")
    print(f"  Match: {max_width == comb(n-2, (n-2)//2)}")

    # H at which bottleneck occurs as fraction of H_max
    print(f"  Bottleneck H / H_max = {bottleneck_H} / {H_max} = {bottleneck_H/H_max:.3f}")

    # ---- METRIC 4: Sink-to-sink distance ----
    if len(sinks) >= 2:
        sink_dist = dist[sinks[0]][sinks[1]]
        print(f"\n  METRIC 4: SINK-TO-SINK DISTANCE = {sink_dist}")
        print(f"  Sinks: {[(s, classes[s]['H']) for s in sinks]}")
    else:
        print(f"\n  METRIC 4: Only 1 sink, no distance.")

    # ---- METRIC 5: H-convexity ----
    print(f"\n  METRIC 5: H-CONVEXITY (is every H-interval connected?)")
    all_convex = True
    for h1 in sorted(H_levels.keys()):
        for h2 in sorted(H_levels.keys()):
            if h2 <= h1: continue
            # Induced subgraph on classes with h1 <= H <= h2
            interval_classes = [c['id'] for c in classes if h1 <= c['H'] <= h2]
            if len(interval_classes) <= 1: continue
            # Check connectivity via BFS
            visited = {interval_classes[0]}
            queue = [interval_classes[0]]
            while queue:
                v = queue.pop(0)
                for w in interval_classes:
                    if w not in visited and adj[v][w]:
                        visited.add(w)
                        queue.append(w)
            if len(visited) < len(interval_classes):
                all_convex = False
                if n <= 5:
                    print(f"    NOT convex: interval [{h1}, {h2}] has {len(visited)}/{len(interval_classes)} connected")
                break
        if not all_convex:
            break

    print(f"  H-convex: {all_convex}")

    # ---- METRIC 6: Degree-H correlation ----
    degrees = [int(sum(adj[c['id']])) for c in classes]
    Hs = [c['H'] for c in classes]
    if np.std(degrees) > 0 and np.std(Hs) > 0:
        corr = np.corrcoef(degrees, Hs)[0,1]
    else:
        corr = 0
    print(f"\n  METRIC 6: DEGREE-H CORRELATION = {corr:.4f}")
    print(f"  {'Positive' if corr > 0 else 'Negative'}: high-H classes have {'more' if corr > 0 else 'fewer'} neighbors")

    # ---- METRIC 7: Degree assortativity ----
    # Are high-degree nodes connected to high-degree nodes?
    deg_pairs = []
    for i in range(N):
        for j in range(i+1, N):
            if adj[i][j]:
                deg_pairs.append((degrees[i], degrees[j]))
    if deg_pairs:
        d1s = [p[0] for p in deg_pairs]
        d2s = [p[1] for p in deg_pairs]
        if np.std(d1s) > 0 and np.std(d2s) > 0:
            assort = np.corrcoef(d1s, d2s)[0,1]
        else:
            assort = 0
        print(f"\n  METRIC 7: DEGREE ASSORTATIVITY = {assort:.4f}")
        print(f"  {'Assortative' if assort > 0 else 'Disassortative'}: hubs connect to {'hubs' if assort > 0 else 'periphery'}")

    # ---- SUMMARY ----
    print(f"\n  SUMMARY for n={n}:")
    print(f"    V={N}, E={sum(adj.sum(axis=1))//2}, diameter={int(dist.max())}")
    print(f"    Bottleneck: H={bottleneck_H} (width={max_width}), at {bottleneck_H/H_max:.0%} of H_max")
    print(f"    Sink distance: {dist[sinks[0]][sinks[1]] if len(sinks) >= 2 else 'N/A'}")
    print(f"    H-convex: {all_convex}")
    print(f"    Degree-H corr: {corr:.3f}")
    print(f"    Entropy range: [{min(entropies):.2f}, {max(entropies):.2f}]")
