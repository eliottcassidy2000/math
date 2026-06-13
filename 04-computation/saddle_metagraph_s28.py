#!/usr/bin/env python3
"""
The H-saddle in the metagraph G_n.
opus-2026-04-03-S28

In the tiling hypercube, H has a saddle structure:
  - Maximum H at certain (hyp, per) corners
  - Minimum H at the transitive tournament (0,0)
  - The saddle creates a "ridge" of high-H tilings

In the metagraph G_n (on isomorphism classes):
  - Each class has a single H value
  - The wiggly layer (d=1) connects tilings that differ in one tile
  - The metagraph collapses the tiling hypercube by the S_n action

QUESTION: Does the saddle create a VALLEY or RIDGE structure in G_n?
Are high-H classes connected to each other, or separated by low-H classes?

Also: the apex tile effect suggests that classes with the apex bit "set"
(i.e., where vertex n → vertex 1 is reversed) tend to have higher H.
Does this create a BIPARTITION of the metagraph?
"""

from collections import defaultdict
from itertools import permutations

def tournament_from_bits(n, bits_list):
    T = [[0]*n for _ in range(n)]
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    for i in range(n, 1, -1):
        T[i-1][i-2] = 1
    for idx, (x, y) in enumerate(tiles):
        if idx < len(bits_list):
            if bits_list[idx] == 0:
                T[x-1][y-1] = 1
            else:
                T[y-1][x-1] = 1
    return T

def count_hp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def canonical_form(T, n):
    """Get canonical form of tournament (min over all relabelings)."""
    min_form = None
    for perm in permutations(range(n)):
        form = []
        for i in range(n):
            for j in range(i+1, n):
                form.append(T[perm[i]][perm[j]])
        form = tuple(form)
        if min_form is None or form < min_form:
            min_form = form
    return min_form

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

for n in [5, 6]:
    m = (n-1)*(n-2)//2
    tiles = get_tiles(n)

    print(f"\n{'='*70}")
    print(f"METAGRAPH SADDLE STRUCTURE: n={n}, m={m}")
    print(f"{'='*70}")

    # Compute all tournaments and their iso classes
    class_data = defaultdict(lambda: {'H': 0, 'count': 0, 'apex_dist': []})

    apex_idx = max(range(m), key=lambda i: tiles[i][0] - tiles[i][1])

    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        H = count_hp(T, n)
        canon = canonical_form(T, n)

        class_data[canon]['H'] = H
        class_data[canon]['count'] += 1
        class_data[canon]['apex_dist'].append(b[apex_idx])

    classes = sorted(class_data.keys(), key=lambda c: class_data[c]['H'])

    print(f"\nIso classes: {len(classes)}")
    print(f"\nClass distribution by H:")
    H_to_classes = defaultdict(list)
    for c in classes:
        H_to_classes[class_data[c]['H']].append(c)

    for H_val in sorted(H_to_classes):
        cls = H_to_classes[H_val]
        count = len(cls)
        total_tilings = sum(class_data[c]['count'] for c in cls)
        # Average apex fraction
        apex_fracs = [sum(class_data[c]['apex_dist'])/len(class_data[c]['apex_dist'])
                      for c in cls]
        avg_apex = sum(apex_fracs) / len(apex_fracs) if apex_fracs else 0
        print(f"  H={H_val:>3}: {count:>3} classes, {total_tilings:>6} tilings, "
              f"avg apex frac={avg_apex:.3f}")

    # Build the wiggly graph on iso classes
    # Two classes are connected if some tiling in class A is a wiggly neighbor
    # of some tiling in class B
    print(f"\n--- WIGGLY CONNECTIONS IN THE METAGRAPH ---")

    # Map each tiling to its class
    tiling_to_class = {}
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        canon = canonical_form(T, n)
        tiling_to_class[bits] = canon

    # Find wiggly edges between classes
    class_edges = set()
    for bits in range(2**m):
        c1 = tiling_to_class[bits]
        for j in range(m):
            bits2 = bits ^ (1 << j)
            c2 = tiling_to_class[bits2]
            if c1 != c2:
                edge = tuple(sorted([c1, c2]))
                class_edges.add(edge)

    # H-difference across wiggly edges
    print(f"\nWiggly edges: {len(class_edges)}")

    delta_H_dist = defaultdict(int)
    for c1, c2 in class_edges:
        delta = abs(class_data[c1]['H'] - class_data[c2]['H'])
        delta_H_dist[delta] += 1

    print(f"\n|delta H| distribution across wiggly edges:")
    for delta in sorted(delta_H_dist):
        count = delta_H_dist[delta]
        pct = count / len(class_edges) * 100
        print(f"  |ΔH|={delta:>3}: {count:>4} edges ({pct:.1f}%)")

    # Which edges connect high-H classes to high-H classes?
    max_H = max(class_data[c]['H'] for c in classes)
    high_threshold = max_H * 0.7

    high_high = 0
    high_low = 0
    low_low = 0
    for c1, c2 in class_edges:
        h1, h2 = class_data[c1]['H'], class_data[c2]['H']
        if h1 >= high_threshold and h2 >= high_threshold:
            high_high += 1
        elif h1 >= high_threshold or h2 >= high_threshold:
            high_low += 1
        else:
            low_low += 1

    total = high_high + high_low + low_low
    print(f"\nEdge type distribution (threshold H >= {high_threshold:.0f}):")
    print(f"  High-High: {high_high} ({high_high/total*100:.1f}%)")
    print(f"  High-Low:  {high_low} ({high_low/total*100:.1f}%)")
    print(f"  Low-Low:   {low_low} ({low_low/total*100:.1f}%)")

    # Is the metagraph connected on high-H classes?
    high_classes = [c for c in classes if class_data[c]['H'] >= high_threshold]
    print(f"\nHigh-H classes: {len(high_classes)} (out of {len(classes)})")

    # Build adjacency for high classes
    high_set = set(high_classes)
    high_adj = defaultdict(set)
    for c1, c2 in class_edges:
        if c1 in high_set and c2 in high_set:
            high_adj[c1].add(c2)
            high_adj[c2].add(c1)

    # BFS to find components
    visited = set()
    components = []
    for c in high_classes:
        if c not in visited:
            component = set()
            queue = [c]
            while queue:
                node = queue.pop(0)
                if node in visited:
                    continue
                visited.add(node)
                component.add(node)
                for neighbor in high_adj.get(node, []):
                    if neighbor not in visited:
                        queue.append(neighbor)
            components.append(component)

    print(f"Connected components among high-H classes: {len(components)}")
    for i, comp in enumerate(components):
        H_vals = [class_data[c]['H'] for c in comp]
        print(f"  Component {i}: {len(comp)} classes, H range [{min(H_vals)}, {max(H_vals)}]")

for n in [5, 6]:
    m = (n-1)*(n-2)//2
    print(f"\n--- APEX BIT AND METAGRAPH PARTITION (n={n}) ---")
    # The apex bit divides tilings into two halves.
    # In each iso class, what fraction of tilings have apex=1?

    tiles = get_tiles(n)
    apex_idx = max(range(m), key=lambda i: tiles[i][0] - tiles[i][1])

    class_apex = defaultdict(list)
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        canon = canonical_form(T, n)
        class_apex[canon].append(b[apex_idx])

    # Show apex fraction per class
    print(f"{'H':>5} {'|class|':>8} {'apex_frac':>10} {'type':>10}")
    for canon in sorted(class_apex, key=lambda c: class_data[c]['H']):
        H = class_data[canon]['H']
        size = len(class_apex[canon])
        frac = sum(class_apex[canon]) / size
        # Type: pure-0, pure-1, or mixed
        if frac == 0: t = "apex=0"
        elif frac == 1: t = "apex=1"
        else: t = "mixed"
        print(f"  {H:>3} {size:>8} {frac:>10.3f} {t:>10}")
