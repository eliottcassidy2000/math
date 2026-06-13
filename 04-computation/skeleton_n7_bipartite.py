#!/usr/bin/env python3
"""
Parse the n=7 skeleton output and check bipartiteness + NSC sidedness.

Uses the already-computed data from skeleton_n7_fast.py output.

kind-pasteur-2026-03-06-S25h
"""
import re
from collections import defaultdict, deque

# Parse the output file
output_file = r"C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-math\tasks\bit8a5ma1.output"

with open(output_file, 'r') as f:
    content = f.read()

# Extract blue edges
edge_pattern = re.compile(r'(\d+)\(scores=.*?\)\s*<-\[(\d+)\]->\s*(\d+)\(scores=')
edges = []
for m in edge_pattern.finditer(content):
    i, w, j = int(m.group(1)), int(m.group(2)), int(m.group(3))
    edges.append((i, j, w))

print(f"Parsed {len(edges)} blue edges")

# Build adjacency
adj = defaultdict(set)
for i, j, w in edges:
    adj[i].add(j)
    adj[j].add(i)

vertices = set()
for i, j, w in edges:
    vertices.add(i)
    vertices.add(j)

print(f"Skeleton: {len(vertices)} vertices, {len(edges)} edges")

# Check bipartiteness
color = {}
is_bip = True
for start in vertices:
    if start in color:
        continue
    color[start] = 0
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if v not in color:
                color[v] = 1 - color[u]
                queue.append(v)
            elif color[v] == color[u]:
                is_bip = False
                print(f"  ODD CYCLE detected: edge {u}-{v}, both color {color[u]}")

print(f"\nBIPARTITE: {is_bip}")

if is_bip:
    side_A = sorted(v for v in vertices if color[v] == 0)
    side_B = sorted(v for v in vertices if color[v] == 1)
    print(f"Side A: {len(side_A)} classes")
    print(f"Side B: {len(side_B)} classes")

    # Degree distribution per side
    deg_A = sorted(len(adj[v]) for v in side_A)
    deg_B = sorted(len(adj[v]) for v in side_B)
    print(f"Side A degrees: min={min(deg_A)}, max={max(deg_A)}, sum={sum(deg_A)}")
    print(f"Side B degrees: min={min(deg_B)}, max={max(deg_B)}, sum={sum(deg_B)}")
    print(f"  (sum should be equal since each edge has one endpoint on each side)")

    # Parse NSC pair info from output
    nsc_pattern = re.compile(r'Pair \((\d+),(\d+)\) scores=.*?: SC_i=(\{[^}]*\}|set\(\)), SC_p=(\{[^}]*\}|set\(\))')
    nsc_pairs = []
    for m_match in nsc_pattern.finditer(content):
        i, p = int(m_match.group(1)), int(m_match.group(2))
        sc_i_str = m_match.group(3)
        sc_p_str = m_match.group(4)
        sc_i = set(eval(sc_i_str)) if sc_i_str != 'set()' else set()
        sc_p = set(eval(sc_p_str)) if sc_p_str != 'set()' else set()
        nsc_pairs.append((i, p, sc_i, sc_p))

    print(f"\nParsed {len(nsc_pairs)} NSC pairs")

    # For each NSC pair: do their SC targets land on the same side or different sides?
    opposite_sided = 0
    same_sided = 0
    mixed = 0
    no_sc = 0

    for i, partner, sc_i, sc_p in nsc_pairs:
        if not sc_i and not sc_p:
            no_sc += 1
            continue

        # Which sides do the SC targets of class i land on?
        sides_i = set()
        for sc_class in sc_i:
            if sc_class in color:
                sides_i.add(color[sc_class])

        # Which sides do the SC targets of partner land on?
        sides_p = set()
        for sc_class in sc_p:
            if sc_class in color:
                sides_p.add(color[sc_class])

        if not sides_i and not sides_p:
            no_sc += 1
        elif not sides_i or not sides_p:
            # Only one has SC targets
            mixed += 1
        elif sides_i == {0} and sides_p == {1}:
            opposite_sided += 1
        elif sides_i == {1} and sides_p == {0}:
            opposite_sided += 1
        elif sides_i == sides_p and len(sides_i) == 1:
            same_sided += 1
        else:
            mixed += 1

    print(f"\nNSC PAIR SIDEDNESS (by SC flip targets):")
    print(f"  Opposite-sided: {opposite_sided}")
    print(f"  Same-sided: {same_sided}")
    print(f"  Mixed/ambiguous: {mixed}")
    print(f"  No SC targets: {no_sc}")
else:
    # Find shortest odd cycle
    print("\nSkeleton is NOT bipartite")
    # Still interesting — what's the girth?
    # BFS to find shortest cycle
    min_cycle = float('inf')
    for start in list(vertices)[:20]:
        dist = {start: 0}
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    queue.append(v)
                elif dist[v] >= dist[u]:  # back edge
                    cycle_len = dist[u] + dist[v] + 1
                    if cycle_len % 2 == 1 and cycle_len < min_cycle:
                        min_cycle = cycle_len
    if min_cycle < float('inf'):
        print(f"  Shortest odd cycle length: {min_cycle}")

print("\nDONE")
