"""
metagraph_structure.py — opus-2026-04-04-S7

Deep structural analysis of the metagraph using multilinear insights.

KEY QUESTIONS:
1. Can the metagraph be expressed as a quotient of the hypercube?
2. What is the relationship between H-values and metagraph distance?
3. Does the adjacency matrix decompose along the perpendicular axes?
4. Can we find a TRANSFER MATRIX representation?
"""
import numpy as np
from collections import defaultdict
from explicit_metagraph import build_metagraph, make_tiles, tiling_to_adj, count_hp

def h_gradient_analysis(G):
    """For each metagraph edge, compute ΔH = |H(C₁) - H(C₂)|.
    This reveals the H-landscape structure."""
    n = G['n']
    class_data = G['class_data']
    num_classes = G['num_classes']

    print(f"\n{'='*60}")
    print(f"H-GRADIENT ANALYSIS n={n}")
    print(f"{'='*60}")

    # Build weighted adjacency with ΔH
    delta_h_dist = defaultdict(int)
    for c1 in range(num_classes):
        for c2 in G['adj'].get(c1, set()):
            if c2 > c1:
                dh = abs(class_data[c1]['H'] - class_data[c2]['H'])
                delta_h_dist[dh] += 1

    print(f"ΔH distribution across edges:")
    for dh in sorted(delta_h_dist.keys()):
        print(f"  ΔH={dh}: {delta_h_dist[dh]} edges")

    # H-levels and connectivity within/between levels
    h_levels = defaultdict(list)
    for cid, data in class_data.items():
        h_levels[data['H']].append(cid)

    print(f"\nH-level structure:")
    for h in sorted(h_levels.keys()):
        classes = h_levels[h]
        # Count edges within this level
        intra = 0
        inter_up = 0
        inter_down = 0
        for c in classes:
            for c2 in G['adj'].get(c, set()):
                h2 = class_data[c2]['H']
                if h2 == h and c2 > c: intra += 1
                elif h2 > h: inter_up += 1
                elif h2 < h: inter_down += 1
        print(f"  H={h}: {len(classes)} classes, "
              f"intra={intra}, up={inter_up}, down={inter_down}")

def metagraph_distance_analysis(G):
    """Compute shortest path distances in the metagraph and relate to H."""
    n = G['n']
    num_classes = G['num_classes']
    class_data = G['class_data']
    adj = G['adj']

    print(f"\n{'='*60}")
    print(f"METAGRAPH DISTANCE ANALYSIS n={n}")
    print(f"{'='*60}")

    # BFS from each class
    dist = np.full((num_classes, num_classes), -1, dtype=int)
    for start in range(num_classes):
        dist[start][start] = 0
        queue = [start]
        front = 0
        while front < len(queue):
            c = queue[front]; front += 1
            for c2 in adj.get(c, set()):
                if dist[start][c2] == -1:
                    dist[start][c2] = dist[start][c] + 1
                    queue.append(c2)

    diameter = dist.max()
    print(f"  Diameter: {diameter}")
    print(f"  Distance distribution:")
    for d in range(diameter + 1):
        count = (dist == d).sum() // 2  # symmetric
        print(f"    d={d}: {count} pairs")

    # Correlation between H-distance and graph distance
    print(f"\n  H-distance vs graph distance:")
    for d in range(1, diameter + 1):
        h_diffs = []
        for i in range(num_classes):
            for j in range(i+1, num_classes):
                if dist[i][j] == d:
                    h_diffs.append(abs(class_data[i]['H'] - class_data[j]['H']))
        if h_diffs:
            print(f"    d={d}: mean |ΔH|={np.mean(h_diffs):.1f}, "
                  f"max |ΔH|={max(h_diffs)}, min |ΔH|={min(h_diffs)}")

    # Distance from transitive (class 0, H=1) to all others
    print(f"\n  Distance from transitive (H=1):")
    for cid in range(num_classes):
        d = dist[0][cid]
        print(f"    Class {cid} (H={class_data[cid]['H']}): d={d}")

def tile_flip_decomposition(G):
    """Decompose metagraph edges by TILE TYPE.
    Each edge comes from flipping a specific tile (x,y).
    The tile's skip determines the "energy" of the flip.
    """
    n = G['n']
    m = G['m']
    tiles = make_tiles(n)
    class_data = G['class_data']
    class_map = G['class_map']

    print(f"\n{'='*60}")
    print(f"TILE-FLIP DECOMPOSITION n={n}")
    print(f"{'='*60}")

    # For each tile, count how many cross-class edges it creates
    tile_edges = defaultdict(set)  # tile_idx → set of (c1, c2) pairs

    for bits in range(1 << m):
        c1 = class_map[bits]
        for tile_idx in range(m):
            flipped = bits ^ (1 << tile_idx)
            c2 = class_map[flipped]
            if c1 != c2:
                edge = (min(c1,c2), max(c1,c2))
                tile_edges[tile_idx].add(edge)

    # Group by skip
    skip_edges = defaultdict(set)
    for tile_idx in range(m):
        skip = tiles[tile_idx][0] - tiles[tile_idx][1]
        skip_edges[skip] |= tile_edges[tile_idx]

    print(f"Edges by skip value:")
    all_edges = set()
    for skip in sorted(skip_edges.keys()):
        edges = skip_edges[skip]
        all_edges |= edges
        print(f"  skip={skip}: {len(edges)} distinct edges")

    total_edges = G['num_edges']
    print(f"\nTotal unique edges: {total_edges}")
    print(f"Union of all skip-edges: {len(all_edges)}")

    # KEY: which skip values COVER all edges?
    # Can we reach all metagraph edges using only skip-2 flips?
    for skip in sorted(skip_edges.keys()):
        coverage = len(skip_edges[skip]) / total_edges * 100
        print(f"  skip={skip} coverage: {coverage:.0f}% ({len(skip_edges[skip])}/{total_edges})")

    # Cumulative coverage
    cumulative = set()
    for skip in sorted(skip_edges.keys()):
        cumulative |= skip_edges[skip]
        coverage = len(cumulative) / total_edges * 100
        print(f"  skip≤{skip} cumulative: {coverage:.0f}%")

    # The "WIGGLY LAYER" decomposition: skip-s edges form a subgraph
    # What are the spectral properties of each layer?
    for skip in sorted(skip_edges.keys()):
        nc = G['num_classes']
        A_skip = np.zeros((nc, nc))
        for c1, c2 in skip_edges[skip]:
            A_skip[c1][c2] = 1
            A_skip[c2][c1] = 1
        eigvals = sorted(np.linalg.eigvalsh(A_skip), reverse=True)
        nonzero = sum(1 for e in eigvals if abs(e) > 0.01)
        print(f"\n  Skip-{skip} layer: {len(skip_edges[skip])} edges, "
              f"rank={nonzero}, λ_max={eigvals[0]:.2f}")

def complement_structure(G):
    """Analyze the complement map on the metagraph.
    Complement: t → 1-t (flip all tiles). Maps class C to class C^c.
    """
    n = G['n']
    m = G['m']
    tiles = make_tiles(n)
    class_data = G['class_data']
    class_map = G['class_map']
    num_classes = G['num_classes']

    print(f"\n{'='*60}")
    print(f"COMPLEMENT STRUCTURE n={n}")
    print(f"{'='*60}")

    # For each class, find its complement class
    comp_map = {}
    for cid in range(num_classes):
        rep = class_data[cid]['members'][0]
        comp = ((1 << m) - 1) ^ rep
        comp_class = class_map[comp]
        comp_map[cid] = comp_class

    # Fixed points = self-complementary classes
    fixed = [c for c in range(num_classes) if comp_map[c] == c]
    orbits_2 = set()
    for c in range(num_classes):
        if comp_map[c] != c:
            pair = (min(c, comp_map[c]), max(c, comp_map[c]))
            orbits_2.add(pair)

    print(f"  SC (fixed) classes: {len(fixed)} — {fixed}")
    print(f"  Complement pairs: {len(orbits_2)}")
    print(f"  Total: {len(fixed) + 2*len(orbits_2)} = {num_classes}")

    # H-values under complement: H(T) = H(T^c) for SC, but H(T) ≠ H(T^c) in general
    # Actually in the TILING model, complement flips all tile bits.
    # The TOURNAMENT complement T^op reverses all arcs, which is NOT the same as
    # the tiling complement (which reverses tile arcs but keeps base path fixed).
    # The tiling complement changes only non-base-path arcs.
    # So: tiling complement ≠ tournament complement.

    # Let's check: does the tiling complement preserve H?
    h_preserved = 0
    h_not_preserved = 0
    for cid in range(num_classes):
        cc = comp_map[cid]
        if class_data[cid]['H'] == class_data[cc]['H']:
            h_preserved += 1
        else:
            h_not_preserved += 1

    print(f"\n  H preserved by tiling complement: {h_preserved}/{num_classes}")
    print(f"  H NOT preserved: {h_not_preserved}/{num_classes}")

    # Show the complement pairs with different H
    for c in range(num_classes):
        cc = comp_map[c]
        if cc > c:  # Only show each pair once
            h1, h2 = class_data[c]['H'], class_data[cc]['H']
            if h1 != h2:
                print(f"    Class {c} (H={h1}) ↔ Class {cc} (H={h2}), ΔH={abs(h1-h2)}")

if __name__ == "__main__":
    G5 = build_metagraph(5)
    h_gradient_analysis(G5)
    metagraph_distance_analysis(G5)
    tile_flip_decomposition(G5)
    complement_structure(G5)
