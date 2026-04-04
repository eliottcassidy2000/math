"""
explicit_metagraph.py — opus-2026-04-04-S7

Build the EXPLICIT metagraph G_n from the multilinear polynomial.

The metagraph has:
  Vertices = isomorphism classes of n-tournaments
  Edges = classes connected by a single arc flip

KEY IDEA: Use the multilinear polynomial H(t) to compute EVERYTHING
about the metagraph, including:
  1. H-value per class (from any representative tiling)
  2. Class adjacency (from tile flips between representatives)
  3. Self-complementary structure (from the complement tiling)
  4. The n→n+2 transition map between metagraphs

NEW APPROACH: Instead of nauty, use a fast canonical form based on
the tournament's INVARIANT VECTOR (score sequence + cycle counts + H-value).
For small n, this is sufficient to distinguish all classes.
"""
import numpy as np
from collections import defaultdict, Counter
from itertools import permutations

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_adj(n, tiles, bits):
    """Convert tiling to adjacency matrix."""
    T = np.zeros((n, n), dtype=np.int8)
    for k in range(n-1):
        T[k+1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1
        if bits & (1 << i):
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def count_hp(T, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    total = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp.get((mask, v), 0)
            if cnt == 0: continue
            if mask == full: total += cnt; continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return total

def canonical_form(T, n):
    """Compute canonical form of tournament T.
    Returns sorted adjacency encoding that's invariant under vertex relabeling.
    Uses brute-force over S_n (feasible for n ≤ 8)."""
    best = None
    for perm in permutations(range(n)):
        # Apply permutation: new T[perm[i]][perm[j]] = T[i][j]
        encoding = []
        for i in range(n):
            for j in range(i+1, n):
                encoding.append(T[perm[i]][perm[j]])
        enc_tuple = tuple(encoding)
        if best is None or enc_tuple < best:
            best = enc_tuple
    return best

def fast_invariant(T, n):
    """Fast invariant for pre-grouping (not always canonical but fast)."""
    scores = tuple(sorted(int(T[i].sum()) for i in range(n)))
    # Count 3-cycles
    c3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (T[a][b] and T[b][c] and T[c][a]) or \
                   (T[a][c] and T[c][b] and T[b][a]):
                    c3 += 1
    return (scores, c3)

def build_metagraph(n):
    """Build the complete metagraph at n."""
    tiles = make_tiles(n)
    m = len(tiles)

    print(f"Building metagraph at n={n}, m={m} tiles, 2^m={1<<m} tilings")

    # Phase 1: Compute H and fast invariant for all tilings
    tilings = {}
    for bits in range(1 << m):
        T = tiling_to_adj(n, tiles, bits)
        h = count_hp(T, n)
        inv = fast_invariant(T, n)
        tilings[bits] = {'T': T, 'H': h, 'inv': inv}

    # Phase 2: Group by fast invariant, then canonicalize within groups
    by_inv = defaultdict(list)
    for bits, data in tilings.items():
        by_inv[data['inv']].append(bits)

    print(f"  Fast invariant groups: {len(by_inv)}")

    # Phase 3: Canonical form within each group
    class_map = {}  # bits → class_id
    class_data = {}  # class_id → info
    class_id = 0

    for inv, members in sorted(by_inv.items()):
        # Within each fast-invariant group, canonicalize
        sub_classes = defaultdict(list)
        for bits in members:
            T = tilings[bits]['T']
            canon = canonical_form(T, n)
            sub_classes[canon].append(bits)

        for canon, members_in_class in sub_classes.items():
            h_val = tilings[members_in_class[0]]['H']
            scores = inv[0]
            c3 = inv[1]

            # Check self-complementary: does complement have same canonical form?
            comp_bits = ((1 << m) - 1) ^ members_in_class[0]
            comp_canon = canonical_form(tilings[comp_bits]['T'], n)
            is_sc = (comp_canon == canon)

            class_data[class_id] = {
                'canon': canon,
                'H': h_val,
                'scores': scores,
                'c3': c3,
                'size': len(members_in_class),
                'is_sc': is_sc,
                'members': members_in_class,
            }
            for bits in members_in_class:
                class_map[bits] = class_id
            class_id += 1

    num_classes = class_id
    print(f"  Isomorphism classes: {num_classes}")

    # Phase 4: Build adjacency matrix of the metagraph
    # Two classes are adjacent iff some representative tilings differ by 1 tile
    adj = defaultdict(set)
    edge_count = defaultdict(int)

    for bits in range(1 << m):
        c1 = class_map[bits]
        for tile in range(m):
            flipped = bits ^ (1 << tile)
            c2 = class_map[flipped]
            if c1 != c2:
                adj[c1].add(c2)
                edge_pair = (min(c1,c2), max(c1,c2))
                edge_count[edge_pair] += 1

    num_edges = len(edge_count)
    print(f"  Metagraph edges: {num_edges}")

    # Phase 5: Spectral analysis
    A = np.zeros((num_classes, num_classes))
    for (c1, c2), count in edge_count.items():
        A[c1][c2] = 1
        A[c2][c1] = 1

    eigvals = sorted(np.linalg.eigvalsh(A), reverse=True)

    print(f"\n  Metagraph adjacency spectrum:")
    print(f"    Eigenvalues: {[round(e,3) for e in eigvals]}")

    # H-value distribution
    h_values = sorted(set(d['H'] for d in class_data.values()))
    print(f"\n  H values: {h_values}")
    print(f"  H range: [{min(h_values)}, {max(h_values)}]")

    # SC/NS classification
    sc_count = sum(1 for d in class_data.values() if d['is_sc'])
    ns_count = num_classes - sc_count
    print(f"\n  Self-complementary classes: {sc_count}")
    print(f"  Non-self-complementary classes: {ns_count}")

    # Degree distribution
    degrees = [len(adj[c]) for c in range(num_classes)]
    print(f"\n  Degree distribution: min={min(degrees)}, max={max(degrees)}, "
          f"mean={np.mean(degrees):.1f}")

    # THE KEY: Express the metagraph in terms of the multilinear structure
    # Each class has a representative tiling. The multilinear polynomial
    # evaluated at that tiling gives H.
    # The CLASS-LEVEL adjacency depends on which tile flips change the class.

    # For each class, count "self-loops" (tile flips preserving the class)
    # vs "cross-edges" (changing the class)
    print(f"\n  Per-class structure:")
    for cid in range(min(num_classes, 15)):
        d = class_data[cid]
        self_loops = 0
        cross = 0
        for bits in d['members'][:1]:  # Just check one representative
            for tile in range(m):
                flipped = bits ^ (1 << tile)
                if class_map[flipped] == cid:
                    self_loops += 1
                else:
                    cross += 1
        print(f"    Class {cid}: H={d['H']}, size={d['size']}, SC={d['is_sc']}, "
              f"scores={d['scores']}, c3={d['c3']}, "
              f"self-loops={self_loops}/{m}, cross={cross}/{m}")

    return {
        'n': n, 'm': m, 'num_classes': num_classes, 'num_edges': num_edges,
        'class_data': class_data, 'class_map': class_map,
        'adj': dict(adj), 'eigvals': eigvals, 'A': A,
    }

def n_to_n2_transition(G_n, G_np2, n):
    """Compute the transition map from G_n classes to G_{n+2} classes.
    Each n-class extends to multiple (n+2)-classes via boundary configurations.
    """
    print(f"\n{'='*60}")
    print(f"n→n+2 TRANSITION: {n} → {n+2}")
    print(f"{'='*60}")

    tiles_n = make_tiles(n)
    tiles_np2 = make_tiles(n+2)
    m_n = len(tiles_n)
    m_np2 = len(tiles_np2)

    # Tile mapping: n-tournament tile (x,y) → (n+2)-tournament tile (x+1, y+1)
    tile_idx_np2 = {t: i for i, t in enumerate(tiles_np2)}
    inner_map = {}
    for i, (x, y) in enumerate(tiles_n):
        inner_map[i] = tile_idx_np2[(x+1, y+1)]

    boundary_idx = [i for i in range(m_np2) if i not in set(inner_map.values())]

    # For each n-class, find which (n+2)-classes it maps to
    transition = defaultdict(lambda: defaultdict(int))  # n_class → {n2_class: count}

    for n_bits in range(1 << m_n):
        n_class = G_n['class_map'][n_bits]

        # Try a sample of boundary configs
        for bnd_sample in range(min(64, 1 << len(boundary_idx))):
            np2_bits = 0
            for i_n, i_np2 in inner_map.items():
                if n_bits & (1 << i_n):
                    np2_bits |= (1 << i_np2)
            for j, b_idx in enumerate(boundary_idx):
                if bnd_sample & (1 << j):
                    np2_bits |= (1 << b_idx)

            if np2_bits in G_np2['class_map']:
                np2_class = G_np2['class_map'][np2_bits]
                transition[n_class][np2_class] += 1

    print(f"  Classes at n={n}: {G_n['num_classes']}")
    print(f"  Classes at n+2={n+2}: {G_np2['num_classes']}")

    # Transition summary
    for n_class in sorted(transition.keys()):
        targets = transition[n_class]
        n_data = G_n['class_data'][n_class]
        print(f"  Class {n_class} (H={n_data['H']}, SC={n_data['is_sc']}) → "
              f"{len(targets)} target classes, "
              f"H range [{min(G_np2['class_data'][c]['H'] for c in targets)}, "
              f"{max(G_np2['class_data'][c]['H'] for c in targets)}]")

if __name__ == "__main__":
    G5 = build_metagraph(5)
    print()
    G6 = build_metagraph(6)

    n_to_n2_transition(G5, G6, 4)  # 4→6
