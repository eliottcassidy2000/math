"""
Deep geometric analysis of tilings on the triangular grid.

Focus: What patterns in the tiling predict the isomorphism class?

Key ideas to test:
1. "Row profiles" and "column profiles" - project tiling onto rows/columns
2. Diagonal symmetry (σ-fixed tilings = self-converse tournaments)
3. Local neighborhoods on the grid - do adjacent tile patterns predict class?
4. "Border" tiles vs "interior" tiles - which matter more?
5. Score sequence from tiling perspective
"""
from itertools import permutations
from collections import defaultdict, Counter
import sys

def build_data(n):
    """Build tiling-to-tournament mapping."""
    # Vertices labeled n, n-1, ..., 1 (the Hamiltonian path order)
    verts = list(range(n, 0, -1))

    # Grid tiles: (r, c) with r>=1, c>=1, r+c <= n-1
    # Tile (r,c) corresponds to the arc between vertices at positions
    # that are r apart and c apart in the base path
    tiles = []
    tile_to_arc = {}
    for r in range(1, n):
        for c in range(1, n):
            if r + c <= n - 1:
                tiles.append((r, c))

    m = len(tiles)
    tile_idx = {t: i for i, t in enumerate(tiles)}

    # Map tiles to arcs: tile (r,c) on strip r+c corresponds to
    # arc between vertices verts[c-1] and verts[c-1+r] ... need to be careful
    # Actually: the standard mapping is:
    # Position i in base path has vertex verts[i]
    # Tile (r,c) corresponds to the arc between verts[c-1] and verts[c-1+r]
    # Wait - let me derive from the definition.

    # From definitions.md: tiling model
    # The base Hamiltonian path visits vertices in order v_1, v_2, ..., v_n
    # A tile at position (r,c) where r+c = k represents the arc between
    # v_c and v_{c+r} = v_k
    # Actually, let me just use the arc-based approach.

    # Simpler: enumerate all arcs not on the base path.
    # Base path: verts[0]->verts[1]->...->verts[n-1]
    # For each pair (i,j) with j > i+1, the arc between verts[i] and verts[j]
    # is a "non-path" arc. Its grid position is (r,c) = (j-i-1, i+1)?
    # Or use the standard: arc from vertex at position i to vertex at position j
    # (j > i+1) maps to tile with r = j-i-1, c = i+1... hmm.

    # Let me just use a direct approach:
    # For n vertices, base path is 0,1,...,n-1
    # Non-path arcs: (i,j) for 0 <= i < j <= n-1, j > i+1
    # These are m = C(n-1, 2) arcs
    # Tile index: enumerate in some order

    arcs = []
    for i in range(n):
        for j in range(i+2, n):
            arcs.append((i, j))

    assert len(arcs) == (n-1)*(n-2)//2

    # Grid position: arc (i,j) -> tile (r,c) where r = j-i-1 (gap-1), c = i+1
    # So strip = r+c = j-i-1+i+1 = j. Hmm, that means strip = j.
    # Actually let me just work with arcs directly and label them by (i,j).

    results = []

    for bits in range(1 << len(arcs)):
        # Build tournament
        adj = [[False]*n for _ in range(n)]
        # Base path arcs: i -> i+1
        for i in range(n-1):
            adj[i][i+1] = True
        # Non-path arcs
        for k, (i, j) in enumerate(arcs):
            if bits & (1 << k):
                adj[j][i] = True  # bit=1 means backward arc j->i
            else:
                adj[i][j] = True  # bit=0 means forward arc i->j

        # Score sequence (sorted)
        scores = sorted([sum(adj[v]) for v in range(n)])

        # Canonical form (brute force)
        canon = None
        for perm in permutations(range(n)):
            mat = tuple(
                tuple(adj[perm[a]][perm[b]] for b in range(n))
                for a in range(n)
            )
            if canon is None or mat < canon:
                canon = mat

        # Tiling properties
        weight = bin(bits).count('1')

        # Row profile: for each row i, count backward arcs from i
        row_back = [0]*n
        for k, (i, j) in enumerate(arcs):
            if bits & (1 << k):
                row_back[j] += 1  # j->i is backward, so row j has a back arc

        # Column profile: for each column j, count backward arcs to j
        col_back = [0]*n
        for k, (i, j) in enumerate(arcs):
            if bits & (1 << k):
                col_back[i] += 1  # j->i, so column i receives a back arc

        results.append({
            'bits': bits,
            'scores': tuple(scores),
            'canon': canon,
            'weight': weight,
            'row_back': tuple(row_back),
            'col_back': tuple(col_back),
        })

    return results, arcs

def analyze_class_geometry(n):
    """What geometric features of tilings distinguish isomorphism classes?"""
    print(f"\n{'='*70}")
    print(f"n = {n}, m = {(n-1)*(n-2)//2}")
    print(f"{'='*70}")

    results, arcs = build_data(n)
    m = len(arcs)

    # Group by canon
    classes = defaultdict(list)
    for r in results:
        classes[r['canon']].append(r)

    print(f"Total tilings: {len(results)}")
    print(f"Isomorphism classes: {len(classes)}")

    # Q1: Can score sequence alone distinguish classes?
    score_to_classes = defaultdict(set)
    for canon, members in classes.items():
        score = members[0]['scores']
        score_to_classes[score].add(canon)

    multi_score = {s: cs for s, cs in score_to_classes.items() if len(cs) > 1}
    print(f"\nScore sequences with multiple classes: {len(multi_score)}")
    for s, cs in sorted(multi_score.items(), key=lambda x: -len(x[1]))[:5]:
        print(f"  {s}: {len(cs)} classes")

    # Q2: Weight distribution per class
    print(f"\nWeight distribution per class:")
    for canon, members in sorted(classes.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
        weights = [r['weight'] for r in members]
        wc = Counter(weights)
        score = members[0]['scores']
        print(f"  size={len(members):4d}, score={score}, weights={dict(sorted(wc.items()))}")

    # Q3: Bit-position entropy - which arc positions vary most within a class?
    print(f"\nBit-position analysis (which arcs are most fixed within a class):")
    # For each class, for each bit position, what fraction is 1?
    bit_variance = [0.0] * m
    for canon, members in classes.items():
        if len(members) < 2:
            continue
        for k in range(m):
            ones = sum(1 for r in members if r['bits'] & (1 << k))
            p = ones / len(members)
            bit_variance[k] += p * (1-p)  # variance contribution

    # Normalize
    nc = sum(1 for c in classes.values() if len(c) >= 2)
    if nc > 0:
        bit_variance = [v/nc for v in bit_variance]

    # Show arcs sorted by variance (low variance = highly predictive)
    arc_var = sorted(enumerate(bit_variance), key=lambda x: x[1])
    print(f"  Most fixed arcs (low variance = same within class):")
    for idx, var in arc_var[:5]:
        i, j = arcs[idx]
        print(f"    arc ({i},{j}), gap={j-i}, variance={var:.4f}")
    print(f"  Most variable arcs (high variance = varies within class):")
    for idx, var in arc_var[-5:]:
        i, j = arcs[idx]
        print(f"    arc ({i},{j}), gap={j-i}, variance={var:.4f}")

    # Q4: Transpose symmetry on the grid
    # σ(r,c) = (c,r) swaps tile at (r,c) with tile at (c,r)
    # In arc terms: arc (i,j) with gap g = j-i maps to...
    # The transpose of a tournament reverses all arcs.
    # Self-converse tournaments: T isomorphic to T^op

    # Check: how many classes are self-converse?
    self_converse = 0
    for canon, members in classes.items():
        # Build T^op canonical
        n_v = len(canon)
        op_adj = tuple(tuple(canon[b][a] for b in range(n_v)) for a in range(n_v))
        op_canon = None
        for perm in permutations(range(n_v)):
            mat = tuple(tuple(op_adj[perm[a]][perm[b]] for b in range(n_v)) for a in range(n_v))
            if op_canon is None or mat < op_canon:
                op_canon = mat
        if op_canon == canon:
            self_converse += 1

    print(f"\nSelf-converse classes: {self_converse}/{len(classes)}")

    # Q5: For each class, look at the "tiling centroid" in {0,1}^m
    # Average bit pattern - which bits tend to be 1?
    print(f"\nClass centroids (avg bit pattern) for largest classes:")
    for canon, members in sorted(classes.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
        centroid = []
        for k in range(m):
            ones = sum(1 for r in members if r['bits'] & (1 << k))
            centroid.append(ones / len(members))
        score = members[0]['scores']
        print(f"  size={len(members)}, score={score}")
        print(f"    centroid: [{', '.join(f'{c:.2f}' for c in centroid)}]")

    # Q6: Hamming distance structure
    # Within each class: what's the diameter and average Hamming distance?
    print(f"\nHamming distance within classes (top 10 by size):")
    for canon, members in sorted(classes.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
        if len(members) == 1:
            continue
        dists = []
        for i in range(min(len(members), 100)):
            for j in range(i+1, min(len(members), 100)):
                d = bin(members[i]['bits'] ^ members[j]['bits']).count('1')
                dists.append(d)
        avg_d = sum(dists)/len(dists)
        max_d = max(dists)
        min_d = min(dists)
        score = members[0]['scores']
        print(f"  size={len(members):4d}, score={score}, "
              f"Hamming: min={min_d}, avg={avg_d:.1f}, max={max_d}")

    return classes, arcs

def analyze_adjacency_on_grid(n):
    """Look at how adjacent tiles on the grid relate within an isomorphism class.

    Two tiles are adjacent on the grid if they share an edge.
    On the triangular grid, tile (r,c) is adjacent to:
    - (r-1,c), (r+1,c), (r,c-1), (r,c+1), (r-1,c+1), (r+1,c-1)
    ... but only within the valid range.
    """
    print(f"\n--- Grid adjacency patterns for n={n} ---")

    results, arcs = build_data(n)
    m = len(arcs)

    # Build grid adjacency
    arc_to_grid = {}
    for k, (i, j) in enumerate(arcs):
        r = j - i - 1  # gap minus 1
        c = i + 1       # position + 1
        arc_to_grid[k] = (r, c)

    # For each pair of adjacent arcs on the grid, compute correlation
    # (do they tend to have the same orientation?)
    print(f"  Arc correlations (same orientation probability):")

    classes = defaultdict(list)
    for r in results:
        classes[r['canon']].append(r['bits'])

    # Compute pairwise bit correlation across ALL tilings
    pair_corr = {}
    for k1 in range(m):
        for k2 in range(k1+1, m):
            same = sum(1 for bits in [r['bits'] for r in results]
                      if ((bits >> k1) & 1) == ((bits >> k2) & 1))
            pair_corr[(k1,k2)] = same / len(results)

    # Group by grid distance
    from math import sqrt
    dist_corr = defaultdict(list)
    for (k1, k2), corr in pair_corr.items():
        r1, c1 = arc_to_grid[k1]
        r2, c2 = arc_to_grid[k2]
        # Use Euclidean distance on grid
        d = abs(r1-r2) + abs(c1-c2)  # Manhattan distance
        dist_corr[d].append(corr)

    print(f"  Avg correlation by Manhattan distance on grid:")
    for d in sorted(dist_corr.keys()):
        vals = dist_corr[d]
        avg = sum(vals)/len(vals)
        print(f"    dist={d}: avg_same_orientation={avg:.4f} ({len(vals)} pairs)")

def find_class_distinguishers(n):
    """For classes with the same score sequence, find what distinguishes them."""
    print(f"\n--- Distinguishing same-score classes for n={n} ---")

    results, arcs = build_data(n)
    m = len(arcs)

    classes = defaultdict(list)
    for r in results:
        classes[r['canon']].append(r)

    # Group classes by score sequence
    score_groups = defaultdict(list)
    for canon, members in classes.items():
        score = members[0]['scores']
        score_groups[score].append((canon, members))

    for score, class_list in sorted(score_groups.items()):
        if len(class_list) < 2:
            continue

        print(f"\n  Score {score}: {len(class_list)} classes")

        for ci, (canon, members) in enumerate(class_list):
            # Compute 3-cycle count
            adj = list(list(row) for row in canon)
            nv = len(adj)
            c3 = 0
            for a in range(nv):
                for b in range(a+1, nv):
                    for c in range(b+1, nv):
                        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
                           (adj[a][c] and adj[c][b] and adj[b][a]):
                            c3 += 1

            # Weight distribution
            weights = Counter(r['weight'] for r in members)

            print(f"    Class {ci}: size={len(members)}, c3={c3}, "
                  f"weights={dict(sorted(weights.items()))}")

if __name__ == '__main__':
    for n in range(3, 7):
        analyze_class_geometry(n)

    print("\n" + "="*70)
    print("ADJACENCY ANALYSIS")
    print("="*70)
    for n in range(3, 7):
        analyze_adjacency_on_grid(n)

    print("\n" + "="*70)
    print("SAME-SCORE CLASS DISTINGUISHERS")
    print("="*70)
    for n in range(3, 7):
        find_class_distinguishers(n)
