#!/usr/bin/env python3
"""
Comprehensive structural analysis of tournament tiling isomorphism classes.

Reproduces the HTML Tournament Tiling Explorer's logic in Python to extract
and verify structural facts about:
- Isomorphism classes and their sizes
- Blue/black self-loops (flip = inversion preserves class)
- Grid symmetry (sigma-fixed tilings)
- Transpose pairing (T vs T^op)
- Parent-child relationships across n values
- The "perpendicular" structure the user observed

Instance: kind-pasteur-2026-03-05-S10
"""

from itertools import permutations
from collections import defaultdict

def build_tournament_data(n):
    """Build all tournament data for a given n, matching the HTML tool's logic."""
    verts = list(range(n, 0, -1))  # [n, n-1, ..., 1]

    # Tiles: triangular grid positions
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))

    m = len(tiles)  # number of tiles = C(n-1, 2)

    # Tile index lookup
    tile_idx = {(x,y): i for i, (x,y) in enumerate(tiles)}

    # Transpose map: (x,y) -> (n-y+1, n-x+1)
    trans_map = [tile_idx[(n-y+1, n-x+1)] for (x,y) in tiles]

    def is_grid_sym(bits):
        """Check if a tiling is sigma-fixed (grid-symmetric)."""
        for i in range(m):
            if trans_map[i] != i and bits[i] != bits[trans_map[i]]:
                return False
        return True

    def bits_to_adj(bits):
        """Convert tiling bits to adjacency matrix."""
        A = [[0]*n for _ in range(n)]
        # Fixed Hamiltonian path: n -> n-1 -> ... -> 1 (indices 0->1->...->n-1)
        for k in range(n-1):
            A[k][k+1] = 1
        for i, (xL, yL) in enumerate(tiles):
            xi = verts.index(xL)
            yi = verts.index(yL)
            if bits[i] == 0:
                A[xi][yi] = 1
            else:
                A[yi][xi] = 1
        return tuple(tuple(row) for row in A)

    def canonicalize(A):
        """Canonical form of adjacency matrix under vertex permutation."""
        best = None
        for p in permutations(range(n)):
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    def adj_scores(A):
        """Out-degree sequence (sorted descending)."""
        scores = [sum(A[i]) for i in range(n)]
        return tuple(sorted(scores, reverse=True))

    def aut_count(A):
        """Count automorphisms of tournament A."""
        count = 0
        for p in permutations(range(n)):
            is_aut = True
            for i in range(n):
                for j in range(n):
                    if A[p[i]][p[j]] != A[i][j]:
                        is_aut = False
                        break
                if not is_aut:
                    break
            if is_aut:
                count += 1
        return count

    def vertex_orbits(A):
        """Count vertex orbits under automorphism group."""
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        def unite(x, y):
            parent[find(x)] = find(y)

        for p in permutations(range(n)):
            is_aut = True
            for i in range(n):
                for j in range(n):
                    if A[p[i]][p[j]] != A[i][j]:
                        is_aut = False
                        break
                if not is_aut:
                    break
            if is_aut:
                for i in range(n):
                    unite(i, p[i])

        return len(set(find(i) for i in range(n)))

    def is_self_converse(A):
        """Check if T is isomorphic to T^op."""
        A_op = tuple(tuple(A[j][i] for j in range(n)) for i in range(n))
        return canonicalize(A) == canonicalize(A_op)

    # Build all tilings
    tilings = []
    for mask in range(1 << m):
        bits = tuple((mask >> k) & 1 for k in range(m))
        adj = bits_to_adj(bits)
        canonical = canonicalize(adj)
        scores = adj_scores(adj)
        grid_sym = is_grid_sym(bits)
        flip_mask = mask ^ ((1 << m) - 1)  # invert all bits
        trans_bits = [0] * m
        for i in range(m):
            trans_bits[trans_map[i]] = bits[i]
        trans_mask = sum(b << k for k, b in enumerate(trans_bits))

        tilings.append({
            'bits': bits,
            'adj': adj,
            'canonical': canonical,
            'scores': scores,
            'mask': mask,
            'grid_sym': grid_sym,
            'flip_mask': flip_mask,
            'trans_mask': trans_mask,
        })

    # Group by isomorphism class
    groups = defaultdict(list)
    for t in tilings:
        groups[t['canonical']].append(t)

    canon_sigs = sorted(groups.keys())
    for ci, sig in enumerate(canon_sigs):
        for t in groups[sig]:
            t['class_index'] = ci

    # Map masks to tilings
    mask_to_tiling = {t['mask']: t for t in tilings}

    # Transpose class targets
    for t in tilings:
        t['trans_class'] = mask_to_tiling[t['trans_mask']]['class_index']

    class_transpose_target = []
    for ci, sig in enumerate(canon_sigs):
        targets = set(t['trans_class'] for t in groups[sig])
        class_transpose_target.append(targets)

    return {
        'n': n,
        'm': m,
        'tiles': tiles,
        'trans_map': trans_map,
        'tilings': tilings,
        'groups': groups,
        'canon_sigs': canon_sigs,
        'class_transpose_target': class_transpose_target,
        'mask_to_tiling': mask_to_tiling,
    }


def analyze_structure(data):
    """Analyze the full structure for a given n."""
    n = data['n']
    m = data['m']
    groups = data['groups']
    canon_sigs = data['canon_sigs']
    class_transpose_target = data['class_transpose_target']
    mask_to_tiling = data['mask_to_tiling']
    num_classes = len(canon_sigs)

    print(f"\n{'='*70}")
    print(f"n = {n}:  m = {m} tiles,  2^m = {1<<m} tilings,  {num_classes} isomorphism classes")
    print(f"{'='*70}")

    # Class details
    print(f"\n--- Class details ---")
    for ci, sig in enumerate(canon_sigs):
        members = groups[sig]
        rep = members[0]
        size = len(members)
        scores = rep['scores']
        grid_sym_count = sum(1 for t in members if t['grid_sym'])

        # Flip analysis: where does each member's flip land?
        flip_classes = set()
        for t in members:
            flip_t = mask_to_tiling[t['flip_mask']]
            flip_classes.add(flip_t['class_index'])

        # Transpose analysis: where does each member's transpose land?
        trans_classes = class_transpose_target[ci]

        # Self-converse check
        sc = len(trans_classes) == 1 and ci in trans_classes

        # Automorphism count
        aut = 1
        for p in permutations(range(n)):
            is_aut = True
            A = rep['adj']
            for i in range(n):
                for j in range(n):
                    if A[p[i]][p[j]] != A[i][j]:
                        is_aut = False
                        break
                if not is_aut:
                    break
            if is_aut:
                aut += 1
        aut -= 1  # we counted the identity separately
        # Actually redo properly
        A = rep['adj']
        aut = sum(1 for p in permutations(range(n))
                  if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)))

        # Vertex orbits
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        def unite(x, y):
            parent[find(x)] = find(y)
        for p in permutations(range(n)):
            is_aut = all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n))
            if is_aut:
                for i in range(n):
                    unite(i, p[i])
        v_orbits = len(set(find(i) for i in range(n)))

        # Classify the flip type
        flip_self = ci in flip_classes
        flip_type = ""
        if flip_self:
            # Check blue vs black: are both endpoints grid-symmetric?
            has_blue_self = False
            has_black_self = False
            for t in members:
                flip_t = mask_to_tiling[t['flip_mask']]
                if flip_t['class_index'] == ci:
                    if t['grid_sym'] and flip_t['grid_sym']:
                        has_blue_self = True
                    else:
                        has_black_self = True
            if has_blue_self and has_black_self:
                flip_type = "BLUE+BLACK self"
            elif has_blue_self:
                flip_type = "BLUE self"
            elif has_black_self:
                flip_type = "BLACK self"

        trans_type = "self-converse" if sc else f"paired with {trans_classes - {ci}}"

        print(f"  #{ci}: size={size}, scores={scores}, |Aut|={aut}, "
              f"v-orbits={v_orbits}, grid-sym={grid_sym_count}/{size}, "
              f"flip->{flip_classes} {'(' + flip_type + ')' if flip_type else ''}, "
              f"trans: {trans_type}")

    # Flip pair analysis (blue/black lines between classes)
    print(f"\n--- Flip pair analysis ---")
    flip_pairs = defaultdict(lambda: {'blue': 0, 'black': 0})
    for t in data['tilings']:
        flip_t = mask_to_tiling[t['flip_mask']]
        ci_a = t['class_index']
        ci_b = flip_t['class_index']
        key = (min(ci_a, ci_b), max(ci_a, ci_b))
        if t['grid_sym'] and flip_t['grid_sym']:
            flip_pairs[key]['blue'] += 1
        else:
            flip_pairs[key]['black'] += 1

    # Each pair is counted twice (once from each direction), except self-loops
    for key in sorted(flip_pairs.keys()):
        blue = flip_pairs[key]['blue']
        black = flip_pairs[key]['black']
        selfloop = " (SELF)" if key[0] == key[1] else ""
        if key[0] == key[1]:
            blue //= 2
            black //= 2
        else:
            blue //= 2
            black //= 2
        if blue > 0:
            print(f"  BLUE  #{key[0]} <-> #{key[1]}: {blue} pairs{selfloop}")
        if black > 0:
            print(f"  BLACK #{key[0]} <-> #{key[1]}: {black} pairs{selfloop}")

    # Transpose pair analysis (red lines)
    print(f"\n--- Transpose pair analysis ---")
    seen = set()
    for ci in range(num_classes):
        targets = class_transpose_target[ci]
        for tgt in targets:
            if tgt != ci:
                pair = (min(ci, tgt), max(ci, tgt))
                if pair not in seen:
                    seen.add(pair)
                    print(f"  RED #{pair[0]} <-> #{pair[1]}")
            else:
                if ci not in seen:
                    seen.add(ci)
                    print(f"  Self-converse: #{ci}")

    # Special classes
    print(f"\n--- Special classes ---")
    # Class 0: transitive (all tiles off)
    t_trans = mask_to_tiling[0]
    print(f"  Transitive tournament: class #{t_trans['class_index']}, "
          f"size {len(groups[canon_sigs[t_trans['class_index']]])}")

    # Full tiling (all tiles on)
    full_mask = (1 << m) - 1
    t_full = mask_to_tiling[full_mask]
    full_ci = t_full['class_index']
    full_size = len(groups[canon_sigs[full_ci]])
    print(f"  Full tiling: class #{full_ci}, size {full_size}")
    print(f"    2^(n-2)+1 = {2**(n-2)+1} (user conjecture for full class size)")

    # H(T) computation for representatives
    print(f"\n--- Hamiltonian path counts ---")
    for ci, sig in enumerate(canon_sigs):
        rep = groups[sig][0]
        A = rep['adj']
        h = count_ham_paths(A, n)
        print(f"  #{ci}: H(T) = {h}, scores = {rep['scores']}")

    return data


def count_ham_paths(A, n):
    """Count Hamiltonian paths by brute force."""
    count = 0
    for p in permutations(range(n)):
        valid = True
        for k in range(n-1):
            if A[p[k]][p[k+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def analyze_parent_child(data_n, data_np):
    """Analyze parent-child relationships between n and n+1."""
    n = data_n['n']
    np = data_np['n']
    print(f"\n{'='*70}")
    print(f"Parent-child: n={n} -> n={np}")
    print(f"{'='*70}")

    # For each class at n+1, remove the top vertex (index 0, label np)
    # and find which class at n it belongs to
    for ci, sig in enumerate(data_np['canon_sigs']):
        rep = data_np['groups'][sig][0]
        A = rep['adj']
        # Sub-tournament on vertices 1..np-1 (removing vertex 0 = label np)
        sub_A = tuple(tuple(A[i][j] for j in range(1, np)) for i in range(1, np))
        sub_can = None
        best = None
        for p in permutations(range(n)):
            s = tuple(sub_A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        # Find which class at n this matches
        parent_ci = -1
        for pci, psig in enumerate(data_n['canon_sigs']):
            if psig == best:
                parent_ci = pci
                break

        print(f"  n={np} class #{ci} (size={len(data_np['groups'][sig])}) "
              f"-> parent class #{parent_ci} at n={n}")


if __name__ == "__main__":
    results = {}
    for n in range(3, 7):
        print(f"\nBuilding data for n={n}...")
        data = build_tournament_data(n)
        results[n] = analyze_structure(data)

    # Parent-child analysis
    for n in range(3, 6):
        if n in results and n+1 in results:
            analyze_parent_child(results[n], results[n+1])

    print(f"\n{'='*70}")
    print("SUMMARY OF STRUCTURAL FACTS")
    print(f"{'='*70}")
    for n in range(3, 7):
        d = results[n]
        nc = len(d['canon_sigs'])
        m = d['m']
        print(f"  n={n}: m={m}, 2^m={1<<m}, classes={nc}")
