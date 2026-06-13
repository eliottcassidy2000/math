#!/usr/bin/env python3
"""
Deep investigation of tournament isomorphism classes and their tiling structure.

Focus areas:
1. How symmetric tilings (grid-symmetric, partially symmetric) group into iso classes
2. "Symmetry kernels" — which structural features persist as n grows
3. The automorphism group's action on tilings within each class
4. Parent-child inheritance of symmetry properties across n values
5. Which iso classes contain grid-symmetric tilings and why

Instance: opus-2026-03-06-S7
"""

from itertools import permutations, combinations
from collections import defaultdict, Counter
import sys

def build_tournament_data(n, verbose=False):
    """Build all tournament data for given n."""
    verts = list(range(n, 0, -1))  # [n, n-1, ..., 1]

    # Tiles: triangular grid
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)
    tile_idx = {(x, y): i for i, (x, y) in enumerate(tiles)}

    # Transpose map: (x,y) -> (n-y+1, n-x+1)
    trans_map = []
    for (x, y) in tiles:
        nx, ny = n - y + 1, n - x + 1
        trans_map.append(tile_idx[(nx, ny)])

    def is_grid_sym(bits):
        for i in range(m):
            if bits[i] != bits[trans_map[i]]:
                return False
        return True

    def bits_to_adj(bits):
        A = [[0]*n for _ in range(n)]
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
        best = None
        for p in permutations(range(n)):
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    def count_ham_paths(A):
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

    def get_automorphisms(A):
        auts = []
        for p in permutations(range(n)):
            if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)):
                auts.append(p)
        return auts

    def get_anti_automorphisms(A):
        """Permutations p with A[p[i]][p[j]] = A[j][i] (= 1-A[i][j] for i!=j)."""
        anti_auts = []
        for p in permutations(range(n)):
            if all(A[p[i]][p[j]] == A[j][i] for i in range(n) for j in range(n)):
                anti_auts.append(p)
        return anti_auts

    # Build all tilings
    tilings = []
    for mask in range(1 << m):
        bits = tuple((mask >> k) & 1 for k in range(m))
        adj = bits_to_adj(bits)
        canonical = canonicalize(adj)
        grid_sym = is_grid_sym(bits)
        flip_mask = mask ^ ((1 << m) - 1)

        trans_bits = [0] * m
        for i in range(m):
            trans_bits[trans_map[i]] = bits[i]
        trans_mask = sum(b << k for k, b in enumerate(trans_bits))

        tilings.append({
            'bits': bits,
            'adj': adj,
            'canonical': canonical,
            'mask': mask,
            'grid_sym': grid_sym,
            'flip_mask': flip_mask,
            'trans_mask': trans_mask,
            'hamming_weight': sum(bits),
        })

    # Group by isomorphism class
    groups = defaultdict(list)
    for t in tilings:
        groups[t['canonical']].append(t)

    canon_sigs = sorted(groups.keys())
    for ci, sig in enumerate(canon_sigs):
        for t in groups[sig]:
            t['class_index'] = ci

    mask_to_tiling = {t['mask']: t for t in tilings}

    # Compute class-level data
    class_data = []
    for ci, sig in enumerate(canon_sigs):
        members = groups[sig]
        rep_adj = members[0]['adj']
        H = count_ham_paths(rep_adj)
        auts = get_automorphisms(rep_adj)
        anti_auts = get_anti_automorphisms(rep_adj)
        scores = tuple(sorted([sum(rep_adj[i]) for i in range(n)], reverse=True))

        grid_sym_members = [t for t in members if t['grid_sym']]

        # Flip targets
        flip_classes = set()
        for t in members:
            flip_t = mask_to_tiling[t['flip_mask']]
            flip_classes.add(flip_t['class_index'])
        is_self_flip = ci in flip_classes

        # Transpose targets
        trans_classes = set()
        for t in members:
            trans_t = mask_to_tiling[t['trans_mask']]
            trans_classes.add(trans_t['class_index'])
        is_self_converse = len(trans_classes) == 1 and ci in trans_classes

        # Hamming weight distribution within class
        hw_dist = Counter(t['hamming_weight'] for t in members)

        class_data.append({
            'ci': ci,
            'size': len(members),
            'H': H,
            'aut_size': len(auts),
            'anti_aut_size': len(anti_auts),
            'scores': scores,
            'grid_sym_count': len(grid_sym_members),
            'is_self_flip': is_self_flip,
            'is_self_converse': is_self_converse,
            'flip_classes': flip_classes,
            'trans_classes': trans_classes,
            'hw_dist': hw_dist,
            'members': members,
            'grid_sym_members': grid_sym_members,
            'auts': auts,
            'anti_auts': anti_auts,
            'adj': rep_adj,
        })

    return {
        'n': n, 'm': m, 'tiles': tiles, 'trans_map': trans_map,
        'tilings': tilings, 'groups': groups, 'canon_sigs': canon_sigs,
        'class_data': class_data, 'mask_to_tiling': mask_to_tiling,
    }


def analyze_symmetry_grouping(data):
    """How do grid-symmetric tilings distribute across isomorphism classes?"""
    n = data['n']
    m = data['m']
    cd = data['class_data']

    print(f"\n{'='*70}")
    print(f"SYMMETRY GROUPING ANALYSIS: n={n}, m={m}")
    print(f"{'='*70}")

    total_gs = sum(c['grid_sym_count'] for c in cd)
    gs_classes = [(c['ci'], c['grid_sym_count'], c['size'], c['H'], c['aut_size'],
                   c['anti_aut_size'], c['is_self_converse'], c['scores'])
                  for c in cd if c['grid_sym_count'] > 0]

    print(f"\nTotal grid-symmetric tilings: {total_gs} (= 2^floor((n-1)^2/4) = 2^{(n-1)**2//4} = {2**((n-1)**2//4)})")
    print(f"Classes containing grid-symmetric tilings: {len(gs_classes)} / {len(cd)}")
    print(f"\nClass | GS tilings | Class size | H(T) | |Aut| | |AntiAut| | SC? | Scores")
    print(f"{'-'*90}")
    for ci, gs, sz, H, aut, aaut, sc, scores in gs_classes:
        ratio = gs / sz if sz > 0 else 0
        print(f"  #{ci:3d} | {gs:10d} | {sz:10d} | {H:5d} | {aut:5d} | {aaut:8d} | {'Y' if sc else 'N':3s} | {scores}")
        print(f"        | GS/size = {ratio:.4f} = {gs}/{sz}")

    # Key insight: GS tilings / class size = |AntiAut| / n! ?
    print(f"\n--- Grid-symmetric fraction analysis ---")
    nfact = 1
    for i in range(1, n+1):
        nfact *= i
    for c in cd:
        if c['grid_sym_count'] > 0:
            gs = c['grid_sym_count']
            sz = c['size']
            aaut = c['anti_aut_size']
            # GS tilings in a class = tilings fixed by sigma = |{t : sigma(t) in same class AND sigma(t) = t}|
            # But sigma maps a tiling to its transpose, which gives T^op
            # So GS tilings correspond to self-converse presentations
            print(f"  #{c['ci']}: GS={gs}, size={sz}, |Aut|={c['aut_size']}, |AntiAut|={aaut}, "
                  f"GS/size={gs/sz:.4f}, |AntiAut|/|Aut|={aaut/c['aut_size'] if c['aut_size'] else 'inf'}")


def analyze_symmetry_kernel(data_list):
    """Track which symmetry properties persist across n values."""
    print(f"\n{'='*70}")
    print(f"SYMMETRY KERNEL ANALYSIS: tracking features across n")
    print(f"{'='*70}")

    for data in data_list:
        n = data['n']
        cd = data['class_data']

        sc_classes = [c for c in cd if c['is_self_converse']]
        sf_classes = [c for c in cd if c['is_self_flip']]
        both = [c for c in cd if c['is_self_converse'] and c['is_self_flip']]
        high_aut = [c for c in cd if c['aut_size'] > 1]
        has_anti = [c for c in cd if c['anti_aut_size'] > 0]

        print(f"\nn={n}: {len(cd)} classes")
        print(f"  Self-converse (T ~ T^op): {len(sc_classes)}")
        print(f"  Self-flip (flip within class): {len(sf_classes)}")
        print(f"  Both SC + SF: {len(both)}")
        print(f"  |Aut| > 1: {len(high_aut)}")
        print(f"  Has anti-automorphisms: {len(has_anti)}")

        # List the SC+SF classes (these are the "maximally symmetric")
        if both:
            print(f"  SC+SF classes:")
            for c in both:
                print(f"    #{c['ci']}: H={c['H']}, |Aut|={c['aut_size']}, |AntiAut|={c['anti_aut_size']}, "
                      f"scores={c['scores']}, GS={c['grid_sym_count']}/{c['size']}")


def analyze_aut_action_on_tilings(data):
    """How does the automorphism group act on tilings within each class?"""
    n = data['n']
    m = data['m']
    tiles = data['tiles']
    cd = data['class_data']
    verts = list(range(n, 0, -1))

    print(f"\n{'='*70}")
    print(f"AUTOMORPHISM ACTION ON TILINGS: n={n}")
    print(f"{'='*70}")

    tile_idx = {(x, y): i for i, (x, y) in enumerate(tiles)}

    for c in cd:
        if c['aut_size'] <= 1:
            continue  # skip trivial automorphism groups

        ci = c['ci']
        auts = c['auts']
        members = c['members']

        print(f"\n--- Class #{ci}: H={c['H']}, |Aut|={c['aut_size']}, size={c['size']}, scores={c['scores']} ---")

        # How does each automorphism permute the tiles?
        # A tile (x, y) under permutation p becomes (p(x), p(y)) — but we need
        # to track both the tile position AND whether the bit flips
        for aut in auts:
            if all(aut[i] == i for i in range(n)):
                continue  # skip identity

            # Map each tile under this automorphism
            aut_label = tuple(verts[aut[verts.index(v)]] if v in verts else v for v in range(1, n+1))
            # Actually: aut permutes vertex indices 0..n-1
            # Tile (x,y) with x,y in verts connects vertex verts.index(x) to verts.index(y)
            # Under aut, this becomes (aut[verts.index(x)], aut[verts.index(y)])
            # = tile connecting verts[aut[verts.index(x)]] to verts[aut[verts.index(y)]]

            fixed_tiles = 0
            flipped_tiles = 0
            for i, (x, y) in enumerate(tiles):
                xi, yi = verts.index(x), verts.index(y)
                new_xi, new_yi = aut[xi], aut[yi]
                new_x, new_y = verts[new_xi], verts[new_yi]
                # The tile connecting new_x, new_y: need max, min
                if new_x > new_y:
                    # same direction
                    if (new_x, new_y) in tile_idx:
                        target = tile_idx[(new_x, new_y)]
                        if target == i:
                            fixed_tiles += 1
                elif new_y > new_x:
                    flipped_tiles += 1

            print(f"  Aut {tuple(aut)}: fixes {fixed_tiles}/{m} tiles, flips orientation on {flipped_tiles}")

        # Orbits of tilings under Aut
        # The class has size = H(T)/|Aut(T)|, meaning |Aut| tilings map to each other
        # Let's verify this by computing actual orbits
        tiling_orbits = []
        seen_masks = set()
        for t in members:
            if t['mask'] in seen_masks:
                continue
            orbit = set()
            # For each automorphism, apply it to this tiling
            for aut in auts:
                new_bits = apply_aut_to_tiling(t['bits'], aut, tiles, verts, tile_idx)
                if new_bits is not None:
                    new_mask = sum(b << k for k, b in enumerate(new_bits))
                    orbit.add(new_mask)
            for om in orbit:
                seen_masks.add(om)
            tiling_orbits.append(orbit)

        orbit_sizes = Counter(len(o) for o in tiling_orbits)
        print(f"  Tiling orbits: {len(tiling_orbits)} orbits, size distribution: {dict(orbit_sizes)}")
        print(f"  Expected orbit size = |Aut| = {c['aut_size']} (Burnside: {c['size']}/{c['aut_size']} = {c['size']/c['aut_size']:.1f} orbits)")

        # Grid-symmetric tilings in this class — are they in their own orbits?
        gs_masks = {t['mask'] for t in c['grid_sym_members']}
        gs_orbits = [o for o in tiling_orbits if o & gs_masks]
        print(f"  Grid-symmetric tilings: {len(gs_masks)}, in {len(gs_orbits)} orbits")


def apply_aut_to_tiling(bits, aut, tiles, verts, tile_idx):
    """Apply automorphism to tiling bits, returning new bits or None if invalid."""
    m = len(tiles)
    new_bits = [0] * m
    for i, (x, y) in enumerate(tiles):
        xi, yi = verts.index(x), verts.index(y)
        new_xi, new_yi = aut[xi], aut[yi]
        new_x, new_y = verts[new_xi], verts[new_yi]

        if new_x > new_y + 1:
            # This is a non-path arc, should be a tile
            if (new_x, new_y) in tile_idx:
                target_tile = tile_idx[(new_x, new_y)]
                new_bits[target_tile] = bits[i]
            else:
                return None
        elif new_y > new_x + 1:
            # Reversed direction
            if (new_y, new_x) in tile_idx:
                target_tile = tile_idx[(new_y, new_x)]
                new_bits[target_tile] = 1 - bits[i]
            else:
                return None
        else:
            # Adjacent vertices (path edge) — this arc is fixed
            # The aut should map path edges to path edges or tiles
            # If new_x = new_y + 1, this is a path edge (always forward)
            # bit doesn't matter for path edges, skip
            pass

    return tuple(new_bits)


def analyze_parent_child_symmetry(data_n, data_np):
    """Track how symmetry properties pass from parent to child classes."""
    n = data_n['n']
    np_val = data_np['n']
    cd_n = data_n['class_data']
    cd_np = data_np['class_data']

    print(f"\n{'='*70}")
    print(f"PARENT-CHILD SYMMETRY INHERITANCE: n={n} -> n={np_val}")
    print(f"{'='*70}")

    # For each class at n+1, find its parent at n (by deleting the highest vertex)
    def canonicalize_adj(A, sz):
        best = None
        for p in permutations(range(sz)):
            s = tuple(A[p[i]][p[j]] for i in range(sz) for j in range(sz))
            if best is None or s < best:
                best = s
        return best

    parent_map = {}
    for c_np in cd_np:
        A = c_np['adj']
        # Delete vertex 0 (highest label)
        sub_A = tuple(tuple(A[i][j] for j in range(1, np_val)) for i in range(1, np_val))
        sub_can = canonicalize_adj(sub_A, n)
        for c_n in cd_n:
            if data_n['canon_sigs'][c_n['ci']] == sub_can:
                parent_map[c_np['ci']] = c_n['ci']
                break

    # Group children by parent
    children_of = defaultdict(list)
    for child_ci, parent_ci in parent_map.items():
        children_of[parent_ci].append(child_ci)

    for parent_ci in sorted(children_of.keys()):
        p = cd_n[parent_ci]
        children = children_of[parent_ci]
        print(f"\n  Parent #{parent_ci} (H={p['H']}, |Aut|={p['aut_size']}, SC={'Y' if p['is_self_converse'] else 'N'}, "
              f"SF={'Y' if p['is_self_flip'] else 'N'}, scores={p['scores']})")

        for child_ci in children:
            ch = cd_np[child_ci]
            # Track which properties are inherited
            inherited = []
            if p['is_self_converse'] and ch['is_self_converse']:
                inherited.append("SC")
            if p['is_self_flip'] and ch['is_self_flip']:
                inherited.append("SF")
            if p['aut_size'] > 1 and ch['aut_size'] > 1:
                inherited.append("|Aut|>1")

            lost = []
            if p['is_self_converse'] and not ch['is_self_converse']:
                lost.append("SC")
            if p['is_self_flip'] and not ch['is_self_flip']:
                lost.append("SF")
            if p['aut_size'] > 1 and ch['aut_size'] <= 1:
                lost.append("|Aut|>1")

            gained = []
            if not p['is_self_converse'] and ch['is_self_converse']:
                gained.append("SC")
            if not p['is_self_flip'] and ch['is_self_flip']:
                gained.append("SF")

            notes = ""
            if inherited:
                notes += f" KEPT:{','.join(inherited)}"
            if lost:
                notes += f" LOST:{','.join(lost)}"
            if gained:
                notes += f" GAINED:{','.join(gained)}"

            print(f"    -> Child #{child_ci} (H={ch['H']}, |Aut|={ch['aut_size']}, "
                  f"SC={'Y' if ch['is_self_converse'] else 'N'}, SF={'Y' if ch['is_self_flip'] else 'N'}, "
                  f"scores={ch['scores']}){notes}")


def analyze_hamming_weight_structure(data):
    """How are tilings distributed by Hamming weight within each class?"""
    n = data['n']
    m = data['m']
    cd = data['class_data']

    print(f"\n{'='*70}")
    print(f"HAMMING WEIGHT STRUCTURE: n={n}")
    print(f"{'='*70}")

    print(f"\nClass | H(T) | Size | HW range | HW median | Symmetric? | Profile")
    print(f"{'-'*90}")

    for c in cd:
        hw = c['hw_dist']
        min_hw = min(hw.keys())
        max_hw = max(hw.keys())
        # Check if the HW distribution is symmetric around m/2
        is_sym = all(hw.get(k, 0) == hw.get(m - k, 0) for k in hw)
        # Compute median
        total = sum(hw.values())
        cumsum = 0
        median_hw = 0
        for k in sorted(hw.keys()):
            cumsum += hw[k]
            if cumsum >= total / 2:
                median_hw = k
                break

        profile = " ".join(f"{k}:{hw[k]}" for k in sorted(hw.keys()))
        if len(profile) > 40:
            profile = profile[:37] + "..."

        print(f"  #{c['ci']:3d} | {c['H']:5d} | {c['size']:5d} | [{min_hw:2d},{max_hw:2d}] | "
              f"{median_hw:5.1f} | {'Y' if is_sym else 'N':10s} | {profile}")


def analyze_grid_symmetric_detail(data):
    """Detailed analysis of grid-symmetric tilings within each class."""
    n = data['n']
    m = data['m']
    cd = data['class_data']
    trans_map = data['trans_map']

    print(f"\n{'='*70}")
    print(f"GRID-SYMMETRIC TILING DETAIL: n={n}")
    print(f"{'='*70}")

    # Count fixed points and orbits of the transpose map
    fixed_tiles = sum(1 for i in range(m) if trans_map[i] == i)
    orbit_pairs = sum(1 for i in range(m) if trans_map[i] > i)
    print(f"Transpose map on tiles: {fixed_tiles} fixed, {orbit_pairs} orbit pairs")
    print(f"Grid-symmetric tilings = 2^({fixed_tiles} + {orbit_pairs}) = 2^{fixed_tiles + orbit_pairs}")

    # For each class with GS tilings, show which "free" bits distinguish them
    for c in cd:
        if c['grid_sym_count'] == 0:
            continue

        gs_members = c['grid_sym_members']
        print(f"\n--- Class #{c['ci']}: {c['grid_sym_count']} GS tilings, H={c['H']}, |Aut|={c['aut_size']}, SC={c['is_self_converse']} ---")

        # Show the GS tilings' bit patterns on the free positions
        free_positions = [i for i in range(m) if trans_map[i] >= i]  # fixed + one from each pair
        print(f"  Free positions ({len(free_positions)}): {free_positions}")

        for t in gs_members[:8]:  # show at most 8
            free_bits = tuple(t['bits'][i] for i in free_positions)
            print(f"    mask={t['mask']:0{m}b} free_bits={free_bits} hw={t['hamming_weight']}")

        if len(gs_members) > 8:
            print(f"    ... ({len(gs_members)} total)")

        # What fraction of all possible GS patterns land in this class?
        total_gs = 2 ** len(free_positions)
        print(f"  {c['grid_sym_count']}/{total_gs} = {c['grid_sym_count']/total_gs:.4f} of all GS tilings")


def analyze_class_persistence(all_data):
    """Track specific tournament structures that persist across n values."""
    print(f"\n{'='*70}")
    print(f"CLASS PERSISTENCE / SYMMETRY KERNEL ANALYSIS")
    print(f"{'='*70}")

    n_values = sorted(all_data.keys())

    # For each n, identify the "most symmetric" classes
    for n_val in n_values:
        cd = all_data[n_val]['class_data']
        print(f"\nn={n_val}: Top symmetric classes (SC + large |Aut|)")
        top = sorted(cd, key=lambda c: (c['is_self_converse'], c['aut_size'], c['anti_aut_size']), reverse=True)[:5]
        for c in top:
            print(f"  #{c['ci']}: H={c['H']}, |Aut|={c['aut_size']}, |AntiAut|={c['anti_aut_size']}, "
                  f"SC={c['is_self_converse']}, SF={c['is_self_flip']}, scores={c['scores']}, "
                  f"GS={c['grid_sym_count']}/{c['size']}")

    # Track the transitive and "full tiling" classes
    print(f"\n--- Transitive tournament class across n ---")
    for n_val in n_values:
        cd = all_data[n_val]['class_data']
        # Class 0 is always transitive
        c = cd[0]
        print(f"  n={n_val}: #{c['ci']}, H={c['H']}, size={c['size']}, |Aut|={c['aut_size']}")

    # Track the cyclic tournament (when n is odd prime)
    print(f"\n--- Maximum H classes across n ---")
    for n_val in n_values:
        cd = all_data[n_val]['class_data']
        max_h = max(c['H'] for c in cd)
        max_classes = [c for c in cd if c['H'] == max_h]
        for c in max_classes:
            print(f"  n={n_val}: #{c['ci']}, H={max_h}, |Aut|={c['aut_size']}, |AntiAut|={c['anti_aut_size']}, "
                  f"SC={c['is_self_converse']}, scores={c['scores']}, size={c['size']}, GS={c['grid_sym_count']}")


def analyze_flip_class_structure(data):
    """Analyze the flip equivalence classes: which iso classes are connected by flips?"""
    n = data['n']
    cd = data['class_data']

    print(f"\n{'='*70}")
    print(f"FLIP CLASS STRUCTURE: n={n}")
    print(f"{'='*70}")

    # Build flip graph: classes connected if some member flips into the other
    flip_edges = set()
    for c in cd:
        for target in c['flip_classes']:
            edge = (min(c['ci'], target), max(c['ci'], target))
            flip_edges.add(edge)

    # Connected components via BFS
    adj = defaultdict(set)
    for a, b in flip_edges:
        adj[a].add(b)
        adj[b].add(a)

    visited = set()
    components = []
    for ci in range(len(cd)):
        if ci in visited:
            continue
        comp = set()
        queue = [ci]
        while queue:
            node = queue.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            for nb in adj[node]:
                if nb not in visited:
                    queue.append(nb)
        components.append(comp)

    print(f"Flip graph: {len(flip_edges)} edges, {len(components)} connected components")
    for i, comp in enumerate(sorted(components, key=lambda c: -len(c))):
        comp_classes = [cd[ci] for ci in sorted(comp)]
        h_values = [c['H'] for c in comp_classes]
        sc_count = sum(1 for c in comp_classes if c['is_self_converse'])
        sf_count = sum(1 for c in comp_classes if c['is_self_flip'])
        total_tilings = sum(c['size'] for c in comp_classes)

        print(f"\n  Component {i}: {len(comp)} classes, {total_tilings} tilings")
        print(f"    H values: {sorted(h_values)}")
        print(f"    SC: {sc_count}/{len(comp)}, SF: {sf_count}/{len(comp)}")
        for c in comp_classes[:10]:
            targets = sorted(c['flip_classes'])
            print(f"    #{c['ci']}: H={c['H']}, flip->{targets}, SC={c['is_self_converse']}, "
                  f"|Aut|={c['aut_size']}")
        if len(comp_classes) > 10:
            print(f"    ... ({len(comp_classes)} total classes)")


if __name__ == '__main__':
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    all_data = {}
    for n in range(3, max_n + 1):
        print(f"\n{'#'*70}")
        print(f"# Building data for n={n}...")
        print(f"{'#'*70}")
        data = build_tournament_data(n)
        all_data[n] = data

        # Per-n analyses
        analyze_symmetry_grouping(data)
        analyze_hamming_weight_structure(data)
        analyze_grid_symmetric_detail(data)
        analyze_flip_class_structure(data)

    # Cross-n analyses
    analyze_symmetry_kernel(list(all_data.values()))
    analyze_class_persistence(all_data)

    # Parent-child analysis
    for n in range(3, max_n):
        if n in all_data and n + 1 in all_data:
            analyze_parent_child_symmetry(all_data[n], all_data[n + 1])

    # Automorphism action (only for small n, expensive)
    for n in range(3, min(max_n + 1, 7)):
        if n in all_data:
            analyze_aut_action_on_tilings(all_data[n])
