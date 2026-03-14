"""
tiling_minimal_state.py -- kind-pasteur-2026-03-14-S78
Find the MINIMAL STATE that makes strip-by-strip transition deterministic.

From S78: (parent_CLASS, strip_bits) is NOT deterministic.
From first principles: (parent_TILING, strip_bits) IS deterministic.

QUESTION: What is the minimal information about the parent tiling
that determines the child class, given the strip bits?

Candidates:
1. Full parent tiling (trivially deterministic)
2. Parent class + something
3. Parent tiling's "top face" = how the last-added vertex connects to others
4. Parent tiling's score sequence of the "active" vertices
5. The adjacency pattern of the top k vertices

The "active boundary" of the tiling at step k is the information needed
to determine how future strips attach. This is like the TRANSFER MATRIX
approach — the state is the "boundary condition."

For the pin grid, adding vertex n means adding arcs (n, 1), (n, 2), ..., (n, n-2).
These arcs determine how vertex n interacts with vertices 1..n-2.
But the EFFECT on the tournament class also depends on how vertices 1..n-1
interact with each other — which is the full parent tiling.

The minimal state might be: the out-neighborhoods of vertex n-1 in the
parent tournament (since vertex n-1 is the "boundary" between old and new).

Let's test: is the (parent_class, score(n-1)) pair deterministic?
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

def build_tiling_data(n):
    tiles = []
    for b in range(1, n-1):
        for a in range(b+2, n+1):
            tiles.append((a, b))
    m = len(tiles)
    tile_idx = {t: i for i, t in enumerate(tiles)}

    trans_map = []
    for a, b in tiles:
        a2, b2 = n+1-b, n+1-a
        if a2 < b2: a2, b2 = b2, a2
        trans_map.append(tile_idx.get((a2, b2), -1))

    def bits_to_adj(mask):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for idx, (a, b) in enumerate(tiles):
            a0, b0 = a-1, b-1
            if (mask >> idx) & 1: A[b0][a0] = 1
            else: A[a0][b0] = 1
        return A

    perms = list(permutations(range(n)))
    def canonicalize(A):
        best = None
        for p in perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    mask_data = {}
    for mask in range(1 << m):
        A = bits_to_adj(mask)
        canon = canonicalize(A)
        mask_data[mask] = {'adj': A, 'canon': canon}

    groups = defaultdict(list)
    for mask, d in mask_data.items():
        groups[d['canon']].append(mask)

    class_list = sorted(groups.keys())
    class_idx = {c: i for i, c in enumerate(class_list)}

    return {'n': n, 'm': m, 'tiles': tiles, 'mask_data': mask_data,
            'groups': groups, 'class_list': class_list, 'class_idx': class_idx}

def main():
    print("=" * 70)
    print("MINIMAL STATE FOR DETERMINISTIC STRIP TRANSITIONS")
    print("kind-pasteur-2026-03-14-S78")
    print("=" * 70)

    for n in [4, 5]:
        data_n = build_tiling_data(n)
        data_p = build_tiling_data(n-1)
        m_n = data_n['m']
        m_p = data_p['m']

        # New strip tiles: arcs involving vertex n (1-indexed)
        new_strip = [i for i, (a, b) in enumerate(data_n['tiles']) if a == n]
        strip_size = len(new_strip)
        old_tiles = [i for i in range(m_n) if i not in new_strip]

        perms_p = list(permutations(range(n-1)))
        def canon_p(A):
            best = None
            for p in perms_p:
                s = ''.join(str(A[p[i]][p[j]]) for i in range(n-1) for j in range(n-1))
                if best is None or s < best: best = s
            return best

        print(f"\n{'='*70}")
        print(f"n = {n}: testing minimal state candidates")
        print(f"{'='*70}")

        # Build transition data
        transitions = []  # (parent_mask, strip_bits, child_class)
        for mask_n in range(1 << m_n):
            A = data_n['mask_data'][mask_n]['adj']
            # Parent sub-adjacency
            sub_A = [[A[i][j] for j in range(1, n)] for i in range(1, n)]
            parent_canon = canon_p(sub_A)
            ci_p = data_p['class_idx'].get(parent_canon, -1)

            # Parent mask (the bits for non-new-strip tiles)
            parent_mask = 0
            for i in old_tiles:
                if (mask_n >> i) & 1:
                    parent_mask |= (1 << old_tiles.index(i))

            strip_config = tuple((mask_n >> idx) & 1 for idx in new_strip)
            ci_n = data_n['class_idx'][data_n['mask_data'][mask_n]['canon']]

            # Score sequence of parent
            parent_scores = tuple(sorted(sum(row) for row in sub_A))

            # Score of vertex n-1 (0-indexed: vertex n-2) in the FULL tournament
            score_top = sum(A[n-2])  # out-degree of vertex n-2 in n-tournament
            # But this includes the new strip arcs! Use parent score instead.
            parent_score_top = sum(sub_A[n-2])  # out-degree of vertex n-2 in (n-1)-tournament

            # The "boundary state": how vertex n-1 connects to 1..n-2 in the parent
            # This is just the last row of the parent adjacency
            boundary = tuple(sub_A[n-2][j] for j in range(n-1))

            transitions.append({
                'parent_mask': parent_mask,
                'parent_class': ci_p,
                'strip': strip_config,
                'child_class': ci_n,
                'parent_scores': parent_scores,
                'parent_score_top': parent_score_top,
                'boundary': boundary,
            })

        # Test candidate minimal states
        candidates = {
            'parent_class': lambda t: t['parent_class'],
            'parent_class + strip': lambda t: (t['parent_class'], t['strip']),
            'parent_class + score_top + strip': lambda t: (t['parent_class'], t['parent_score_top'], t['strip']),
            'parent_scores + strip': lambda t: (t['parent_scores'], t['strip']),
            'parent_mask + strip': lambda t: (t['parent_mask'], t['strip']),
            'boundary + strip': lambda t: (t['boundary'], t['strip']),
        }

        for name, key_fn in candidates.items():
            # Group by key, check if child_class is unique
            groups = defaultdict(set)
            for t in transitions:
                key = key_fn(t)
                groups[key].add(t['child_class'])

            deterministic = all(len(v) == 1 for v in groups.values())
            max_ambiguity = max(len(v) for v in groups.values())
            n_ambiguous = sum(1 for v in groups.values() if len(v) > 1)

            print(f"  {name:40s}: deterministic={deterministic}, "
                  f"max_ambiguity={max_ambiguity}, n_ambiguous={n_ambiguous}/{len(groups)}")

        # The TILING-level state (parent_mask + strip) should be deterministic
        # because it fully determines the child tiling
        # But class(child) depends on the full child tiling, which IS determined
        # by (parent_tiling, strip_bits). So parent_mask + strip is deterministic.

        # KEY QUESTION: What is the MINIMUM state between class and full tiling?
        # Try: parent_class + boundary
        boundary_test = defaultdict(set)
        for t in transitions:
            key = (t['parent_class'], t['boundary'], t['strip'])
            boundary_test[key].add(t['child_class'])

        det_boundary = all(len(v) == 1 for v in boundary_test.values())
        print(f"\n  parent_class + boundary + strip: deterministic={det_boundary}")

        if det_boundary:
            print(f"    THE BOUNDARY ROW IS THE MINIMAL STATE!")
            print(f"    The transition depends on parent CLASS + how the top vertex")
            print(f"    connects to the rest, plus the new strip configuration.")

    # ========================================
    # PART 2: What does deterministic boundary mean?
    # ========================================
    print(f"\n{'='*70}")
    print("PART 2: THE BOUNDARY STATE — WHAT IT TELLS US")
    print("  The boundary = last row of parent adjacency matrix")
    print("  = how vertex n-1 beats/loses to vertices 1..n-2")
    print("  Combined with parent class, this determines child class")
    print(f"{'='*70}")

    n = 5
    data_n = build_tiling_data(n)
    data_p = build_tiling_data(n-1)
    m_n = data_n['m']

    new_strip = [i for i, (a, b) in enumerate(data_n['tiles']) if a == n]
    old_tiles = [i for i in range(m_n) if i not in new_strip]

    perms_p = list(permutations(range(n-1)))
    def canon_p(A):
        best = None
        for p in perms_p:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(n-1) for j in range(n-1))
            if best is None or s < best: best = s
        return best

    # For each parent class, enumerate the possible boundaries
    class_boundaries = defaultdict(set)
    for mask_n in range(1 << m_n):
        A = data_n['mask_data'][mask_n]['adj']
        sub_A = [[A[i][j] for j in range(1, n)] for i in range(1, n)]
        parent_canon = canon_p(sub_A)
        ci_p = data_p['class_idx'].get(parent_canon, -1)
        boundary = tuple(sub_A[n-2][j] for j in range(n-1))
        class_boundaries[ci_p].add(boundary)

    print(f"\n  n={n}: Boundaries per parent class:")
    for ci_p in sorted(class_boundaries.keys()):
        boundaries = class_boundaries[ci_p]
        print(f"    class {ci_p}: {len(boundaries)} distinct boundaries")
        for b in sorted(boundaries):
            print(f"      {b} (out-degree from top = {sum(b)})")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
