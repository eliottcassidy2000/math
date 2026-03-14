"""
tiling_monadic_structure.py -- kind-pasteur-2026-03-14-S78
Monadic structure of the tiling model: strip-by-strip composition.

THE MONAD:
- Object: isomorphism class of tournament
- Unit: the unique class at n=2 (single arc)
- Bind: extend tiling by adding vertex n (= adding strip Str(n-1))
  A tiling at n-1 with class C, when extended by the choices in the
  new strip, produces a tiling at n with some class C'.
  The "bind" maps (C, strip_bits) -> C'.

This gives a TRANSITION SYSTEM:
  T(n-1 -> n)[C][s] = C'
where C is a class at n-1, s is the strip choice (k bits), C' is class at n.

QUESTIONS:
1. Is the transition deterministic? (C, s) always -> unique C'?
   YES — the tiling uniquely determines the tournament.
2. How many classes at n come from each class at n-1?
3. Is the transition "regular" — same branching for each class?
4. Can we compute the class transition matrix?
5. What is the monadic MULTIPLICATION (join)?
   Extend from n-2 to n in one step vs. two steps (n-2 -> n-1 -> n).

ALSO: The relationship between tiling structure and isomorphism class.
- Each class has a "fiber" of tilings over it
- The fiber size = |Aut(T)| * (orbit stuff)
- Blue/black lines connect fibers via flip

FAST FORMULAS to derive:
- Class sizes: |C| = n! / |Aut(T)| (orbit-stabilizer)
- GS count in class: related to automorphisms preserving the transpose
- H value: from OCF
- Total tilings: 2^m = sum of class sizes
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def build_tiling_data(n):
    """Build complete tiling data for size n."""
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

    def is_gs(mask):
        for i in range(m):
            j = trans_map[i]
            if j != i and ((mask >> i) & 1) != ((mask >> j) & 1):
                return False
        return True

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

    def count_ham(A):
        dp = {}
        for v in range(n): dp[(1 << v, v)] = 1
        for ms in range(2, n+1):
            for mask in range(1 << n):
                if bin(mask).count('1') != ms: continue
                for v in range(n):
                    if not (mask & (1 << v)): continue
                    pm = mask ^ (1 << v)
                    t = sum(dp.get((pm, u), 0) for u in range(n) if (pm & (1 << u)) and A[u][v])
                    if t: dp[(mask, v)] = t
        return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

    # Build everything
    mask_data = {}
    for mask in range(1 << m):
        A = bits_to_adj(mask)
        canon = canonicalize(A)
        gs = is_gs(mask)
        mask_data[mask] = {'adj': A, 'canon': canon, 'gs': gs}

    groups = defaultdict(list)
    for mask, d in mask_data.items():
        groups[d['canon']].append(mask)

    class_list = sorted(groups.keys())
    class_idx = {c: i for i, c in enumerate(class_list)}

    # Compute H for each class
    class_H = {}
    class_scores = {}
    class_c3 = {}
    for canon in class_list:
        rep_mask = groups[canon][0]
        A = mask_data[rep_mask]['adj']
        class_H[canon] = count_ham(A)
        class_scores[canon] = tuple(sorted(sum(row) for row in A))
        c3 = sum(1 for i in range(n) for j in range(n) for k in range(n)
                 if A[i][j] and A[j][k] and A[k][i]) // 3
        class_c3[canon] = c3

    return {
        'n': n, 'm': m, 'tiles': tiles, 'trans_map': trans_map,
        'mask_data': mask_data, 'groups': groups,
        'class_list': class_list, 'class_idx': class_idx,
        'class_H': class_H, 'class_scores': class_scores, 'class_c3': class_c3,
    }

def main():
    print("=" * 70)
    print("MONADIC STRUCTURE OF THE TILING MODEL")
    print("kind-pasteur-2026-03-14-S78")
    print("=" * 70)

    # Build data for n = 3, 4, 5, 6
    all_data = {}
    for n in [3, 4, 5, 6]:
        print(f"\n  Building n={n}...", end=" ", flush=True)
        all_data[n] = build_tiling_data(n)
        print(f"{len(all_data[n]['class_list'])} classes, {1 << all_data[n]['m']} tilings")

    # ========================================
    # PART 1: PARENT-CHILD CLASS TRANSITIONS
    # ========================================
    print(f"\n{'='*70}")
    print("PART 1: PARENT-CHILD CLASS TRANSITIONS (MONADIC BIND)")
    print("  For each class at n, which class at n-1 does it come from?")
    print("  'Come from' = delete vertex n (the top vertex) from the tournament")
    print(f"{'='*70}")

    for n in [4, 5, 6]:
        data_n = all_data[n]
        data_p = all_data[n-1]

        perms_p = list(permutations(range(n-1)))
        def canon_p(A):
            best = None
            for p in perms_p:
                s = ''.join(str(A[p[i]][p[j]]) for i in range(n-1) for j in range(n-1))
                if best is None or s < best: best = s
            return best

        # For each class at n, find its parent class at n-1
        parent_map = {}  # class_at_n -> class_at_n-1
        children_map = defaultdict(list)  # class_at_n-1 -> [classes_at_n]

        for canon in data_n['class_list']:
            rep_mask = data_n['groups'][canon][0]
            A = data_n['mask_data'][rep_mask]['adj']
            # Delete vertex 0 (which is vertex n in 1-indexed = top vertex)
            # Sub-adjacency on vertices 1..n-1
            sub_A = [[A[i][j] for j in range(1, n)] for i in range(1, n)]
            parent_canon = canon_p(sub_A)

            ci_n = data_n['class_idx'][canon]
            ci_p = data_p['class_idx'].get(parent_canon, -1)
            parent_map[ci_n] = ci_p
            children_map[ci_p].append(ci_n)

        print(f"\n  n={n} -> n-1={n-1}:")
        print(f"  {'Parent(n-1)':>12} {'H_p':>5} {'#Children':>10} {'Children H values':>30}")
        for ci_p in sorted(children_map.keys()):
            canon_p_sig = data_p['class_list'][ci_p]
            H_p = data_p['class_H'][canon_p_sig]
            kids = children_map[ci_p]
            kid_Hs = sorted([data_n['class_H'][data_n['class_list'][k]] for k in kids])
            print(f"  {ci_p:12d} {H_p:5d} {len(kids):10d} {kid_Hs}")

        # Transition matrix: how many TILINGS in child class come from parent?
        # Actually: each tiling at n is an extension of a tiling at n-1.
        # The "new" arcs are those in the top strip Str(n-1).
        # Strip Str(n-1) has n-2 tiles: arcs (n, 1), (n, 2), ..., (n, n-2)
        # (In 1-indexed: arcs from vertex n to vertices 1..n-2)

        # Count transitions per tiling
        print(f"\n  Transition counts (tilings in child from parent class):")
        trans_counts = defaultdict(lambda: defaultdict(int))

        for mask_n in range(1 << data_n['m']):
            A = data_n['mask_data'][mask_n]['adj']
            # Parent tiling: restrict to vertices 1..n-1 (0-indexed: 1..n-1)
            sub_A = [[A[i][j] for j in range(1, n)] for i in range(1, n)]
            parent_canon = canon_p(sub_A)
            ci_p = data_p['class_idx'].get(parent_canon, -1)

            ci_n = data_n['class_idx'][data_n['mask_data'][mask_n]['canon']]
            trans_counts[ci_p][ci_n] += 1

        for ci_p in sorted(trans_counts.keys()):
            canon_p_sig = data_p['class_list'][ci_p]
            H_p = data_p['class_H'][canon_p_sig]
            total = sum(trans_counts[ci_p].values())
            for ci_n in sorted(trans_counts[ci_p].keys()):
                count = trans_counts[ci_p][ci_n]
                canon_n_sig = data_n['class_list'][ci_n]
                H_n = data_n['class_H'][canon_n_sig]
                print(f"    parent {ci_p}(H={H_p}) -> child {ci_n}(H={H_n}): {count} tilings")

    # ========================================
    # PART 2: CLASS SIZE FORMULA
    # ========================================
    print(f"\n{'='*70}")
    print("PART 2: CLASS SIZE = n! / |Aut(T)| (ORBIT-STABILIZER)")
    print(f"{'='*70}")

    for n in [3, 4, 5, 6]:
        data = all_data[n]
        print(f"\n  n={n}:")
        total_check = 0
        for canon in data['class_list']:
            size = len(data['groups'][canon])
            H = data['class_H'][canon]
            # |Aut| = n! / size
            aut_size = math.factorial(n) // size
            ci = data['class_idx'][canon]
            total_check += size
            print(f"    class {ci:3d}: size={size:4d}, |Aut|={aut_size:3d}, H={H:3d}, "
                  f"check n!/|Aut|={math.factorial(n)//aut_size}")
        print(f"  Total: {total_check} = 2^{data['m']} = {1 << data['m']}? {total_check == 1 << data['m']}")

    # ========================================
    # PART 3: GS COUNT PER CLASS
    # ========================================
    print(f"\n{'='*70}")
    print("PART 3: GS TILINGS PER CLASS — WHAT DETERMINES THE COUNT?")
    print(f"{'='*70}")

    for n in [4, 5, 6]:
        data = all_data[n]
        print(f"\n  n={n}:")
        print(f"  {'CI':>4} {'Size':>5} {'GS':>3} {'nGS':>4} {'|Aut|':>5} {'H':>4} {'SC':>3} {'GS/Size':>8} {'GS*|Aut|/n!':>12}")

        for canon in data['class_list']:
            members = data['groups'][canon]
            size = len(members)
            gs_count = sum(1 for mask in members if data['mask_data'][mask]['gs'])
            ngs = size - gs_count
            aut_size = math.factorial(n) // size
            H = data['class_H'][canon]
            ci = data['class_idx'][canon]

            # Is this class self-converse?
            rep_A = data['mask_data'][members[0]]['adj']
            # Transpose: vertex v -> n-1-v (0-indexed)
            trans_A = [[rep_A[n-1-i][n-1-j] for j in range(n)] for i in range(n)]
            trans_canon = None
            for p in permutations(range(n)):
                s = ''.join(str(trans_A[p[i]][p[j]]) for i in range(n) for j in range(n))
                if trans_canon is None or s < trans_canon: trans_canon = s
            is_sc = (trans_canon == canon)

            gs_ratio = gs_count / size if size > 0 else 0
            gs_aut_ratio = gs_count * aut_size / math.factorial(n) if gs_count > 0 else 0

            sc_str = 'SC' if is_sc else ''
            print(f"  {ci:4d} {size:5d} {gs_count:3d} {ngs:4d} {aut_size:5d} {H:4d} {sc_str:>3s} "
                  f"{gs_ratio:8.4f} {gs_aut_ratio:12.4f}")

    # ========================================
    # PART 4: MONADIC COMPOSITION — STRIP DECOMPOSITION
    # ========================================
    print(f"\n{'='*70}")
    print("PART 4: STRIP-BY-STRIP DECOMPOSITION")
    print("  The tiling is built strip by strip: Str(2), Str(3), ..., Str(n-1)")
    print("  Strip Str(k) has k-1 tiles")
    print("  Total tiles: 1 + 2 + ... + (n-2) = C(n-1, 2)")
    print(f"{'='*70}")

    n = 5
    data = all_data[n]
    tiles = data['tiles']
    m = data['m']

    # Organize tiles by strip
    strips = defaultdict(list)
    for idx, (a, b) in enumerate(tiles):
        strip_num = a  # Strip k contains tiles (k, 1), (k, 2), ..., (k, k-2)
        # Actually strip Str(k) has positions (r,c) with r+c = k
        # In terms of arcs: r = a-b-1, c = b, so r+c = a-1
        # So strip number = a-1 (0-indexed from strip 2)
        strips[a].append(idx)

    print(f"\n  n={n}: Strips:")
    for k in sorted(strips.keys()):
        tile_list = [(tiles[i][0], tiles[i][1]) for i in strips[k]]
        print(f"    Vertex {k} (strip Str({k-1})): {len(strips[k])} tiles = {tile_list}")

    # For each tiling, extract strip-by-strip configuration
    print(f"\n  Strip-by-strip state transitions (first 20 tilings):")
    for mask in range(min(20, 1 << m)):
        A = data['mask_data'][mask]['adj']
        ci = data['class_idx'][data['mask_data'][mask]['canon']]
        H = data['class_H'][data['class_list'][ci]]

        strip_bits = {}
        for k in sorted(strips.keys()):
            bits = tuple((mask >> idx) & 1 for idx in strips[k])
            strip_bits[k] = bits

        print(f"    mask={mask:06b}: strips={dict(strip_bits)}, class={ci}, H={H}")

    # ========================================
    # PART 5: FAST FORMULAS SUMMARY
    # ========================================
    print(f"\n{'='*70}")
    print("PART 5: FAST FORMULAS — WHAT WE CAN BE CERTAIN OF")
    print(f"{'='*70}")

    formulas = [
        ("Total tilings", "2^{C(n-1,2)}", "exact"),
        ("GS tilings", "2^{(C(n-1,2) + floor((n-1)/2)) / 2}", "proved"),
        ("GS fixed points", "floor((n-1)/2)", "proved"),
        ("GS paired positions", "(C(n-1,2) - floor((n-1)/2)) / 2", "proved"),
        ("GS weight enumerator", "(1+z)^f * (1+z^2)^p", "proved (S76)"),
        ("Class size", "n! / |Aut(T)|", "orbit-stabilizer theorem"),
        ("Sum of class sizes", "2^{C(n-1,2)}", "partition of tilings"),
        ("Flip of GS is GS", "always", "proved (complement of equal = equal)"),
        ("Blue skeleton bipartite (odd n)", "by t3 parity", "proved (THM-060)"),
        ("Blue line weights", "always even", "verified n=3..6"),
        ("Blueself at odd n", "0", "proved (THM-023)"),
        ("H = I(Omega, 2)", "OCF formula", "proved (Grinberg-Stanley)"),
        ("deg(H) in arc variables", "2*floor((n-1)/2)", "proved (S72, path reversal)"),
        ("Fourier level-2 magnitude", "(n-2)!/2^{n-2}", "proved (S75)"),
    ]

    for name, formula, status in formulas:
        print(f"  {name:40s}: {formula:35s} [{status}]")

    # ========================================
    # PART 6: THE MONADIC BIND — COMPUTING TRANSITIONS
    # ========================================
    print(f"\n{'='*70}")
    print("PART 6: THE MONADIC BIND — TRANSITION MATRIX")
    print("  M[parent_class, strip_config] = child_class")
    print("  This is the CORE of the monadic structure")
    print(f"{'='*70}")

    for n in [4, 5]:
        data_n = all_data[n]
        data_p = all_data[n-1]

        # The new strip when going from n-1 to n:
        # Arcs from vertex n to vertices 1..n-2 (not n-1, that's backbone)
        # In 0-indexed: arcs from vertex n-1 to vertices 0..n-3
        # These are tiles (n, b) for b = 1..n-2
        new_strip = [i for i, (a, b) in enumerate(data_n['tiles']) if a == n]
        strip_size = len(new_strip)

        print(f"\n  n={n}: new strip has {strip_size} tiles (arcs from vertex {n})")

        # Build transition matrix
        perms_p = list(permutations(range(n-1)))
        def canon_p(A):
            best = None
            for p in perms_p:
                s = ''.join(str(A[p[i]][p[j]]) for i in range(n-1) for j in range(n-1))
                if best is None or s < best: best = s
            return best

        # For each parent tiling × strip config → child class
        transitions = defaultdict(Counter)  # (parent_class, strip_bits) -> child_class_counts

        for mask_n in range(1 << data_n['m']):
            # Extract parent: sub-tournament on vertices 1..n-1 (0-indexed: 1..n-1)
            A = data_n['mask_data'][mask_n]['adj']
            sub_A = [[A[i][j] for j in range(1, n)] for i in range(1, n)]
            parent_canon = canon_p(sub_A)
            ci_p = data_p['class_idx'].get(parent_canon, -1)

            # Extract strip configuration
            strip_config = tuple((mask_n >> idx) & 1 for idx in new_strip)

            # Child class
            ci_n = data_n['class_idx'][data_n['mask_data'][mask_n]['canon']]

            transitions[(ci_p, strip_config)][ci_n] += 1

        # Display as transition table
        n_parent_classes = len(data_p['class_list'])
        n_strip_configs = 1 << strip_size

        print(f"\n  Transition table: {n_parent_classes} parent classes x {n_strip_configs} strip configs")
        print(f"  {'Parent':>7} {'Strip':>10} -> {'Child classes (ci: count)':>40}")

        for ci_p in range(n_parent_classes):
            H_p = data_p['class_H'][data_p['class_list'][ci_p]]
            for s in range(n_strip_configs):
                strip_bits = tuple((s >> i) & 1 for i in range(strip_size))
                key = (ci_p, strip_bits)
                if key in transitions:
                    child_dist = transitions[key]
                    child_str = ', '.join(f"{ci}(H={data_n['class_H'][data_n['class_list'][ci]]}): {c}"
                                          for ci, c in sorted(child_dist.items()))
                    print(f"  {ci_p:4d}(H={H_p:2d}) {str(strip_bits):>10s} -> {child_str}")

        # KEY: Is the child class UNIQUE for each (parent, strip)?
        # I.e., is the transition deterministic?
        deterministic = all(len(v) == 1 for v in transitions.values())
        print(f"\n  Transition deterministic? {deterministic}")
        if not deterministic:
            for key, dist in transitions.items():
                if len(dist) > 1:
                    print(f"    NON-DETERMINISTIC: parent={key[0]}, strip={key[1]} -> {dict(dist)}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
