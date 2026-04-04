"""
recursive_preservation.py — opus-2026-04-04-S3

DISCOVERY: Multilinear coefficients c_S at n are preserved at n+1 for
"old" tile subsets (tiles not involving the new vertex n).

Test: does this extend to ALL degrees? If so, the multilinear polynomial
has a PERFECTLY RECURSIVE STRUCTURE — every coefficient involving only
old tiles is unchanged when n increases.

Also test: the disjoint pair formula in detail.
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def test_recursive_preservation():
    """Test if ALL multilinear coefficients (all degrees) are preserved n -> n+1."""
    print("=" * 70)
    print("FULL RECURSIVE PRESERVATION TEST")
    print("=" * 70)

    for n_small, n_big in [(4, 5), (5, 6), (6, 7)]:
        print(f"\n--- n={n_small} -> n={n_big} ---")

        tiles_s = make_tiles(n_small)
        tiles_b = make_tiles(n_big)
        m_s = len(tiles_s)
        m_b = len(tiles_b)

        H_s, _, _ = compute_all_H(n_small)
        H_b, _, _ = compute_all_H(n_big)
        coeffs_s = mobius_inversion(H_s, m_s)
        coeffs_b = mobius_inversion(H_b, m_b)

        # Map: big tile index -> small tile index (for old tiles only)
        tile_map = {}  # big_index -> small_index
        for j, t in enumerate(tiles_b):
            if t[0] == n_big:  # new tiles have x = n_big
                continue
            for i, t_s in enumerate(tiles_s):
                if t_s == t:
                    tile_map[j] = i
                    break

        # For each coefficient in the small polynomial, find the corresponding big coefficient
        preserved = defaultdict(int)  # degree -> count preserved
        changed = defaultdict(int)    # degree -> count changed
        examples = []

        for S_s, c_s in coeffs_s.items():
            # Convert S_s (bitmask in small indexing) to S_b (bitmask in big indexing)
            bits_s = [i for i in range(m_s) if S_s & (1 << i)]
            # Find the corresponding big indices
            inverse_map = {v: k for k, v in tile_map.items()}
            bits_b = [inverse_map[i] for i in bits_s if i in inverse_map]
            if len(bits_b) != len(bits_s):
                continue  # some tiles not found (shouldn't happen)

            S_b = 0
            for b in bits_b:
                S_b |= (1 << b)

            c_b = coeffs_b.get(S_b, 0)
            deg = len(bits_s)

            if c_b == c_s:
                preserved[deg] += 1
            else:
                changed[deg] += 1
                if len(examples) < 5:
                    tile_set = [tiles_s[i] for i in bits_s]
                    examples.append((tile_set, deg, c_s, c_b))

        print(f"  By degree:")
        max_deg = max(list(preserved.keys()) + list(changed.keys()))
        total_preserved = 0
        total_changed = 0
        for d in range(max_deg + 1):
            p = preserved.get(d, 0)
            c = changed.get(d, 0)
            total_preserved += p
            total_changed += c
            if p + c > 0:
                print(f"    d={d}: preserved={p}, changed={c} "
                      f"{'ALL PRESERVED' if c == 0 else 'SOME CHANGED'}")

        print(f"  TOTAL: preserved={total_preserved}, changed={total_changed}")
        if total_changed == 0:
            print(f"  *** PERFECT RECURSIVE PRESERVATION at ALL degrees! ***")
        else:
            print(f"  Examples of changes:")
            for ts, d, c_s, c_b in examples:
                print(f"    {ts} (deg {d}): small={c_s}, big={c_b}, diff={c_b-c_s}")

def investigate_disjoint_formula_deep():
    """Deep investigation of the disjoint pair coefficient formula."""
    print(f"\n{'='*70}")
    print(f"DISJOINT PAIR FORMULA — DEEP ANALYSIS")
    print(f"{'='*70}")

    for n in [6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        print(f"\n  n={n}:")

        disjoint_data = []
        for S, c in coeffs.items():
            if bin(S).count('1') != 2: continue
            bits = [i for i in range(m) if S & (1 << i)]
            t1, t2 = tiles[bits[0]], tiles[bits[1]]
            x1, y1 = t1; x2, y2 = t2
            # Check disjoint (no shared vertex of any kind)
            if {x1,y1} & {x2,y2}: continue

            s1, s2 = x1-y1, x2-y2
            overlap = min(x1,x2) - max(y1,y2)

            # Staircase: row = x-y-1, col = y-1
            r1, c1 = x1-y1-1, y1-1
            r2, c2 = x2-y2-1, y2-1

            # Key: the interval of vertices touched by each tile
            # Tile (x,y) touches vertices y, y+1, ..., x (in 1-indexed)
            # Actually in the base path, tile (x,y) has endpoints x-1 and y-1 (0-indexed).
            # Its "range" on the path is [y, x] (1-indexed) = all vertices from y to x.
            range1 = set(range(y1, x1+1))
            range2 = set(range(y2, x2+1))
            range_overlap = len(range1 & range2)

            # Number of internal vertices shared
            # Internal vertices of tile (x,y): {y+1, ..., x-1}
            internal1 = set(range(y1+1, x1))
            internal2 = set(range(y2+1, x2))
            internal_overlap = len(internal1 & internal2)

            disjoint_data.append({
                't1': t1, 't2': t2, 's1': s1, 's2': s2,
                'c': c, 'overlap': overlap,
                'range_overlap': range_overlap,
                'internal_overlap': internal_overlap,
                'r1': r1, 'c1': c1, 'r2': r2, 'c2': c2,
                'ratio': Fraction(c, 2**(s1+s2-2))
            })

        # Group by (overlap, internal_overlap)
        by_ov = defaultdict(list)
        for d in disjoint_data:
            key = (d['overlap'], d['internal_overlap'])
            by_ov[key].append(d)

        print(f"  Grouped by (overlap, internal_overlap):")
        for key in sorted(by_ov.keys()):
            entries = by_ov[key]
            ratios = [d['ratio'] for d in entries]
            distinct = sorted(set(ratios))
            print(f"    ov={key[0]}, int_ov={key[1]}: {len(entries)} pairs, "
                  f"ratios={distinct}")

        # Group by range_overlap
        by_range = defaultdict(list)
        for d in disjoint_data:
            by_range[d['range_overlap']].append(d)

        print(f"\n  Grouped by range_overlap:")
        for ro in sorted(by_range.keys()):
            entries = by_range[ro]
            ratios = sorted(set(d['ratio'] for d in entries))
            print(f"    range_overlap={ro}: {len(entries)} pairs, ratios={ratios}")

        # Try: ratio depends on number of "shared intermediate vertices"?
        # For disjoint tiles (x1,y1) and (x2,y2) with say y1 < y2 < x1 < x2:
        # The tiles overlap in the middle, sharing vertices y2,...,x1.
        # Number of shared vertices = x1 - y2 + 1 when y2 <= x1.
        # (= overlap + 1 when overlap >= 0)
        print(f"\n  Attempt: ratio = f(range_overlap - 2)?")
        for d in sorted(disjoint_data, key=lambda x: x['range_overlap']):
            rov = d['range_overlap']
            print(f"    {d['t1']}-{d['t2']} s=({d['s1']},{d['s2']}) "
                  f"range_ov={rov}: ratio={d['ratio']}")

def investigate_cubic_recursive():
    """Test cubic coefficient preservation and structure."""
    print(f"\n{'='*70}")
    print(f"CUBIC COEFFICIENT RECURSIVE STRUCTURE")
    print(f"{'='*70}")

    for n in [5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        cubics = {S: c for S, c in coeffs.items() if bin(S).count('1') == 3}
        print(f"\n  n={n}: {len(cubics)} cubic coefficients")

        # Group by tile geometry
        for S, c in sorted(cubics.items(), key=lambda x: abs(x[1]), reverse=True)[:15]:
            bits = [i for i in range(m) if S & (1 << i)]
            tile_set = [tiles[i] for i in bits]
            skips = [t[0]-t[1] for t in tile_set]

            # Classify by shared vertices
            all_verts = set()
            for x, y in tile_set:
                all_verts.add(x)
                all_verts.add(y)

            # Pairwise relationships
            rels = []
            for i in range(3):
                for j in range(i+1, 3):
                    t1, t2 = tile_set[i], tile_set[j]
                    v1 = {t1[0], t1[1]}
                    v2 = {t2[0], t2[1]}
                    shared = v1 & v2
                    if not shared:
                        rels.append("disj")
                    elif t1[0] == t2[0] or t1[1] == t2[1]:
                        rels.append("same")
                    else:
                        rels.append("cross")

            print(f"    {tile_set} skips={skips}: c₃={c}, rels={rels}, |verts|={len(all_verts)}")

def investigate_H_symmetry():
    """Check symmetries of H(t): complement, reflection, etc."""
    print(f"\n{'='*70}")
    print(f"H POLYNOMIAL SYMMETRY ANALYSIS")
    print(f"{'='*70}")

    for n in [5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)

        # Complement symmetry: T_comp has all tiles flipped
        # T_comp(t) = T(1-t1, ..., 1-tm) → tiling = (2^m - 1) XOR t
        full = (1 << m) - 1
        complement_match = all(H[t] == H[t ^ full] for t in range(1 << m))
        print(f"\n  n={n}: H(t) = H(complement(t))? {complement_match}")

        # Grid reflection: (x,y) -> (n+1-y, n+1-x)
        # This maps tile (x,y) to tile (n+1-y, n+1-x)
        reflection_map = {}
        for i, (x, y) in enumerate(tiles):
            rx, ry = n+1-y, n+1-x
            for j, (x2, y2) in enumerate(tiles):
                if x2 == rx and y2 == ry:
                    reflection_map[i] = j
                    break

        if len(reflection_map) == m:
            print(f"  Grid reflection is a permutation of tiles")
            # Does H(t) = H(reflected(t))?
            def reflect_tiling(t):
                new_t = 0
                for i in range(m):
                    if t & (1 << i):
                        new_t |= (1 << reflection_map[i])
                return new_t

            reflection_match = all(H[t] == H[reflect_tiling(t)] for t in range(1 << m))
            print(f"  H(t) = H(reflected(t))? {reflection_match}")
        else:
            print(f"  Grid reflection: only {len(reflection_map)}/{m} tiles map to valid tiles")

        # Coefficients: does the multilinear polynomial have reflection symmetry?
        coeffs = mobius_inversion(H, m)
        if len(reflection_map) == m:
            coeff_reflect_match = 0
            coeff_total = 0
            for S, c in coeffs.items():
                if S == 0: continue
                # Reflected S
                bits = [i for i in range(m) if S & (1 << i)]
                reflected_bits = [reflection_map[i] for i in bits]
                S_ref = 0
                for b in reflected_bits:
                    S_ref |= (1 << b)
                c_ref = coeffs.get(S_ref, 0)
                coeff_total += 1
                if c == c_ref:
                    coeff_reflect_match += 1

            print(f"  Coefficient reflection symmetry: {coeff_reflect_match}/{coeff_total}")

def investigate_degree_sum_pattern():
    """Investigate the degree-sum pattern (sums of coefficients by degree)."""
    print(f"\n{'='*70}")
    print(f"DEGREE-SUM PATTERN")
    print(f"{'='*70}")

    for n in [3, 4, 5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        if m == 0:
            continue
        H, _, _ = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        by_deg = defaultdict(int)
        for S, c in coeffs.items():
            by_deg[bin(S).count('1')] += c

        max_d = max(by_deg.keys()) if by_deg else 0
        sums = [by_deg.get(d, 0) for d in range(max_d + 1)]

        H_at_1 = sum(sums)  # = H(1,...,1)
        H_at_0 = sums[0]    # = H(0,...,0) = 1

        print(f"  n={n}: degree sums = {sums}")
        print(f"    H(0)={H_at_0}, H(1)={H_at_1}")
        print(f"    H(0)+H(1) = {H_at_0 + H_at_1}")

        # Alternate evaluations
        # H(1,...,1) - H(0,...,0) = sum_{d>=1} S_d
        # = H(all backward) - H(all forward) = H(anti-transitive) - 1
        print(f"    H(anti-transitive) = {H[(1<<m)-1]}")

        # The generating function of degree sums is H(p,...,p) = sum S_d p^d
        # We know S_0 = 1. Can we predict S_d from the structure?

        # S_1 = sum of linear coefficients = sum_{tiles k} 2^(skip(k)-1)
        expected_S1 = sum(2**(t[0]-t[1]-1) for t in tiles)
        print(f"    S_1 = {sums[1] if len(sums) > 1 else 'N/A'}, expected = {expected_S1}")

if __name__ == "__main__":
    test_recursive_preservation()
    investigate_disjoint_formula_deep()
    investigate_cubic_recursive()
    investigate_H_symmetry()
    investigate_degree_sum_pattern()
