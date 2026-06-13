#!/usr/bin/env python3
"""
Connection between blueself classes and SC H-maximizers.

HYPOTHESIS: At even n, the blueself classes are exactly the SC classes
that achieve maximum H within their score sequence class.

From THM-022:
  - n=4: 1 blueself class with H=5 (global maximum)
  - n=6: 2 blueself classes with H=37 and H=45

From T091 (SC maximizer):
  - n=4: SC max within each score class
  - n=6: SC max = 45 (score (3,3,3,2,2,2)), SC max = 37 (score (4,3,3,2,2,1))

QUESTION: Are the blueself classes EXACTLY the SC maximizers?

Also investigates: what makes blueself tournaments special?
The grid-symmetry of the tiling + self-flip = specific structural constraint.

Author: opus-2026-03-06-S17
"""

import sys
import os
import itertools
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, hamiltonian_path_count,
                             opposite_tournament, find_odd_cycles, conflict_graph)


def canonical_form(T):
    n = len(T)
    best = None
    for perm in itertools.permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best


def is_self_converse(T):
    n = len(T)
    Top = [[T[j][i] for j in range(n)] for i in range(n)]
    return canonical_form(T) == canonical_form(Top)


def score_sequence(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)], reverse=True))


def indep_poly(adj):
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


# Tiling model functions
def build_tiles(n):
    tiles = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            tiles.append((x, y))
    return tiles


def build_transpose_map(n, tiles):
    tile_idx = {(x, y): i for i, (x, y) in enumerate(tiles)}
    return [tile_idx[(n - y + 1, n - x + 1)] for i, (x, y) in enumerate(tiles)]


def is_grid_symmetric(bits, trans_map):
    for i in range(len(bits)):
        if trans_map[i] != i and bits[i] != bits[trans_map[i]]:
            return False
    return True


def are_isomorphic(A1, A2, n):
    for p in itertools.permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if A1[p[i]][p[j]] != A2[i][j]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False


def bits_to_adj(n, tiles, bits):
    verts = list(range(n, 0, -1))
    vert_idx = {v: i for i, v in enumerate(verts)}
    A = [[0] * n for _ in range(n)]
    for k in range(n - 1):
        A[k][k + 1] = 1
    for i, (x, y) in enumerate(tiles):
        xi, yi = vert_idx[x], vert_idx[y]
        if bits[i] == 0:
            A[xi][yi] = 1
        else:
            A[yi][xi] = 1
    return A


# ============================================================
# Analyze blueself classes vs SC maximizers at n=4,6
# ============================================================
for n in [4, 6]:
    print(f"\n{'='*70}")
    print(f"  BLUESELF vs SC MAXIMIZER at n={n}")
    print(f"{'='*70}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)
    total = 1 << m

    # Find all blueself tilings
    blueself_masks = []
    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        if not is_grid_symmetric(bits, trans_map):
            continue
        # Check if self-flip
        flip_bits = tuple(1 - b for b in bits)
        A = bits_to_adj(n, tiles, bits)
        A_f = bits_to_adj(n, tiles, flip_bits)
        if are_isomorphic(A, A_f, n):
            blueself_masks.append(mask)

    print(f"  Blueself tilings found: {len(blueself_masks)}")

    # Group blueself tilings by isomorphism class
    blueself_classes = {}
    for mask in blueself_masks:
        bits = tuple((mask >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        cf = canonical_form(A)
        if cf not in blueself_classes:
            H = hamiltonian_path_count(A)
            ss = score_sequence(A)
            cycles = find_odd_cycles(A)
            cg = conflict_graph(cycles) if cycles else []
            ip = indep_poly(cg)
            blueself_classes[cf] = {
                'H': H, 'scores': ss, 'ip': ip,
                'tilings': [], 'is_SC': is_self_converse(A)
            }
        blueself_classes[cf]['tilings'].append(mask)

    print(f"  Blueself isomorphism classes: {len(blueself_classes)}")
    for cf, info in blueself_classes.items():
        print(f"    H={info['H']}, scores={info['scores']}, I.P.={info['ip']}, "
              f"SC={info['is_SC']}, #tilings={len(info['tilings'])}")

    # Now find ALL SC maximizers within score classes
    all_classes = {}
    for bits_int in range(total):
        T = tournament_from_bits(n, bits_int)
        cf = canonical_form(T)
        if cf not in all_classes:
            H = hamiltonian_path_count(T)
            ss = score_sequence(T)
            sc = is_self_converse(T)
            all_classes[cf] = {'H': H, 'scores': ss, 'is_SC': sc}

    # Group by score sequence
    by_scores = defaultdict(list)
    for cf, info in all_classes.items():
        by_scores[info['scores']].append(info)

    print(f"\n  Score sequence analysis:")
    for ss, classes in sorted(by_scores.items()):
        sc_classes = [c for c in classes if c['is_SC']]
        nsc_classes = [c for c in classes if not c['is_SC']]
        max_H_sc = max(c['H'] for c in sc_classes) if sc_classes else -1
        max_H_nsc = max(c['H'] for c in nsc_classes) if nsc_classes else -1
        max_H_overall = max(c['H'] for c in classes)

        # Check if any blueself class has this score
        blueself_in_score = [cf for cf, info in blueself_classes.items()
                            if info['scores'] == ss]

        marker = ""
        if blueself_in_score:
            blueself_H = [blueself_classes[cf]['H'] for cf in blueself_in_score]
            if max(blueself_H) == max_H_sc:
                marker = " <<< BLUESELF = SC MAX"
            else:
                marker = f" <<< BLUESELF H={blueself_H} ≠ SC MAX {max_H_sc}"

        print(f"    scores={ss}: {len(sc_classes)} SC, {len(nsc_classes)} NSC | "
              f"max H: SC={max_H_sc}, NSC={max_H_nsc}, overall={max_H_overall}{marker}")

    # Check: are ALL blueself classes SC maximizers within their score class?
    all_match = True
    for cf, info in blueself_classes.items():
        ss = info['scores']
        sc_max = max(c['H'] for c in by_scores[ss] if c['is_SC'])
        if info['H'] != sc_max:
            all_match = False
            print(f"  MISMATCH: blueself H={info['H']} but SC max H={sc_max} for {ss}")

    if all_match and blueself_classes:
        print(f"\n  *** CONFIRMED: Blueself classes ARE exactly the SC maximizers ***")
    elif not blueself_classes:
        print(f"\n  (No blueself classes at n={n})")
    else:
        print(f"\n  *** HYPOTHESIS FAILS: blueself ≠ SC maximizer ***")

    # Are SC maximizers always blueself? Check the converse
    sc_max_classes = set()
    for ss, classes in by_scores.items():
        sc_classes = [c for c in classes if c['is_SC']]
        if sc_classes:
            max_sc = max(c['H'] for c in sc_classes)
            for c in sc_classes:
                if c['H'] == max_sc:
                    # Find the canonical form for this class
                    for cf, info in all_classes.items():
                        if info['H'] == c['H'] and info['scores'] == c['scores'] and info['is_SC']:
                            sc_max_classes.add(cf)

    blueself_cfs = set(blueself_classes.keys())
    print(f"\n  SC maximizer classes: {len(sc_max_classes)}")
    print(f"  Blueself classes: {len(blueself_cfs)}")
    extra_sc_max = sc_max_classes - blueself_cfs
    if extra_sc_max:
        print(f"  SC maximizers that are NOT blueself: {len(extra_sc_max)}")
        for cf in extra_sc_max:
            info = all_classes[cf]
            print(f"    H={info['H']}, scores={info['scores']}")
    else:
        print(f"  ALL SC maximizers are blueself (converse also holds)")

print(f"\n{'='*70}")
print("SUMMARY")
print("="*70)
print("""
The blueself property combines two constraints:
1. Grid-symmetry: tiling invariant under (x,y) -> (n+1-y, n+1-x)
2. Self-flip: tiling isomorphic to its bit complement

The SC maximizer property is:
3. Self-converse: T ≅ T^op
4. Maximum H within the score sequence class

QUESTION: Does (1)+(2) imply (3)+(4)?
- (2) => isomorphic to flip, which is related but not identical to T^op
- Grid-symmetry imposes score structure that may force near-regularity
- Near-regularity -> high H -> SC max candidate

This connection between the TILING model's grid symmetry and the
TOURNAMENT model's self-converse H-maximization is a deep structural bridge.
""")
