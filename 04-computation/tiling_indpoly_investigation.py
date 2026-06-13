"""
Deep investigation of the relationship between:
  - The tiling model (flip, grid-symmetry, transpose) from tournament-tiling-explorer.html
  - The independence polynomial I(Omega(T), x) of the odd-cycle conflict graph

Key definitions from the HTML visualization:
  - Base path: n -> n-1 -> ... -> 1 (fixed Hamiltonian path)
  - Tiles: non-path edges (a,b) with a >= b+2; m = C(n-1,2) tiles
  - Tiling: binary vector in {0,1}^m; bit 0 = forward arc a->b, bit 1 = backward b->a
  - Flip: complement all bits (T -> "opposite direction on non-path edges")
  - Grid symmetry (transpose): reflect tile grid across y=x; corresponds to T -> T^op (reverse all arcs)
    but adjusted for the fixed path structure
  - Isomorphism class: equivalence under vertex permutations

Derived concepts:
  - "blueself": a tiling T whose FLIP PARTNER lands in the SAME isomorphism class,
    AND T is grid-symmetric (self-converse in the tiling sense)
  - "blackself": flip partner in same class, but T is NOT grid-symmetric
  - "blue-cross": flip partner in DIFFERENT class, both grid-symmetric
  - "black-cross": flip partner in different class, at least one not grid-symmetric

Author: opus-2026-03-06-S16
"""

import itertools
from collections import defaultdict, Counter
import numpy as np


# ─────────────────────────────────────────────────────────────
# CORE: Tournament construction from tiling bits
# ─────────────────────────────────────────────────────────────

def build_tiles(n):
    """Build tile list matching the HTML: (a,b) with a >= b+2, ordered as in the HTML."""
    tiles = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            tiles.append((x, y))
    return tiles


def build_transpose_map(n, tiles):
    """Build the grid-transpose map: tile (x,y) -> tile (n-y+1, n-x+1)."""
    tile_idx = {(x, y): i for i, (x, y) in enumerate(tiles)}
    trans_map = []
    for i, (x, y) in enumerate(tiles):
        tx, ty = n - y + 1, n - x + 1
        trans_map.append(tile_idx[(tx, ty)])
    return trans_map


def bits_to_adj(n, tiles, bits):
    """Convert tiling bits to adjacency matrix. VERTS = [n, n-1, ..., 1]."""
    verts = list(range(n, 0, -1))  # [n, n-1, ..., 1]
    vert_idx = {v: i for i, v in enumerate(verts)}
    A = [[0] * n for _ in range(n)]
    # Path edges: k -> k+1 in vertex-index terms (n -> n-1 -> ... -> 1)
    for k in range(n - 1):
        A[k][k + 1] = 1  # verts[k] -> verts[k+1]
    # Tile edges
    for i, (x, y) in enumerate(tiles):
        xi, yi = vert_idx[x], vert_idx[y]
        if bits[i] == 0:
            A[xi][yi] = 1  # x -> y (forward)
        else:
            A[yi][xi] = 1  # y -> x (backward)
    return A


def adj_to_tuple(A):
    return tuple(tuple(row) for row in A)


def is_grid_symmetric(bits, trans_map):
    """Check if tiling is symmetric under the grid transpose."""
    for i in range(len(bits)):
        if trans_map[i] != i and bits[i] != bits[trans_map[i]]:
            return False
    return True


def flip_bits(bits):
    """Flip all bits (complement)."""
    return tuple(1 - b for b in bits)


def transpose_bits(bits, trans_map):
    """Apply grid transpose to bits."""
    out = [0] * len(bits)
    for i in range(len(bits)):
        out[trans_map[i]] = bits[i]
    return tuple(out)


# ─────────────────────────────────────────────────────────────
# CANONICALIZATION (isomorphism classes)
# ─────────────────────────────────────────────────────────────

def canonicalize(A, n):
    """Canonical form of adjacency matrix under vertex permutations."""
    perms = list(itertools.permutations(range(n)))
    best = None
    for p in perms:
        s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best


# ─────────────────────────────────────────────────────────────
# ODD CYCLES AND INDEPENDENCE POLYNOMIAL
# ─────────────────────────────────────────────────────────────

def find_directed_cycles(A, n, length):
    """Find all directed cycles of given length. Returns list of frozensets of vertices."""
    cycles = set()
    for subset in itertools.combinations(range(n), length):
        min_v = subset[0]
        rest = list(subset[1:])
        for perm in itertools.permutations(rest):
            path = [min_v] + list(perm)
            valid = True
            for k in range(length):
                if A[path[k]][path[(k + 1) % length]] != 1:
                    valid = False
                    break
            if valid:
                cycles.add(frozenset(subset))
    return list(cycles)


def find_all_odd_cycles(A, n):
    """Find all directed odd cycles in tournament."""
    all_cycles = []
    for length in range(3, n + 1, 2):
        all_cycles.extend(find_directed_cycles(A, n, length))
    return all_cycles


def independence_polynomial(cycles):
    """Compute I(Omega, x) coefficients [alpha_0, alpha_1, ..., alpha_d].
    Omega: conflict graph where cycles adjacent iff sharing a vertex.
    Independent sets: pairwise vertex-disjoint cycle collections."""
    m = len(cycles)
    if m == 0:
        return [1]

    # Build conflict matrix
    conflict = [[False] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if cycles[i] & cycles[j]:
                conflict[i][j] = conflict[j][i] = True

    # Count independent sets by size (up to size 4 for safety)
    alpha = [1, m, 0, 0, 0]

    # alpha_2
    for i in range(m):
        for j in range(i + 1, m):
            if not conflict[i][j]:
                alpha[2] += 1

    # alpha_3
    for i in range(m):
        for j in range(i + 1, m):
            if conflict[i][j]:
                continue
            for k in range(j + 1, m):
                if not conflict[i][k] and not conflict[j][k]:
                    alpha[3] += 1

    # alpha_4
    for i in range(m):
        for j in range(i + 1, m):
            if conflict[i][j]:
                continue
            for k in range(j + 1, m):
                if conflict[i][k] or conflict[j][k]:
                    continue
                for l in range(k + 1, m):
                    if not conflict[i][l] and not conflict[j][l] and not conflict[k][l]:
                        alpha[4] += 1

    # Trim trailing zeros
    while len(alpha) > 1 and alpha[-1] == 0:
        alpha.pop()
    return alpha


def eval_poly(coeffs, x):
    return sum(c * x ** k for k, c in enumerate(coeffs))


def poly_roots_all_real(coeffs):
    """Check if polynomial has all real roots."""
    if len(coeffs) <= 2:
        return True
    poly_c = list(reversed(coeffs))
    roots = np.roots(poly_c)
    return all(abs(r.imag) < 1e-8 for r in roots)


def score_sequence(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))


# ─────────────────────────────────────────────────────────────
# AUTOMORPHISM GROUP
# ─────────────────────────────────────────────────────────────

def automorphism_count(A, n):
    """Count automorphisms of tournament with adjacency matrix A."""
    count = 0
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)):
            count += 1
    return count


# ─────────────────────────────────────────────────────────────
# MAIN INVESTIGATION
# ─────────────────────────────────────────────────────────────

def investigate(n, verbose=True):
    print(f"\n{'=' * 70}")
    print(f"  TILING MODEL ↔ INDEPENDENCE POLYNOMIAL INVESTIGATION: n = {n}")
    print(f"{'=' * 70}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)
    total = 1 << m

    print(f"  Tiles: {m}, Total tilings: {total}")

    # Enumerate all tilings
    data = []
    canon_groups = defaultdict(list)

    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        canon = canonicalize(A, n)
        gs = is_grid_symmetric(bits, trans_map)
        flip_b = flip_bits(bits)
        flip_mask = sum(b << k for k, b in enumerate(flip_b))
        trans_b = transpose_bits(bits, trans_map)
        trans_mask = sum(b << k for k, b in enumerate(trans_b))

        cycles = find_all_odd_cycles(A, n)
        ip = independence_polynomial(cycles)
        H = eval_poly(ip, 2)
        scores = score_sequence(A, n)

        entry = {
            'mask': mask,
            'bits': bits,
            'canon': canon,
            'grid_sym': gs,
            'flip_mask': flip_mask,
            'trans_mask': trans_mask,
            'cycles': cycles,
            'indpoly': ip,
            'H': H,
            'scores': scores,
            'n_3cycles': len([c for c in cycles if len(c) == 3]),
            'n_5cycles': len([c for c in cycles if len(c) == 5]),
            'n_7cycles': len([c for c in cycles if len(c) == 7]),
        }
        data.append(entry)
        canon_groups[canon].append(entry)

    # Assign class indices
    canon_sigs = sorted(canon_groups.keys())
    ci_map = {sig: i for i, sig in enumerate(canon_sigs)}
    for entry in data:
        entry['class_idx'] = ci_map[entry['canon']]

    # For each entry, find flip partner's class
    mask_to_entry = {e['mask']: e for e in data}
    for entry in data:
        flip_entry = mask_to_entry[entry['flip_mask']]
        entry['flip_class'] = flip_entry['class_idx']
        trans_entry = mask_to_entry[entry['trans_mask']]
        entry['trans_class'] = trans_entry['class_idx']

    num_classes = len(canon_sigs)
    print(f"  Isomorphism classes: {num_classes}")

    # ─────────────────────────────────────────────────────────
    # CLASSIFY TILING TYPES
    # ─────────────────────────────────────────────────────────

    # For each tiling, determine its type:
    # self_flip: flip partner in same class
    # cross_flip: flip partner in different class
    # grid_sym: is grid-symmetric
    # Combinations:
    #   blueself:  self_flip AND grid_sym
    #   blackself: self_flip AND NOT grid_sym
    #   bluecross: cross_flip AND grid_sym (both endpoints)
    #   blackcross: cross_flip AND NOT grid_sym (at least one)

    type_counts = Counter()
    for e in data:
        self_flip = (e['class_idx'] == e['flip_class'])
        if self_flip:
            if e['grid_sym']:
                e['type'] = 'blueself'
            else:
                e['type'] = 'blackself'
        else:
            flip_e = mask_to_entry[e['flip_mask']]
            if e['grid_sym'] and flip_e['grid_sym']:
                e['type'] = 'bluecross'
            else:
                e['type'] = 'blackcross'
        type_counts[e['type']] += 1

    print(f"\n  Tiling type distribution:")
    for t in ['blueself', 'blackself', 'bluecross', 'blackcross']:
        print(f"    {t:12s}: {type_counts[t]:5d}  ({100*type_counts[t]/total:.1f}%)")

    # ─────────────────────────────────────────────────────────
    # PER-CLASS ANALYSIS
    # ─────────────────────────────────────────────────────────

    print(f"\n  Per-class analysis:")
    print(f"  {'CI':>3s} {'Size':>5s} {'Scores':>20s} {'H':>5s} {'I.P.':>20s} {'#GS':>4s} {'#BS':>4s} {'#BlS':>4s} {'#BkS':>4s} {'FlipClass':>10s} {'TransClass':>10s}")

    class_data = []
    for ci, sig in enumerate(canon_sigs):
        members = canon_groups[sig]
        size = len(members)
        rep = members[0]
        H = rep['H']
        ip = rep['indpoly']
        scores = rep['scores']

        n_gs = sum(1 for e in members if e['grid_sym'])
        n_bs = sum(1 for e in members if e['type'] == 'blueself')
        n_bls = sum(1 for e in members if e['type'] == 'blueself')
        n_bks = sum(1 for e in members if e['type'] == 'blackself')

        # Flip class (should be consistent within class)
        flip_classes = set(e['flip_class'] for e in members)
        trans_classes = set(e['trans_class'] for e in members)

        cd = {
            'ci': ci, 'size': size, 'H': H, 'ip': ip, 'scores': scores,
            'n_gs': n_gs, 'n_blueself': n_bs, 'n_blackself': n_bks,
            'flip_classes': flip_classes, 'trans_classes': trans_classes,
            'members': members,
            'self_flip': ci in flip_classes,
            'self_transpose': ci in trans_classes,
        }
        class_data.append(cd)

        fc_str = ','.join(str(x) for x in sorted(flip_classes))
        tc_str = ','.join(str(x) for x in sorted(trans_classes))
        print(f"  {ci:3d} {size:5d} {str(scores):>20s} {H:5.0f} {str(ip):>20s} {n_gs:4d} {n_bs:4d} {n_bls:4d} {n_bks:4d} {fc_str:>10s} {tc_str:>10s}")

    # ─────────────────────────────────────────────────────────
    # THEOREM CANDIDATES
    # ─────────────────────────────────────────────────────────

    print(f"\n{'─' * 70}")
    print("  THEOREM CANDIDATES")
    print(f"{'─' * 70}")

    # (1) Are blueself and blackself mutually exclusive within a CLASS?
    print("\n  (1) Blueself vs blackself mutual exclusivity within classes:")
    both_count = 0
    for cd in class_data:
        has_blueself = cd['n_blueself'] > 0
        has_blackself = cd['n_blackself'] > 0
        if has_blueself and has_blackself:
            both_count += 1
            print(f"      Class {cd['ci']}: has BOTH blueself ({cd['n_blueself']}) and blackself ({cd['n_blackself']})")
    if both_count == 0:
        print(f"      *** THEOREM: At n={n}, no class has both blueself and blackself tilings.")
        print(f"      This means blueself and blackself are mutually exclusive AT THE CLASS LEVEL.")
    else:
        print(f"      {both_count} classes have both types.")

    # (1b) Are blueself and blackself mutually exclusive at the TILING level?
    # By definition they are: a tiling is either grid-sym or not, so a self-flip tiling
    # is either blueself OR blackself, never both. But the question is about classes.
    print(f"\n      Note: At the individual tiling level, blueself and blackself are")
    print(f"      TRIVIALLY mutually exclusive (a tiling is either grid-symmetric or not).")
    print(f"      The non-trivial question is whether a CLASS can contain both types.")

    # (2) Does the flip operation preserve I(Omega(T), x)?
    print(f"\n  (2) Does flipping all tiles preserve I(Omega(T), x)?")
    flip_preserves_ip = True
    flip_preserves_H = True
    for e in data:
        flip_e = mask_to_entry[e['flip_mask']]
        if e['indpoly'] != flip_e['indpoly']:
            flip_preserves_ip = False
        if e['H'] != flip_e['H']:
            flip_preserves_H = False
    print(f"      Preserves I(Omega, x) coefficients: {flip_preserves_ip}")
    print(f"      Preserves H = I(Omega, 2): {flip_preserves_H}")

    if not flip_preserves_ip:
        # Analyze how flip changes I.P.
        ip_changes = Counter()
        for e in data:
            flip_e = mask_to_entry[e['flip_mask']]
            change = (tuple(e['indpoly']), tuple(flip_e['indpoly']))
            ip_changes[change] += 1
        print(f"      Distinct (I.P., flip-I.P.) pairs: {len(ip_changes)}")
        if len(ip_changes) <= 20:
            for (ip1, ip2), count in sorted(ip_changes.items()):
                if ip1 != ip2:
                    print(f"        {list(ip1)} <-> {list(ip2)}: {count} tilings")

    # (3) Does transpose preserve I(Omega(T), x)?
    print(f"\n  (3) Does transpose (T -> T^op) preserve I(Omega(T), x)?")
    trans_preserves_ip = True
    for e in data:
        trans_e = mask_to_entry[e['trans_mask']]
        if e['indpoly'] != trans_e['indpoly']:
            trans_preserves_ip = False
            break
    print(f"      Preserves I(Omega, x) coefficients: {trans_preserves_ip}")

    # (4) Relationship between grid-symmetry (self-converse) and I.P.
    print(f"\n  (4) Independence polynomial conditioned on grid-symmetry (self-converse):")
    gs_ips = Counter()
    ngs_ips = Counter()
    for e in data:
        ip_tuple = tuple(e['indpoly'])
        if e['grid_sym']:
            gs_ips[ip_tuple] += 1
        else:
            ngs_ips[ip_tuple] += 1

    print(f"      Grid-symmetric tilings: {sum(gs_ips.values())}")
    print(f"        Distinct I.P.s: {len(gs_ips)}")
    for ip, cnt in sorted(gs_ips.items(), key=lambda x: x[0]):
        print(f"          {list(ip)} (H={eval_poly(ip,2):.0f}): {cnt}")

    print(f"      Non-grid-symmetric tilings: {sum(ngs_ips.values())}")
    print(f"        Distinct I.P.s: {len(ngs_ips)}")
    for ip, cnt in sorted(ngs_ips.items(), key=lambda x: x[0]):
        print(f"          {list(ip)} (H={eval_poly(ip,2):.0f}): {cnt}")

    # (5) Self-flip classes: characterize by I.P.
    print(f"\n  (5) Self-flip classes (flip maps class to itself):")
    self_flip_classes = [cd for cd in class_data if cd['self_flip']]
    cross_flip_classes = [cd for cd in class_data if not cd['self_flip']]
    print(f"      Self-flip classes: {len(self_flip_classes)}")
    for cd in self_flip_classes:
        print(f"        Class {cd['ci']}: size={cd['size']}, H={cd['H']}, I.P.={cd['ip']}, scores={cd['scores']}, "
              f"#blueself={cd['n_blueself']}, #blackself={cd['n_blackself']}, #gs={cd['n_gs']}")
    print(f"      Cross-flip classes: {len(cross_flip_classes)}")

    # (5b) Do self-flip classes have special I.P. structure?
    sf_H = [cd['H'] for cd in self_flip_classes]
    cf_H = [cd['H'] for cd in cross_flip_classes]
    print(f"      Self-flip H values: {sorted(set(sf_H))}")
    print(f"      Cross-flip H values: {sorted(set(cf_H))}")

    # (6) Self-transpose classes
    print(f"\n  (6) Self-transpose classes (transpose maps class to itself):")
    self_trans = [cd for cd in class_data if cd['self_transpose']]
    cross_trans = [cd for cd in class_data if not cd['self_transpose']]
    print(f"      Self-transpose: {len(self_trans)}")
    print(f"      Cross-transpose: {len(cross_trans)}")

    # (7) Relationship: self-flip AND self-transpose
    print(f"\n  (7) Combined: self-flip AND self-transpose:")
    both_self = [cd for cd in class_data if cd['self_flip'] and cd['self_transpose']]
    flip_only = [cd for cd in class_data if cd['self_flip'] and not cd['self_transpose']]
    trans_only = [cd for cd in class_data if not cd['self_flip'] and cd['self_transpose']]
    neither = [cd for cd in class_data if not cd['self_flip'] and not cd['self_transpose']]
    print(f"      Both self-flip and self-transpose: {len(both_self)}")
    print(f"      Self-flip only: {len(flip_only)}")
    print(f"      Self-transpose only: {len(trans_only)}")
    print(f"      Neither: {len(neither)}")

    # (8) H values vs tiling type
    print(f"\n  (8) H value statistics by tiling type:")
    for t in ['blueself', 'blackself', 'bluecross', 'blackcross']:
        H_vals = [e['H'] for e in data if e['type'] == t]
        if H_vals:
            print(f"      {t:12s}: count={len(H_vals)}, H range=[{min(H_vals):.0f}, {max(H_vals):.0f}], "
                  f"mean={np.mean(H_vals):.2f}, distinct H={len(set(H_vals))}")

    # (9) Are blueself tilings always self-converse tournaments?
    print(f"\n  (9) Blueself tilings as self-converse tournaments:")
    blueself_tilings = [e for e in data if e['type'] == 'blueself']
    if blueself_tilings:
        # Grid-symmetric means transpose(bits) = bits, which means T^op has the
        # same tiling bits under the grid transform. Check if T is actually
        # isomorphic to T^op (self-converse).
        sc_count = 0
        for e in blueself_tilings:
            trans_e = mask_to_entry[e['trans_mask']]
            if e['canon'] == trans_e['canon']:
                sc_count += 1
        print(f"      Blueself tilings that are self-converse: {sc_count}/{len(blueself_tilings)}")
        # Note: grid-symmetric does NOT mean T = T^op exactly; it means the
        # tiling bits are symmetric. T^op has ALL arcs reversed, but the path
        # arcs are also reversed, so grid-sym =/= self-converse in general.

    # (10) Disjoint cycle pair count by type
    print(f"\n  (10) Average disjoint cycle pairs (alpha_2) by tiling type:")
    for t in ['blueself', 'blackself', 'bluecross', 'blackcross']:
        a2_vals = [e['indpoly'][2] if len(e['indpoly']) > 2 else 0 for e in data if e['type'] == t]
        if a2_vals:
            print(f"      {t:12s}: mean alpha_2={np.mean(a2_vals):.3f}, "
                  f"max={max(a2_vals)}, zero%={100*a2_vals.count(0)/len(a2_vals):.1f}%")

    # (11) Does I.P. determine the class? (i.e., is I.P. a complete invariant?)
    print(f"\n  (11) Is I(Omega, x) a complete isomorphism invariant?")
    ip_to_classes = defaultdict(set)
    for cd in class_data:
        ip_to_classes[tuple(cd['ip'])].add(cd['ci'])
    collisions = {ip: classes for ip, classes in ip_to_classes.items() if len(classes) > 1}
    if collisions:
        print(f"      NO — {len(collisions)} I.P. values shared by multiple classes:")
        for ip, classes in sorted(collisions.items()):
            class_info = [(ci, class_data[ci]['size'], class_data[ci]['scores']) for ci in sorted(classes)]
            print(f"        I.P.={list(ip)}: classes {[x[0] for x in class_info]}, "
                  f"scores={[x[2] for x in class_info]}")
    else:
        print(f"      YES — I.P. distinguishes all {num_classes} classes at n={n}.")

    # (12) Cycle counts by tiling type
    print(f"\n  (12) Cycle count statistics by tiling type:")
    for t in ['blueself', 'blackself', 'bluecross', 'blackcross']:
        entries = [e for e in data if e['type'] == t]
        if entries:
            c3 = [e['n_3cycles'] for e in entries]
            c5 = [e['n_5cycles'] for e in entries]
            print(f"      {t:12s}: mean c3={np.mean(c3):.2f}, mean c5={np.mean(c5):.2f}")

    # (13) KEY: Does flip preserve the NUMBER of odd cycles?
    print(f"\n  (13) Does flip preserve the total number of odd cycles?")
    flip_preserves_ncycles = True
    flip_preserves_c3 = True
    for e in data:
        flip_e = mask_to_entry[e['flip_mask']]
        if len(e['cycles']) != len(flip_e['cycles']):
            flip_preserves_ncycles = False
        if e['n_3cycles'] != flip_e['n_3cycles']:
            flip_preserves_c3 = False
    print(f"      Preserves total #cycles: {flip_preserves_ncycles}")
    print(f"      Preserves #3-cycles: {flip_preserves_c3}")

    # (14) Flip operation on the tournament level
    print(f"\n  (14) What does flip DO to the tournament?")
    # Flip inverts all non-path edges. Path edges stay: n->n-1->...->1
    # The flip of T has the same Hamiltonian path but all other edges reversed.
    # This is related to T^op but NOT the same (path edges are preserved, not reversed).
    # Let's check: is flip(T) isomorphic to T^op?
    flip_is_top = 0
    flip_not_top = 0
    for e in data:
        flip_e = mask_to_entry[e['flip_mask']]
        trans_e = mask_to_entry[e['trans_mask']]
        if flip_e['canon'] == trans_e['canon']:
            flip_is_top += 1
        else:
            flip_not_top += 1
    print(f"      flip(T) isomorphic to T^op: {flip_is_top}/{total}")
    print(f"      flip(T) NOT isomorphic to T^op: {flip_not_top}/{total}")

    # (15) Flip AND transpose commutation
    print(f"\n  (15) Do flip and transpose commute?")
    commute_count = 0
    for e in data:
        # flip(transpose(T)) vs transpose(flip(T))
        ft = flip_bits(transpose_bits(e['bits'], trans_map))
        tf = transpose_bits(flip_bits(e['bits']), trans_map)
        if ft == tf:
            commute_count += 1
    print(f"      Commute: {commute_count}/{total}")

    return data, class_data, canon_groups


if __name__ == '__main__':
    for n in [4, 5]:
        investigate(n)

    # n=6 is expensive (2^10 = 1024 tilings, but canonicalization is O(n!) per tiling)
    # Let's try it
    print("\n\nAttempting n=6 (1024 tilings, canonicalization expensive)...")
    investigate(6)
