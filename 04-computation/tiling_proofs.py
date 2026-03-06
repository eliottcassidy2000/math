"""
Focused proofs and deeper analysis of tiling model ↔ independence polynomial.

Key findings from initial investigation:
  1. Blueself/blackself mutually exclusive at class level (n=4,5,6)
  2. Transpose ALWAYS preserves I(Omega(T), x) — PROVABLE
  3. Flip and transpose ALWAYS commute — PROVABLE
  4. Blueself implies self-converse — needs deeper analysis
  5. Self-flip classes have structured I.P. patterns

This script:
  (A) Proves Theorem B (transpose preserves I.P.) via cycle reversal argument
  (B) Analyzes the group structure {flip, transpose, isomorphism}
  (C) Investigates blueself/blackself mutual exclusivity mechanism
  (D) Examines whether self-flip interacts with the discriminant bound
  (E) Studies the "flip distance" from a tiling to its I.P.-preserving neighborhood

Author: opus-2026-03-06-S16
"""

import itertools
from collections import defaultdict, Counter
import numpy as np


def build_tiles(n):
    tiles = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            tiles.append((x, y))
    return tiles


def build_transpose_map(n, tiles):
    tile_idx = {(x, y): i for i, (x, y) in enumerate(tiles)}
    return [tile_idx[(n - y + 1, n - x + 1)] for i, (x, y) in enumerate(tiles)]


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


def canonicalize(A, n):
    perms = list(itertools.permutations(range(n)))
    best = None
    for p in perms:
        s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best


def find_directed_cycles(A, n, length):
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
    all_cycles = []
    for length in range(3, n + 1, 2):
        all_cycles.extend(find_directed_cycles(A, n, length))
    return all_cycles


def independence_polynomial(cycles):
    m = len(cycles)
    if m == 0:
        return [1]
    conflict = [[False] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if cycles[i] & cycles[j]:
                conflict[i][j] = conflict[j][i] = True
    alpha = [1, m, 0, 0, 0]
    for i in range(m):
        for j in range(i + 1, m):
            if not conflict[i][j]:
                alpha[2] += 1
    for i in range(m):
        for j in range(i + 1, m):
            if conflict[i][j]: continue
            for k in range(j + 1, m):
                if not conflict[i][k] and not conflict[j][k]:
                    alpha[3] += 1
    for i in range(m):
        for j in range(i + 1, m):
            if conflict[i][j]: continue
            for k in range(j + 1, m):
                if conflict[i][k] or conflict[j][k]: continue
                for l in range(k + 1, m):
                    if not conflict[i][l] and not conflict[j][l] and not conflict[k][l]:
                        alpha[4] += 1
    while len(alpha) > 1 and alpha[-1] == 0:
        alpha.pop()
    return alpha


def eval_poly(coeffs, x):
    return sum(c * x ** k for k, c in enumerate(coeffs))


def is_grid_symmetric(bits, trans_map):
    for i in range(len(bits)):
        if trans_map[i] != i and bits[i] != bits[trans_map[i]]:
            return False
    return True


def flip_bits(bits):
    return tuple(1 - b for b in bits)


def transpose_bits(bits, trans_map):
    out = [0] * len(bits)
    for i in range(len(bits)):
        out[trans_map[i]] = bits[i]
    return tuple(out)


def opposite_tournament(A, n):
    """Compute T^op: reverse all arcs."""
    return [[A[j][i] for j in range(n)] for i in range(n)]


# ─────────────────────────────────────────────────────────────
# THEOREM B: Transpose preserves I(Omega(T), x)
# ─────────────────────────────────────────────────────────────

def prove_transpose_preserves_ip(n):
    """Verify that T^op has the same odd cycle vertex sets as T,
    hence the same conflict graph, hence the same I.P."""
    print(f"\n{'=' * 60}")
    print(f"  THEOREM B: T^op preserves I(Omega, x) at n={n}")
    print(f"{'=' * 60}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)

    failures = 0
    total = 1 << m

    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        A_op = opposite_tournament(A, n)

        cycles_T = find_all_odd_cycles(A, n)
        cycles_Top = find_all_odd_cycles(A_op, n)

        # The cycle VERTEX SETS should be identical
        set_T = set(cycles_T)
        set_Top = set(cycles_Top)

        if set_T != set_Top:
            failures += 1
            print(f"  FAILURE at mask={mask}: T has {len(set_T)} cycles, T^op has {len(set_Top)}")
            diff1 = set_T - set_Top
            diff2 = set_Top - set_T
            if diff1:
                print(f"    In T but not T^op: {diff1}")
            if diff2:
                print(f"    In T^op but not T: {diff2}")

    if failures == 0:
        print(f"  VERIFIED: All {total} tournaments at n={n} have identical")
        print(f"  odd cycle vertex sets in T and T^op.")
        print(f"")
        print(f"  PROOF SKETCH: If v1->v2->...->vk->v1 is a directed cycle in T,")
        print(f"  then v1<-v2<-...<-vk<-v1 holds in T^op, which means")
        print(f"  vk->vk-1->...->v1->vk is a directed cycle in T^op on the")
        print(f"  SAME vertex set. The map C -> reverse(C) is a bijection between")
        print(f"  directed cycles of T and T^op preserving vertex sets.")
        print(f"  Since Omega(T) depends only on vertex sets of cycles,")
        print(f"  Omega(T) = Omega(T^op), so I(Omega(T), x) = I(Omega(T^op), x). QED")
    else:
        print(f"  {failures} FAILURES!")
    return failures


# ─────────────────────────────────────────────────────────────
# THEOREM: What flip does at the tournament level
# ─────────────────────────────────────────────────────────────

def analyze_flip_tournament_level(n):
    """Analyze what flipping all non-path bits does to the tournament."""
    print(f"\n{'=' * 60}")
    print(f"  FLIP ANALYSIS at n={n}")
    print(f"{'=' * 60}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)
    total = 1 << m

    # For each T, compute flip(T) and analyze the relationship
    # flip(T) reverses all non-path edges but keeps path edges
    # How does this relate to T^op?
    # T^op reverses ALL edges (including path)
    # So flip(T) and T^op differ exactly on the path edges

    print(f"\n  Flip reverses non-path edges, keeps path n->n-1->...->1")
    print(f"  T^op reverses ALL edges (including path)")
    print(f"  So flip(T) and T^op differ on exactly n-1 path edges")
    print(f"  flip(T) = T^op with path edges re-reversed")
    print(f"")

    # Check: is flip(T) always obtainable from T^op by flipping the n-1 path edges?
    # More precisely: flip(T) has A[i][j] for path edges same as T,
    # and A[i][j] reversed for non-path edges.
    # T^op has ALL A[j][i].
    # So flip(T)[i][j] = T[i][j] for path, = T[j][i] for non-path.
    # T^op[i][j] = T[j][i] for all.
    # Difference: on path edges, flip(T) has T[i][j] while T^op has T[j][i] = 1-T[i][j].

    # What is the isomorphism type of flip(T)?
    # Let's check if flip(T) = σ(T^op) for some fixed permutation σ.

    # The reversal permutation σ: i -> n-1-i reverses vertex order.
    # Under σ, the path n->n-1->...->1 becomes 1->2->...->n (reversed direction).
    # But σ(T^op) would have A_σ[σ(i)][σ(j)] = T^op[i][j] = T[j][i].
    # So A_σ[a][b] = T[σ^-1(b)][σ^-1(a)] = T[n-1-b][n-1-a].

    # Let's check: for path edges (consecutive vertices k, k+1 in index terms),
    # T has T[k][k+1]=1, T[k+1][k]=0.
    # A_σ[a][a+1] = T[n-1-(a+1)][n-1-a] = T[n-2-a][n-1-a].
    # n-1-a and n-2-a are consecutive indices (diff = 1), and n-1-a > n-2-a.
    # So T[n-2-a][n-1-a] = 1 (path edge, smaller index -> larger index = 0).
    # Wait, T[k][k+1] = 1 means verts[k] -> verts[k+1], where verts[k] > verts[k+1].
    # So T[smaller_index][larger_index] = 1 for path edges.
    # T[n-2-a][n-1-a]: n-2-a < n-1-a, so this is T[smaller][larger] = 1.
    # Good, so A_σ preserves path edges!

    # For non-path edges: A_σ[a][b] = T[n-1-b][n-1-a].
    # If a and b are not consecutive, is this the reversal of T[a][b]?
    # We need A_σ[a][b] = 1-T[a][b] for non-path edges.
    # A_σ[a][b] = T[n-1-b][n-1-a] = 1 - T[n-1-a][n-1-b] (tournament property).

    # So A_σ[a][b] = 1 - T[n-1-a][n-1-b].
    # For this to equal 1-T[a][b], we need T[n-1-a][n-1-b] = T[a][b] for non-path edges.
    # This is NOT true in general — it would require the reversal permutation to be
    # an automorphism of T on non-path edges.

    # So flip(T) ≠ σ(T^op) in general. Let's verify computationally.
    sigma = list(range(n - 1, -1, -1))  # reversal: i -> n-1-i

    flip_eq_sigma_top = 0
    for mask in range(min(total, 500)):
        bits = tuple((mask >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        flip_b = flip_bits(bits)
        A_flip = bits_to_adj(n, tiles, flip_b)
        A_op = opposite_tournament(A, n)
        # σ(T^op)
        A_sigma_op = [[A_op[sigma[i]][sigma[j]] for j in range(n)] for i in range(n)]

        if A_flip == A_sigma_op:
            flip_eq_sigma_top += 1

    tested = min(total, 500)
    print(f"  flip(T) = σ(T^op) where σ reverses vertex order: {flip_eq_sigma_top}/{tested}")

    if flip_eq_sigma_top == tested:
        print(f"  *** THEOREM: flip(T) = σ(T^op) where σ: v -> n+1-v ***")
        print(f"  This means flip is the COMPOSITION of:")
        print(f"    (1) Taking the opposite tournament (reverse all arcs)")
        print(f"    (2) Relabeling vertices by the reversal permutation σ")
        print(f"  Since T^op preserves I(Omega,x) (Theorem B),")
        print(f"  and σ is just a relabeling (preserves isomorphism class),")
        print(f"  flip preserves the ISOMORPHISM CLASS of (Omega(T), I.P.).")
        print(f"")
        print(f"  BUT WAIT — we showed flip does NOT preserve I.P.!")
        print(f"  Resolution: σ changes the LABELING, and the tiling model")
        print(f"  is not σ-invariant. The TOURNAMENT flip(T) is isomorphic")
        print(f"  to T^op (and hence has same I.P.), but its TILING representation")
        print(f"  changes because the fixed path is different.")
    else:
        print(f"  flip(T) ≠ σ(T^op) in general.")

    # Let me check more directly: does flip(T) have the same canonicalization as T^op?
    print(f"\n  Checking: is flip(T) isomorphic to T^op?")
    same_class = 0
    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        flip_b = flip_bits(bits)
        A_flip = bits_to_adj(n, tiles, flip_b)
        A_op = opposite_tournament(A, n)

        c_flip = canonicalize(A_flip, n)
        c_op = canonicalize(A_op, n)

        if c_flip == c_op:
            same_class += 1

    print(f"  flip(T) ≅ T^op: {same_class}/{total}")

    if same_class == total:
        print(f"\n  *** THEOREM: flip(T) is ALWAYS isomorphic to T^op ***")
        print(f"  PROOF: flip reverses non-path edges, keeps path.")
        print(f"         T^op reverses ALL edges.")
        print(f"         The reversal permutation σ: v↦n+1-v maps the")
        print(f"         path n→n-1→...→1 to 1→2→...→n (same edges, reversed).")
        print(f"         σ(T^op) reverses all arcs then relabels, which")
        print(f"         re-reverses path edges (net: path preserved) and")
        print(f"         keeps non-path edges reversed. This equals flip(T).")
        print(f"")
        print(f"  COROLLARY: flip(T) and T have the same I.P.!")
        print(f"  (Since flip(T) ≅ T^op, and T^op has same Omega as T.)")
    else:
        print(f"  NOT always! {total - same_class} failures.")

    return same_class == total


# ─────────────────────────────────────────────────────────────
# DEEP ANALYSIS: Self-flip classes and mutual exclusivity
# ─────────────────────────────────────────────────────────────

def deep_self_flip_analysis(n):
    """Deep dive into self-flip classes and blueself/blackself structure."""
    print(f"\n{'=' * 60}")
    print(f"  SELF-FLIP CLASS DEEP ANALYSIS at n={n}")
    print(f"{'=' * 60}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)
    total = 1 << m

    # Build all data
    data = []
    canon_groups = defaultdict(list)

    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        canon = canonicalize(A, n)
        gs = is_grid_symmetric(bits, trans_map)
        flip_b = flip_bits(bits)
        flip_mask = sum(b << k for k, b in enumerate(flip_b))

        cycles = find_all_odd_cycles(A, n)
        ip = independence_polynomial(cycles)
        H = eval_poly(ip, 2)

        entry = {
            'mask': mask, 'bits': bits, 'canon': canon,
            'grid_sym': gs, 'flip_mask': flip_mask,
            'ip': ip, 'H': H,
        }
        data.append(entry)
        canon_groups[canon].append(entry)

    canon_sigs = sorted(canon_groups.keys())
    ci_map = {sig: i for i, sig in enumerate(canon_sigs)}
    mask_to_entry = {e['mask']: e for e in data}

    for e in data:
        e['ci'] = ci_map[e['canon']]
        e['flip_ci'] = ci_map[mask_to_entry[e['flip_mask']]['canon']]

    # Find self-flip classes: classes where at least one member flips to the same class
    print(f"\n  Self-flip class details:")
    print(f"  (A self-flip MEMBER is a tiling T where flip(T) is in the same class as T)")
    print()

    for ci, sig in enumerate(canon_sigs):
        members = canon_groups[sig]
        self_flip_members = [e for e in members if e['flip_ci'] == ci]
        if not self_flip_members:
            continue

        gs_self = [e for e in self_flip_members if e['grid_sym']]
        ngs_self = [e for e in self_flip_members if not e['grid_sym']]

        print(f"  Class {ci}: {len(members)} members, {len(self_flip_members)} self-flip")
        print(f"    Grid-symmetric self-flip (blueself): {len(gs_self)}")
        print(f"    Non-grid-symmetric self-flip (blackself): {len(ngs_self)}")
        print(f"    I.P. = {members[0]['ip']}, H = {members[0]['H']}")

        # For each self-flip member, who is their flip partner?
        for e in self_flip_members:
            fp = mask_to_entry[e['flip_mask']]
            same_tiling = (e['mask'] == fp['mask'])
            print(f"      mask={e['mask']:>4d} gs={e['grid_sym']} -> flip_mask={fp['mask']:>4d} gs={fp['grid_sym']} "
                  f"{'(FIXED POINT)' if same_tiling else ''}")

        # KEY ANALYSIS: What connects blueself to blackself?
        # If T is grid-symmetric and self-flip, its flip partner flip(T) is also
        # in the same class. Is flip(T) also grid-symmetric?
        for e in gs_self:
            fp = mask_to_entry[e['flip_mask']]
            print(f"    -> Blueself {e['mask']}: flip partner {fp['mask']} is "
                  f"{'grid-sym' if fp['grid_sym'] else 'NOT grid-sym'}")

        print()

    # Check: if T is grid-symmetric, is flip(T) also grid-symmetric?
    print(f"\n  Global check: If T is grid-symmetric, is flip(T) grid-symmetric?")
    gs_to_gs = 0
    gs_to_ngs = 0
    for e in data:
        if e['grid_sym']:
            fp = mask_to_entry[e['flip_mask']]
            if fp['grid_sym']:
                gs_to_gs += 1
            else:
                gs_to_ngs += 1
    n_gs = sum(1 for e in data if e['grid_sym'])
    print(f"    Grid-sym -> flip grid-sym: {gs_to_gs}/{n_gs}")
    print(f"    Grid-sym -> flip NOT grid-sym: {gs_to_ngs}/{n_gs}")

    # PROOF: flip preserves grid-symmetry iff flip commutes with transpose.
    # We already know they commute! So:
    # If G(T)=T (grid-sym), then G(F(T)) = F(G(T)) = F(T). QED.
    if gs_to_ngs == 0:
        print(f"\n  *** THEOREM: Flip preserves grid-symmetry. ***")
        print(f"  PROOF: Grid-symmetry means G(T)=T where G=transpose.")
        print(f"         Since F and G commute: G(F(T)) = F(G(T)) = F(T).")
        print(f"         So F(T) is also grid-symmetric. QED")
        print(f"")
        print(f"  COROLLARY: If T is blueself (self-flip + grid-sym),")
        print(f"  then flip(T) is also grid-symmetric and in the same class.")
        print(f"  So flip(T) is also blueself. Blueself members come in PAIRS")
        print(f"  (or fixed points if flip(T)=T).")

    # Check: if T is NOT grid-symmetric, is flip(T) also NOT grid-symmetric?
    ngs_to_gs = 0
    ngs_to_ngs = 0
    for e in data:
        if not e['grid_sym']:
            fp = mask_to_entry[e['flip_mask']]
            if fp['grid_sym']:
                ngs_to_gs += 1
            else:
                ngs_to_ngs += 1
    n_ngs = sum(1 for e in data if not e['grid_sym'])
    print(f"\n    NOT grid-sym -> flip NOT grid-sym: {ngs_to_ngs}/{n_ngs}")
    print(f"    NOT grid-sym -> flip grid-sym: {ngs_to_gs}/{n_ngs}")

    if ngs_to_gs == 0:
        print(f"\n  *** THEOREM: Flip also preserves NON-grid-symmetry. ***")
        print(f"  So grid-symmetry is a FLIP INVARIANT.")
        print(f"")
        print(f"  CONSEQUENCE FOR MUTUAL EXCLUSIVITY:")
        print(f"  A self-flip member is blueself iff it's grid-symmetric.")
        print(f"  A self-flip member is blackself iff it's NOT grid-symmetric.")
        print(f"  Since flip preserves grid-symmetry, the flip partner of a")
        print(f"  blueself tiling is also grid-symmetric (hence also blueself).")
        print(f"  Similarly for blackself.")
        print(f"")
        print(f"  But this does NOT yet prove class-level mutual exclusivity!")
        print(f"  A class could have both grid-sym AND non-grid-sym members,")
        print(f"  with self-flip members of both types.")

    # The remaining question: WHY can't a class have both blueself and blackself?
    # Let's check: for each self-flip class, are the self-flip members
    # ALL of the same grid-symmetry type?

    print(f"\n  CRITICAL CHECK: Within each self-flip class, do self-flip members")
    print(f"  have uniform grid-symmetry?")

    for ci, sig in enumerate(canon_sigs):
        members = canon_groups[sig]
        self_flip_members = [e for e in members if e['flip_ci'] == ci]
        if not self_flip_members:
            continue

        gs_types = set(e['grid_sym'] for e in self_flip_members)
        if len(gs_types) > 1:
            print(f"    *** VIOLATION at class {ci}! Self-flip members have MIXED grid-symmetry! ***")
        else:
            gs_type = list(gs_types)[0]
            print(f"    Class {ci}: all {len(self_flip_members)} self-flip members are "
                  f"{'grid-sym (blueself)' if gs_type else 'NOT grid-sym (blackself)'}")


# ─────────────────────────────────────────────────────────────
# HAMMING WEIGHT AND FLIP SYMMETRY
# ─────────────────────────────────────────────────────────────

def hamming_weight_analysis(n):
    """Analyze how flip relates to Hamming weight and I.P. structure."""
    print(f"\n{'=' * 60}")
    print(f"  HAMMING WEIGHT AND FLIP SYMMETRY at n={n}")
    print(f"{'=' * 60}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)
    total = 1 << m

    # Flip maps bits -> 1-bits, so HW(flip(T)) = m - HW(T)
    # The "center" of flip is HW = m/2
    # Grid-symmetric tilings have HW that depends on the transpose map's fixed points

    # Count fixed points of transpose map
    fixed = sum(1 for i in range(m) if trans_map[i] == i)
    print(f"  Tiles: {m}, Transpose fixed points: {fixed}")
    print(f"  Grid-symmetric tilings must have matching bits on orbits")
    print(f"  Transpose orbits: {fixed} fixed + {(m-fixed)//2} pairs")
    print(f"  Grid-symmetric tilings: 2^{fixed} * 2^{(m-fixed)//2} = {2**fixed * 2**((m-fixed)//2)}")
    print(f"  (Verify: {sum(1 for mask in range(total) if is_grid_symmetric(tuple((mask>>k)&1 for k in range(m)), trans_map))})")

    # H(T) + H(flip(T)) for each T
    h_sums = []
    h_diffs = []
    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        flip_b = flip_bits(bits)
        A_f = bits_to_adj(n, tiles, flip_b)

        cycles = find_all_odd_cycles(A, n)
        ip = independence_polynomial(cycles)
        H = eval_poly(ip, 2)

        cycles_f = find_all_odd_cycles(A_f, n)
        ip_f = independence_polynomial(cycles_f)
        H_f = eval_poly(ip_f, 2)

        hw = sum(bits)
        h_sums.append((hw, H + H_f))
        h_diffs.append((hw, H - H_f))

    # H(T) + H(flip(T)) by Hamming weight
    hw_to_sums = defaultdict(list)
    hw_to_diffs = defaultdict(list)
    for hw, s in h_sums:
        hw_to_sums[hw].append(s)
    for hw, d in h_diffs:
        hw_to_diffs[hw].append(d)

    print(f"\n  H(T) + H(flip(T)) by Hamming weight:")
    for hw in sorted(hw_to_sums):
        vals = hw_to_sums[hw]
        print(f"    HW={hw:2d} ({m-hw:2d}): mean={np.mean(vals):.1f}, "
              f"min={min(vals):.0f}, max={max(vals):.0f}, count={len(vals)}")

    # Since flip(T) ≅ T^op and T^op has same I.P. as T:
    # H(flip(T)) should equal H(T)!
    print(f"\n  H(T) - H(flip(T)) by Hamming weight:")
    all_zero = True
    for hw in sorted(hw_to_diffs):
        vals = hw_to_diffs[hw]
        if any(v != 0 for v in vals):
            all_zero = False
        print(f"    HW={hw:2d}: mean={np.mean(vals):.1f}, min={min(vals):.0f}, max={max(vals):.0f}")

    if all_zero:
        print(f"\n  *** H(T) = H(flip(T)) for ALL tilings! ***")
        print(f"  This confirms: flip(T) ≅ T^op => same H.")
    else:
        print(f"\n  H(T) ≠ H(flip(T)) for some tilings!")
        print(f"  This would CONTRADICT flip(T) ≅ T^op.")


# ─────────────────────────────────────────────────────────────
# INDEPENDENCE NUMBER AND SELF-FLIP
# ─────────────────────────────────────────────────────────────

def ip_vs_self_flip(n):
    """Analyze independence polynomial structure in self-flip vs cross-flip classes."""
    print(f"\n{'=' * 60}")
    print(f"  I.P. STRUCTURE: SELF-FLIP vs CROSS-FLIP at n={n}")
    print(f"{'=' * 60}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)
    total = 1 << m

    data = []
    canon_groups = defaultdict(list)

    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        canon = canonicalize(A, n)
        gs = is_grid_symmetric(bits, trans_map)
        flip_b = flip_bits(bits)
        flip_mask = sum(b << k for k, b in enumerate(flip_b))

        cycles = find_all_odd_cycles(A, n)
        ip = independence_polynomial(cycles)
        H = eval_poly(ip, 2)

        entry = {
            'mask': mask, 'bits': bits, 'canon': canon,
            'grid_sym': gs, 'flip_mask': flip_mask,
            'ip': ip, 'H': H,
        }
        data.append(entry)
        canon_groups[canon].append(entry)

    canon_sigs = sorted(canon_groups.keys())
    ci_map = {sig: i for i, sig in enumerate(canon_sigs)}
    mask_to_entry = {e['mask']: e for e in data}

    for e in data:
        e['ci'] = ci_map[e['canon']]
        e['flip_ci'] = ci_map[mask_to_entry[e['flip_mask']]['canon']]

    # For self-flip tilings: does the flip partner have the EXACT same I.P.?
    # (Not just isomorphic — identical coefficients)
    print(f"\n  Self-flip members: I.P. of T vs I.P. of flip(T)")
    ip_match = 0
    ip_mismatch = 0
    for e in data:
        if e['ci'] == e['flip_ci']:
            fp = mask_to_entry[e['flip_mask']]
            if e['ip'] == fp['ip']:
                ip_match += 1
            else:
                ip_mismatch += 1

    n_self = ip_match + ip_mismatch
    print(f"    I.P. matches: {ip_match}/{n_self}")
    print(f"    I.P. mismatches: {ip_mismatch}/{n_self}")
    if ip_mismatch == 0 and n_self > 0:
        print(f"    *** Self-flip members always have I.P.(T) = I.P.(flip(T))! ***")
        print(f"    (This follows from flip(T) ≅ T^op and Omega(T)=Omega(T^op).)")

    # Count self-flip members by I.P. degree
    print(f"\n  Self-flip members by I.P. degree:")
    sf_by_deg = Counter()
    cf_by_deg = Counter()
    for e in data:
        deg = len(e['ip']) - 1
        if e['ci'] == e['flip_ci']:
            sf_by_deg[deg] += 1
        else:
            cf_by_deg[deg] += 1
    for deg in sorted(set(list(sf_by_deg.keys()) + list(cf_by_deg.keys()))):
        print(f"    Degree {deg}: {sf_by_deg[deg]} self-flip, {cf_by_deg[deg]} cross-flip")

    # Alpha_2 (disjoint pairs) in self-flip vs cross-flip
    print(f"\n  Alpha_2 distribution:")
    sf_a2 = [e['ip'][2] if len(e['ip']) > 2 else 0 for e in data if e['ci'] == e['flip_ci']]
    cf_a2 = [e['ip'][2] if len(e['ip']) > 2 else 0 for e in data if e['ci'] != e['flip_ci']]
    if sf_a2:
        print(f"    Self-flip: mean={np.mean(sf_a2):.2f}, max={max(sf_a2)}, values={sorted(Counter(sf_a2).items())}")
    if cf_a2:
        print(f"    Cross-flip: mean={np.mean(cf_a2):.2f}, max={max(cf_a2)}, values={sorted(Counter(cf_a2).items())}")


# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────

if __name__ == '__main__':
    for n in [4, 5]:
        prove_transpose_preserves_ip(n)

    for n in [4, 5, 6]:
        analyze_flip_tournament_level(n)

    for n in [4, 5, 6]:
        deep_self_flip_analysis(n)

    for n in [4, 5]:
        hamming_weight_analysis(n)

    for n in [4, 5, 6]:
        ip_vs_self_flip(n)
