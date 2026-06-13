"""
Q-009: Prove OCF at n=7 via exhaustive arc-flip identity verification.

Strategy: For i=0, j=1 (WLOG), verify that for ALL 2^20 arc configurations
of the remaining arcs, the swap involution identity holds:

  U_{T'} - U_T = 2*sum(s_x * H(B_x)) + 2*(gained_5 - lost_5) + 2*(gained_7 - lost_7)

Combined with base case (transitive tournament) and arc-flip reachability,
this proves H(T) = I(Omega(T), 2) for all n=7 tournaments.

Optimizations:
- Precompute path shapes as bitmask conditions
- Precompute H(B_x) lookup tables (64 entries each)
- Use numpy vectorized operations over all 2^20 assignments

Author: opus-2026-03-05-S4
"""

import numpy as np
from itertools import permutations, combinations
import time


def build_arc_encoding():
    """Build mapping from arc (a,b) to (bit_position, is_positive).

    T[a][b] = 1 iff (assignment >> bit_pos) & 1 == (1 if is_positive else 0).

    Variables (20 total):
    - bits 0-4: T[x][0] for x in {2,3,4,5,6}
    - bits 5-9: T[1][x] for x in {2,3,4,5,6}
    - bits 10-19: T[a][b] for a<b in {2,3,4,5,6}
    """
    others = [2, 3, 4, 5, 6]
    arc_bit = {}
    bit_idx = 0

    for x in others:
        arc_bit[(x, 0)] = (bit_idx, True)
        arc_bit[(0, x)] = (bit_idx, False)
        bit_idx += 1

    for x in others:
        arc_bit[(1, x)] = (bit_idx, True)
        arc_bit[(x, 1)] = (bit_idx, False)
        bit_idx += 1

    for idx_a, a in enumerate(others):
        for b in others[idx_a+1:]:
            arc_bit[(a, b)] = (bit_idx, True)
            arc_bit[(b, a)] = (bit_idx, False)
            bit_idx += 1

    assert bit_idx == 20
    return arc_bit


def get_arc_masks(a, b, arc_bit, for_T=True):
    """Get (ones_bit, zeros_bit) for requiring T[a][b]=1.

    Returns (bit_position, need_one) or (None, always_value) for fixed arcs.
    """
    if for_T:
        if a == 0 and b == 1: return None, True
        if a == 1 and b == 0: return None, False
    else:
        if a == 0 and b == 1: return None, False
        if a == 1 and b == 0: return None, True

    bit_pos, positive = arc_bit[(a, b)]
    return bit_pos, positive


def precompute_shapes(n, arc_bit, for_T=True):
    """Precompute all path shapes as bitmask conditions.

    Returns list of (ones_mask, zeros_mask, block_conditions).
    block_conditions = list of (bit_pos, blocks_when_one) for blocking checks.
    """
    i_v, j_v = (0, 1) if for_T else (1, 0)
    others = [2, 3, 4, 5, 6]
    shapes = []

    for pos in range(n - 1):
        remaining_positions = [k for k in range(n) if k != pos and k != pos + 1]
        for perm in permutations(others):
            path = [None] * n
            path[pos] = i_v
            path[pos + 1] = j_v
            for idx, rp in enumerate(remaining_positions):
                path[rp] = perm[idx]

            ones_mask = 0
            zeros_mask = 0
            impossible = False

            for k in range(n - 1):
                a, b = path[k], path[k+1]
                bit_pos, need = get_arc_masks(a, b, arc_bit, for_T)
                if bit_pos is None:
                    if not need:
                        impossible = True
                        break
                else:
                    if need:
                        ones_mask |= (1 << bit_pos)
                    else:
                        zeros_mask |= (1 << bit_pos)

            if impossible or (ones_mask & zeros_mask):
                continue

            # Blocking conditions
            # For T-paths: pred blocks if T[pred][j]=0, succ blocks if T[i][succ]=0
            # For T'-paths: pred blocks if T'[pred][i]=0, succ blocks if T'[j][succ]=0
            block_conds = []

            if for_T:
                if pos > 0:
                    pred = path[pos - 1]
                    # blocks if T[pred][1] = 0
                    bit_pos, positive = arc_bit[(pred, 1)]
                    # T[pred][1]=0 means: bit doesn't match positive
                    # blocks_when_one = not positive (if bit=1 and positive=False, T=0, blocks)
                    block_conds.append((bit_pos, not positive))
                if pos + 1 < n - 1:
                    succ = path[pos + 2]
                    # blocks if T[0][succ] = 0
                    bit_pos, positive = arc_bit[(0, succ)]
                    block_conds.append((bit_pos, not positive))
            else:
                # T' paths: j=1 before i=0
                if pos > 0:
                    pred = path[pos - 1]
                    # blocks if T'[pred][0] = 0
                    # T'[pred][0] = T[pred][0] (not involving i-j arc directly)
                    bit_pos, positive = arc_bit[(pred, 0)]
                    block_conds.append((bit_pos, not positive))
                if pos + 1 < n - 1:
                    succ = path[pos + 2]
                    # blocks if T'[1][succ] = 0
                    bit_pos, positive = arc_bit[(1, succ)]
                    block_conds.append((bit_pos, not positive))

            shapes.append((ones_mask, zeros_mask, block_conds))

    return shapes


def compute_unmatched_vectorized(shapes, num_vars=20):
    """Compute U (unmatched count) for all 2^num_vars assignments using numpy."""
    N = 1 << num_vars
    assignments = np.arange(N, dtype=np.int32)
    unmatched = np.zeros(N, dtype=np.int32)

    for ones_mask, zeros_mask, block_conds in shapes:
        check_mask = ones_mask | zeros_mask
        valid = (assignments & check_mask) == ones_mask

        if not block_conds:
            continue  # no blocking possible, these paths are always matched

        blocked = np.zeros(N, dtype=bool)
        for bit_pos, blocks_when_one in block_conds:
            bit_val = (assignments >> bit_pos) & 1
            if blocks_when_one:
                blocked |= (bit_val == 1)
            else:
                blocked |= (bit_val == 0)

        unmatched += (valid & blocked).astype(np.int32)

    return unmatched


def precompute_h4_lookup(verts, arc_bit):
    """Precompute H(T[verts]) for all 2^6 configurations of arcs among verts.

    verts is a list of 4 vertices from {2,3,4,5,6}.
    The 6 arcs among these 4 vertices correspond to specific bit positions.
    Returns (sub_bits, lookup) where:
    - sub_bits: list of (bit_pos, pair_index) for extracting sub-mask
    - lookup: array of 64 values, indexed by sub-mask
    """
    pairs = []
    for idx_a, a in enumerate(verts):
        for b in verts[idx_a+1:]:
            pairs.append((a, b))
    assert len(pairs) == 6

    # Map each pair to its bit position in the global assignment
    sub_bits = []
    for a, b in pairs:
        bit_pos, positive = arc_bit[(a, b)]
        sub_bits.append((bit_pos, positive))

    # For each of 64 sub-assignments, compute H
    lookup = np.zeros(64, dtype=np.int32)
    for sub_mask in range(64):
        # Build adjacency
        T = {}
        for idx, (a, b) in enumerate(pairs):
            if (sub_mask >> idx) & 1:
                T[(a, b)] = 1
                T[(b, a)] = 0
            else:
                T[(a, b)] = 0
                T[(b, a)] = 1

        # Count Hamiltonian paths
        h = 0
        for perm in permutations(verts):
            valid = True
            for k in range(3):
                if T.get((perm[k], perm[k+1]), 0) != 1:
                    valid = False
                    break
            if valid:
                h += 1
        lookup[sub_mask] = h

    return sub_bits, lookup


def extract_sub_mask(assignments, sub_bits):
    """Extract sub-mask for the 6 among-others arcs from global assignments."""
    result = np.zeros_like(assignments)
    for idx, (bit_pos, positive) in enumerate(sub_bits):
        bit_val = (assignments >> bit_pos) & 1
        if not positive:
            bit_val = 1 - bit_val
        result |= (bit_val << idx)
    return result


def compute_formula_vectorized(arc_bit, num_vars=20):
    """Compute delta_I formula for all 2^num_vars assignments.

    delta_I = 2*sum(s_x * H(B_x)) + 2*(gained_5 - lost_5) + 2*(gained_7 - lost_7)
    Convention: H(T') - H(T)
    """
    N = 1 << num_vars
    others = [2, 3, 4, 5, 6]
    assignments = np.arange(N, dtype=np.int32)

    # Compute s_x for each x in others
    # s_x = 1 - T[x][0] - T[1][x]
    s = {}
    for x in others:
        bit_x0, pos_x0 = arc_bit[(x, 0)]
        bit_1x, pos_1x = arc_bit[(1, x)]
        T_x0 = (assignments >> bit_x0) & 1
        if not pos_x0:
            T_x0 = 1 - T_x0
        T_1x = (assignments >> bit_1x) & 1
        if not pos_1x:
            T_1x = 1 - T_1x
        s[x] = (1 - T_x0 - T_1x).astype(np.int32)

    # Compute H(B_x) for each x using lookup tables
    formula_sum = np.zeros(N, dtype=np.int32)
    for x in others:
        B_x = [v for v in others if v != x]
        sub_bits, lookup = precompute_h4_lookup(B_x, arc_bit)
        sub_mask = extract_sub_mask(assignments, sub_bits)
        h_bx = lookup[sub_mask]
        formula_sum += s[x] * h_bx

    # 5-cycle counts
    # Lost 5-cycle: (0,1,v1,v2,v3,0) needs T[1][v1], T[v1][v2], T[v2][v3], T[v3][0]
    # Gained 5-cycle: (1,0,v1,v2,v3,1) needs T[0][v1], T[v1][v2], T[v2][v3], T[v3][1]
    lost_5 = np.zeros(N, dtype=np.int32)
    gained_5 = np.zeros(N, dtype=np.int32)

    for subset in combinations(others, 3):
        for perm in permutations(subset):
            v1, v2, v3 = perm

            # Lost: check T[1][v1], T[v1][v2], T[v2][v3], T[v3][0]
            arcs_lost = [(1, v1), (v1, v2), (v2, v3), (v3, 0)]
            valid_lost = np.ones(N, dtype=bool)
            for a, b in arcs_lost:
                bit_pos, positive = arc_bit[(a, b)]
                bit_val = (assignments >> bit_pos) & 1
                if positive:
                    valid_lost &= (bit_val == 1)
                else:
                    valid_lost &= (bit_val == 0)
            lost_5 += valid_lost.astype(np.int32)

            # Gained: check T[0][v1], T[v1][v2], T[v2][v3], T[v3][1]
            arcs_gained = [(0, v1), (v1, v2), (v2, v3), (v3, 1)]
            valid_gained = np.ones(N, dtype=bool)
            for a, b in arcs_gained:
                bit_pos, positive = arc_bit[(a, b)]
                bit_val = (assignments >> bit_pos) & 1
                if positive:
                    valid_gained &= (bit_val == 1)
                else:
                    valid_gained &= (bit_val == 0)
            gained_5 += valid_gained.astype(np.int32)

    # 7-cycle counts
    # Lost 7-cycle: (0,1,v1,v2,v3,v4,v5,0) needs T[1][v1], ..., T[v5][0]
    # Gained 7-cycle: (1,0,v1,v2,v3,v4,v5,1) needs T[0][v1], ..., T[v5][1]
    lost_7 = np.zeros(N, dtype=np.int32)
    gained_7 = np.zeros(N, dtype=np.int32)

    for perm in permutations(others):  # 120 permutations
        v1, v2, v3, v4, v5 = perm

        # Lost: T[1][v1], T[v1][v2], T[v2][v3], T[v3][v4], T[v4][v5], T[v5][0]
        arcs_lost = [(1, v1), (v1, v2), (v2, v3), (v3, v4), (v4, v5), (v5, 0)]
        valid_lost = np.ones(N, dtype=bool)
        for a, b in arcs_lost:
            bit_pos, positive = arc_bit[(a, b)]
            bit_val = (assignments >> bit_pos) & 1
            if positive:
                valid_lost &= (bit_val == 1)
            else:
                valid_lost &= (bit_val == 0)
        lost_7 += valid_lost.astype(np.int32)

        # Gained: T[0][v1], T[v1][v2], T[v2][v3], T[v3][v4], T[v4][v5], T[v5][1]
        arcs_gained = [(0, v1), (v1, v2), (v2, v3), (v3, v4), (v4, v5), (v5, 1)]
        valid_gained = np.ones(N, dtype=bool)
        for a, b in arcs_gained:
            bit_pos, positive = arc_bit[(a, b)]
            bit_val = (assignments >> bit_pos) & 1
            if positive:
                valid_gained &= (bit_val == 1)
            else:
                valid_gained &= (bit_val == 0)
        gained_7 += valid_gained.astype(np.int32)

    delta_I = 2 * formula_sum + 2 * (gained_5 - lost_5) + 2 * (gained_7 - lost_7)
    return delta_I


if __name__ == "__main__":
    print("=== n=7 Exhaustive Proof of OCF ===")
    print(f"Checking 2^20 = {1<<20} arc configurations...\n")

    t0 = time.time()
    arc_bit = build_arc_encoding()

    # Precompute shapes
    print("Precomputing T-path shapes...", end=" ", flush=True)
    t1 = time.time()
    shapes_T = precompute_shapes(7, arc_bit, for_T=True)
    print(f"{len(shapes_T)} shapes in {time.time()-t1:.1f}s")

    print("Precomputing T'-path shapes...", end=" ", flush=True)
    t1 = time.time()
    shapes_Tp = precompute_shapes(7, arc_bit, for_T=False)
    print(f"{len(shapes_Tp)} shapes in {time.time()-t1:.1f}s")

    # Compute U_T for all assignments
    print("\nComputing U_T (unmatched T-paths)...", flush=True)
    t1 = time.time()
    U_T = compute_unmatched_vectorized(shapes_T)
    print(f"  Done in {time.time()-t1:.1f}s")

    print("Computing U_T' (unmatched T'-paths)...", flush=True)
    t1 = time.time()
    U_Tp = compute_unmatched_vectorized(shapes_Tp)
    print(f"  Done in {time.time()-t1:.1f}s")

    # Compute formula
    print("\nComputing delta_I formula...", flush=True)
    t1 = time.time()
    delta_I = compute_formula_vectorized(arc_bit)
    print(f"  Done in {time.time()-t1:.1f}s")

    # Compare
    delta_H = U_Tp - U_T
    matches = np.sum(delta_H == delta_I)
    total = 1 << 20

    print(f"\n{'='*60}")
    print(f"Results: {matches}/{total} match")

    if matches == total:
        print("PROVED: delta_H = delta_I for ALL n=7 arc configurations.")
        print("Combined with base case, this proves OCF (H(T) = I(Omega(T), 2)) at n=7.")
        print("This extends the proof frontier from n<=6 to n<=7.")
    else:
        failures = np.where(delta_H != delta_I)[0]
        print(f"FAILED: {total - matches} counterexamples found.")
        for idx in failures[:5]:
            print(f"  mask={idx}: delta_H={delta_H[idx]}, delta_I={delta_I[idx]}, "
                  f"diff={delta_H[idx]-delta_I[idx]}")

    print(f"\nTotal time: {time.time()-t0:.1f}s")
