"""
Q-009: Prove OCF at n=8 via exhaustive arc-flip identity verification.

At n=8, the simplified formula from n<=7 FAILS due to VD 3-5 pairs.
We use the full A-clique formula:
  delta_I = 2 * [sum_{C' gained} H(comp(C')) - sum_{C lost} H(comp(C))]

27 arc variables, 2^27 = 134,217,728 configurations.
Processed in 128 chunks of 2^20 = 1,048,576 assignments each.

Author: opus-2026-03-05-S4
"""

import numpy as np
from itertools import permutations, combinations
import time


def build_arc_encoding(others):
    """Build mapping from arc (a,b) to (bit_position, is_positive)."""
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
    return arc_bit, bit_idx


def precompute_shapes(n, others, arc_bit, for_T=True):
    """Precompute path shapes as bitmask conditions."""
    i_v, j_v = (0, 1) if for_T else (1, 0)
    shapes = []
    for pos in range(n - 1):
        remaining = [k for k in range(n) if k != pos and k != pos + 1]
        for perm in permutations(others):
            path = [None] * n
            path[pos] = i_v
            path[pos + 1] = j_v
            for idx, rp in enumerate(remaining):
                path[rp] = perm[idx]

            ones_mask = 0
            zeros_mask = 0
            impossible = False
            for k in range(n - 1):
                a, b = path[k], path[k+1]
                if for_T:
                    if a == 0 and b == 1: continue
                    if a == 1 and b == 0: impossible = True; break
                else:
                    if a == 1 and b == 0: continue
                    if a == 0 and b == 1: impossible = True; break
                bit_pos, positive = arc_bit[(a, b)]
                if positive: ones_mask |= (1 << bit_pos)
                else: zeros_mask |= (1 << bit_pos)

            if impossible or (ones_mask & zeros_mask):
                continue

            block_conds = []
            if for_T:
                if pos > 0:
                    pred = path[pos - 1]
                    bp, p = arc_bit[(pred, 1)]
                    block_conds.append((bp, not p))
                if pos + 1 < n - 1:
                    succ = path[pos + 2]
                    bp, p = arc_bit[(0, succ)]
                    block_conds.append((bp, not p))
            else:
                if pos > 0:
                    pred = path[pos - 1]
                    bp, p = arc_bit[(pred, 0)]
                    block_conds.append((bp, not p))
                if pos + 1 < n - 1:
                    succ = path[pos + 2]
                    bp, p = arc_bit[(1, succ)]
                    block_conds.append((bp, not p))
            shapes.append((ones_mask, zeros_mask, block_conds))
    return shapes


def compute_unmatched_chunk(shapes, assignments):
    """Compute unmatched count for a chunk of assignments."""
    N = len(assignments)
    unmatched = np.zeros(N, dtype=np.int32)
    for ones_mask, zeros_mask, block_conds in shapes:
        if not block_conds:
            continue
        check_mask = ones_mask | zeros_mask
        valid = (assignments & check_mask) == ones_mask
        blocked = np.zeros(N, dtype=bool)
        for bit_pos, blocks_when_one in block_conds:
            bit_val = (assignments >> bit_pos) & 1
            if blocks_when_one:
                blocked |= (bit_val == 1)
            else:
                blocked |= (bit_val == 0)
        unmatched += (valid & blocked).astype(np.int32)
    return unmatched


def precompute_h_lookup(verts, arc_bit):
    """Precompute H(T[verts]) lookup table indexed by sub-tournament config."""
    pairs = []
    for idx_a, a in enumerate(verts):
        for b in verts[idx_a+1:]:
            pairs.append((a, b))
    num_pairs = len(pairs)

    sub_bits = []
    for a, b in pairs:
        bp, pos = arc_bit[(a, b)]
        sub_bits.append((bp, pos))

    lookup = np.zeros(1 << num_pairs, dtype=np.int32)
    for sub_mask in range(1 << num_pairs):
        T = {}
        for idx, (a, b) in enumerate(pairs):
            if (sub_mask >> idx) & 1:
                T[(a, b)] = 1; T[(b, a)] = 0
            else:
                T[(a, b)] = 0; T[(b, a)] = 1
        h = 0
        for perm in permutations(verts):
            valid = True
            for k in range(len(verts) - 1):
                if T.get((perm[k], perm[k+1]), 0) != 1:
                    valid = False; break
            if valid: h += 1
        lookup[sub_mask] = h
    return sub_bits, lookup


def extract_sub_mask(assignments, sub_bits):
    """Extract sub-mask for lookup."""
    result = np.zeros(len(assignments), dtype=np.int32)
    for idx, (bit_pos, positive) in enumerate(sub_bits):
        bit_val = (assignments >> bit_pos) & 1
        if not positive:
            bit_val = 1 - bit_val
        result |= (bit_val.astype(np.int32) << idx)
    return result


def precompute_cycle_shapes(others, arc_bit):
    """Precompute cycle checking patterns for all odd lengths.

    For each (L, subset, permutation), store:
    - (ones_mask, zeros_mask) for lost cycle (j -> perm -> i)
    - (ones_mask, zeros_mask) for gained cycle (i -> perm -> j)
    - complement vertices for H lookup
    """
    cycle_patterns = []  # list of (lost_ones, lost_zeros, gained_ones, gained_zeros, complement)

    for L in range(3, len(others) + 3, 2):
        num_inter = L - 2
        if num_inter > len(others):
            break
        for subset in combinations(others, num_inter):
            complement = tuple(x for x in others if x not in subset)

            for perm in permutations(subset):
                # Lost: j=1 -> perm[0] -> ... -> perm[-1] -> i=0
                lost_ones = 0
                lost_zeros = 0
                lost_ok = True
                arcs = [(1, perm[0])]
                for k in range(len(perm) - 1):
                    arcs.append((perm[k], perm[k+1]))
                arcs.append((perm[-1], 0))
                for a, b in arcs:
                    bp, positive = arc_bit[(a, b)]
                    if positive: lost_ones |= (1 << bp)
                    else: lost_zeros |= (1 << bp)
                if lost_ones & lost_zeros:
                    lost_ok = False

                # Gained: i=0 -> perm[0] -> ... -> perm[-1] -> j=1
                gained_ones = 0
                gained_zeros = 0
                gained_ok = True
                arcs = [(0, perm[0])]
                for k in range(len(perm) - 1):
                    arcs.append((perm[k], perm[k+1]))
                arcs.append((perm[-1], 1))
                for a, b in arcs:
                    bp, positive = arc_bit[(a, b)]
                    if positive: gained_ones |= (1 << bp)
                    else: gained_zeros |= (1 << bp)
                if gained_ones & gained_zeros:
                    gained_ok = False

                cycle_patterns.append((
                    lost_ones if lost_ok else None,
                    lost_zeros if lost_ok else None,
                    gained_ones if gained_ok else None,
                    gained_zeros if gained_ok else None,
                    complement
                ))

    return cycle_patterns


def compute_delta_I_chunk(assignments, cycle_patterns, complement_lookups):
    """Compute delta_I for a chunk of assignments."""
    N = len(assignments)
    delta_I = np.zeros(N, dtype=np.int32)

    for lost_ones, lost_zeros, gained_ones, gained_zeros, complement in cycle_patterns:
        # Get H(complement) for these assignments
        if complement in complement_lookups:
            sub_bits, lookup = complement_lookups[complement]
            sub_mask = extract_sub_mask(assignments, sub_bits)
            H_comp = lookup[sub_mask]
        else:
            H_comp = np.ones(N, dtype=np.int32)

        contrib = np.zeros(N, dtype=np.int32)

        if gained_ones is not None:
            check = gained_ones | gained_zeros
            valid = (assignments & check) == gained_ones
            contrib += valid.astype(np.int32)

        if lost_ones is not None:
            check = lost_ones | lost_zeros
            valid = (assignments & check) == lost_ones
            contrib -= valid.astype(np.int32)

        delta_I += 2 * contrib * H_comp

    return delta_I


def prove_n(n):
    """Prove OCF at n via exhaustive arc-flip identity verification."""
    others = list(range(2, n))
    arc_bit, num_vars = build_arc_encoding(others)

    print(f"=== n={n} Exhaustive Proof of OCF ===")
    print(f"Variables: {num_vars}, Configurations: 2^{num_vars} = {1 << num_vars:,}")

    CHUNK_BITS = min(20, num_vars)
    CHUNK_SIZE = 1 << CHUNK_BITS
    NUM_CHUNKS = 1 << (num_vars - CHUNK_BITS)

    print(f"Chunks: {NUM_CHUNKS} x {CHUNK_SIZE:,}\n")

    # Precompute shapes
    t1 = time.time()
    shapes_T = precompute_shapes(n, others, arc_bit, for_T=True)
    shapes_Tp = precompute_shapes(n, others, arc_bit, for_T=False)
    print(f"Path shapes: T={len(shapes_T)}, T'={len(shapes_Tp)} ({time.time()-t1:.1f}s)")

    # Precompute cycle patterns
    t1 = time.time()
    cycle_patterns = precompute_cycle_shapes(others, arc_bit)
    print(f"Cycle patterns: {len(cycle_patterns)} ({time.time()-t1:.1f}s)")

    # Precompute complement H lookups
    t1 = time.time()
    complements = set(cp[4] for cp in cycle_patterns)
    complement_lookups = {}
    for comp in complements:
        if len(comp) >= 2:
            complement_lookups[comp] = precompute_h_lookup(list(comp), arc_bit)
        # len 0 or 1: H = 1, handled as default
    print(f"Complement lookups: {len(complement_lookups)} ({time.time()-t1:.1f}s)")

    # Main verification loop
    total_ok = 0
    total_checked = 0
    t0 = time.time()

    for chunk_idx in range(NUM_CHUNKS):
        base = chunk_idx * CHUNK_SIZE
        assignments = np.arange(base, base + CHUNK_SIZE, dtype=np.int32)

        U_T = compute_unmatched_chunk(shapes_T, assignments)
        U_Tp = compute_unmatched_chunk(shapes_Tp, assignments)
        delta_H = U_Tp - U_T

        delta_I = compute_delta_I_chunk(assignments, cycle_patterns, complement_lookups)

        matches = int(np.sum(delta_H == delta_I))
        total_ok += matches
        total_checked += CHUNK_SIZE

        elapsed = time.time() - t0
        rate = total_checked / elapsed if elapsed > 0 else 0
        remaining = (1 << num_vars) - total_checked
        eta = remaining / rate if rate > 0 else 0

        if chunk_idx % max(1, NUM_CHUNKS // 20) == 0 or matches != CHUNK_SIZE:
            print(f"  Chunk {chunk_idx+1}/{NUM_CHUNKS}: {matches}/{CHUNK_SIZE} "
                  f"[{elapsed:.0f}s elapsed, ETA {eta:.0f}s]")

        if matches != CHUNK_SIZE:
            failures = np.where(delta_H != delta_I)[0]
            for idx in failures[:5]:
                print(f"    FAIL at {base+idx}: dH={delta_H[idx]}, dI={delta_I[idx]}")
            print("ABORTING — formula incorrect at this n.")
            return False

    print(f"\n{'='*60}")
    print(f"Results: {total_ok}/{total_checked} match")
    total_time = time.time() - t0
    print(f"Total time: {total_time:.1f}s")

    if total_ok == (1 << num_vars):
        print(f"PROVED: delta_H = delta_I for ALL n={n} arc configurations.")
        print(f"Combined with base case, this proves OCF (H(T) = I(Omega(T), 2)) at n={n}.")
        return True
    return False


if __name__ == "__main__":
    t_total = time.time()

    # First verify n=7 (should be fast, serves as correctness check)
    print("Running n=7 as sanity check...")
    ok7 = prove_n(7)
    print()

    if ok7:
        print("n=7 passed. Now attempting n=8...\n")
        prove_n(8)

    print(f"\nGrand total time: {time.time()-t_total:.1f}s")
