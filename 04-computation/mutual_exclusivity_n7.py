"""
Test blueself/blackself mutual exclusivity at n=7.

Strategy: n=7 has 2^15 = 32768 tilings. Full canonicalization (7! per tiling)
is expensive. We use score-sequence filtering: two tournaments can only be
isomorphic if they have the same score sequence. This dramatically reduces
the number of full canonical checks needed.

We also investigate:
  - Whether mutual exclusivity holds or fails at n=7
  - The structure of the flip operation on tournament classes
  - Counting patterns for self-flip at larger n

Author: opus-2026-03-06-S16
"""

import itertools
from collections import defaultdict, Counter
import random
import time


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


def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))


def adj_signature(A, n):
    """Fast isomorphism-invariant signature: sorted score sequence + sorted degree-sequence of
    the 'wins-against-higher-scorer' matrix."""
    scores = [sum(A[i]) for i in range(n)]
    # Create a refined signature using neighbor score multisets
    sig = []
    for i in range(n):
        out_scores = tuple(sorted([scores[j] for j in range(n) if A[i][j] == 1]))
        in_scores = tuple(sorted([scores[j] for j in range(n) if A[j][i] == 1]))
        sig.append((scores[i], out_scores, in_scores))
    return tuple(sorted(sig))


def canonicalize(A, n):
    perms = list(itertools.permutations(range(n)))
    best = None
    for p in perms:
        s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best


def is_grid_symmetric(bits, trans_map):
    for i in range(len(bits)):
        if trans_map[i] != i and bits[i] != bits[trans_map[i]]:
            return False
    return True


def flip_bits(bits):
    return tuple(1 - b for b in bits)


def are_isomorphic(A1, A2, n):
    """Check if two tournaments are isomorphic by trying all permutations."""
    # First quick check: score sequences must match
    if score_seq(A1, n) != score_seq(A2, n):
        return False
    # Refined check
    if adj_signature(A1, n) != adj_signature(A2, n):
        return False
    # Full check
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


def investigate_n7():
    n = 7
    print(f"{'=' * 60}")
    print(f"  MUTUAL EXCLUSIVITY TEST AT n = {n}")
    print(f"{'=' * 60}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)
    total = 1 << m
    print(f"  Tiles: {m}, Total tilings: {total}")

    # Step 1: Find all self-flip members
    # A tiling T is self-flip iff T and flip(T) are isomorphic.
    # Use score-sequence filtering to avoid expensive canonical checks.

    t0 = time.time()
    self_flip_members = []
    checked = 0
    skipped = 0

    # We only need to check masks < total/2 (each pair {T, flip(T)} only needs
    # checking once, and the flip of mask m is (2^M - 1) - m)
    all_ones = (1 << m) - 1

    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        flip_mask = mask ^ all_ones
        if flip_mask == mask:
            # Self-flip at BIT level — impossible since m > 0 and all bits
            # can't equal their complements
            continue

        # Only process each pair once
        if mask > flip_mask:
            continue

        A = bits_to_adj(n, tiles, bits)
        flip_b = flip_bits(bits)
        A_f = bits_to_adj(n, tiles, flip_b)

        if score_seq(A, n) != score_seq(A_f, n):
            skipped += 1
            continue

        if adj_signature(A, n) != adj_signature(A_f, n):
            skipped += 1
            continue

        checked += 1
        if are_isomorphic(A, A_f, n):
            gs = is_grid_symmetric(bits, trans_map)
            gs_f = is_grid_symmetric(flip_b, trans_map)
            self_flip_members.append({
                'mask': mask, 'flip_mask': flip_mask,
                'gs': gs, 'gs_flip': gs_f,
                'scores': score_seq(A, n),
            })

        if mask % 5000 == 0 and mask > 0:
            elapsed = time.time() - t0
            rate = mask / elapsed
            eta = (total // 2 - mask) / rate
            print(f"    {mask}/{total//2} ({100*mask*2/total:.0f}%) "
                  f"found {len(self_flip_members)} self-flip pairs, "
                  f"checked {checked}, skipped {skipped}, "
                  f"ETA {eta:.0f}s", flush=True)

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.1f}s")
    print(f"  Self-flip pairs found: {len(self_flip_members)}")
    print(f"  Full isomorphism checks: {checked}")
    print(f"  Skipped (score mismatch): {skipped}")

    # Step 2: Check mutual exclusivity
    # For each self-flip pair, both T and flip(T) are self-flip members.
    # Group by isomorphism class (use canonical form for the self-flip members only).

    print(f"\n  Canonicalizing self-flip members...")
    class_map = {}  # canonical -> list of (mask, gs) pairs
    for sf in self_flip_members:
        A = bits_to_adj(n, tiles, tuple((sf['mask'] >> k) & 1 for k in range(m)))
        canon = canonicalize(A, n)

        if canon not in class_map:
            class_map[canon] = []
        class_map[canon].append({'mask': sf['mask'], 'gs': sf['gs']})
        class_map[canon].append({'mask': sf['flip_mask'], 'gs': sf['gs_flip']})

    # Deduplicate within classes
    for canon in class_map:
        seen = set()
        unique = []
        for entry in class_map[canon]:
            if entry['mask'] not in seen:
                seen.add(entry['mask'])
                unique.append(entry)
        class_map[canon] = unique

    print(f"  Self-flip classes: {len(class_map)}")

    # Check mutual exclusivity
    violations = 0
    for canon, members in class_map.items():
        gs_types = set(m['gs'] for m in members)
        if len(gs_types) > 1:
            violations += 1
            n_gs = sum(1 for m in members if m['gs'])
            n_ngs = sum(1 for m in members if not m['gs'])
            print(f"  *** VIOLATION: class has {n_gs} grid-sym + {n_ngs} non-grid-sym self-flip members ***")
            print(f"      Scores: {score_seq(bits_to_adj(n, tiles, tuple((members[0]['mask']>>k)&1 for k in range(m))), n)}")

    if violations == 0:
        print(f"\n  *** THEOREM HOLDS at n={n}: Blueself and blackself are mutually exclusive ***")
        print(f"  at the class level. {len(class_map)} self-flip classes checked.")
    else:
        print(f"\n  THEOREM FAILS at n={n}: {violations} violations found!")

    # Statistics
    print(f"\n  Self-flip class statistics:")
    n_blueself_classes = sum(1 for members in class_map.values()
                            if all(m['gs'] for m in members))
    n_blackself_classes = sum(1 for members in class_map.values()
                             if all(not m['gs'] for m in members))
    print(f"    Blueself classes (all grid-sym): {n_blueself_classes}")
    print(f"    Blackself classes (all non-grid-sym): {n_blackself_classes}")

    sizes = Counter(len(members) for members in class_map.values())
    print(f"    Self-flip members per class: {dict(sorted(sizes.items()))}")

    # Score sequences of self-flip classes
    print(f"\n  Self-flip classes by score sequence:")
    for canon, members in sorted(class_map.items(),
                                  key=lambda x: x[1][0].get('gs', False)):
        bits = tuple((members[0]['mask'] >> k) & 1 for k in range(m))
        A = bits_to_adj(n, tiles, bits)
        sc = score_seq(A, n)
        gs_type = 'blueself' if members[0]['gs'] else 'blackself'
        print(f"    {gs_type:10s}: scores={sc}, {len(members)} self-flip members")


if __name__ == '__main__':
    investigate_n7()
