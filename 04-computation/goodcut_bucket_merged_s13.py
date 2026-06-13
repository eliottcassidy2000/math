#!/usr/bin/env python3
"""
goodcut_bucket_merged_s13.py

opus-2026-05-29-S13

Formalisation companion for TournamentH7.GoodCuts.lean.

The Lean module proves that the good-cut bucket d=1 is impossible and that
grid reflection preserves d. This script asks how the bucket index interacts
with the merged tournament quotient G_n/Z_2:

  tiling hypercube Q_m  --quotient-->  merged tournament classes
          |
          +-- good-cut bucket d = #cuts hit by upward tiles

Exact for n=3..6.
"""

from collections import Counter, defaultdict
from itertools import permutations

from projection_defect_waggly_layers_s1 import (
    hamiltonian_paths,
    merged_tournament_canon,
    tile_pairs,
    tiling_to_tournament,
)


def good_cut_set(bits, tiles):
    cuts = set()
    for k, (hi, lo) in enumerate(tiles):
        if bits & (1 << k):
            cuts.update(range(lo + 1, hi + 1))
    return frozenset(cuts)


def reflect_bits(bits, n, tiles):
    tile_index = {tile: k for k, tile in enumerate(tiles)}
    out = 0
    for k, (hi, lo) in enumerate(tiles):
        r_hi = n - 1 - lo
        r_lo = n - 1 - hi
        rk = tile_index[(r_hi, r_lo)]
        if bits & (1 << k):
            out |= 1 << rk
    return out


def build_merged_classes(n, tiles):
    perms = list(permutations(range(n)))
    canon_to_id = {}
    cls = {}
    H_by_class = {}
    for bits in range(1 << len(tiles)):
        A = tiling_to_tournament(bits, n, tiles)
        canon = merged_tournament_canon(A, perms)
        if canon not in canon_to_id:
            cid = len(canon_to_id)
            canon_to_id[canon] = cid
            H_by_class[cid] = hamiltonian_paths(A)
        cls[bits] = canon_to_id[canon]
    return cls, H_by_class


def transition_matrix(n, tiles):
    matrix = Counter()
    m = len(tiles)
    for bits in range(1 << m):
        d0 = len(good_cut_set(bits, tiles))
        for k in range(m):
            other = bits ^ (1 << k)
            if other <= bits:
                continue
            d1 = len(good_cut_set(other, tiles))
            matrix[(d0, d1)] += 1
            matrix[(d1, d0)] += 1
    return matrix


def analyze_n(n):
    tiles = tile_pairs(n)
    m = len(tiles)
    N = 1 << m
    cls, H_by_class = build_merged_classes(n, tiles)

    bucket_counts = Counter()
    cutset_counts = Counter()
    class_buckets = defaultdict(Counter)
    bucket_classes = defaultdict(set)
    reflection_failures = []

    for bits in range(N):
        g = good_cut_set(bits, tiles)
        d = len(g)
        bucket_counts[d] += 1
        cutset_counts[g] += 1
        cid = cls[bits]
        class_buckets[cid][d] += 1
        bucket_classes[d].add(cid)

        rb = reflect_bits(bits, n, tiles)
        if len(good_cut_set(rb, tiles)) != d:
            reflection_failures.append((bits, rb, d, len(good_cut_set(rb, tiles))))

    pure_classes = sum(1 for c in class_buckets.values() if len(c) == 1)
    mixed_classes = len(class_buckets) - pure_classes
    max_span = max(max(c) - min(c) for c in class_buckets.values())
    widest = [
        (cid, min(c), max(c), H_by_class[cid], dict(sorted(c.items())))
        for cid, c in class_buckets.items()
        if max(c) - min(c) == max_span
    ]

    print("\n" + "=" * 96)
    print(f"n={n}, m={m}, tilings={N}, merged classes={len(class_buckets)}")
    print(f"bucket counts d -> tilings: {dict(sorted(bucket_counts.items()))}")
    print(f"distinct good-cut sets: {len(cutset_counts)} of {2 ** (n - 1)} possible cut subsets")
    print(f"reflection bucket failures: {len(reflection_failures)}")
    print(f"class bucket purity: pure={pure_classes}, mixed={mixed_classes}, max bucket span={max_span}")
    print("classes touched by bucket:")
    for d in sorted(bucket_counts):
        print(f"  d={d}: {len(bucket_classes[d])} merged classes")

    print("widest merged classes by good-cut bucket span:")
    for cid, dmin, dmax, H, profile in widest[:10]:
        print(f"  class={cid:3d} H={H:3d} d=[{dmin},{dmax}] profile={profile}")

    print("single-tile bucket transition counts (directed, symmetrized):")
    matrix = transition_matrix(n, tiles)
    for (a, b), count in sorted(matrix.items()):
        if count:
            print(f"  {a}->{b}: {count}")


def main():
    print("GOOD-CUT BUCKETS VS MERGED TOURNAMENT CLASSES")
    print("opus-2026-05-29-S13")
    print("bit 1/upward tile contributes interval {lo+1,...,hi} of good cuts")
    for n in range(3, 7):
        analyze_n(n)


if __name__ == "__main__":
    main()
