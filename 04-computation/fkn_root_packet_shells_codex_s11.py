#!/usr/bin/env python3
"""
fkn_root_packet_shells_codex_s11.py

codex-2026-06-15-S11 follow-on to T823 / the FKN tiling-cube thread.

The goal is local, not full-cube: analyze the radius-1 and radius-2 Hamming
shells around the transitive ground state in the fixed-base-path tiling model.

Main checks:
  1. One-flip shell:
     - score defect is exactly e_y - e_x;
     - c3 = x-y-1;
     - H = 2^(x-y-1) + 1;
     - all one-flip tilings are pairwise non-isomorphic.
  2. Two-flip shell:
     - c3 is controlled exactly by interval incidence:
         same-end  -> g1 + g2 - 1
         cross-end -> g1 + g2 + 1
         disjoint  -> g1 + g2
       where g = x-y-1 is the free gap;
     - the quadratic H-defect q = H - (1 + 2^g1 + 2^g2) has the expected sign
       from THM-301/302 on same-end / cross-end pairs, and is nonnegative on all
       tested disjoint pairs.

Tournament Analysis note:
  The default tournament-analysis wrapper is not clean here. The objects under
  study already ARE tournaments, and the meaningful local structure is the
  interval/root packet geometry of the Hamming shells, not a derived pairwise
  dominance relation on some quotient. We therefore report shell and packet
  statistics directly.
"""

from collections import Counter
from itertools import combinations, permutations
from math import comb


def tiles(n):
    return [
        (x, y)
        for y in range(1, n - 1)
        for x in range(n, y + 1, -1)
        if x - y >= 2
    ]


def gap(tile):
    x, y = tile
    return x - y - 1


def relation(t1, t2):
    x1, y1 = t1
    x2, y2 = t2
    if x1 == x2 or y1 == y2:
        return "same-end"
    if x1 == y2 or y1 == x2:
        return "cross-end"
    return "disjoint"


def adj_from_active(n, active, all_tiles):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i > j:
                A[i][j] = 1
    active = set(active)
    for idx, (x, y) in enumerate(all_tiles):
        if idx in active:
            A[x - 1][y - 1] = 0
            A[y - 1][x - 1] = 1
    return A


def scores(A):
    return [sum(row) for row in A]


def count_c3(A):
    n = len(A)
    out = 0
    for i, j, k in combinations(range(n), 3):
        e = A[i][j] + A[j][k] + A[k][i]
        f = A[j][i] + A[k][j] + A[i][k]
        if e == 3 or f == 3:
            out += 1
    return out


def count_H(A):
    n = len(A)
    out_masks = [0] * n
    for i in range(n):
        mask = 0
        for j in range(n):
            if A[i][j]:
                mask |= 1 << j
        out_masks[i] = mask
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            paths = row[last]
            if not paths:
                continue
            nxt = out_masks[last] & ~mask
            while nxt:
                bit = nxt & -nxt
                j = bit.bit_length() - 1
                dp[mask | bit][j] += paths
                nxt ^= bit
    return sum(dp[full])


def canon_key(A):
    n = len(A)
    best = None
    for perm in permutations(range(n)):
        key = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or key < best:
            best = key
    return best


def shell_one(n):
    all_tiles = tiles(n)
    transitive_scores = list(range(n))
    iso_keys = set()
    score_keys = set()
    height_data = Counter()
    for idx, tile in enumerate(all_tiles):
        A = adj_from_active(n, [idx], all_tiles)
        sc = scores(A)
        c3 = count_c3(A)
        H = count_H(A)
        g = gap(tile)
        x, y = tile
        defect = [sc[v] - transitive_scores[v] for v in range(n)]
        expected_defect = [0] * n
        expected_defect[y - 1] = 1
        expected_defect[x - 1] = -1
        if defect != expected_defect:
            raise AssertionError(("score defect", n, tile, defect, expected_defect))
        if c3 != g:
            raise AssertionError(("shell1 c3", n, tile, c3, g))
        if H != 2**g + 1:
            raise AssertionError(("shell1 H", n, tile, H, 2**g + 1))
        iso_keys.add(canon_key(A))
        score_keys.add(tuple(sorted(sc)))
        height_data[g] += 1
    print(f"n={n} shell1: m={len(all_tiles)} one-flip atoms")
    print(f"  unique iso classes  = {len(iso_keys)} (expect m)")
    print(f"  unique score multisets = {len(score_keys)} (expect m)")
    print("  height summary (gap g -> count, c3, H):")
    for g in sorted(height_data):
        print(f"    g={g}: count={height_data[g]}, c3={g}, H={2**g + 1}")
    if len(iso_keys) != len(all_tiles):
        raise AssertionError(("shell1 iso distinctness", n, len(iso_keys), len(all_tiles)))
    if len(score_keys) != len(all_tiles):
        raise AssertionError(("shell1 score distinctness", n, len(score_keys), len(all_tiles)))


def expected_shell2_c3(t1, t2):
    g1 = gap(t1)
    g2 = gap(t2)
    kind = relation(t1, t2)
    if kind == "same-end":
        return g1 + g2 - 1
    if kind == "cross-end":
        return g1 + g2 + 1
    return g1 + g2


def shell_two(n):
    all_tiles = tiles(n)
    relation_counts = Counter()
    q_signs = Counter()
    q_ranges = {}
    exact = 0
    total = 0
    for i, j in combinations(range(len(all_tiles)), 2):
        t1 = all_tiles[i]
        t2 = all_tiles[j]
        kind = relation(t1, t2)
        A = adj_from_active(n, [i, j], all_tiles)
        c3 = count_c3(A)
        H = count_H(A)
        g1 = gap(t1)
        g2 = gap(t2)
        expected = expected_shell2_c3(t1, t2)
        if c3 != expected:
            raise AssertionError(("shell2 c3", n, t1, t2, c3, expected))
        q = H - (1 + 2**g1 + 2**g2)
        relation_counts[kind] += 1
        if q < 0:
            q_signs[(kind, "neg")] += 1
        elif q > 0:
            q_signs[(kind, "pos")] += 1
        else:
            q_signs[(kind, "zero")] += 1
        if kind not in q_ranges:
            q_ranges[kind] = [q, q]
        else:
            q_ranges[kind][0] = min(q_ranges[kind][0], q)
            q_ranges[kind][1] = max(q_ranges[kind][1], q)
        exact += 1
        total += 1
    print(f"n={n} shell2: {total} two-flip packets")
    print(f"  c3 packet law exact on all {exact} pairs")
    for kind in ("same-end", "cross-end", "disjoint"):
        if kind in q_ranges:
            lo, hi = q_ranges[kind]
            q_range = f"[{lo},{hi}]"
        else:
            q_range = "[]"
        print(
            f"  {kind:9s}: count={relation_counts[kind]:3d}, "
            f"q-signs="
            f"(neg={q_signs[(kind, 'neg')]}, "
            f"zero={q_signs[(kind, 'zero')]}, "
            f"pos={q_signs[(kind, 'pos')]}), "
            f"q-range={q_range}"
        )
    if q_signs[("same-end", "neg")] != relation_counts["same-end"]:
        raise AssertionError(("same-end sign", n, q_signs, relation_counts))
    if q_signs[("cross-end", "pos")] != relation_counts["cross-end"]:
        raise AssertionError(("cross-end sign", n, q_signs, relation_counts))
    if q_signs[("disjoint", "neg")] != 0:
        raise AssertionError(("disjoint sign", n, q_signs))


def main():
    print("=" * 72)
    print("FKN ROOT PACKET SHELLS AROUND THE TRANSITIVE TILING")
    print("=" * 72)
    print(
        "Free tiles are the non-simple interval roots of the fixed-path chamber; "
        "shell1 gives the dictator atoms, shell2 their first packet interactions."
    )
    for n in range(4, 9):
        print("\n" + "-" * 72)
        print(f"n={n}  (m=C(n-1,2)={comb(n - 1, 2)} free tiles)")
        print("-" * 72)
        shell_one(n)
        shell_two(n)
    print("\nDONE.")


if __name__ == "__main__":
    main()
