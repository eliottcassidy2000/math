#!/usr/bin/env python3
"""
fkn_near_transitive_tiling_codex.py

FKN-style stability probe for the staircase tiling model.

The transitive tiling is the all-zero corner of the tiling cube. We study:

1. Shells by Hamming distance from that corner.
2. The degree-1 Walsh profile of H on the tiling cube.
3. The top Möbius coefficient on the vertex-deletion lattice
   mu_f(T) = sum_{U subset V} (-1)^(n-|U|) f(T[U]),
   interpreted as the irreducible n-body interaction of f.

This is not a proof of an FKN theorem for tournaments. It is a diagnostic:
does the first Fourier layer dominate near the transitive corner, and where
does genuine interacting structure first appear?
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations
from math import comb


def tiles_for_n(n: int) -> list[tuple[int, int]]:
    tiles = []
    for b in range(1, n - 1):
        for a in range(b + 2, n + 1):
            tiles.append((a, b))
    return tiles


def bits_to_adj(bits: int, n: int, tiles: list[tuple[int, int]]) -> list[list[int]]:
    adj = [[0] * n for _ in range(n)]
    for i in range(1, n):
        adj[i][i - 1] = 1
    for idx, (a, b) in enumerate(tiles):
        a0 = a - 1
        b0 = b - 1
        if (bits >> idx) & 1:
            adj[b0][a0] = 1
        else:
            adj[a0][b0] = 1
    return adj


def flat_mask(adj: list[list[int]]) -> tuple[int, ...]:
    return tuple(
        sum(adj[i][j] << j for j in range(len(adj))) for i in range(len(adj))
    )


def induced_submask(mask_rows: tuple[int, ...], verts: tuple[int, ...]) -> tuple[int, ...]:
    pos = {v: i for i, v in enumerate(verts)}
    out = [0] * len(verts)
    for old_i in verts:
        new_i = pos[old_i]
        row_mask = 0
        row = mask_rows[old_i]
        for old_j in verts:
            if (row >> old_j) & 1:
                row_mask |= 1 << pos[old_j]
        out[new_i] = row_mask
    return tuple(out)


@lru_cache(maxsize=None)
def hamiltonian_paths(mask_rows: tuple[int, ...]) -> int:
    n = len(mask_rows)
    if n <= 1:
        return 1
    full = 1 << n
    dp = [[0] * n for _ in range(full)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(full):
        for last in range(n):
            cur = dp[mask][last]
            if cur == 0:
                continue
            remaining = (~mask) & (full - 1)
            nxt = remaining
            while nxt:
                bit = nxt & -nxt
                v = bit.bit_length() - 1
                if (mask_rows[last] >> v) & 1:
                    dp[mask | bit][v] += cur
                nxt ^= bit
    return sum(dp[full - 1])


def count_c3(mask_rows: tuple[int, ...]) -> int:
    n = len(mask_rows)
    total = 0
    for i, j, k in combinations(range(n), 3):
        s = ((mask_rows[i] >> j) & 1) + ((mask_rows[j] >> k) & 1) + ((mask_rows[k] >> i) & 1)
        if s in (0, 3):
            total += 1
    return total


def score_sequence(mask_rows: tuple[int, ...]) -> tuple[int, ...]:
    n = len(mask_rows)
    return tuple(sorted(mask_rows[i].bit_count() for i in range(n)))


def canonical_label(mask_rows: tuple[int, ...]) -> tuple[int, ...]:
    n = len(mask_rows)
    best = None
    for perm in permutations(range(n)):
        candidate = []
        for i in perm:
            row_mask = 0
            for new_j, old_j in enumerate(perm):
                if (mask_rows[i] >> old_j) & 1:
                    row_mask |= 1 << new_j
            candidate.append(row_mask)
        candidate_t = tuple(candidate)
        if best is None or candidate_t < best:
            best = candidate_t
    return best


def all_subset_values(mask_rows: tuple[int, ...]) -> tuple[dict[int, int], dict[int, int]]:
    n = len(mask_rows)
    h_by_subset: dict[int, int] = {}
    c3_by_subset: dict[int, int] = {}
    for subset_mask in range(1 << n):
        verts = tuple(v for v in range(n) if (subset_mask >> v) & 1)
        sub = induced_submask(mask_rows, verts)
        h_by_subset[subset_mask] = hamiltonian_paths(sub)
        c3_by_subset[subset_mask] = count_c3(sub)
    return h_by_subset, c3_by_subset


def mobius_top(values: dict[int, int], n: int) -> int:
    total = 0
    for subset_mask, value in values.items():
        sign = -1 if (n - subset_mask.bit_count()) % 2 else 1
        total += sign * value
    return total


def walsh_degree1_coeffs(values: list[int], m: int) -> list[float]:
    denom = float(1 << m)
    coeffs = []
    for i in range(m):
        acc = 0
        for bits, value in enumerate(values):
            if (bits >> i) & 1:
                acc -= value
            else:
                acc += value
        coeffs.append(acc / denom)
    return coeffs


def walsh_energies(values: list[int], coeffs1: list[float]) -> tuple[float, float]:
    denom = float(len(values))
    mean = sum(values) / denom
    total = sum((v - mean) ** 2 for v in values) / denom
    level1 = sum(c * c for c in coeffs1)
    return level1, total


def tile_strip(tile: tuple[int, int]) -> int:
    a, b = tile
    return a - b - 1


def analyze_n(n: int) -> list[str]:
    tiles = tiles_for_n(n)
    m = len(tiles)
    all_h = [0] * (1 << m)
    shell_rows: dict[int, list[dict[str, object]]] = defaultdict(list)
    single_flip_records: list[tuple[int, tuple[int, int], int, int, int, int]] = []
    canon_counts: Counter[tuple[int, ...]] = Counter()
    canon_by_shell: dict[int, Counter[tuple[int, ...]]] = defaultdict(Counter)

    for bits in range(1 << m):
        adj = bits_to_adj(bits, n, tiles)
        mask_rows = flat_mask(adj)
        h_val = hamiltonian_paths(mask_rows)
        c3_val = count_c3(mask_rows)
        score = score_sequence(mask_rows)
        all_h[bits] = h_val
        shell = bits.bit_count()
        h_by_subset, c3_by_subset = all_subset_values(mask_rows)
        mu_h = mobius_top(h_by_subset, n)
        mu_c3 = mobius_top(c3_by_subset, n)
        record = {
            "bits": bits,
            "h": h_val,
            "c3": c3_val,
            "score": score,
            "mu_h": mu_h,
            "mu_c3": mu_c3,
        }
        shell_rows[shell].append(record)
        if n <= 6:
            canon = canonical_label(mask_rows)
            canon_counts[canon] += 1
            canon_by_shell[shell][canon] += 1
        if shell == 1:
            idx = bits.bit_length() - 1
            tile = tiles[idx]
            single_flip_records.append(
                (idx, tile, tile_strip(tile), h_val, c3_val, mu_h)
            )

    coeffs1 = walsh_degree1_coeffs(all_h, m)
    level1, total_var = walsh_energies(all_h, coeffs1)

    lines = []
    lines.append("=" * 78)
    lines.append(f"n={n}  staircase_bits={m}  tilings={1 << m}")
    lines.append("=" * 78)
    if n <= 6:
        lines.append(f"isomorphism classes among gauge-fixed tilings: {len(canon_counts)}")
    else:
        lines.append("isomorphism classes among gauge-fixed tilings: skipped at n=7")
    lines.append(
        "degree-1 Walsh energy ratio for H on the tiling cube: "
        f"{level1:.6f} / {total_var:.6f} = {level1 / total_var:.6f}"
        if total_var
        else "degree-1 Walsh energy ratio for H on the tiling cube: 0 / 0"
    )

    coord_data = sorted(
        enumerate(coeffs1),
        key=lambda item: abs(item[1]),
        reverse=True,
    )
    lines.append("top degree-1 coordinates for H:")
    for idx, coeff in coord_data[: min(6, len(coord_data))]:
        tile = tiles[idx]
        lines.append(
            f"  tile {tile} strip={tile_strip(tile)} coeff={coeff:.6f}"
        )

    strip_coeffs: dict[int, list[float]] = defaultdict(list)
    for idx, coeff in enumerate(coeffs1):
        strip_coeffs[tile_strip(tiles[idx])].append(coeff)
    lines.append("degree-1 coefficients grouped by strip:")
    for strip in sorted(strip_coeffs):
        vals = strip_coeffs[strip]
        avg = sum(vals) / len(vals)
        span = max(vals) - min(vals)
        lines.append(
            f"  strip {strip}: count={len(vals)} avg={avg:.6f} min={min(vals):.6f} "
            f"max={max(vals):.6f} span={span:.6f}"
        )

    lines.append("shell summary (distance from transitive corner):")
    for shell in range(m + 1):
        rows = shell_rows[shell]
        h_vals = sorted({row['h'] for row in rows})
        mu_h_vals = sorted({row['mu_h'] for row in rows})
        mu_c3_vals = sorted({row['mu_c3'] for row in rows})
        top_examples = Counter((row["h"], row["c3"], row["score"]) for row in rows).most_common(5)
        lines.append(
            f"  d={shell}: count={len(rows)} distinct_H={len(h_vals)} "
            f"H_range=[{min(h_vals)},{max(h_vals)}] "
            f"mu_H={mu_h_vals[:4]}{'...' if len(mu_h_vals) > 4 else ''} "
            f"mu_c3={mu_c3_vals[:4]}{'...' if len(mu_c3_vals) > 4 else ''}"
        )
        for (h_val, c3_val, score), ct in top_examples:
            lines.append(
                f"    signature H={h_val} c3={c3_val} score={score} occurs {ct} times"
            )
        if n <= 6:
            lines.append(
                f"    iso classes represented in shell: {len(canon_by_shell[shell])}"
            )

    if single_flip_records:
        lines.append("single-flip shell d=1:")
        for idx, tile, strip, h_val, c3_val, mu_h in sorted(
            single_flip_records,
            key=lambda row: (row[2], row[1][1], row[1][0]),
            reverse=True,
        ):
            lines.append(
                f"  tile {tile} strip={strip} -> H={h_val} c3={c3_val} mu_H={mu_h}"
            )
        best = max(single_flip_records, key=lambda row: row[3])
        lines.append(
            "best single flip by H: "
            f"tile {best[1]} strip={best[2]} gives H={best[3]}, c3={best[4]}, mu_H={best[5]}"
        )

    transitive_h = all_h[0]
    lines.append(f"transitive corner: H={transitive_h}, c3=0, mu_H={mobius_top(all_subset_values(flat_mask(bits_to_adj(0, n, tiles)))[0], n)}")
    lines.append("")
    return lines


def main() -> None:
    print("FKN-style near-transitive tiling probe")
    print("Board = free transitive corner; tournament = interacting perturbation.")
    print("Metrics: shell distance, degree-1 Walsh layer of H, top Möbius interaction.")
    print("")
    for n in range(4, 8):
        for line in analyze_n(n):
            print(line)
    print("Interpretation:")
    print("  1. If the FKN analogy were literally dictator-like, degree-1 energy would dominate.")
    print("  2. What the data actually asks is subtler: does one coordinate define the")
    print("     leading escape direction from the transitive corner?")
    print("  3. The candidate dictator-shadow is the extreme tile (n,1), i.e. reversing")
    print("     the min-max arc. The shell d=1 table shows whether that coordinate is unique.")
    print("  4. mu_H measures the irreducible n-vertex interaction after subtracting every")
    print("     lower-order subtournament shadow. Nonzero mu_H is the first clean sign that")
    print("     the tournament is more than a free assembly of its smaller faces.")


if __name__ == "__main__":
    main()
