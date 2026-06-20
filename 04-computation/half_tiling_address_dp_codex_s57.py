#!/usr/bin/env python3
"""S57 half-address DP scaffold for grid-symmetric tilings.

This is the first HYP-2689 test after THM-551.  It enumerates the
grid-symmetric/self-converse subcube directly from half-address representatives
instead of enumerating all full tilings and filtering by reflection symmetry.

This is intentionally scoped: one bit per mirror orbit is valid only on this
grid-symmetric subcube.  The full complement quotient needs unordered pair
states on nonfixed mirror orbits.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product


Tile = tuple[int, int]
Address = tuple[int, int]


def tiles(n: int) -> list[Tile]:
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]


def addr(tile: Tile) -> Address:
    a, b = tile
    return (a, a + b - 1)


def reflect_tile(n: int, tile: Tile) -> Tile:
    a, b = tile
    return (n + 1 - b, n + 1 - a)


def half_reps(n: int) -> list[Tile]:
    reps = [t for t in tiles(n) if addr(t)[1] <= n]
    return sorted(reps, key=lambda t: (addr(t)[1], addr(t)[0]))


def crossing_layers(n: int) -> dict[int, int]:
    return dict(sorted(Counter(addr(t)[1] for t in half_reps(n)).items()))


def full_bits_from_half(n: int, rep_bits: tuple[int, ...]) -> dict[Tile, int]:
    reps = half_reps(n)
    bits: dict[Tile, int] = {}
    for tile, bit in zip(reps, rep_bits):
        mirror = reflect_tile(n, tile)
        bits[tile] = bit
        bits[mirror] = bit
    assert len(bits) == len(tiles(n))
    return bits


def tournament_adj(n: int, bits: dict[Tile, int]) -> list[list[int]]:
    adj = [[0] * (n + 1) for _ in range(n + 1)]
    for k in range(n, 1, -1):
        adj[k][k - 1] = 1
    for a, b in tiles(n):
        bit = bits[(a, b)]
        if bit == 0:
            adj[a][b] = 1
        else:
            adj[b][a] = 1
    return adj


def c3_count(adj: list[list[int]], n: int) -> int:
    total = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            for k in range(j + 1, n + 1):
                outs = (
                    adj[i][j] + adj[i][k],
                    adj[j][i] + adj[j][k],
                    adj[k][i] + adj[k][j],
                )
                if outs == (1, 1, 1):
                    total += 1
    return total


def hamiltonian_paths(adj: list[list[int]], n: int) -> int:
    full = (1 << n) - 1
    dp: list[defaultdict[int, int]] = [defaultdict(int) for _ in range(1 << n)]
    for v in range(1, n + 1):
        dp[1 << (v - 1)][v] = 1
    for mask in range(1 << n):
        if not dp[mask]:
            continue
        for last, count in list(dp[mask].items()):
            for nxt in range(1, n + 1):
                bit = 1 << (nxt - 1)
                if mask & bit:
                    continue
                if adj[last][nxt]:
                    dp[mask | bit][nxt] += count
    return sum(dp[full].values())


def score_multiset(adj: list[list[int]], n: int) -> tuple[int, ...]:
    return tuple(sorted(sum(adj[i][j] for j in range(1, n + 1)) for i in range(1, n + 1)))


def analyze_n(n: int) -> dict[str, object]:
    reps = half_reps(n)
    h_spec: Counter[int] = Counter()
    c3_spec: Counter[int] = Counter()
    score_spec: Counter[tuple[int, ...]] = Counter()
    leaders: list[tuple[int, int, tuple[int, ...], tuple[int, ...]]] = []
    for rep_bits in product((0, 1), repeat=len(reps)):
        bits = full_bits_from_half(n, rep_bits)
        adj = tournament_adj(n, bits)
        c3 = c3_count(adj, n)
        h = hamiltonian_paths(adj, n)
        scores = score_multiset(adj, n)
        c3_spec[c3] += 1
        h_spec[h] += 1
        score_spec[scores] += 1
        if not leaders or h > leaders[0][0]:
            leaders = [(h, c3, scores, rep_bits)]
        elif h == leaders[0][0] and len(leaders) < 5:
            leaders.append((h, c3, scores, rep_bits))
    return {
        "n": n,
        "half_coords": len(reps),
        "assignments": 2 ** len(reps),
        "layers": crossing_layers(n),
        "h_spec": h_spec,
        "c3_spec": c3_spec,
        "score_classes": len(score_spec),
        "leaders": leaders,
    }


def main() -> None:
    print("S57 half-tiling address DP scaffold")
    print("scope: grid-symmetric/self-converse subcube, one bit per mirror orbit")
    print("warning: full complement quotient requires unordered pair states")
    print()
    print(" n half_coords assignments crossing_layers Hmax H_spectrum c3_range score_multisets")
    for n in range(3, 9):
        row = analyze_n(n)
        h_values = row["h_spec"]
        c3_values = row["c3_spec"]
        hmax = max(h_values)
        c3_range = (min(c3_values), max(c3_values))
        print(
            f"{n:2d} {row['half_coords']:11d} {row['assignments']:11d} "
            f"{row['layers']} {hmax:4d} {len(h_values):10d} {c3_range!s:>9} "
            f"{row['score_classes']:15d}"
        )
        for h, c3, scores, rep_bits in row["leaders"][:3]:
            print(f"    leader H={h} c3={c3} scores={scores} half_bits={''.join(map(str, rep_bits))}")
    print()
    print("SYNTHESIS")
    print("  This directly enumerates 2^h_n grid-symmetric tilings by crossing layers.")
    print("  It avoids the full 2^F_n filter when the target is the self-converse half cube.")
    print("  It does not yet solve the full complement quotient; that needs unordered")
    print("  pair states for nonfixed mirror-cell orbits.")


if __name__ == "__main__":
    main()
