#!/usr/bin/env python3
"""S217: diagonal-layer flow model for tournament half-tilings.

The prompt proposes a recursive tiling model where growing a tournament from
n to n+1 adds a diagonal layer of binary outputs, and consecutive layers of
sizes k and k+1 are connected by k^2+k lines.  This script treats those lines
as pairwise XOR/equality observations between tile bits.

Main algebra: the line system is the GF(2) coboundary/cut-space presentation
delta(x)(uv)=x_u+x_v on consecutive layer graphs.  K_{k,k+1} has k(k+1)
lines but XOR rank only 2k.  The k(k-1) missing dimensions are exactly the
4-cycle parity redundancies.  Globally, all adjacent layer lines have rank
#tiles - 1: the only invisible bit is global complement.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, permutations
from math import comb


KNOWN_U = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}


def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {pair: idx for idx, pair in enumerate(combinations(range(n), 2))}


PAIR_INDEX: dict[int, dict[tuple[int, int], int]] = {}


def idx(n: int) -> dict[tuple[int, int], int]:
    if n not in PAIR_INDEX:
        PAIR_INDEX[n] = pair_index(n)
    return PAIR_INDEX[n]


def edge(mask: int, n: int, i: int, j: int) -> int:
    if i < j:
        return (mask >> idx(n)[(i, j)]) & 1
    return 1 - ((mask >> idx(n)[(j, i)]) & 1)


def relabel(mask: int, n: int, perm: tuple[int, ...]) -> int:
    """Return mask in new labels where new vertex a is old vertex perm[a]."""
    out = 0
    out_idx = idx(n)
    for a, b in combinations(range(n), 2):
        if edge(mask, n, perm[a], perm[b]):
            out |= 1 << out_idx[(a, b)]
    return out


def canonical(mask: int, n: int) -> int:
    return min(relabel(mask, n, perm) for perm in permutations(range(n)))


def automorphism_count(mask: int, n: int) -> int:
    return sum(1 for perm in permutations(range(n)) if relabel(mask, n, perm) == mask)


def hamiltonian_path_count(mask: int, n: int) -> int:
    total = 0
    for perm in permutations(range(n)):
        if all(edge(mask, n, perm[i], perm[i + 1]) for i in range(n - 1)):
            total += 1
    return total


def fixed_path_tournament(n: int, tile_bits: int) -> int:
    """Tournament with base path 0->1->...->n-1 and nonpath tiles free."""
    mask = 0
    out_idx = idx(n)
    bit = 0
    for i, j in combinations(range(n), 2):
        if j == i + 1:
            value = 1
        else:
            value = (tile_bits >> bit) & 1
            bit += 1
        if value:
            mask |= 1 << out_idx[(i, j)]
    return mask


def fixed_path_fibers(n: int) -> dict[int, int]:
    tile_count = comb(n - 1, 2)
    fibers: Counter[int] = Counter()
    for bits in range(1 << tile_count):
        fibers[canonical(fixed_path_tournament(n, bits), n)] += 1
    return dict(fibers)


def layer_flow_rows(layer_sizes: list[int]) -> tuple[int, int, int]:
    """Return line count, GF(2) rank, and redundancy for adjacent K layers."""
    tile_count = sum(layer_sizes)
    line_count = sum(a * b for a, b in zip(layer_sizes, layer_sizes[1:]))
    nonempty = [size for size in layer_sizes if size]
    if len(nonempty) <= 1:
        rank = 0
    else:
        # Consecutive complete bipartite links make one connected component.
        rank = tile_count - 1
    return line_count, rank, line_count - rank


def local_layer_table(max_k: int) -> list[dict[str, int]]:
    rows = []
    for k in range(1, max_k + 1):
        lines = k * (k + 1)
        rank = 2 * k
        rows.append(
            {
                "k": k,
                "lines": lines,
                "xor_rank": rank,
                "four_cycle_redundancy": lines - rank,
                "rectangle_generators": k * (k - 1),
            }
        )
    return rows


def global_flow_table(max_n: int) -> list[dict[str, int]]:
    rows = []
    for n in range(2, max_n + 1):
        full_sizes = list(range(1, n))
        fixed_sizes = list(range(0, n - 1))
        full_lines, full_rank, full_red = layer_flow_rows(full_sizes)
        fixed_lines, fixed_rank, fixed_red = layer_flow_rows(fixed_sizes)
        fixed_tiles = comb(n - 1, 2)
        path_reflection_fixed_bits = ((n - 1) * (n - 1)) // 4
        path_reflection_merged = (
            2**fixed_tiles + 2**path_reflection_fixed_bits
        ) // 2
        rows.append(
            {
                "n": n,
                "U_n": KNOWN_U.get(n, -1),
                "full_tiles": comb(n, 2),
                "full_lines": full_lines,
                "full_rank": full_rank,
                "full_redundancy": full_red,
                "fixed_path_tiles": fixed_tiles,
                "fixed_path_tilings": 2**fixed_tiles,
                "fixed_lines": fixed_lines,
                "fixed_rank": fixed_rank,
                "fixed_redundancy": fixed_red,
                "path_reflection_fixed_bits": path_reflection_fixed_bits,
                "path_reflection_merged_tilings": path_reflection_merged,
            }
        )
    return rows


def print_local_table() -> None:
    print("1. LOCAL DIAGONAL-LAYER FLOW K_{k,k+1}")
    print("   k  lines=k(k+1)  xor_rank=2k  redundancy  rectangle_generators")
    for row in local_layer_table(8):
        print(
            f"   {row['k']:1d}  {row['lines']:12d}  {row['xor_rank']:11d}"
            f"  {row['four_cycle_redundancy']:10d}  {row['rectangle_generators']:20d}"
        )
    print(
        "   Redundancy k(k-1) is the cycle space of K_{k,k+1}: every 4-cycle says\n"
        "   L(a,b)+L(a,b')+L(a',b)+L(a',b')=0 over GF(2)."
    )
    print(
        "   A spanning-tree basis of 2k lines determines all others by\n"
        "   L(a,b)=L(a,b0)+L(a0,b)+L(a0,b0)."
    )


def print_global_table() -> None:
    print()
    print("2. GLOBAL ADJACENT-LAYER FLOW")
    print("   Full layers have sizes 1,2,...,n-1.  Fixed-path layers drop the path arcs.")
    print(
        "   n  U(n)  full_tiles  full_lines  full_rank  full_red  "
        "fixed_tiles  fixed_lines  fixed_rank  fixed_red  fixed_cube  Z2_path"
    )
    for row in global_flow_table(8):
        print(
            f"   {row['n']:1d}  {row['U_n']:4d}  {row['full_tiles']:10d}"
            f"  {row['full_lines']:10d}  {row['full_rank']:9d}"
            f"  {row['full_redundancy']:8d}  {row['fixed_path_tiles']:11d}"
            f"  {row['fixed_lines']:11d}  {row['fixed_rank']:10d}"
            f"  {row['fixed_redundancy']:9d}  {row['fixed_path_tilings']:10d}"
            f"  {row['path_reflection_merged_tilings']:7d}"
        )
    print(
        "   Full redundancy equals 2*C(n-1,3)+C(n-2,2): local rectangle cycles\n"
        "   plus extra hourglass cycles linking adjacent layer bridges."
    )


def print_fixed_path_fibers(max_n: int = 6) -> None:
    print()
    print("3. FIXED-PATH HALF-TILINGS SURJECT TO A000568 CLASSES")
    print("   A fixed Hamiltonian path gives 2^C(n-1,2) smaller half-tilings.")
    print("   The fiber over an isomorphism class is H(T)/|Aut(T)|.")
    for n in range(3, max_n + 1):
        fibers = fixed_path_fibers(n)
        checks = []
        for canon, count in fibers.items():
            h = hamiltonian_path_count(canon, n)
            aut = automorphism_count(canon, n)
            checks.append(count == h // aut and h % aut == 0)
        hist = dict(sorted(Counter(fibers.values()).items()))
        print(
            f"   n={n}: fixed_tilings={2 ** comb(n - 1, 2)} "
            f"iso_classes={len(fibers)} A000568={KNOWN_U[n]} "
            f"fiber_hist={hist} H_over_Aut_check={all(checks)}"
        )


def print_reading() -> None:
    print()
    print("READING")
    print(
        "  A labelled tournament is a binary upper-triangular half-tiling modulo "
        "the full S_n relabelling action.  A fixed Hamiltonian path gives the "
        "smaller cube with C(n-1,2) free nonpath tiles; it surjects onto A000568 "
        "classes with uneven fiber H(T)/|Aut(T)|."
    )
    print(
        "  The proposed k^2+k inter-layer lines are algebraically redundant in a "
        "controlled way: they are a GF(2) coboundary/cut-space presentation.  "
        "K_{k,k+1} has k(k+1) pair lines but only 2k XOR rank; the k(k-1) "
        "surplus lines are 4-cycle duplication laws.  Globally the adjacent-layer "
        "flow has rank #tiles-1, so it reconstructs all tile bits up to global "
        "complement."
    )
    print(
        "  Path reversal plus converse is therefore a useful diagonal half-tiling "
        "quotient, but it is not A000568 by itself.  The Z2 fixed-path counts "
        "already exceed U(n) at n=4,5,6.  A000568 is the orbit count after the full "
        "relabeling/path-presentation redundancy is also divided out."
    )
    print(
        "  LRC translation: the line redundancies should be used as packet sidecars, "
        "not scalars.  Rectangle/4-cycle relations are local duplication laws; "
        "hourglass cycles are bridge-to-bridge duplication laws; Hamiltonian-path "
        "fibers are global duplication laws; diagonal fixed sets are "
        "path-reversal/converse duplication laws."
    )


def main() -> None:
    print("=" * 80)
    print("S217: tournament diagonal-layer flow and half-tiling quotient laws")
    print("=" * 80)
    print_local_table()
    print_global_table()
    print_fixed_path_fibers()
    print_reading()


if __name__ == "__main__":
    main()
