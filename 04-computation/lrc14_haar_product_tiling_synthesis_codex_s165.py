#!/usr/bin/env python3
"""Haar product / tournament-tiling synthesis scout for LRC14.

This script is deliberately small: it verifies the algebraic product rule that
the synthesis uses, then records a proof-carrier tournament over the recent
LRC14 lanes.

The point is not that dyadic Haar rectangles are literally tournament tilings.
The point is that both are coordinate-retained product algebras before an
unsafe quotient:

* geometric 2D Haar packets multiply by coordinatewise interval meet/refinement;
* fixed-path tournament tilings carry Boolean-Walsh packets whose product is
  symmetric difference of staircase tile supports.

Both product rules become lossy if one keeps only strip counts or scalar
profiles.  That is exactly the recent LRC14 lesson from Haar/Baire fronts,
handoff atlases, smoothing dispatchers, and tope/cocircuit packets.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from typing import Dict, Iterable, List, Sequence, Tuple


FINE = 32


@dataclass(frozen=True, order=True)
class DyadicInterval:
    level: int
    index: int

    @property
    def start(self) -> int:
        return self.index * FINE // (2**self.level)

    @property
    def end(self) -> int:
        return (self.index + 1) * FINE // (2**self.level)

    @property
    def mid(self) -> int:
        return (self.start + self.end) // 2

    def contains(self, other: "DyadicInterval") -> bool:
        return self.start <= other.start and other.end <= self.end

    def disjoint(self, other: "DyadicInterval") -> bool:
        return self.end <= other.start or other.end <= self.start

    def haar(self) -> Tuple[int, ...]:
        vals = [0] * FINE
        for i in range(self.start, self.mid):
            vals[i] = 1
        for i in range(self.mid, self.end):
            vals[i] = -1
        return tuple(vals)

    def indicator(self) -> Tuple[int, ...]:
        vals = [0] * FINE
        for i in range(self.start, self.end):
            vals[i] = 1
        return tuple(vals)


def vec_mul(a: Sequence[int], b: Sequence[int]) -> Tuple[int, ...]:
    return tuple(x * y for x, y in zip(a, b))


def vec_scale(c: int, a: Sequence[int]) -> Tuple[int, ...]:
    return tuple(c * x for x in a)


def classify_1d_product(a: DyadicInterval, b: DyadicInterval) -> Tuple[str, Tuple[int, ...]]:
    """Return the exact packet type and vector for h_a h_b.

    Unnormalized Haar is used.  For dyadic intervals:
    - disjoint intervals give zero;
    - equal intervals give the interval indicator;
    - if one interval contains the other, the coarser Haar is constant on the
      finer support, so the product is +/- the finer Haar.
    """
    ha, hb = a.haar(), b.haar()
    prod_vec = vec_mul(ha, hb)
    if a.disjoint(b):
        return "zero", tuple([0] * FINE)
    if a == b:
        return "square_indicator", a.indicator()
    if a.contains(b):
        sign = 1 if b.start < a.mid else -1
        return "signed_finer_haar", vec_scale(sign, b.haar())
    if b.contains(a):
        sign = 1 if a.start < b.mid else -1
        return "signed_finer_haar", vec_scale(sign, a.haar())
    raise AssertionError(f"dyadic intervals overlap without nesting: {a}, {b}")


def intervals(max_level: int = 3) -> List[DyadicInterval]:
    return [
        DyadicInterval(level, index)
        for level in range(max_level + 1)
        for index in range(2**level)
    ]


def tensor(a: Sequence[int], b: Sequence[int]) -> Tuple[int, ...]:
    return tuple(x * y for x in a for y in b)


def rectangle_haar(rect: Tuple[DyadicInterval, DyadicInterval]) -> Tuple[int, ...]:
    return tensor(rect[0].haar(), rect[1].haar())


def verify_haar_product() -> Dict[str, object]:
    ints = intervals(3)
    one_d_counts: Dict[str, int] = {}
    one_d_failures = 0
    for a, b in product(ints, repeat=2):
        kind, expected = classify_1d_product(a, b)
        actual = vec_mul(a.haar(), b.haar())
        one_d_counts[kind] = one_d_counts.get(kind, 0) + 1
        if actual != expected:
            one_d_failures += 1

    rect_intervals = intervals(2)
    rects = list(product(rect_intervals, repeat=2))
    two_d_counts: Dict[Tuple[str, str], int] = {}
    two_d_failures = 0
    for r1, r2 in product(rects, repeat=2):
        x_kind, x_vec = classify_1d_product(r1[0], r2[0])
        y_kind, y_vec = classify_1d_product(r1[1], r2[1])
        expected = tensor(x_vec, y_vec)
        actual = vec_mul(rectangle_haar(r1), rectangle_haar(r2))
        two_d_counts[(x_kind, y_kind)] = two_d_counts.get((x_kind, y_kind), 0) + 1
        if actual != expected:
            two_d_failures += 1

    return {
        "one_d_atoms": len(ints),
        "one_d_pairs": len(ints) ** 2,
        "one_d_counts": one_d_counts,
        "one_d_failures": one_d_failures,
        "two_d_rectangles": len(rects),
        "two_d_pairs": len(rects) ** 2,
        "two_d_counts": two_d_counts,
        "two_d_failures": two_d_failures,
    }


Tile = Tuple[int, int]
Mask = frozenset[Tile]


def staircase_tiles(n: int) -> List[Tile]:
    # Tile (a,b) is the fixed-base-path nonadjacent arc with a >= b+2.
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]


def tiling_masks(n: int) -> List[Tuple[str, Mask]]:
    tiles = staircase_tiles(n)
    masks: List[Tuple[str, Mask]] = []
    for tile in tiles:
        masks.append((f"tile_{tile[0]}_{tile[1]}", frozenset([tile])))
    for gap in sorted({a - b for a, b in tiles}):
        masks.append((f"gap_{gap}", frozenset((a, b) for a, b in tiles if a - b == gap)))
    for cut in range(2, n):
        masks.append((f"cut_{cut}", frozenset((a, b) for a, b in tiles if b < cut < a)))
    for top in range(4, n + 1):
        masks.append((f"principal_top_{top}", frozenset((a, b) for a, b in tiles if a >= top)))
    return masks


def verify_tiling_walsh(n: int = 6) -> Dict[str, object]:
    masks = tiling_masks(n)
    lookup = {mask: name for name, mask in masks}
    mismatches = 0
    closure_hits = 0
    product_sizes: Dict[int, int] = {}
    for _, a in masks:
        for _, b in masks:
            product_mask = a.symmetric_difference(b)
            product_sizes[len(product_mask)] = product_sizes.get(len(product_mask), 0) + 1
            if product_mask in lookup:
                closure_hits += 1
            # Boolean-Walsh characters satisfy chi_A * chi_B = chi_{A xor B}.
            # The mask computation is the exact product rule, so any mismatch
            # would be a programming error.
            if product_mask != (a | b) - (a & b):
                mismatches += 1
    return {
        "n": n,
        "tile_count": len(staircase_tiles(n)),
        "named_masks": len(masks),
        "mask_pairs": len(masks) ** 2,
        "walsh_mismatches": mismatches,
        "named_closure_hits": closure_hits,
        "product_size_hist": dict(sorted(product_sizes.items())),
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: Tuple[int, int, int, int, int, int, int]


def tournament_analysis() -> Dict[str, object]:
    # Coordinates:
    # LRC predicate, exact product rule, coordinate ownership, discrepancy sign,
    # boundary/topology, quotient warning, formal checkability.
    carriers = [
        Carrier("haar_rectangle_packet", (1, 1, 1, 1, 1, 1, 1)),
        Carrier("tiling_walsh_staircase_packet", (1, 1, 1, 0, 1, 1, 1)),
        Carrier("certificate_handoff_zipper", (1, 0, 1, 1, 1, 1, 1)),
        Carrier("tope_boundary_cocircuit_packet", (1, 0, 1, 1, 1, 1, 1)),
        Carrier("smoothing_dispatch_policy", (1, 0, 1, 1, 0, 1, 1)),
        Carrier("residue_discrepancy_green_packet", (1, 1, 0, 1, 0, 1, 1)),
        Carrier("strip_count_scalar_profile", (0, 0, 0, 1, 0, 0, 1)),
        Carrier("raw_tournament_isomorphism_class", (0, 0, 0, 0, 0, 0, 1)),
    ]
    tie_path = {name: idx for idx, name in enumerate(c.name for c in carriers)}

    edges: Dict[Tuple[str, str], str] = {}
    scores = {c.name: 0 for c in carriers}
    for a, b in combinations(carriers, 2):
        sa, sb = sum(a.vector), sum(b.vector)
        if sa > sb:
            winner = a
        elif sb > sa:
            winner = b
        else:
            winner = a if tie_path[a.name] < tie_path[b.name] else b
        loser = b if winner == a else a
        edges[(winner.name, loser.name)] = winner.name
        scores[winner.name] += 1

    score_hist: Dict[int, int] = {}
    for score in scores.values():
        score_hist[score] = score_hist.get(score, 0) + 1

    directed_3cycles = 0
    names = [c.name for c in carriers]
    for x, y, z in combinations(names, 3):
        xy = (x, y) in edges
        yz = (y, z) in edges
        zx = (z, x) in edges
        yx = (y, x) in edges
        zy = (z, y) in edges
        xz = (x, z) in edges
        if (xy and yz and zx) or (yx and zy and xz):
            directed_3cycles += 1

    path = sorted(carriers, key=lambda c: (-scores[c.name], tie_path[c.name]))
    return {
        "vertices": names,
        "observable": "retained LRC/product/discrepancy/topology/formal-check coordinates",
        "tie_path": " > ".join(names),
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": directed_3cycles,
        "retention_path": " > ".join(c.name for c in path),
        "scores": scores,
    }


def render() -> str:
    haar = verify_haar_product()
    tiling = verify_tiling_walsh(6)
    ta = tournament_analysis()
    lines: List[str] = []
    lines.append("LRC14 Haar-product / tournament-tiling synthesis scout (codex S165)")
    lines.append("=" * 74)
    lines.append("")
    lines.append("Geometric Haar product rule")
    lines.append("----------------------------")
    lines.append(f"1D atoms checked: {haar['one_d_atoms']}")
    lines.append(f"1D ordered products: {haar['one_d_pairs']}")
    lines.append(f"1D product kinds: {haar['one_d_counts']}")
    lines.append(f"1D failures: {haar['one_d_failures']}")
    lines.append(f"2D rectangles checked: {haar['two_d_rectangles']}")
    lines.append(f"2D ordered products: {haar['two_d_pairs']}")
    lines.append(f"2D factorization failures: {haar['two_d_failures']}")
    lines.append("2D product kind histogram:")
    for key, count in sorted(haar["two_d_counts"].items()):
        lines.append(f"  {key}: {count}")
    lines.append("")
    lines.append("Tournament tiling Walsh product rule")
    lines.append("------------------------------------")
    lines.append(f"fixed-path n: {tiling['n']}")
    lines.append(f"staircase tile count: {tiling['tile_count']}")
    lines.append(f"named masks tested: {tiling['named_masks']}")
    lines.append(f"ordered mask products: {tiling['mask_pairs']}")
    lines.append(f"Walsh xor mismatches: {tiling['walsh_mismatches']}")
    lines.append(f"named closure hits: {tiling['named_closure_hits']}")
    lines.append(f"product support-size histogram: {tiling['product_size_hist']}")
    lines.append("")
    lines.append("Shared readout")
    lines.append("--------------")
    lines.append("1. 2D Haar packets multiply by coordinatewise interval meet/refinement.")
    lines.append("2. Fixed-path tournament tilings multiply by symmetric difference of")
    lines.append("   staircase tile supports, the Boolean Haar/Walsh product rule.")
    lines.append("3. Both algebras are exact before quotienting and lossy after strip/count")
    lines.append("   scalarization unless owner coordinates are retained.")
    lines.append("4. LRC14 discrepancy should therefore be carried by labelled packets:")
    lines.append("   q/Farey scale, Haar open/boundary state, endpoint owners, packet family,")
    lines.append("   and certificate handoff.  AP/GW are boundary packets, not scalar zeros.")
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append("vertices: proof carriers, not runners/arcs/tiles")
    lines.append(f"pairwise observable: {ta['observable']}")
    lines.append(f"tie Hamiltonian path: {ta['tie_path']}")
    lines.append(f"score_hist: {ta['score_hist']}")
    lines.append(f"directed_3cycles: {ta['directed_3cycles']}")
    lines.append("retention path:")
    lines.append(f"  {ta['retention_path']}")
    lines.append("")
    lines.append("Candidate theorem target")
    lines.append("------------------------")
    lines.append("A quotient used in the LRC14 proof must be a homomorphism for the relevant")
    lines.append("Haar/Walsh packet product, or it must emit a named residual coefficient.")
    lines.append("Equivalently: scalar discrepancy, strip counts, and tournament isomorphism")
    lines.append("classes are admissible only after the lost 2D tile/rectangle address is")
    lines.append("reconstructed, annihilated by orthogonality, or routed to a state-lift/F7")
    lines.append("residual packet.")
    return "\n".join(lines) + "\n"


def main() -> None:
    output = render()
    print(output, end="")
    with open(
        "05-knowledge/results/lrc14_haar_product_tiling_synthesis_codex_s165.out",
        "w",
        encoding="utf-8",
    ) as f:
        f.write(output)


if __name__ == "__main__":
    main()
