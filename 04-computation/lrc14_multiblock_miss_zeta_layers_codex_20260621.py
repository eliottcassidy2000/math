#!/usr/bin/env python3
"""Miss-zeta layer scout for the LRC14 multi-block carrier-product gap.

The current HYP-2714 residual is the moderate-span multi-block carrier error.
Route E already gives the decorrelated carrier product.  This script rewrites
that product in miss-zeta coordinates:

    1[all inner sectors hit] = sum_{R subset {1..6}} (-1)^|R| 1[R missed].

For a split row, it compares the actual miss-zeta profile with the shared-slow-x
independent-carrier product and reports which residual sizes pay the difference.
The goal is to identify the proof currency before reaching for a blunt
multi-dimensional Erdos-Turan constant.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, gcd


FULL = (1 << 6) - 1

CAP = {
    8: F(2243, 5880),
    9: F(1979, 4004),
    10: F(55, 91),
    11: F(66, 91),
    12: F(6, 7),
}


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def mask_size(mask: int) -> int:
    return mask.bit_count()


def mask_to_tuple(mask: int) -> tuple[int, ...]:
    return tuple(i + 1 for i in range(6) if mask & (1 << i))


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for x in row:
        g = gcd(g, x)
    return g == 1


@lru_cache(maxsize=None)
def row_breakpoints(row: tuple[int, ...]) -> tuple[F, ...]:
    bps = {F(0), F(1)}
    for e in row:
        if e == 0:
            continue
        for a in range(7 * abs(e) + 1):
            bps.add(F(a, 7 * abs(e)))
    return tuple(sorted(x for x in bps if 0 <= x <= 1))


@lru_cache(maxsize=None)
def actual_hit_profile(row: tuple[int, ...]) -> tuple[F, ...]:
    """Exact hit-mask law of the actual row over x."""
    row = tuple(sorted(set(row)))
    mass = [F(0) for _ in range(64)]
    xs = row_breakpoints(row)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = 0
        for e in row:
            sec = int(((e * mid) % 1) * 7)
            if 1 <= sec <= 6:
                hit |= 1 << (sec - 1)
        mass[hit] += hi - lo
    assert sum(mass, F(0)) == 1
    return tuple(mass)


def zeta_from_hit_profile(hit_profile: tuple[F, ...]) -> tuple[F, ...]:
    """z[R] = Pr(R is contained in the missed set)."""
    z = [F(0) for _ in range(64)]
    for hit, mass in enumerate(hit_profile):
        missed = FULL ^ hit
        sub = missed
        while True:
            z[sub] += mass
            if sub == 0:
                break
            sub = (sub - 1) & missed
    return tuple(z)


def cover_from_zeta(z: tuple[F, ...]) -> F:
    return sum(((-1) ** mask_size(mask)) * z[mask] for mask in range(64))


def krawtchouk(n: int, h: int, j: int) -> int:
    """Binary Krawtchouk value K_h(j) on the n-cube."""
    return sum(
        (-1) ** t * comb(j, t) * comb(n - j, h - t)
        for t in range(max(0, h - (n - j)), min(h, j) + 1)
    )


def krawtchouk_shadow(delta: tuple[F, ...]) -> tuple[tuple[F, ...], tuple[F, ...]]:
    """Radial weight enumerator and its Krawtchouk/MacWilliams shadow."""
    weight = [F(0) for _ in range(7)]
    for mask, value in enumerate(delta):
        weight[mask_size(mask)] += value
    shadow = [
        sum(weight[j] * krawtchouk(6, h, j) for j in range(7))
        for h in range(7)
    ]
    return tuple(weight), tuple(shadow)


@lru_cache(maxsize=None)
def block_x_breakpoints(block: tuple[int, ...]) -> tuple[F, ...]:
    """Internal slow-x breakpoints for a coherent block shape."""
    bps = {F(0), F(1)}
    for d in block:
        if d == 0:
            continue
        for a in range(7 * abs(d) + 1):
            bps.add(F(a, 7 * abs(d)))
    return tuple(sorted(x for x in bps if 0 <= x <= 1))


@lru_cache(maxsize=None)
def carrier_hit_law_at_x(block: tuple[int, ...], x: F) -> tuple[F, ...]:
    """Law over a uniform carrier theta of the inner sectors hit by theta+block*x."""
    base = tuple((d * x) % 1 for d in block)
    cuts = {F(0), F(1)}
    for b in base:
        for sec in range(7):
            cuts.add((F(sec, 7) - b) % 1)
    cut_list = sorted(cuts)
    mass = [F(0) for _ in range(64)]
    for lo, hi in zip(cut_list, cut_list[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = 0
        for b in base:
            sec = int(((b + mid) % 1) * 7)
            if 1 <= sec <= 6:
                hit |= 1 << (sec - 1)
        mass[hit] += hi - lo
    assert sum(mass, F(0)) == 1
    return tuple(mass)


@lru_cache(maxsize=None)
def carrier_miss_zeta_at_x(block: tuple[int, ...], x: F) -> tuple[F, ...]:
    return zeta_from_hit_profile(carrier_hit_law_at_x(block, x))


def product_miss_zeta_shared_x(blocks: tuple[tuple[int, ...], ...]) -> tuple[F, ...]:
    """Shared-slow-x independent-carrier product miss-zeta profile."""
    cuts = {F(0), F(1)}
    for block in blocks:
        cuts.update(block_x_breakpoints(block))
    xs = sorted(cuts)
    out = [F(0) for _ in range(64)]
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        zetas = [carrier_miss_zeta_at_x(block, mid) for block in blocks]
        for mask in range(64):
            v = F(1)
            for z in zetas:
                v *= z[mask]
            out[mask] += (hi - lo) * v
    return tuple(out)


def row_from_blocks(offset_blocks: tuple[tuple[int, tuple[int, ...]], ...]) -> tuple[int, ...]:
    row = [0]
    for offset, block in offset_blocks:
        row.extend(offset + d for d in block)
    return tuple(sorted(set(row)))


def layer_report(name: str, offset_blocks: tuple[tuple[int, tuple[int, ...]], ...]) -> dict[str, object]:
    blocks = tuple(block for _offset, block in offset_blocks)
    row = row_from_blocks(offset_blocks)
    actual_z = zeta_from_hit_profile(actual_hit_profile(row))
    product_z = product_miss_zeta_shared_x(blocks)
    actual = cover_from_zeta(actual_z)
    product = cover_from_zeta(product_z)
    assert actual == sum(actual_hit_profile(row)[hit] for hit in range(64) if hit == FULL)

    by_size = defaultdict(F)
    unsigned = defaultdict(F)
    delta = tuple(product_z[mask] - actual_z[mask] for mask in range(64))
    radial_weight, shadow = krawtchouk_shadow(delta)
    worst = (F(0), 0)
    for mask in range(64):
        term = ((-1) ** mask_size(mask)) * delta[mask]
        by_size[mask_size(mask)] += term
        unsigned[mask_size(mask)] += abs(term)
        if abs(term) > worst[0]:
            worst = (abs(term), mask)
    assert shadow[6] == product - actual

    cap = CAP.get(len(row))
    print("=" * 92)
    print(name)
    print(f"  row={row}")
    print(f"  block shapes={blocks}; primitive={primitive(row)}")
    print(f"  actual p0={fmt(actual)}")
    print(f"  product cover={fmt(product)}")
    print(f"  product-actual={fmt(product - actual)}  upper? {product >= actual}")
    if cap is not None:
        print(f"  cap_{len(row)}-actual={fmt(cap - actual)}; cap-product={fmt(cap - product)}")
        ratio = abs(product - actual) / (cap - product) if cap > product else None
        print(f"  abs(product-actual)/(cap-product)={fmt(ratio) if ratio is not None else 'n/a'}")
    print("  inclusion-exclusion layers for product-actual:")
    for size in range(7):
        if by_size[size] or unsigned[size]:
            print(
                f"    |R|={size}: signed={fmt(by_size[size])} "
                f"unsigned_L1={fmt(unsigned[size])}"
            )
    radial_l1 = sum(abs(x) for x in radial_weight)
    coord_l1 = sum(abs(x) for x in delta)
    print("  radial discrepancy weights W_j=sum_{|R|=j}(z_product-z_actual):")
    for size, value in enumerate(radial_weight):
        if value:
            print(f"    W_{size}={fmt(value)}")
    print("  Krawtchouk/MacWilliams shadow M_h=sum_j W_j K_h(j):")
    for h, value in enumerate(shadow):
        if value:
            marker = "  <-- Product-actual" if h == 6 else ""
            print(f"    M_{h}={fmt(value)}{marker}")
    if radial_l1:
        print(f"  |M_6|/sum_j|W_j|={fmt(abs(shadow[6]) / radial_l1)}")
    if coord_l1:
        print(f"  |M_6|/coordinate_L1={fmt(abs(shadow[6]) / coord_l1)}")
    print(
        "  largest coordinate term: "
        f"R={mask_to_tuple(worst[1])}, |term|={fmt(worst[0])}"
    )
    print()

    return {
        "name": name,
        "row": row,
        "actual": actual,
        "product": product,
        "gap": product - actual,
        "upper": product >= actual,
        "cap_slack_product": (cap - product) if cap is not None else None,
        "shadow": shadow,
        "radial_l1": radial_l1,
        "coord_l1": coord_l1,
        "worst_size": mask_size(worst[1]),
    }


def tournament(rows: list[dict[str, object]]) -> None:
    print("=" * 92)
    print("TOURNAMENT ANALYSIS")
    print("  vertices: tested multi-block proof obligations")
    print("  pairwise observable: larger product-cover slack, then smaller product-actual error")
    print("  switch/gauge: miss-zeta residual masks before scalar cover")
    scores = [0] * len(rows)
    edges = set()
    for i, a in enumerate(rows):
        for j, b in enumerate(rows):
            if i >= j:
                continue
            ai = (a["cap_slack_product"] or F(0), -abs(a["gap"]), -i)
            bj = (b["cap_slack_product"] or F(0), -abs(b["gap"]), -j)
            if ai >= bj:
                edges.add((i, j))
                scores[i] += 1
            else:
                edges.add((j, i))
                scores[j] += 1
    cycles = 0
    for a, b, c in combinations(range(len(rows)), 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian proof-pressure path:")
    order = sorted(
        range(len(rows)),
        key=lambda i: (rows[i]["cap_slack_product"] or F(0), -abs(rows[i]["gap"])),
        reverse=True,
    )
    for i in order:
        print(
            f"    {rows[i]['name']}: cap-product={rows[i]['cap_slack_product']} "
            f"product-actual={rows[i]['gap']}"
        )
    print()


def character_tournament(rows: list[dict[str, object]]) -> None:
    print("=" * 92)
    print("KRAWTCHOUK CHARACTER TOURNAMENT")
    print("  vertices: residual-character weights h=0..6 in the Boolean cube")
    print("  observable: larger aggregate |M_h| across tested split rows")
    print("  quotient: preserves the cover correction at h=6; destroys sector labels")
    risk = []
    for h in range(7):
        total = sum(abs(row["shadow"][h]) for row in rows)  # type: ignore[index]
        risk.append(total)
    scores = [0] * 7
    edges = set()
    for i, j in combinations(range(7), 2):
        if (risk[i], -i) >= (risk[j], -j):
            edges.add((i, j))
            scores[i] += 1
        else:
            edges.add((j, i))
            scores[j] += 1
    cycles = 0
    for a, b, c in combinations(range(7), 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  aggregate-risk path:")
    for h in sorted(range(7), key=lambda idx: (risk[idx], -idx), reverse=True):
        marker = "  <-- cover-error character" if h == 6 else ""
        print(f"    h={h}: sum|M_h|={fmt(risk[h])}{marker}")
    print()


def main() -> None:
    print("LRC14 multi-block miss-zeta layer scout")
    print("Exact arithmetic: actual row vs shared-slow-x independent carrier product.\n")
    cases = [
        ("two 4-blocks, moderate gap", ((14, (0, 1, 2, 3)), (30, (0, 1, 2, 3)))),
        ("two 4-blocks, wider gap", ((30, (0, 1, 2, 3)), (80, (0, 1, 2, 3)))),
        ("5+3 split", ((20, (0, 1, 2, 3, 4)), (55, (0, 1, 2)))),
        ("3+3+2 split", ((18, (0, 1, 2)), (45, (0, 1, 2)), (90, (0, 1)))),
        (
            "five 2-blocks",
            (
                (15, (0, 1)),
                (30, (0, 1)),
                (46, (0, 1)),
                (63, (0, 1)),
                (81, (0, 1)),
            ),
        ),
        (
            "seven singleton carriers",
            (
                (19, (0,)),
                (31, (0,)),
                (44, (0,)),
                (58, (0,)),
                (73, (0,)),
                (89, (0,)),
                (106, (0,)),
            ),
        ),
    ]
    reports = [layer_report(name, blocks) for name, blocks in cases]
    tournament(reports)
    character_tournament(reports)
    print("SYNTHESIS")
    print("  The carrier-product comparison is naturally a miss-zeta statement:")
    print("    product cover - actual p0 = sum_R (-1)^|R|(z_product(R)-z_actual(R)).")
    print("  Important correction: in this anchor-separated gauge, ProductCover is")
    print("  not an upper bound in the tested rows; actual p0 can exceed the product.")
    print("  The needed lemma is the symmetric budget statement")
    print("    |actual - product| <= cap_k - product,")
    print("  and the tested ratios are small relative to that budget.")
    print("  The useful proof currency is residual-size layers, not raw pointwise")
    print("  residual-coordinate dominance.  This avoids the false small-R cone route")
    print("  while preserving the exact product structure Route E needs.")
    print("  Boolean-cube reading: the same error is the h=6 Krawtchouk/MacWilliams")
    print("  coordinate M_6 of the radial miss-zeta discrepancy.  The proof target")
    print("  should bound this signed top character after routing low-height")
    print("  resonances, not the full residual-coordinate L1 norm.")
    print("  Next target: prove a signed Erdos-Turan/Koksma bound for the zeta layers")
    print("  whose unsigned mass is large, and use the product cap slack for the rest.")


if __name__ == "__main__":
    main()
