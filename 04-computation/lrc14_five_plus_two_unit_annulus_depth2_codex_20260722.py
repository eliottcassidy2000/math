#!/usr/bin/env python3
"""Exact Z/13^3 capacity obstruction for the 5+2 valuation pattern (1,2).

After normalizing the guard to one, write the divided deep coefficients as
13*u and 169*w with u,w units.  On 13^3-torsion, the 169*w tooth is safe
exactly off the zero residue column modulo 13.  For each of the 78 possible
u classes modulo 169 up to sign, this script forms the guard-safe points that
are also safe for both deep teeth.  It then checks all 1014 unit coefficient
classes modulo 2197 up to sign.  Five unit masks cannot cover the surviving
universe even if each is granted the largest individual mask for that u.

All predicates are integer cross-products.  The pairwise observable between a
deep class u and a unit class q is the exact surviving mask size.  It is a
bipartite capacity table, not a binary orientation; a tournament switch or tie
Hamiltonian path would forget the absolute five-mask capacity used by the
proof.  The row-size and margin histograms are the faithful fingerprints.

No floating point arithmetic or optional package is used.  Runtime checks are
ordinary exceptions and remain active under ``python -O``.
"""

from __future__ import annotations

from collections import Counter


P = 13
P2 = P * P
MODULUS = P**3


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def norm_mod(x: int, modulus: int) -> int:
    x %= modulus
    return min(x, modulus - x)


def sign_class_representatives(modulus: int) -> tuple[int, ...]:
    return tuple(
        r
        for r in range(1, (modulus + 1) // 2)
        if r % P != 0
    )


DEEP_CLASSES = sign_class_representatives(P2)
UNIT_CLASSES = sign_class_representatives(MODULUS)

BASE_POINTS = tuple(
    z
    for z in range(MODULUS)
    if z % P != 0 and 7 * norm_mod(z, MODULUS) > MODULUS
)
POSITION = {z: i for i, z in enumerate(BASE_POINTS)}


def mask_from_predicate(predicate) -> int:
    mask = 0
    for z in BASE_POINTS:
        if predicate(z):
            mask |= 1 << POSITION[z]
    return mask


UNIT_MASKS = tuple(
    mask_from_predicate(
        lambda z, q=q: 14 * norm_mod(q * z, MODULUS) < MODULUS
    )
    for q in UNIT_CLASSES
)


def deep_safe_mask(u: int) -> int:
    return mask_from_predicate(
        lambda z: 14 * norm_mod(u * z, P2) > P2
    )


def main() -> None:
    require(MODULUS == 2197, "modulus changed")
    require(len(DEEP_CLASSES) == 78, "deep sign-class count changed")
    require(len(UNIT_CLASSES) == 1014, "unit sign-class count changed")
    require(len(BASE_POINTS) == 1450, "base annulus size changed")

    universe_histogram: Counter[int] = Counter()
    margin_histogram: Counter[int] = Counter()
    worst_margin = None
    worst_rows: list[tuple[int, int, int, tuple[int, ...]]] = []

    for u in DEEP_CLASSES:
        universe = deep_safe_mask(u)
        universe_size = universe.bit_count()
        intersections = tuple((mask & universe).bit_count() for mask in UNIT_MASKS)
        maximum = max(intersections)
        maximizers = tuple(
            q for q, size in zip(UNIT_CLASSES, intersections) if size == maximum
        )
        margin = universe_size - 5 * maximum

        universe_histogram[universe_size] += 1
        margin_histogram[margin] += 1
        require(margin > 0, f"five-mask capacity failed for u={u}")

        row = (u, universe_size, maximum, maximizers)
        if worst_margin is None or margin < worst_margin:
            worst_margin = margin
            worst_rows = [row]
        elif margin == worst_margin:
            worst_rows.append(row)

    require(
        universe_histogram
        == Counter({1210: 2, 1218: 1, 1222: 3, 1226: 21, 1228: 40, 1230: 11}),
        "universe histogram changed",
    )
    require(worst_margin == 30, "worst capacity margin changed")
    require(worst_rows == [(14, 1230, 240, (183,))], "worst row changed")

    print("five-plus-two exact unit-annulus certificate")
    print(
        f"modulus={MODULUS}; deep_classes={len(DEEP_CLASSES)}; "
        f"unit_classes={len(UNIT_CLASSES)}; base_annulus={len(BASE_POINTS)}"
    )
    print("universe_sizes=" + repr(dict(sorted(universe_histogram.items()))))
    print("capacity_margins=" + repr(dict(sorted(margin_histogram.items()))))
    print(f"worst_margin={worst_margin}; worst_rows={worst_rows}")
    print("PASS")


if __name__ == "__main__":
    main()
