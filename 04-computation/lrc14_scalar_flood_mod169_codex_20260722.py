#!/usr/bin/env python3
"""Exact mod-169 obstruction for the two scalar flood tails in THM-2133.

Normalize the guard coefficient to one modulo 169 and restrict to guard-safe
torsion points whose numerator is nonzero modulo 13.  A terminal of exact
13-adic valuation one is safe at every such point.  The remaining five or six
unit terminals therefore have to cover the 110-point universe below.

The program constructs all 78 unit coefficient classes modulo sign.  Their
restrictions give 77 distinct masks; repeated masks cannot enlarge a union.
It then computes the exact maximum unions of five and six masks by a complete
include/exclude search with the certified upper bound obtained by summing the
largest remaining individual marginal gains.

The pairwise observable is exact mask intersection size.  Its faithful carrier
is an undirected weighted intersection graph plus the masks themselves.
Orienting an edge by label or weight would discard zero cuts and union data, so
there is no cover-preserving switch/gauge or tie Hamiltonian path.  Tournament
fingerprints are therefore intentionally not reported; the intersection
histogram and exact maximum-union witnesses are the relevant fingerprints.

No floating point arithmetic or optional package is used.  Runtime checks are
ordinary exceptions and remain active under ``python -O``.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations


P = 13
MODULUS = P * P
GUARD_CUTOFF = MODULUS // 7
TERMINAL_CUTOFF = MODULUS // 14


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def norm(x: int) -> int:
    x %= MODULUS
    return min(x, MODULUS - x)


UNIVERSE = tuple(
    z
    for z in range(MODULUS)
    if z % P != 0 and 7 * norm(z) > MODULUS
)
POSITION = {z: i for i, z in enumerate(UNIVERSE)}


def terminal_mask(r: int) -> int:
    mask = 0
    for z in UNIVERSE:
        if norm(r * z) <= TERMINAL_CUTOFF:
            mask |= 1 << POSITION[z]
    return mask


def coefficient_classes() -> tuple[tuple[int, ...], tuple[int, ...]]:
    labels = tuple(range(1, (MODULUS + 1) // 2))
    labels = tuple(r for r in labels if r % P != 0)
    return labels, tuple(terminal_mask(r) for r in labels)


COEFFICIENT_LABELS, COEFFICIENT_MASKS = coefficient_classes()


def distinct_masks() -> tuple[tuple[int, ...], tuple[int, ...]]:
    labels: list[int] = []
    masks: list[int] = []
    for label, mask in zip(COEFFICIENT_LABELS, COEFFICIENT_MASKS):
        if mask not in masks:
            labels.append(label)
            masks.append(mask)
    return tuple(labels), tuple(masks)


LABELS, MASKS = distinct_masks()


def greedy_union(cardinality: int) -> tuple[int, tuple[int, ...]]:
    mask = 0
    remaining = list(range(len(MASKS)))
    chosen: list[int] = []
    for _ in range(cardinality):
        index = max(
            remaining,
            key=lambda i: ((mask | MASKS[i]).bit_count(), -LABELS[i]),
        )
        remaining.remove(index)
        chosen.append(index)
        mask |= MASKS[index]
    return mask.bit_count(), tuple(chosen)


def exact_union_maximum(cardinality: int) -> tuple[int, tuple[int, ...], int, int]:
    """Return the exact maximum using a valid marginal-gain branch bound."""

    best_size, best_choice = greedy_union(cardinality)
    nodes = 0
    prunes = 0
    ordered = tuple(reversed(range(len(MASKS))))

    def search(pos: int, mask: int, left: int, chosen: tuple[int, ...]) -> None:
        nonlocal best_size, best_choice, nodes, prunes
        nodes += 1
        if left == 0:
            size = mask.bit_count()
            if size > best_size:
                best_size = size
                best_choice = chosen
            return
        if len(ordered) - pos < left:
            return

        current = mask.bit_count()
        gains = sorted(
            (
                (mask | MASKS[ordered[j]]).bit_count() - current
                for j in range(pos, len(ordered))
            ),
            reverse=True,
        )
        # Future masks overlap one another, so the sum of their individual
        # current marginal gains is an upper bound on every continuation.
        if current + sum(gains[:left]) <= best_size:
            prunes += 1
            return

        index = ordered[pos]
        search(pos + 1, mask | MASKS[index], left - 1, chosen + (index,))
        search(pos + 1, mask, left, chosen)

    search(0, 0, cardinality, ())
    return best_size, best_choice, nodes, prunes


def main() -> None:
    require(GUARD_CUTOFF == 24, "guard cutoff changed")
    require(TERMINAL_CUTOFF == 12, "terminal cutoff changed")
    require(len(UNIVERSE) == 110, "universe size changed")
    require(len(COEFFICIENT_LABELS) == 78, "unit sign-class count changed")
    require(len(MASKS) == 77, "restricted mask count changed")

    size_histogram = Counter(mask.bit_count() for mask in MASKS)
    intersection_histogram = Counter(
        (MASKS[i] & MASKS[j]).bit_count()
        for i, j in combinations(range(len(MASKS)), 2)
    )
    require(
        size_histogram == Counter({18: 40, 16: 21, 20: 11, 12: 3, 8: 1, 0: 1}),
        "mask-size histogram changed",
    )

    results = {}
    for cardinality, expected in ((5, 88), (6, 96)):
        maximum, witness, nodes, prunes = exact_union_maximum(cardinality)
        require(maximum == expected, f"maximum for {cardinality} masks changed")
        require(maximum < len(UNIVERSE), "a forbidden cover appeared")
        witness_labels = tuple(sorted(LABELS[i] for i in witness))
        results[cardinality] = (maximum, witness_labels, nodes, prunes)

    print("scalar flood exact mod-169 certificate")
    print(
        f"universe={len(UNIVERSE)}; unit_sign_classes={len(COEFFICIENT_LABELS)}; "
        f"distinct_restricted_masks={len(MASKS)}"
    )
    print("mask_sizes=" + repr(dict(sorted(size_histogram.items()))))
    print("pair_intersections=" + repr(dict(sorted(intersection_histogram.items()))))
    for cardinality in (5, 6):
        maximum, witness, nodes, prunes = results[cardinality]
        print(
            f"max_{cardinality}_union={maximum}; witness={witness}; "
            f"uncovered={len(UNIVERSE)-maximum}; nodes={nodes}; prunes={prunes}"
        )
    print("PASS")


if __name__ == "__main__":
    main()
