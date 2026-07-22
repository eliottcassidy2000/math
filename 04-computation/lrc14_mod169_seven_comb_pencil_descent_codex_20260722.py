#!/usr/bin/env python3
"""Exact mod-169 certificate for THM-2128.

Universe.  At the valuation-one guard H=13 (after multiplication by a
unit), the strict guard-safe 169-torsion points are

    U = {z mod 169 : z mod 13 is not 0,+1,-1}.

For a terminal coefficient r prime to 13, its strict radius-1/14 danger
set on U is

    S_r = {z in U : min(rz mod 169, -rz mod 169) <= 12}.

The integer 12 is exact because 12 < 169/14 < 13.  The script proves that
seven such sets cannot cover U in two independent ways:

1. the small zero-intersection graph certificate used in the paper proof;
2. an exhaustive branch-and-bound maximum-union search over all C(78,7)
   sign classes, with a valid marginal-gain upper bound.

No floating point arithmetic or optional package is used.  Runtime checks
remain active under ``python -O``.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations


P = 13
MODULUS = P * P
TERMINAL_RADIUS_NUMERATOR = 12
SAFE_RESIDUES = frozenset(range(2, 12))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_residue_norm(x: int) -> int:
    x %= MODULUS
    return min(x, MODULUS - x)


UNIVERSE = tuple(z for z in range(MODULUS) if z % P in SAFE_RESIDUES)
POSITION = {z: i for i, z in enumerate(UNIVERSE)}
FULL_MASK = (1 << len(UNIVERSE)) - 1


def terminal_mask(r: int) -> int:
    mask = 0
    for z in UNIVERSE:
        if circle_residue_norm(r * z) <= TERMINAL_RADIUS_NUMERATOR:
            mask |= 1 << POSITION[z]
    return mask


def sign_classes() -> tuple[tuple[int, ...], tuple[int, ...]]:
    masks: list[int] = []
    labels: list[int] = []
    for r in range(1, MODULUS):
        if r % P == 0:
            continue
        mask = terminal_mask(r)
        if mask not in masks:
            masks.append(mask)
            labels.append(r)
    return tuple(masks), tuple(labels)


MASKS, LABELS = sign_classes()
INDEX = {label: i for i, label in enumerate(LABELS)}


def residue_type(label: int) -> int:
    residue = label % P
    return min(residue, P - residue)


def union_size(indices: tuple[int, ...] | list[int]) -> int:
    mask = 0
    for i in indices:
        mask |= MASKS[i]
    return mask.bit_count()


def zero_graph() -> tuple[frozenset[int], ...]:
    neighbors = [set() for _ in MASKS]
    for i, j in combinations(range(len(MASKS)), 2):
        if MASKS[i] & MASKS[j] == 0:
            neighbors[i].add(j)
            neighbors[j].add(i)
    return tuple(frozenset(row) for row in neighbors)


ZERO = zero_graph()


def certificate_row(i: int) -> tuple[int, int, int, int]:
    pair_common_cap = max(
        len(ZERO[i] & ZERO[j]) for j in range(len(MASKS)) if j != i
    )
    others = tuple(j for j in range(len(MASKS)) if j != i)
    triple_common_cap = max(
        len(ZERO[i] & ZERO[j] & ZERO[k]) for j, k in combinations(others, 2)
    )
    six_neighbor_union_cap = max(
        union_size((i,) + row) for row in combinations(sorted(ZERO[i]), 6)
    )
    return (
        len(ZERO[i]),
        pair_common_cap,
        triple_common_cap,
        six_neighbor_union_cap,
    )


EXPECTED_ROWS = {
    1: (7, 4, 3, 110),
    2: (7, 4, 3, 102),
    3: (6, 4, 3, 96),
    4: (6, 4, 3, 96),
    5: (8, 4, 3, 112),
    6: (8, 4, 3, 104),
}


def greedy_seven_union() -> tuple[int, tuple[int, ...]]:
    mask = 0
    remaining = list(range(len(MASKS)))
    chosen: list[int] = []
    for _ in range(7):
        i = max(
            remaining,
            key=lambda j: ((mask | MASKS[j]).bit_count(), -LABELS[j]),
        )
        remaining.remove(i)
        chosen.append(i)
        mask |= MASKS[i]
    return mask.bit_count(), tuple(chosen)


def exact_seven_union_maximum() -> tuple[int, tuple[int, ...], int, int]:
    """Exhaust all seven-subsets with a certified marginal-gain bound."""

    seed_size, seed = greedy_seven_union()
    best_size = seed_size
    best_choice = seed
    nodes = 0
    prunes = 0

    # The fixed order makes the include/exclude tree enumerate every subset
    # exactly once.  Descending index is slightly faster for this instance.
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
        # Adding the largest individual marginal gains ignores their future
        # overlaps and is therefore an upper bound on every continuation.
        if current + sum(gains[:left]) <= best_size:
            prunes += 1
            return

        i = ordered[pos]
        search(pos + 1, mask | MASKS[i], left - 1, chosen + (i,))
        search(pos + 1, mask, left, chosen)

    search(0, 0, 7, ())
    return best_size, best_choice, nodes, prunes


def main() -> None:
    require(len(UNIVERSE) == 130, "guard-safe universe size changed")
    require(len(MASKS) == 78, "unit sign-class count changed")
    require(len(set(MASKS)) == 78, "duplicate sign classes survived")

    columns = {a: tuple(z for z in UNIVERSE if z % P == a) for a in SAFE_RESIDUES}
    require(all(len(column) == 13 for column in columns.values()), "column size")
    for mask in MASKS:
        require(mask.bit_count() == 20, "terminal set does not have size 20")
        for column in columns.values():
            require(
                sum((mask >> POSITION[z]) & 1 for z in column) == 2,
                "terminal set is not a double toothpick in every column",
            )
        for z in UNIVERSE:
            neg = (-z) % MODULUS
            require(
                ((mask >> POSITION[z]) & 1) == ((mask >> POSITION[neg]) & 1),
                "negation symmetry failed",
            )

    intersection_histogram = Counter(
        (MASKS[i] & MASKS[j]).bit_count()
        for i, j in combinations(range(len(MASKS)), 2)
    )
    require(
        intersection_histogram == Counter({0: 273, 2: 1690, 4: 728, 6: 182, 8: 52, 10: 78}),
        "pair-intersection histogram changed",
    )
    require(
        all(size % 2 == 0 for size in intersection_histogram),
        "pair intersections should be even by negation",
    )

    rows = {}
    for kind in range(1, 7):
        representative = INDEX[kind]
        row = certificate_row(representative)
        rows[kind] = row
        require(row == EXPECTED_ROWS[kind], f"certificate row {kind} changed")

    # Multiplication of the torsion variable by a unit congruent to +/-1
    # modulo 13 preserves U.  It is transitive on each residue type, so the
    # six representative rows must also occur at every one of the 78 labels.
    for i, label in enumerate(LABELS):
        require(
            certificate_row(i) == EXPECTED_ROWS[residue_type(label)],
            f"orbit transport failed at label {label}",
        )

    table_cap = max(row[3] for row in rows.values())
    require(table_cap == 112 < len(UNIVERSE), "zero-graph table no longer excludes cover")

    maximum, witness, nodes, prunes = exact_seven_union_maximum()
    witness_labels = tuple(sorted(LABELS[i] for i in witness))
    require(maximum == 116, "independent seven-union maximum changed")
    require(union_size(list(witness)) == maximum, "maximum witness mismatch")
    require(maximum < len(UNIVERSE), "seven translates unexpectedly cover U")

    print("THM-2128 exact mod-169 certificate")
    print(f"universe={len(UNIVERSE)}; sign_classes={len(MASKS)}; set_size=20; columns=10x13")
    print("pair_intersections=" + repr(dict(sorted(intersection_histogram.items()))))
    print("orbit_table: type -> (zero_degree,pair_common,triple_common,six_neighbor_union)")
    for kind in range(1, 7):
        print(f"  {kind} -> {rows[kind]}")
    print(f"zero_graph_cover_cap={table_cap}")
    print(
        f"independent_max_7_union={maximum}; witness={witness_labels}; "
        f"uncovered={len(UNIVERSE)-maximum}"
    )
    print(f"branch_and_bound_nodes={nodes}; prunes={prunes}")
    print("PASS")


if __name__ == "__main__":
    main()
