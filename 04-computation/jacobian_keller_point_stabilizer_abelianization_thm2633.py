#!/usr/bin/env python3
"""Exact finite-group controls for THM-2633.

The geometric proof that a Keller point stabilizer surjects onto the
monodromy abelianization is in the theorem.  This companion exhausts the
small transitive permutation groups used by its degree-three and degree-four
corollaries and checks the D4 derangement character directly.
"""

from __future__ import annotations

import ast
from itertools import permutations
from pathlib import Path


Perm = tuple[int, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def compose(a: Perm, b: Perm) -> Perm:
    """Return a after b."""
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a: Perm) -> Perm:
    ans = [0] * len(a)
    for i, ai in enumerate(a):
        ans[ai] = i
    return tuple(ans)


def generated(generators: list[Perm]) -> set[Perm]:
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        a = frontier.pop()
        for b in generators:
            c = compose(a, b)
            if c not in group:
                group.add(c)
                frontier.append(c)
    return group


def stabilizer(group: set[Perm], point: int = 0) -> set[Perm]:
    return {g for g in group if g[point] == point}


def derived_subgroup(group: set[Perm]) -> set[Perm]:
    commutators = [
        compose(compose(compose(a, b), inverse(a)), inverse(b))
        for a in group
        for b in group
    ]
    return generated(commutators)


def parity(g: Perm) -> int:
    return sum(
        g[i] > g[j]
        for i in range(len(g))
        for j in range(i + 1, len(g))
    ) % 2


def gate_data(group: set[Perm]) -> tuple[int, int, int, bool]:
    """Return |G|, |H|, |[G,G]|, and H[G,G]=G."""
    h = stabilizer(group)
    derived = derived_subgroup(group)
    image = generated(list(h | derived))
    return len(group), len(h), len(derived), image == group


def main() -> None:
    id3 = (0, 1, 2)
    c3 = (1, 2, 0)
    t3 = (0, 2, 1)
    c3_group = generated([c3])
    s3_group = generated([c3, t3])
    require(id3 in c3_group, "C3 identity missing")

    id4 = (0, 1, 2, 3)
    r = (1, 2, 3, 0)
    s = (0, 3, 2, 1)
    z = compose(r, r)
    v = (1, 0, 3, 2)
    c4_group = generated([r])
    v4_group = generated([z, v])
    d4_group = generated([r, s])
    s4_group = set(permutations(range(4)))
    a4_group = {g for g in s4_group if parity(g) == 0}
    require(id4 in c4_group, "C4 identity missing")

    groups = [
        ("C3", c3_group, False),
        ("S3", s3_group, True),
        ("C4", c4_group, False),
        ("V4", v4_group, False),
        ("D4", d4_group, False),
        ("A4", a4_group, True),
        ("S4", s4_group, True),
    ]

    rows: list[tuple[str, int, int, int, bool]] = []
    for name, group, expected in groups:
        order, h_order, derived_order, passes = gate_data(group)
        require(passes == expected, f"{name} point-stabilizer gate changed")
        rows.append((name, order, h_order, derived_order, passes))

    # D4 deck character: in the normal form r^i s^j, the r-parity character
    # is trivial on the point stabilizer and every odd element is a
    # derangement in the four-sheet action.
    coordinates: dict[Perm, tuple[int, int]] = {}
    for i in range(4):
        ri = id4
        for _ in range(i):
            ri = compose(r, ri)
        for j in range(2):
            coordinates[compose(ri, s if j else id4)] = (i, j)
    require(len(coordinates) == 8, "D4 normal form is not unique")

    def chi_deck(g: Perm) -> int:
        return coordinates[g][0] % 2

    h_d4 = stabilizer(d4_group)
    deck_odd = {g for g in d4_group if chi_deck(g)}
    require(all(chi_deck(g) == 0 for g in h_d4),
            "D4 deck character is nontrivial on the point stabilizer")
    require(len(deck_odd) == 4, "D4 deck-odd class count changed")
    require(all(all(g[i] != i for i in range(4)) for g in deck_odd),
            "a D4 deck-odd element gained a fixed sheet")
    fixed_hist: dict[int, int] = {}
    for g in d4_group:
        fixed = sum(g[i] == i for i in range(4))
        fixed_hist[fixed] = fixed_hist.get(fixed, 0) + 1
    require(fixed_hist == {4: 1, 0: 5, 2: 2},
            "D4 fixed-sheet histogram changed")

    # The optimized interpreter must retain every truth-bearing check.
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(assert_count == 0, "companion contains Python assert statements")

    print("THM-2633 affine Keller point-stabilizer abelianization controls")
    print("group_gate_rows=(name,order,point_stabilizer,derived,passes)")
    for row in rows:
        print(row)
    print("degree4_excluded_by_gate=('C4','V4','D4')")
    print("degree4_survivors=('A4','S4')")
    print("D4_deck_character_trivial_on_H=PASS odd_elements=4 all_derangements=PASS")
    print("D4_fixed_sheet_hist=((0,5),(2,2),(4,1))")
    print("S3_positive_control=H_order_2_maps_onto_abelianization_C2")
    print("exact assertions passed")


if __name__ == "__main__":
    main()
