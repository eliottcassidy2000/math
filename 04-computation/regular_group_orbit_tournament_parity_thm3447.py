#!/usr/bin/env python3
"""Exact companion for THM-3447.

Classify invariant tournaments and generalized directed relations on a regular
finite-group orbit through inverse pairs in the difference group.  The script
uses only finite exact enumeration and integer arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
from typing import Callable, Hashable, Iterable


EXPECTED_SEMANTIC_SHA256 = "0672eab70672e58ff4f44298e7547bb44b9843f0f8f47987374fd03895b142c9"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def digest(value: object) -> str:
    return sha256(repr(value).encode("ascii")).hexdigest()


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


Element = Hashable


@dataclass(frozen=True)
class FiniteGroup:
    name: str
    elements: tuple[Element, ...]
    identity: Element
    multiply: Callable[[Element, Element], Element]

    def validate(self) -> dict[Element, Element]:
        elems = set(self.elements)
        require(len(elems) == len(self.elements), (self.name, "duplicate elements"))
        require(self.identity in elems, (self.name, "missing identity"))
        for x in self.elements:
            require(self.multiply(self.identity, x) == x, (self.name, "left identity", x))
            require(self.multiply(x, self.identity) == x, (self.name, "right identity", x))
        for x, y in product(self.elements, repeat=2):
            require(self.multiply(x, y) in elems, (self.name, "closure", x, y))
        for x, y, z in product(self.elements, repeat=3):
            require(
                self.multiply(self.multiply(x, y), z) == self.multiply(x, self.multiply(y, z)),
                (self.name, "associativity", x, y, z),
            )
        inverse: dict[Element, Element] = {}
        for x in self.elements:
            hits = tuple(
                y
                for y in self.elements
                if self.multiply(x, y) == self.identity
                and self.multiply(y, x) == self.identity
            )
            require(len(hits) == 1, (self.name, "inverse", x, hits))
            inverse[x] = hits[0]
        return inverse


def cyclic_group(n: int) -> FiniteGroup:
    return FiniteGroup(f"C{n}", tuple(range(n)), 0, lambda x, y: (x + y) % n)


def vector_group(modulus: int, rank: int, name: str | None = None) -> FiniteGroup:
    elements = tuple(product(range(modulus), repeat=rank))
    return FiniteGroup(
        name or f"(Z/{modulus})^{rank}",
        elements,
        (0,) * rank,
        lambda x, y: tuple((a + b) % modulus for a, b in zip(x, y)),
    )


def symmetric_group_three() -> FiniteGroup:
    elements = tuple(permutations(range(3)))

    def compose(x: tuple[int, ...], y: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(x[y[i]] for i in range(3))

    return FiniteGroup("S3", elements, (0, 1, 2), compose)


def heisenberg_group_three() -> FiniteGroup:
    elements = tuple(product(range(3), repeat=3))

    def multiply(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
        a, b, c = x
        d, e, f = y
        return ((a + d) % 3, (b + e) % 3, (c + f + a * e) % 3)

    return FiniteGroup("Heis(F3)", elements, (0, 0, 0), multiply)


def inverse_packet(
    group: FiniteGroup, inverse: dict[Element, Element]
) -> tuple[tuple[Element, ...], tuple[tuple[Element, Element], ...]]:
    involutions = tuple(
        x for x in group.elements if x != group.identity and inverse[x] == x
    )
    seen: set[Element] = set(involutions)
    pairs: list[tuple[Element, Element]] = []
    for x in group.elements:
        if x == group.identity or x in seen:
            continue
        y = inverse[x]
        require(y != x, (group.name, "unexpected fixed inverse", x))
        pairs.append((x, y))
        seen.add(x)
        seen.add(y)
    require(len(seen) == len(group.elements) - 1, (group.name, "inverse partition", seen))
    return involutions, tuple(pairs)


def difference(
    group: FiniteGroup,
    inverse: dict[Element, Element],
    x: Element,
    y: Element,
) -> Element:
    return group.multiply(inverse[x], y)


def arcs_from_connection_set(
    group: FiniteGroup,
    inverse: dict[Element, Element],
    connection_set: set[Element],
) -> set[tuple[Element, Element]]:
    return {
        (x, y)
        for x in group.elements
        for y in group.elements
        if x != y and difference(group, inverse, x, y) in connection_set
    }


def relation_kind(
    group: FiniteGroup,
    arcs: set[tuple[Element, Element]],
    maximal_partial_arc_count: int,
) -> tuple[bool, bool, bool, bool, bool]:
    tournament = True
    partial = True
    semicomplete = True
    symmetric = True
    for x, y in combinations(group.elements, 2):
        count = int((x, y) in arcs) + int((y, x) in arcs)
        tournament = tournament and count == 1
        partial = partial and count <= 1
        semicomplete = semicomplete and count >= 1
        symmetric = symmetric and count != 1
    maximal_partial = partial and len(arcs) == maximal_partial_arc_count
    return tournament, partial, semicomplete, symmetric, maximal_partial


def check_invariance(group: FiniteGroup, arcs: set[tuple[Element, Element]]) -> None:
    for a, x, y in product(group.elements, repeat=3):
        if x == y:
            continue
        ax = group.multiply(a, x)
        ay = group.multiply(a, y)
        require(
            ((x, y) in arcs) == ((ax, ay) in arcs),
            (group.name, "left-translation invariance", a, x, y),
        )


def forced_involution_matchings(
    group: FiniteGroup,
    involutions: Iterable[Element],
) -> tuple[frozenset[frozenset[Element]], ...]:
    matchings: list[frozenset[frozenset[Element]]] = []
    all_edges: set[frozenset[Element]] = set()
    for h in involutions:
        edges = frozenset(
            frozenset((x, group.multiply(x, h))) for x in group.elements
        )
        require(len(edges) == len(group.elements) // 2, (group.name, "matching", h))
        require(all(len(edge) == 2 for edge in edges), (group.name, "loop", h))
        require(not (all_edges & set(edges)), (group.name, "matching overlap", h))
        all_edges.update(edges)
        matchings.append(edges)
    return tuple(matchings)


def canonical_half_set(pairs: tuple[tuple[Element, Element], ...]) -> set[Element]:
    return {min(pair, key=repr) for pair in pairs}


def find_directed_triangle(
    group: FiniteGroup, arcs: set[tuple[Element, Element]]
) -> tuple[Element, Element, Element] | None:
    for x, y, z in combinations(group.elements, 3):
        for a, b, c in ((x, y, z), (x, z, y)):
            if (a, b) in arcs and (b, c) in arcs and (c, a) in arcs:
                return a, b, c
    return None


def audit_group(group: FiniteGroup, enumerate_connections: bool) -> tuple[object, ...]:
    inverse = group.validate()
    involutions, pairs = inverse_packet(group, inverse)
    order = len(group.elements)
    t = len(involutions)
    pair_count = len(pairs)
    require(2 * pair_count + t == order - 1, (group.name, "orbit invoice"))

    matchings = forced_involution_matchings(group, involutions)
    forced_pairs = sum(len(matching) for matching in matchings)
    require(forced_pairs == order * t // 2, (group.name, "forced-pair invoice"))

    predicted = {
        "generalized": 4**pair_count * 2**t,
        "partial": 3**pair_count,
        "semicomplete": 3**pair_count,
        "symmetric": 2 ** (pair_count + t),
        "maximal_partial": 2**pair_count,
        "tournament": 0 if t else 2**pair_count,
    }
    require(predicted["generalized"] == 2 ** (order - 1), (group.name, predicted))

    observed: dict[str, int] | None = None
    if enumerate_connections:
        nonidentity = tuple(x for x in group.elements if x != group.identity)
        counts = {key: 0 for key in predicted}
        for mask in range(1 << len(nonidentity)):
            connection_set = {
                nonidentity[i]
                for i in range(len(nonidentity))
                if (mask >> i) & 1
            }
            arcs = arcs_from_connection_set(group, inverse, connection_set)
            check_invariance(group, arcs)
            tournament, partial, semicomplete, symmetric, maximal_partial = relation_kind(
                group, arcs, order * pair_count
            )
            counts["generalized"] += 1
            counts["tournament"] += int(tournament)
            counts["partial"] += int(partial)
            counts["semicomplete"] += int(semicomplete)
            counts["symmetric"] += int(symmetric)
            counts["maximal_partial"] += int(maximal_partial)
        require(counts == predicted, (group.name, counts, predicted))
        observed = counts

    triangle: tuple[Element, Element, Element] | None = None
    if t == 0 and order > 1:
        connection_set = canonical_half_set(pairs)
        arcs = arcs_from_connection_set(group, inverse, connection_set)
        require(
            relation_kind(group, arcs, order * pair_count)[0],
            (group.name, "half-set is not tournament"),
        )
        check_invariance(group, arcs)
        triangle = find_directed_triangle(group, arcs)
        require(triangle is not None, (group.name, "nontrivial regular tournament is transitive"))

    return (
        group.name,
        order,
        t,
        pair_count,
        forced_pairs,
        tuple(predicted.items()),
        None if observed is None else tuple(observed.items()),
        triangle,
    )


def hensel_lattice_rows() -> tuple[tuple[object, ...], ...]:
    rows: list[tuple[object, ...]] = []
    for p, n, r in ((2, 1, 2), (2, 2, 1), (2, 2, 2), (2, 3, 3), (3, 1, 2), (3, 2, 1), (5, 1, 2)):
        order = p ** (n * r)
        involutions = 0 if p % 2 else 2**r - 1
        inverse_pairs = (order - 1 - involutions) // 2
        forced_pairs = order * involutions // 2
        tournaments = 0 if involutions else 2**inverse_pairs
        rows.append((p, n, r, order, involutions, inverse_pairs, forced_pairs, tournaments))
    require(rows[0][6] == 6, rows[0])  # F_2^2: every pair is forced.
    require(rows[2][6] == 24, rows[2])  # (Z/4)^2: three forced matchings.
    return tuple(rows)


def main() -> None:
    small_groups = (
        cyclic_group(3),
        cyclic_group(4),
        vector_group(2, 2, "V4"),
        cyclic_group(5),
        symmetric_group_three(),
        cyclic_group(6),
        cyclic_group(7),
        vector_group(2, 3, "F2^3"),
        cyclic_group(9),
    )
    small_rows = tuple(audit_group(group, True) for group in small_groups)

    nonabelian_odd = heisenberg_group_three()
    nonabelian_row = audit_group(nonabelian_odd, False)
    require(nonabelian_row[1:5] == (27, 0, 13, 0), nonabelian_row)
    require(dict(nonabelian_row[5])["tournament"] == 8192, nonabelian_row)

    hensel_rows = hensel_lattice_rows()
    size_rows = tuple(row for row in small_rows if row[1] in (4, 6))
    require(tuple((row[0], row[2], row[4]) for row in size_rows) == (
        ("C4", 1, 2),
        ("V4", 3, 6),
        ("S3", 3, 9),
        ("C6", 1, 3),
    ), size_rows)

    verdict = (
        "regular_invariant_tournament_iff_no_nonidentity_involution_iff_odd_group_order",
        "tournament_count_at_odd_order=2^((|G|-1)/2)",
        "each_involution_forces_a_disjoint_missing_or_bidirected_perfect_matching",
        "all_generalized_relations=4^pairs*2^involutions=2^(|G|-1)",
        "partial_and_semicomplete_counts=3^pairs",
        "symmetric_count=2^(pairs+involutions);maximal_partial_count=2^pairs",
        "nonabelian_scope_PASS_on_S3_and_Heis(F3)",
        "Hensel_p_odd_admits_half_set_gauges;p_2_forces_2^r-1_matchings",
        "first_two_adic_lift_is_XOR_all_pairs_forced",
        "no_LRC_current_or_canonical_orientation_follows",
    )
    semantic_surface = (small_rows, nonabelian_row, hensel_rows, size_rows, verdict)
    semantic_sha = digest(semantic_surface)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, (semantic_sha, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    print("THM-3447 regular group-orbit tournament parity exact companion")
    print("status=EXACT_DETERMINISTIC_COMPANION_FOR_THM3447;finite_checks_do_not_replace_proof")
    print("arithmetic=finite_group_tables;exact_integer_enumeration;no_randomness")
    print(f"small_group_rows={small_rows}")
    print(f"nonabelian_odd_control={nonabelian_row}")
    print(f"hensel_lattice_rows={hensel_rows}")
    print("hensel_columns=(p,depth_gap,rank,orbit_size,involutions,inverse_pairs,forced_nonsingle_pairs,invariant_tournaments)")
    print(f"size4_size6_controls={size_rows}")
    print(f"verdict={verdict}")
    print("scope=one_regular_group_orbit;cross_orbit_relations_need_orbit-bank_sidecars;G-invariance_not_full-automorphism-canonicity")
    print("loss_boundary=half-set_is_an_orientation_gauge;missing/bidirected_states_survive_involutions;no_physical_LRC_transport")
    print(f"semantic_sha256={semantic_sha}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;all_checks_survive_python_O")
    print("commands=python -B 04-computation/regular_group_orbit_tournament_parity_thm3447.py;python -B -O 04-computation/regular_group_orbit_tournament_parity_thm3447.py")
    print("PASS")


if __name__ == "__main__":
    main()
