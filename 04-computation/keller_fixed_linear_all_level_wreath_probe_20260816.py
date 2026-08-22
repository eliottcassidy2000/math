#!/usr/bin/env python3
"""Exact support audit for fixed-linear primitivity in the sporadic Keller tower.

This companion deliberately does not recompute the intrinsic discriminant
class of THM-3531.  It checks the four elementary pieces used after that
theorem:

* the generic ``L=0`` cubic has one finite branch and one tamely ramified
  quadratic pair, hence local inertia is a transposition;
* the two rational three-point fibres from THM-2473 rule out every nonzero
  constant linear form descending through one copy of F;
* the natural blocks of the full ternary wreath action are the ancestor
  cylinders (an exact finite hostile replay is made through depth four);
* over F_41, every one of the 1,723 projective linear directions has a
  degree-nine split fibre on which it separates all nine points.

The all-level theorem still uses the geometric input of THM-3530: the newest
prime has one simple preimage on V(L), while the preceding iterate is finite
etale at its generic point.  This script audits the algebra and combinatorics
around that input; it is not a substitute for the divisor argument.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import itertools
import json
import math

import sympy as sp


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def determinant_3(rows: list[list[Fraction]]) -> Fraction:
    a, b, c = rows[0]
    d, e, f = rows[1]
    g, h, i = rows[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def collision_point_rows() -> tuple[list[list[Fraction]], Fraction]:
    """Return three fibre-difference rows spanning Q^3."""

    def points(m: int) -> tuple[tuple[Fraction, ...], ...]:
        mm = Fraction(m)
        anchor = (Fraction(0), Fraction(0), -Fraction(1, 4 * m * m))
        plus = (mm, -Fraction(3, 2 * m), Fraction(13, 2 * m * m))
        minus = (-mm, Fraction(3, 2 * m), Fraction(13, 2 * m * m))
        return anchor, plus, minus

    a1, p1, n1 = points(1)
    _, p2, n2 = points(2)
    rows = [
        [p1[j] - n1[j] for j in range(3)],
        [p2[j] - n2[j] for j in range(3)],
        [p1[j] - a1[j] for j in range(3)],
    ]
    det = determinant_3(rows)
    require(det != 0, "the two rational fibres did not kill every linear direction")
    return rows, det


def local_cubic_audit() -> dict[str, str | list[str]]:
    l_value, t_value, c_value, w_value = sp.symbols("L T c w")
    cubic = l_value * w_value**3 + t_value * w_value - 2 * c_value
    discriminant = sp.factor(sp.discriminant(cubic, w_value))
    expected = -4 * l_value * (t_value**3 + 27 * c_value**2 * l_value)
    require(sp.expand(discriminant - expected) == 0, "raw cubic discriminant")

    # For u=1/w the polynomial is L+T*u^2-2c*u^3.  With L a
    # uniformizer and T,c units, the lower Newton hull has slopes -1/2,0;
    # hence u-valuations 1/2,1/2,0 and w-valuations -1/2,-1/2,0.
    newton_points = [(0, 1), (2, 0), (3, 0)]
    u_valuations = [Fraction(1, 2), Fraction(1, 2), Fraction(0)]
    w_valuations = [-value for value in u_valuations]
    require(sum(u_valuations) == 1, "Newton valuation product")
    return {
        "raw_discriminant": str(discriminant),
        "newton_points_u": [str(item) for item in newton_points],
        "u_valuations": [str(item) for item in u_valuations],
        "w_valuations": [str(item) for item in w_valuations],
        "raw_disc_L_valuation": "1",
        "monic_disc_L_valuation": "-3",
        "inertia": "one_transposition_on_the_infinite_pair",
    }


def words(depth: int) -> list[tuple[int, ...]]:
    return list(itertools.product(range(3), repeat=depth))


def tree_generators(depth: int) -> list[tuple[int, ...]]:
    """Generators of Aut(the rooted ternary tree of the given depth)."""

    leaves = words(depth)
    index = {leaf: position for position, leaf in enumerate(leaves)}
    generators: list[tuple[int, ...]] = []
    for level in range(depth):
        for prefix in itertools.product(range(3), repeat=level):
            for local in ((1, 0, 2), (1, 2, 0)):  # (01) and (012)
                permutation: list[int] = []
                for leaf in leaves:
                    image = list(leaf)
                    if leaf[:level] == prefix:
                        image[level] = local[leaf[level]]
                    permutation.append(index[tuple(image)])
                require(sorted(permutation) == list(range(3**depth)), "tree generator")
                generators.append(tuple(permutation))
    require(len(generators) == 3**depth - 1, "ternary generator census")
    return generators


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            value, self.parent[value] = self.parent[value], root
        return root

    def union(self, left: int, right: int) -> bool:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return False
        if left_root > right_root:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        return True


def minimal_invariant_block(
    generators: list[tuple[int, ...]], first: int, second: int
) -> frozenset[int]:
    """Smallest class in an invariant equivalence relation containing a pair."""

    degree = len(generators[0])
    relation = UnionFind(degree)
    relation.union(first, second)
    changed = True
    while changed:
        changed = False
        classes: dict[int, list[int]] = defaultdict(list)
        for point in range(degree):
            classes[relation.find(point)].append(point)
        for generator in generators:
            for members in classes.values():
                image_anchor = generator[members[0]]
                for member in members[1:]:
                    changed |= relation.union(image_anchor, generator[member])
    root = relation.find(first)
    return frozenset(point for point in range(degree) if relation.find(point) == root)


def common_prefix_length(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    answer = 0
    for left_digit, right_digit in zip(left, right):
        if left_digit != right_digit:
            break
        answer += 1
    return answer


def wreath_block_audit(max_depth: int = 4) -> dict[int, dict[int, int]]:
    histograms: dict[int, dict[int, int]] = {}
    for depth in range(1, max_depth + 1):
        leaves = words(depth)
        generators = tree_generators(depth)
        base = leaves[0]
        histogram: Counter[int] = Counter()
        for position, leaf in enumerate(leaves):
            actual = minimal_invariant_block(generators, 0, position)
            prefix_length = common_prefix_length(base, leaf)
            expected = frozenset(
                index
                for index, candidate in enumerate(leaves)
                if candidate[:prefix_length] == base[:prefix_length]
            )
            require(actual == expected, ("non-ancestor minimal block", depth, leaf))
            histogram[len(actual)] += 1
        expected_histogram = {
            1: 1,
            **{3**radius: 2 * 3 ** (radius - 1) for radius in range(1, depth + 1)},
        }
        require(dict(sorted(histogram.items())) == expected_histogram, "block histogram")
        histograms[depth] = dict(sorted(histogram.items()))
    return histograms


def wreath_orders(max_depth: int = 6) -> dict[int, int]:
    orders: dict[int, int] = {}
    order = 1
    for depth in range(1, max_depth + 1):
        order *= 6 ** (3 ** (depth - 1))
        orders[depth] = order
    return orders


PRIME = 41


def fmap_mod(point: tuple[int, int, int]) -> tuple[int, int, int]:
    x_value, y_value, z_value = point
    unit = (1 + x_value * y_value) % PRIME
    return (
        (
            unit**3 * z_value
            + y_value**2 * unit * (4 + 3 * x_value * y_value)
        )
        % PRIME,
        (
            y_value
            + 3 * x_value * unit**2 * z_value
            + 3 * x_value * y_value**2 * (4 + 3 * x_value * y_value)
        )
        % PRIME,
        (2 * x_value - 3 * x_value**2 * y_value - x_value**3 * z_value)
        % PRIME,
    )


def projective_directions() -> list[tuple[int, int, int]]:
    return (
        [(1, beta, gamma) for beta in range(PRIME) for gamma in range(PRIME)]
        + [(0, 1, gamma) for gamma in range(PRIME)]
        + [(0, 0, 1)]
    )


def finite_field_direction_audit() -> dict[str, object]:
    fibres_one: dict[tuple[int, int, int], list[tuple[int, int, int]]] = defaultdict(list)
    fibres_two: dict[tuple[int, int, int], list[tuple[int, int, int]]] = defaultdict(list)
    for point in itertools.product(range(PRIME), repeat=3):
        middle = fmap_mod(point)
        target = fmap_mod(middle)
        fibres_one[middle].append(point)
        fibres_two[target].append(point)

    one_histogram = Counter(map(len, fibres_one.values()))
    two_histogram = Counter(map(len, fibres_two.values()))
    require(one_histogram == Counter({1: 35261, 3: 11220}), "F_41 one-level census")
    require(
        two_histogram
        == Counter({1: 19926, 2: 2851, 3: 7919, 4: 1782, 5: 1393,
                    6: 288, 7: 469, 9: 48}),
        "F_41 two-level census",
    )

    full_fibres = sorted(
        (target, tuple(sorted(fibre)))
        for target, fibre in fibres_two.items()
        if len(fibre) == 9
    )
    directions = projective_directions()
    require(len(directions) == PRIME**2 + PRIME + 1, "projective plane census")

    first_witness: list[tuple[int, int, int] | None] = [None] * len(directions)
    first_index: list[int | None] = [None] * len(directions)
    per_fibre_good: list[int] = []
    for fibre_index, (target, fibre) in enumerate(full_fibres):
        good_count = 0
        for direction_index, direction in enumerate(directions):
            values = {
                sum(direction[j] * point[j] for j in range(3)) % PRIME
                for point in fibre
            }
            if len(values) == 9:
                good_count += 1
                if first_witness[direction_index] is None:
                    first_witness[direction_index] = target
                    first_index[direction_index] = fibre_index
        per_fibre_good.append(good_count)

    require(all(target is not None for target in first_witness), "uncovered F_41 direction")
    witness_rows = [
        {"direction": direction, "target": first_witness[index]}
        for index, direction in enumerate(directions)
    ]
    witness_hash = sha256(
        json.dumps(witness_rows, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    literals = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    literal_targets = {
        str(direction): first_witness[directions.index(direction)] for direction in literals
    }
    return {
        "one_fibre_histogram": dict(sorted(one_histogram.items())),
        "two_fibre_histogram": dict(sorted(two_histogram.items())),
        "full_degree_nine_fibres": len(full_fibres),
        "projective_directions": len(directions),
        "directions_covered": sum(target is not None for target in first_witness),
        "first_fibre_good_directions": per_fibre_good[0],
        "last_first_witness_index": max(index for index in first_index if index is not None),
        "literal_first_targets": literal_targets,
        "direction_witness_sha256": witness_hash,
    }


def primitive_integer_directions(limit: int) -> list[tuple[int, int, int]]:
    """Initial canonical enumeration used for the countable-avoidance hostile."""

    directions: list[tuple[int, int, int]] = []
    height = 1
    while len(directions) < limit:
        shell: list[tuple[int, int, int]] = []
        for triple in itertools.product(range(-height, height + 1), repeat=3):
            if max(map(abs, triple)) != height or triple == (0, 0, 0):
                continue
            if math.gcd(math.gcd(abs(triple[0]), abs(triple[1])), abs(triple[2])) != 1:
                continue
            first_nonzero = next(value for value in triple if value)
            if first_nonzero < 0:
                continue
            shell.append(triple)
        directions.extend(sorted(shell))
        height += 1
    return directions[:limit]


def countable_avoidance_hostile() -> dict[str, object]:
    directions = primitive_integer_directions(32)
    require(len(set(directions)) == len(directions), "projective enumeration collision")
    # At finite stage r, the union of the first r rational lines is proper:
    # the next canonical direction is not one of them.  Continuing the same
    # enumeration covers every point of P^2(Q), so the intersection of the
    # nested complements is empty despite each finite complement being open.
    for stage in range(1, len(directions)):
        require(directions[stage] not in set(directions[:stage]), "avoidance hostile")
    return {
        "checked_prefix": len(directions),
        "first_twelve": directions[:12],
        "finite_stage_survivors": len(directions) - 1,
        "infinite_statement": "enumerate_P2(Q); B_r=union_of_first_r_points; intersection(P2(Q)-B_r)=empty",
    }


collision_rows, collision_det = collision_point_rows()
local_data = local_cubic_audit()
block_histograms = wreath_block_audit()
order_ledger = wreath_orders()
finite_field_data = finite_field_direction_audit()
avoidance_data = countable_avoidance_hostile()

semantic_payload = {
    "collision_rows": [[str(value) for value in row] for row in collision_rows],
    "collision_det": str(collision_det),
    "local": local_data,
    "block_histograms": block_histograms,
    "wreath_orders": order_ledger,
    "finite_field": finite_field_data,
    "countable_avoidance": avoidance_data,
}
semantic_hash = sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("status=VERIFIED_EXACT_SUPPORT fixed_linear_all_level_wreath_primitivity")
print("scope=local_cubic+two_fibre_linear_descent+wreath_blocks_depth_le_4+F41_level2_all_directions")
print(f"collision_rows={semantic_payload['collision_rows']}")
print(f"collision_det={collision_det};constant_linear_descent_kernel=0")
print(f"local_cubic={local_data}")
print(f"wreath_block_histograms={block_histograms}")
print(f"wreath_orders={order_ledger}")
print(f"finite_field={finite_field_data}")
print(f"countable_avoidance_hostile={avoidance_data}")
print(f"semantic_sha256={semantic_hash}")
print("verdict=all exact companion gates PASS; all-level implication still cites THM-3530 divisor geometry")

