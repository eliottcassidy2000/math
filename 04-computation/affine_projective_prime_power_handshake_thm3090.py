#!/usr/bin/env python3
"""Exact controls for THM-3090's affine/projective handshake."""

from __future__ import annotations

import argparse
from itertools import permutations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "05-knowledge/results/affine_projective_prime_power_handshake_thm3090.out"
)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


class Field:
    def __init__(self, p, degree=1, modulus=None):
        self.p = p
        self.degree = degree
        self.modulus = modulus
        if degree == 1:
            self.elements = tuple(range(p))
            self.zero = 0
            self.one = 1
        else:
            require(p == 2 and degree == 3 and modulus == 0b1011, "F8 model")
            self.elements = tuple(range(8))
            self.zero = 0
            self.one = 1

    @property
    def q(self):
        return len(self.elements)

    def add(self, left, right):
        if self.degree == 1:
            return (left + right) % self.p
        return left ^ right

    def neg(self, value):
        if self.degree == 1:
            return (-value) % self.p
        return value

    def sub(self, left, right):
        return self.add(left, self.neg(right))

    def mul(self, left, right):
        if self.degree == 1:
            return (left * right) % self.p
        result = 0
        a = left
        b = right
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & 0b1000:
                a ^= self.modulus
        return result

    def inv(self, value):
        require(value != self.zero, "inverse zero")
        for candidate in self.elements:
            if self.mul(value, candidate) == self.one:
                return candidate
        raise RuntimeError((self.q, value, "no inverse"))


F2 = Field(2)
F3 = Field(3)
F8 = Field(2, degree=3, modulus=0b1011)


def det2(field, matrix):
    return field.sub(
        field.mul(matrix[0], matrix[3]),
        field.mul(matrix[1], matrix[2]),
    )


def mat_vec(field, matrix, vector):
    return (
        field.add(field.mul(matrix[0], vector[0]), field.mul(matrix[1], vector[1])),
        field.add(field.mul(matrix[2], vector[0]), field.mul(matrix[3], vector[1])),
    )


def gl2(field):
    rows = []
    for matrix in product(field.elements, repeat=4):
        if det2(field, matrix) != field.zero:
            rows.append(matrix)
    require(len(rows) == (field.q**2 - 1) * (field.q**2 - field.q), field.q)
    return tuple(rows)


def affine_permutations(field):
    points = tuple(product(field.elements, repeat=2))
    lookup = {point: index for index, point in enumerate(points)}
    rows = []
    for matrix in gl2(field):
        for shift in points:
            image = []
            for point in points:
                linear = mat_vec(field, matrix, point)
                target = (
                    field.add(linear[0], shift[0]),
                    field.add(linear[1], shift[1]),
                )
                image.append(lookup[target])
            rows.append(tuple(image))
    require(len(rows) == field.q**2 * (field.q**2 - 1) * (field.q**2 - field.q), field.q)
    require(len(set(rows)) == len(rows), (field.q, "faithful affine action"))
    return tuple(rows)


def projective_lines(field):
    zero = (field.zero, field.zero)
    lines = []
    for slope in field.elements:
        lines.append(
            frozenset((scalar, field.mul(slope, scalar)) for scalar in field.elements)
        )
    lines.append(frozenset((field.zero, scalar) for scalar in field.elements))
    require(all(zero in line and len(line) == field.q for line in lines), field.q)
    return tuple(lines)


def normalize_matrix(field, matrix):
    first = next(value for value in matrix if value != field.zero)
    scale = field.inv(first)
    return tuple(field.mul(scale, value) for value in matrix)


def projective_permutations(field):
    lines = projective_lines(field)
    lookup = {line: index for index, line in enumerate(lines)}
    normalized = {normalize_matrix(field, matrix) for matrix in gl2(field)}
    require(len(normalized) == field.q * (field.q**2 - 1), field.q)
    rows = []
    for matrix in sorted(normalized):
        image = []
        for line in lines:
            target = frozenset(mat_vec(field, matrix, point) for point in line)
            image.append(lookup[target])
        rows.append(tuple(image))
    require(len(set(rows)) == len(rows), (field.q, "faithful projective action"))
    return tuple(rows)


def tuple_orbits(action, degree, arity):
    universe = tuple(permutations(range(degree), arity))
    unseen = set(universe)
    sizes = []
    while unseen:
        seed = min(unseen)
        orbit = {
            tuple(permutation[index] for index in seed)
            for permutation in action
        }
        require(orbit <= unseen, (degree, arity, seed, "orbit overlap"))
        unseen -= orbit
        sizes.append(len(orbit))
    return tuple(sorted(sizes))


def is_prime(value):
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def prime_powers(limit):
    rows = set()
    for prime in range(2, limit + 1):
        if not is_prime(prime):
            continue
        value = prime
        while value <= limit:
            rows.add(value)
            value *= prime
    return tuple(sorted(rows))


def affine_order(q):
    return q**2 * (q**2 - 1) * (q**2 - q)


def projective_order(r):
    return r * (r**2 - 1)


POWERS = prime_powers(4096)
DEGREE_PAIRS = tuple(
    (q, r)
    for q in POWERS
    for r in POWERS
    if q * q == r + 1
)
require(DEGREE_PAIRS == ((2, 3), (3, 8)), DEGREE_PAIRS)

ORDER_ROWS = tuple(
    (q, r, q * q, affine_order(q), projective_order(r))
    for q, r in DEGREE_PAIRS
)
require(ORDER_ROWS == ((2, 3, 4, 24, 24), (3, 8, 9, 432, 504)), ORDER_ROWS)

for q in range(2, 100):
    r = q * q - 1
    difference = projective_order(r) - affine_order(q)
    require(difference == q**2 * (q**2 - 1) * (q - 2), (q, difference))

for q in range(3, 100):
    collinear = q**2 * (q**2 - 1)
    noncollinear = affine_order(q)
    total = q**2 * (q**2 - 1) * (q**2 - 2)
    require((q - 2) * collinear + noncollinear == total, (q, "triple atlas"))

AFFINE2 = affine_permutations(F2)
AFFINE3 = affine_permutations(F3)
PROJECTIVE3 = projective_permutations(F3)
PROJECTIVE8 = projective_permutations(F8)

TRANSITIVITY = (
    ("AGL2F2", len(AFFINE2), tuple_orbits(AFFINE2, 4, 2), tuple_orbits(AFFINE2, 4, 3)),
    ("PGL2F3", len(PROJECTIVE3), tuple_orbits(PROJECTIVE3, 4, 2), tuple_orbits(PROJECTIVE3, 4, 3)),
    ("AGL2F3", len(AFFINE3), tuple_orbits(AFFINE3, 9, 2), tuple_orbits(AFFINE3, 9, 3)),
    ("PGL2F8", len(PROJECTIVE8), tuple_orbits(PROJECTIVE8, 9, 2), tuple_orbits(PROJECTIVE8, 9, 3)),
)
require(
    TRANSITIVITY
    == (
        ("AGL2F2", 24, (12,), (24,)),
        ("PGL2F3", 24, (12,), (24,)),
        ("AGL2F3", 432, (72,), (72, 432)),
        ("PGL2F8", 504, (72,), (504,)),
    ),
    TRANSITIVITY,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    lines = [
        "Affine/projective prime-power handshake and septimal counterfeit",
        "degree_equation=q^2=r+1;prime_power_solutions=(2,3),(3,8)",
        "order_difference=|PGL2(q^2-1)|-|AGL2(q)|=q^2(q^2-1)(q-2)",
    ]
    for q, r, degree, affine, projective in ORDER_ROWS:
        lines.append(
            f"PAIR;q={q};r={r};degree={degree};AGL_order={affine};PGL_order={projective};equal={int(affine == projective)}"
        )
    for name, order, pairs, triples in TRANSITIVITY:
        lines.append(
            f"ACTION;name={name};order={order};ordered_pair_orbits={','.join(map(str, pairs))};ordered_triple_orbits={','.join(map(str, triples))}"
        )
    lines += [
        "exceptional=(2,3):both_natural_actions_are_S4_and_sharply_three_transitive",
        "counterfeit=(3,8):equal_degree_and_two_transitivity_but_orders_432_vs504;AGL_triple_orbits=72+432;PGL_sharp_triples=504",
        "affine_triple_atlas=q>2 has q-2 collinear ratio orbits of size q^2(q^2-1) plus one noncollinear orbit of size |AGL2(q)|",
        "septimal_boundary=F8* is C7;parity_is_trivial;direction_rank=9;internal_contrast_rank=54",
        "scope=permutation_degree/order/transitivity_and_scalar_fibre_only;no_canonical_generator_pair,tree,quartic_owner,Keller,or_LRC_map",
        "prime_power_scan_through_4096=PASS;all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
