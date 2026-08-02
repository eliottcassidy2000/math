#!/usr/bin/env python3
"""Exact distinguished-B deletion law for THM-3110's two Ewens currents."""

import argparse
import hashlib
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "05-knowledge/results/gmc_product_gamma_distinguished_packet_deletion_cubic_scout.out"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


A, B, A2, AB, B2, A3, A2B, AB2, B3 = sp.symbols(
    "A B A2 AB B2 A3 A2B AB2 B3"
)
MOMENT_VARIABLES = (A, B, A2, AB, B2, A3, A2B, AB2, B3)
LINEAR_FORM = {
    A: (1, 0),
    B: (0, 1),
    A2: (2, 0),
    AB: (1, 1),
    B2: (0, 2),
    A3: (3, 0),
    A2B: (2, 1),
    AB2: (1, 2),
    B3: (0, 3),
}
y, z = sp.symbols("y z")
MOMENT_VALUE = {
    (0, 0): 1,
    (1, 0): A,
    (0, 1): B,
    (2, 0): A2,
    (1, 1): AB,
    (0, 2): B2,
    (3, 0): A3,
    (2, 1): A2B,
    (1, 2): AB2,
    (0, 3): B3,
}


def expectation(polynomial):
    return sp.expand(
        sum(
            coefficient * MOMENT_VALUE[(i, j)] / A**i / B**j
            for (i, j), coefficient in sp.Poly(sp.expand(polynomial), y, z).terms()
        )
    )


U = y - 1
V = z - y
g11, g12, g22 = expectation(U**2), expectation(U * V), expectation(V**2)
t111, t112 = expectation(U**3), expectation(U**2 * V)
t122, t222 = expectation(U * V**2), expectation(V**3)
INVARIANTS = (
    3 * t112 * g11 * g22 - t222 * g11**2 - 2 * t111 * g12 * g22,
    3 * t122 * g11 * g22 - 2 * t222 * g12 * g11 - t111 * g22**2,
)


def bank(invariant, denominator):
    numerator, actual = sp.cancel(invariant).as_numer_denom()
    require(sp.expand(actual - denominator) == 0, "denominator")
    answer = []
    for monomial, coefficient in sp.Poly(-sp.expand(numerator), *MOMENT_VARIABLES).terms():
        row = []
        for variable, exponent in zip(MOMENT_VARIABLES, monomial):
            row.extend([LINEAR_FORM[variable]] * exponent)
        answer.append((int(coefficient), tuple(sorted(row, reverse=True))))
    return tuple(answer)


BANKS = (bank(INVARIANTS[0], A**5 * B**3), bank(INVARIANTS[1], A**5 * B**4))
require(tuple(map(len, BANKS)) == (24, 25), "bank census")


def set_partitions(items):
    answer = []

    def recurse(position, blocks):
        if position == len(items):
            answer.append(tuple(tuple(block) for block in blocks))
            return
        item = items[position]
        for index in range(len(blocks)):
            blocks[index].append(item)
            recurse(position + 1, blocks)
            blocks[index].pop()
        blocks.append([item])
        recurse(position + 1, blocks)
        blocks.pop()

    recurse(0, [])
    return tuple(answer)


def canonical(blocks):
    return tuple(
        sorted(
            (tuple(sorted(block)) for block in blocks if block),
            key=lambda block: (block[0], len(block), block),
        )
    )


def colour_type(partition):
    return tuple(
        sorted(
            (
                (sum(item < 5 for item in block), sum(item >= 5 for item in block))
                for block in partition
            ),
            reverse=True,
        )
    )


def zeta_current(bank_rows, b_count):
    universe = tuple(range(5 + b_count))
    partitions = set_partitions(universe)
    coefficient = {tuple(sorted(row, reverse=True)): value for value, row in bank_rows}
    multiplicity = Counter(colour_type(partition) for partition in partitions)
    omega = tuple(
        (partition, Fraction(coefficient[colour_type(partition)], multiplicity[colour_type(partition)]))
        for partition in partitions
        if coefficient.get(colour_type(partition), 0)
    )
    cache = {}
    zeta = defaultdict(Fraction)
    for partition, weight in omega:
        choices = []
        for block in partition:
            if block not in cache:
                cache[block] = set_partitions(block)
            choices.append(cache[block])
        for refinements in product(*choices):
            refined = canonical(sum((list(piece) for piece in refinements), []))
            zeta[refined] += weight
    return partitions, zeta


P3, W3 = zeta_current(BANKS[0], 3)
P4, W4 = zeta_current(BANKS[1], 4)
require(
    Counter(1 if W3[canonical(p)] > 0 else -1 if W3[canonical(p)] < 0 else 0 for p in P3)
    == Counter({0: 3660, 1: 285, -1: 195}),
    "B1 signs",
)
require(
    Counter(1 if W4[canonical(p)] > 0 else -1 if W4[canonical(p)] < 0 else 0 for p in P4)
    == Counter({0: 19527, -1: 900, 1: 720}),
    "B2 signs",
)

STAR = 8


def singleton_lift(partition):
    return canonical((*partition, (STAR,)))


def attached_lift(partition, block_index):
    blocks = [list(block) for block in partition]
    blocks[block_index].append(STAR)
    return canonical(blocks)


for partition in P3:
    require(W4[singleton_lift(partition)] == -W3[canonical(partition)] / 2, partition)

PUSH = {}
for partition in P3:
    PUSH[canonical(partition)] = sum(
        (W4[attached_lift(partition, index)] for index in range(len(partition))),
        Fraction(),
    )

TYPES = (
    (((2, 0), (1, 0), (1, 0), (1, 0), (0, 3)), Fraction(1, 20), 10),
    (((2, 0), (1, 2), (1, 0), (1, 0), (0, 1)), Fraction(-1, 60), 90),
    (((2, 1), (2, 0), (1, 0), (0, 1), (0, 1)), Fraction(1, 60), 90),
    (((3, 0), (2, 0), (0, 1), (0, 1), (0, 1)), Fraction(-1, 20), 10),
)
observed = Counter((colour_type(partition), weight) for partition, weight in PUSH.items() if weight)
expected = Counter({(kind, weight): count for kind, weight, count in TYPES})
require(observed == expected, (observed, expected))
require(all(8 - len(partition) == 3 for partition, weight in PUSH.items() if weight), "rank-three support")
orbit_masses = tuple(weight * count for _kind, weight, count in TYPES)
require(orbit_masses == (Fraction(1, 2), Fraction(-3, 2), Fraction(3, 2), Fraction(-1, 2)), orbit_masses)
require(sum(orbit_masses) == 0, orbit_masses)

record = (
    tuple(
        (len(rows), sum(c for c, _ in rows if c > 0), sum(-c for c, _ in rows if c < 0))
        for rows in BANKS
    ),
    tuple((kind, weight, count) for kind, weight, count in TYPES),
    orbit_masses,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    lines = [
        "Product-Gamma distinguished-B deletion cubic exact scout",
        "banks=24/25;positive_negative_masses=37/37,39/39",
        "singleton_suspension=W4(pi+star)=-W3(pi)/2;checks=4140",
        "attached_pushforward_support=rank3;nonzero_partitions=200;types=4",
    ]
    for index, (kind, weight, count) in enumerate(TYPES):
        lines.append(
            f"TYPE{index};colour={kind};point_weight={weight};count={count};orbit_mass={weight * count}"
        )
    lines.extend(
        [
            "orbit_mass_vector=1/2,-3/2,3/2,-1/2=one_half_times_(1-T)^3",
            "boundary=this_is_a_deletion_pushforward_identity_not_a_positive_decomposition_or_all_degree_closure",
            "record_sha256=" + hashlib.sha256(repr(record).encode()).hexdigest(),
            "all_exact_controls=PASS",
        ]
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
