#!/usr/bin/env python3
"""Exact THM-3245 pointed-divisor median-cube and modular no-go referee."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3174 = ROOT / "04-computation/lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.py"
OUTPUT_3174 = ROOT / "05-knowledge/results/lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.out"
SOURCE_3179 = ROOT / "04-computation/lrc14_j7_k3_z222_composite_terminal_descent_cap221_thm3179.py"
OUTPUT_3179 = ROOT / "05-knowledge/results/lrc14_j7_k3_z222_composite_terminal_descent_cap221_thm3179.out"
SOURCE_3207 = ROOT / "04-computation/lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.py"
OUTPUT_3207 = ROOT / "05-knowledge/results/lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.out"
SOURCE_3218 = ROOT / "04-computation/lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.py"
OUTPUT_3218 = ROOT / "05-knowledge/results/lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.out"
SOURCE_3230 = ROOT / "04-computation/lrc14_j7_k3_z219_common_gcd_three_terminal_descent_cap218_thm3230.py"
OUTPUT_3230 = ROOT / "05-knowledge/results/lrc14_j7_k3_z219_common_gcd_three_terminal_descent_cap218_thm3230.out"
RESIDUAL_3230 = ROOT / "05-knowledge/results/lrc14_j7_k3_z219_residual_masks_thm3245.data"

SOURCE_3174_SHA256 = "3ac95e58861078828671e606e1c97705ac4def2728019ac21bd78eba0b9f1c18"
OUTPUT_3174_SHA256 = "4edc6c7efb64ce1062aa863a5f4c17e4735f42bf7b4ac371c6bca660346abe43"
SEMANTIC_3174_SHA256 = "f6493459d78f6af58b999bd54c8b72e2a8f989c533a1e5a91d5e6bc337352ec7"
SOURCE_3179_SHA256 = "f094bb0fa276e1af2e41c2c7db907ad14477c0230ab5a37a7ad9e61c0a350c27"
OUTPUT_3179_SHA256 = "75801a1dc7ba1086b213fba1630639a54be7c17443742335be07e0db53ccbd0e"
SOURCE_3207_SHA256 = "b3d1f0c451087017c1363d42be8789df78d5eec7db7a05b49dc5ca9e194f2091"
OUTPUT_3207_SHA256 = "1e4b01f6e5bfb179ad5fe6ba786124ff2c9fdf3833eed4917c0b3ce0abb7b76d"
SOURCE_3218_SHA256 = "97c55e0a9acc9b42c7c7a856889bd2b8fbd49854ac39a1fe974a488f63a3cff8"
OUTPUT_3218_SHA256 = "af7fe69ba68579611f171b9b2252f9a9923bdd2785720a6cb1b208a611ee0565"
SEMANTIC_3218_SHA256 = "ed67ff55840d07965a6cf8a478c0fc49cab637511ad3d684b2de8831dfcfcd94"
SOURCE_3230_SHA256 = "0b14b1a86e2a08aae20433ee12d601128090f04658c0934fc6b263fb31d9e723"
OUTPUT_3230_SHA256 = "d8a4e87382da55879ee5f517f79277e4c3913823d5f968931ded22f84ace3bcf"
SEMANTIC_3230_SHA256 = "999e0186c706441074e50ca1a2e0d689a9f19c6fa2c13137abc730eb51f92a49"
RESIDUAL_3230_SHA256 = "91a574490b4df60bc5d666e068ef48d2b3af68129991c0a9f32c1f33d9d81d1a"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

EXPECTED_LAYERS = {
    222: (
        (219, 199, 20),
        "68ea5b7128d3deafd4bb9a14f6d1d09ba867a5f117645a3e1f4f3d03e761de73",
        (((2, True), 1), ((6, False), 20), ((6, True), 198)),
    ),
    221: (
        (90, 83, 7),
        "ab4f5b7fe3b5e1330d6e51dfb150527fa27cac56b755ed5e2ae2522f0f27ceb4",
        (((1, False), 7), ((1, True), 83)),
    ),
    220: (
        (289, 249, 40),
        "bc75fef8aebdc6b350957ea96842b38a24bad2174200646ebaf2817844515799",
        (
            ((4, False), 20),
            ((4, True), 52),
            ((20, False), 20),
            ((20, True), 163),
            ((44, True), 5),
            ((220, True), 29),
        ),
    ),
    219: (
        (16, 15, 1),
        "6a39c83d47b8080fd52e64479b886f1b6a887150f32cfbf30a676fbb0aaa54ba",
        (((3, False), 1), ((3, True), 15)),
    ),
    218: (
        (119, 117, 2),
        "d416b484148e008e9f58273f8533cc9d05a171e5ba4963214d3fc73d7fc20bc7",
        (((2, False), 2), ((2, True), 117)),
    ),
    217: (
        (8, 7, 1),
        "885fa20869364ee4a808324612ad373278464a339bf47e183418c1731af4ea6d",
        (((7, False), 1), ((7, True), 7)),
    ),
    216: (
        (480, 447, 33),
        "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649",
        (
            ((8, True), 19),
            ((18, True), 1),
            ((24, False), 24),
            ((24, True), 111),
            ((36, True), 15),
            ((72, False), 9),
            ((72, True), 301),
        ),
    ),
    215: (
        (13, 13, 0),
        "f883ca2c03c242f6cd3cff65782b613b7fae9b8bcf6ac592b8500a457d543ff9",
        (((1, True), 12), ((5, True), 1)),
    ),
}

Z219_POSITION_CENSUS = (
    ((0,), 92),
    ((0, 1), 45),
    ((0, 2), 24),
    ((1,), 136),
    ((1, 2), 52),
    ((1, 3), 4),
    ((2,), 66),
    ((3,), 5),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3174):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def factor(number):
    require(number >= 1, number)
    residue = number
    out = []
    prime = 2
    while prime * prime <= residue:
        exponent = 0
        while residue % prime == 0:
            residue //= prime
            exponent += 1
        if exponent:
            out.append((prime, exponent))
        prime += 1
    if residue > 1:
        out.append((residue, 1))
    return tuple(out)


def valuation(number, prime):
    exponent = 0
    while number % prime == 0:
        number //= prime
        exponent += 1
    return exponent


def divisors(number):
    data = factor(number)
    return tuple(
        sorted(
            product_value(exponents, data)
            for exponents in product(*(range(power + 1) for _prime, power in data))
        )
    )


def product_value(exponents, data):
    value = 1
    for exponent, (prime, _power) in zip(exponents, data):
        value *= prime**exponent
    return value


def coordinates(value, data):
    return tuple(valuation(value, prime) for prime, _power in data)


def threshold_bits(value, data):
    vector = coordinates(value, data)
    return tuple(
        int(exponent >= height)
        for exponent, (_prime, power) in zip(vector, data)
        for height in range(1, power + 1)
    )


def hamming(left, right):
    return sum(a != b for a, b in zip(left, right))


def lattice_distance(left, right, data):
    return sum(abs(a - b) for a, b in zip(coordinates(left, data), coordinates(right, data)))


def coordinate_median(left, middle, right, data):
    exponents = tuple(
        sorted(values)[1]
        for values in zip(
            coordinates(left, data), coordinates(middle, data), coordinates(right, data)
        )
    )
    return product_value(exponents, data)


def lies_between(left, middle, right, data):
    return lattice_distance(left, right, data) == (
        lattice_distance(left, middle, data) + lattice_distance(middle, right, data)
    )


def audit_lattice(g):
    data = factor(g)
    divs = divisors(g)
    grid = tuple(product(*(range(power + 1) for _prime, power in data)))
    expected_size = 1
    for _prime, power in data:
        expected_size *= power + 1
    require(tuple(sorted(coordinates(value, data) for value in divs)) == tuple(sorted(grid)), (g, "grid"))
    require(len(divs) == expected_size == len(grid), (g, "cardinality"))
    require(len(set(coordinates(value, data) for value in divs)) == len(divs), (g, "bijection"))
    for left in divs:
        require(product_value(coordinates(left, data), data) == left, (g, left, "inverse"))
        for right in divs:
            require(
                coordinates(lcm(left, right), data)
                == tuple(max(a, b) for a, b in zip(coordinates(left, data), coordinates(right, data))),
                (g, left, right, "join"),
            )
            require(
                coordinates(gcd(left, right), data)
                == tuple(min(a, b) for a, b in zip(coordinates(left, data), coordinates(right, data))),
                (g, left, right, "meet"),
            )
            require(
                lattice_distance(left, right, data)
                == hamming(threshold_bits(left, data), threshold_bits(right, data)),
                (g, left, right, "isometry"),
            )
    for triple in product(divs, repeat=3):
        median = coordinate_median(*triple, data)
        candidates = tuple(
            value
            for value in divs
            if lies_between(triple[0], value, triple[1], data)
            and lies_between(triple[1], value, triple[2], data)
            and lies_between(triple[2], value, triple[0], data)
        )
        require(candidates == (median,), (g, triple, candidates, median))
    translation_records = []
    identity = tuple(divs)
    for anchor in divs:
        image = tuple(lcm(anchor, value) for value in divs)
        image_twice = tuple(lcm(anchor, value) for value in image)
        bijective = len(set(image)) == len(divs)
        require(image_twice == image, (g, anchor, "idempotence"))
        require(bijective == (anchor == 1), (g, anchor, "internal unit"))
        require((image == identity) == (anchor == 1), (g, anchor, "identity"))
        require((tuple(image[divs.index(value)] for value in image) == identity) == (anchor == 1), (g, anchor, "C2"))
        third = tuple(image[divs.index(value)] for value in image_twice)
        require((third == identity) == (anchor == 1), (g, anchor, "C3"))
        if anchor > 1:
            require(lcm(anchor, 1) == lcm(anchor, anchor), (g, anchor, "collision"))
        translation_records.append((anchor, len(set(image))))
    for left_anchor, right_anchor in product(divs, repeat=2):
        composite = tuple(
            lcm(left_anchor, lcm(right_anchor, value)) for value in divs
        )
        joined = tuple(lcm(lcm(left_anchor, right_anchor), value) for value in divs)
        require(composite == joined, (g, left_anchor, right_anchor, "closure composition"))
    return (
        g,
        data,
        len(divs),
        sum(power for _prime, power in data),
        tuple(translation_records),
    )


def audit_layers():
    wrapper = load("pointed_divisor_wrapper")
    entry = wrapper.load("pointed_divisor_entry")
    source = entry.load("pointed_divisor_source")
    driver = source.load("pointed_divisor_driver")
    predecessor = driver.load("pointed_divisor_predecessor")
    bridge = predecessor.load("pointed_divisor_bridge")
    base = bridge.load("pointed_divisor_base")
    require(base.thm.ATLAS_SHA256 == ATLAS_SHA256, base.thm.ATLAS_SHA256)
    packets = []
    for level in EXPECTED_LAYERS:
        rows = entry.atlas_rows(base, level)
        expected_census, expected_hash, expected_strata = EXPECTED_LAYERS[level]
        census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
        row_hash = hashlib.sha256(repr(rows).encode()).hexdigest()
        strata = tuple(sorted(Counter((gcd(level, row[1]), row[3]) for row in rows).items()))
        require(census == expected_census, (level, census))
        require(row_hash == expected_hash, (level, row_hash))
        require(strata == expected_strata, (level, strata))
        packets.append((level, census, row_hash, strata))
    return tuple(packets)


def bits6(value):
    require(value in (1, 2, 3, 6), value)
    return (int(value % 2 == 0), int(value % 3 == 0))


def from_bits6(bits):
    return (2 if bits[0] else 1) * (3 if bits[1] else 1)


def audit_v4_and_farey():
    divs = (1, 2, 3, 6)
    for shift in product((0, 1), repeat=2):
        image = tuple(from_bits6(tuple(a ^ b for a, b in zip(bits6(value), shift))) for value in divs)
        require(len(set(image)) == 4, (shift, "V4 permutation"))
        require(tuple(from_bits6(tuple(a ^ b for a, b in zip(bits6(value), shift))) for value in image) == divs, (shift, "V4 involution"))
        join_preserving = all(
            image[divs.index(lcm(left, right))]
            == lcm(image[divs.index(left)], image[divs.index(right)])
            for left in divs
            for right in divs
        )
        require(join_preserving == (shift == (0, 0)), (shift, join_preserving))
    atom_swap = {1: 1, 2: 3, 3: 2, 6: 6}
    require(
        all(
            atom_swap[lcm(left, right)] == lcm(atom_swap[left], atom_swap[right])
            and atom_swap[gcd(left, right)] == gcd(atom_swap[left], atom_swap[right])
            for left in divs
            for right in divs
        ),
        "atom swap",
    )
    require(
        {atom_swap[value] for value in (3, 6)} == {2, 6},
        "occupied ambient edge moved",
    )
    actual_left, actual_right = Fraction(3, 6), Fraction(6, 6)
    atom_left, atom_right = Fraction(2, 6), Fraction(3, 6)
    mediant = Fraction(atom_left.numerator + atom_right.numerator, atom_left.denominator + atom_right.denominator)
    z220_left, z220_right = Fraction(2, 20), Fraction(4, 20)
    det = lambda left, right: abs(
        left.numerator * right.denominator - right.numerator * left.denominator
    )
    require(det(actual_left, actual_right) == 1, "accidental z222 Farey edge")
    require(det(atom_left, atom_right) == 1, "ambient atom Farey edge")
    require(mediant == Fraction(2, 5), mediant)
    require(Fraction(lcm(2, 3), 6) == 1 and Fraction(lcm(2, 3), 6) != mediant, "join is not mediant")
    require(det(z220_left, z220_right) == 5, "z220 non-Farey cover")
    cover_records = []
    for g in (2, 4, 6, 20, 44, 220):
        data = factor(g)
        for q in divisors(g):
            for prime, top_exponent in data:
                if valuation(q, prime) >= top_exponent:
                    continue
                upper = prime * q
                left = Fraction(q, g)
                right = Fraction(upper, g)
                determinant = det(left, right)
                expected = (prime - 1) * g // (prime * q)
                require(determinant == expected, (g, q, prime, determinant, expected))
                require((determinant == 1) == (prime == 2 and upper == g), (g, q, prime, determinant))
                cover_records.append((g, q, upper, determinant))
    return (
        tuple((shift, tuple(from_bits6(tuple(a ^ b for a, b in zip(bits6(value), shift))) for value in divs)) for shift in product((0, 1), repeat=2)),
        tuple(atom_swap[value] for value in divs),
        (actual_left, actual_right, det(actual_left, actual_right)),
        (atom_left, atom_right, mediant),
        (z220_left, z220_right, det(z220_left, z220_right)),
        tuple(cover_records),
    )


def audit_two_state_prime_boundary():
    records = []
    for prime in (2, 3):
        divs = (1, prime)
        identity = divs
        restoration = tuple(lcm(prime, value) for value in divs)
        flip = (prime, 1)
        flip_twice = tuple(flip[divs.index(value)] for value in flip)
        require(restoration == (prime, prime), (prime, restoration))
        require(len(set(restoration)) == 1, (prime, "restoration collision"))
        require(flip_twice == identity, (prime, "external flip order"))
        require(flip[0] != 1, (prime, "external flip moves bottom"))
        join_preserving = all(
            flip[divs.index(lcm(left, right))]
            == lcm(flip[divs.index(left)], flip[divs.index(right)])
            for left in divs
            for right in divs
        )
        require(not join_preserving, (prime, "external flip falsely preserves join"))
        records.append((prime, len(divs), restoration, flip))
    require(len(divisors(3)) == 2, "Div(3) is not a three-state set")
    require(len(divisors(4)) == 3, "Div(4) is not a three-state set")
    return tuple(records)


def audit_z219_sorted_position_graph(output3230):
    require(
        "residual_quotient_counts:((3, 3, 424),)" in output3230,
        "z219 top quotient census",
    )
    require(lf_sha(RESIDUAL_3230) == RESIDUAL_3230_SHA256, "z219 residual data")
    thm3230 = load("thm3245_supplier_thm3230", SOURCE_3230)
    records = {}
    for line in RESIDUAL_3230.read_text(encoding="utf-8").splitlines():
        if not line.startswith("RECORD\t"):
            continue
        fields = line.split("\t", 3)
        require(len(fields) == 4, line)
        index = int(fields[1])
        body = ast.literal_eval(fields[2])
        bank = ast.literal_eval(fields[3])
        require(index not in records, (index, "duplicate residual record"))
        records[index] = (body, bank)
    require(len(records) == len(thm3230.RESIDUAL) == 10, len(records))
    require({body for body, _bank in records.values()} == set(thm3230.RESIDUAL), "residual bodies")

    position_counts = Counter()
    quotient_counts = Counter()
    for _index, (body, bank) in sorted(records.items()):
        expected_count, expected_sha = thm3230.RESIDUAL[body]
        require(len(bank) == expected_count, (body, len(bank), expected_count))
        require(hashlib.sha256(repr(bank).encode()).hexdigest() == expected_sha, (body, "bank hash"))
        ruler = 14 * lcm(*body)
        require(gcd(219, ruler) == 3, (body, ruler))
        first_d = ruler // 3
        top = valuation(ruler, 3)
        for ds in bank:
            require(lcm(*ds) == ruler, (body, ds, "full ruler"))
            quotient_counts[(3, lcm(*ds) // first_d)] += 1
            top_positions = tuple(index for index, value in enumerate(ds) if valuation(value, 3) == top)
            require(top_positions, (body, ds, "missing top position"))
            position_counts[top_positions] += 1
    require(tuple(sorted(quotient_counts.items())) == (((3, 3), 424),), quotient_counts)
    require(tuple(sorted(position_counts.items())) == Z219_POSITION_CENSUS, position_counts)

    weighted = dict(position_counts)
    support = frozenset(weighted)
    vertices = frozenset(item for subset in support for item in subset)
    edges = frozenset(subset for subset in support if len(subset) == 2)
    require(vertices == frozenset(range(4)), vertices)
    require(
        edges == frozenset(((0, 1), (0, 2), (1, 2), (1, 3))),
        edges,
    )

    def image(subset, permutation):
        return tuple(sorted(permutation[item] for item in subset))

    unweighted = []
    weighted_aut = []
    for permutation in permutations(range(4)):
        if frozenset(image(subset, permutation) for subset in support) == support:
            unweighted.append(permutation)
        if all(weighted.get(image(subset, permutation)) == count for subset, count in weighted.items()):
            weighted_aut.append(permutation)
    require(tuple(unweighted) == ((0, 1, 2, 3), (2, 1, 0, 3)), unweighted)
    require(tuple(weighted_aut) == ((0, 1, 2, 3),), weighted_aut)

    reduced_support = frozenset(subset for subset in support if 3 not in subset)
    reduced_unweighted = []
    reduced_weighted = []
    for permutation in permutations(range(3)):
        if frozenset(image(subset, permutation) for subset in reduced_support) == reduced_support:
            reduced_unweighted.append(permutation)
        if all(
            weighted.get(image(subset, permutation)) == weighted[subset]
            for subset in reduced_support
        ):
            reduced_weighted.append(permutation)
    require(len(reduced_unweighted) == 6, reduced_unweighted)
    require(tuple(reduced_weighted) == ((0, 1, 2),), reduced_weighted)
    require(weighted[(3,)] + weighted[(1, 3)] == 9, "real slot-3 mass")
    return (
        tuple(sorted(position_counts.items())),
        tuple(sorted(edges)),
        tuple(unweighted),
        tuple(weighted_aut),
        tuple(reduced_unweighted),
        tuple(reduced_weighted),
    )


def top_signature(value, ruler):
    return tuple(valuation(value, prime) == valuation(ruler, prime) for prime in (2, 3))


def audit_supplier_and_cross_level():
    body = (1, 9, 10, 11, 12, 14)
    ruler = 14 * lcm(*body)
    separate = (2, 9702, 21560, 32340)
    coupled = (72, 5390, 9702, 32340)
    require(ruler == 194040, ruler)
    for ds in (separate, coupled):
        require(lcm(*ds) == ruler, ds)
        require(lcm(*ds) // (ruler // gcd(222, ruler)) == 6, ds)
    separate_signatures = tuple(top_signature(value, ruler) for value in separate)
    coupled_signatures = tuple(top_signature(value, ruler) for value in coupled)
    require(any(a for a, _b in separate_signatures) and any(b for _a, b in separate_signatures), separate_signatures)
    require(not any(a and b for a, b in separate_signatures), separate_signatures)
    require(any(a and b for a, b in coupled_signatures), coupled_signatures)

    shared_body = (1, 6, 8, 10, 11, 14)
    shared_ruler = 14 * lcm(*shared_body)
    shared_ds = (4620, 6468, 21560, 129360)
    require(shared_ruler == lcm(*shared_ds) == 129360, shared_ruler)
    q222 = lcm(*shared_ds) // (shared_ruler // gcd(222, shared_ruler))
    q221 = lcm(*shared_ds) // (shared_ruler // gcd(221, shared_ruler))
    require((q222, q221) == (6, 1), (q222, q221))
    return (
        (body, ruler, separate, separate_signatures, coupled, coupled_signatures),
        (shared_body, shared_ruler, shared_ds, q222, q221),
    )


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def permutation_power(value, exponent):
    result = tuple(range(len(value)))
    for _ in range(exponent):
        result = compose(value, result)
    return result


def audit_legal_sidecar_control():
    states = tuple(permutations((0, 1, 2)))
    index = {state: position for position, state in enumerate(states)}
    transposition = (1, 0, 2)
    three_cycle = (1, 2, 0)

    def left_action(group_element):
        return tuple(index[compose(group_element, state)] for state in states)

    involution = left_action(transposition)
    order_three = left_action(three_cycle)
    identity = tuple(range(len(states)))
    require(permutation_power(involution, 2) == identity, "legal S^2")
    require(permutation_power(order_three, 3) == identity, "legal R^3")
    require(involution != identity and order_three != identity, "legal exact orders")
    require(compose(involution, order_three) != compose(order_three, involution), "legal noncommutation")
    require(len(set(involution)) == len(states) and len(set(order_three)) == len(states), "legal bijections")
    return involution, order_three


def audit_external_lattice_action_hostile():
    divs = divisors(30)
    primes = (2, 3, 5)

    def bits(value):
        return tuple(valuation(value, prime) for prime in primes)

    def value(vector):
        return product_value(vector, tuple((prime, 1) for prime in primes))

    swap = tuple(value((bits(q)[1], bits(q)[0], bits(q)[2])) for q in divs)
    cycle = tuple(value((bits(q)[2], bits(q)[0], bits(q)[1])) for q in divs)
    identity = tuple(divs)

    def map_compose(left, right):
        return tuple(left[divs.index(item)] for item in right)

    require(map_compose(swap, swap) == identity, "external lattice S^2")
    require(map_compose(cycle, map_compose(cycle, cycle)) == identity, "external lattice R^3")
    require(map_compose(swap, cycle) != map_compose(cycle, swap), "external lattice noncommutation")
    for mapping in (swap, cycle):
        require(
            all(
                mapping[divs.index(lcm(left, right))]
                == lcm(mapping[divs.index(left)], mapping[divs.index(right)])
                and mapping[divs.index(gcd(left, right))]
                == gcd(mapping[divs.index(left)], mapping[divs.index(right)])
                for left in divs
                for right in divs
            ),
            "external lattice automorphism",
        )
    return swap, cycle


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    dependency_paths = (
        (SOURCE_3174, SOURCE_3174_SHA256),
        (OUTPUT_3174, OUTPUT_3174_SHA256),
        (SOURCE_3179, SOURCE_3179_SHA256),
        (OUTPUT_3179, OUTPUT_3179_SHA256),
        (SOURCE_3207, SOURCE_3207_SHA256),
        (OUTPUT_3207, OUTPUT_3207_SHA256),
        (SOURCE_3218, SOURCE_3218_SHA256),
        (OUTPUT_3218, OUTPUT_3218_SHA256),
        (SOURCE_3230, SOURCE_3230_SHA256),
        (OUTPUT_3230, OUTPUT_3230_SHA256),
        (RESIDUAL_3230, RESIDUAL_3230_SHA256),
    )
    for path, expected in dependency_paths:
        require(lf_sha(path) == expected, (path, lf_sha(path), expected))
    require(
        f"semantic_sha256={SEMANTIC_3174_SHA256}" in OUTPUT_3174.read_text(),
        "THM3174 semantic",
    )
    output3179 = OUTPUT_3179.read_text()
    output3207 = OUTPUT_3207.read_text()
    output3218 = OUTPUT_3218.read_text()
    output3230 = OUTPUT_3230.read_text()
    require("residual_quotients:(((6, 3), 38), ((6, 6), 1012))" in output3179, "z222 quotient census")
    require("coupling:(((3, False), 38), ((6, False), 700), ((6, True), 312))" in output3179, "z222 coupling census")
    require("fixed_modulus=gcd(221,L):1;first_d_equals_L:1" in output3207, "z221 coprime reset")
    require(f"semantic_sha256={SEMANTIC_3218_SHA256}" in output3218, "THM3218 semantic")
    require(f"semantic_sha256={SEMANTIC_3230_SHA256}" in output3230, "THM3230 semantic")

    lattice_packet = tuple(
        audit_lattice(g) for g in (1, 2, 3, 4, 5, 6, 7, 8, 18, 20, 24, 36, 44, 72, 220)
    )
    layer_packet = audit_layers()
    v4_farey_packet = audit_v4_and_farey()
    prime_boundary_packet = audit_two_state_prime_boundary()
    z219_position_packet = audit_z219_sorted_position_graph(output3230)
    supplier_packet = audit_supplier_and_cross_level()
    sidecar_packet = audit_legal_sidecar_control()
    external_action_packet = audit_external_lattice_action_hostile()
    packet = (
        "lrc-pointed-divisor-join-lattice-modular-no-go-v4",
        tuple(expected for _path, expected in dependency_paths),
        ATLAS_SHA256,
        lattice_packet,
        layer_packet,
        v4_farey_packet,
        prime_boundary_packet,
        z219_position_packet,
        supplier_packet,
        sidecar_packet,
        external_action_packet,
    )
    semantic = hashlib.sha256(repr(packet).encode()).hexdigest()

    lines = []
    emit = lines.append
    emit("LRC THM3245 pointed-divisor median-cube and modular no-go")
    emit(f"dependencies=THM3174:{SOURCE_3174_SHA256}/{OUTPUT_3174_SHA256};THM3179:{SOURCE_3179_SHA256}/{OUTPUT_3179_SHA256};THM3207:{SOURCE_3207_SHA256}/{OUTPUT_3207_SHA256};THM3218:{SOURCE_3218_SHA256}/{OUTPUT_3218_SHA256}/{SEMANTIC_3218_SHA256};THM3230:{SOURCE_3230_SHA256}/{OUTPUT_3230_SHA256}/{SEMANTIC_3230_SHA256};z219_residual_data:{RESIDUAL_3230_SHA256};atlas:{ATLAS_SHA256}")
    for g, data, size, cube_dimension, translations in lattice_packet:
        emit(f"lattice=g:{g};factor:{data};divisors:{size};threshold_dimension:{cube_dimension};median:1;internal_join_restorations:{translations};commuting_idempotent_band:1;nontrivial_units:0")
    for level, census, row_hash, strata in layer_packet:
        emit(f"layer=z:{level};census:{census};row_sha256:{row_hash};gcd_wall_strata:{strata}")
    emit("z222_residual=q3:38;q6_separate:700;q6_coupled:312;occupied_vertices:(3,6);induced_ambient_Hasse_edge:1;no_transition_asserted:1")
    emit("z221_residual=quotient_singleton:1;coprime_reset:1")
    emit("z219_residual=g:3;q:3;q_over_g:1;occupied_quotient_support:(3,);C3_action:0")
    emit("z218_metadata=g:2;rows:119;ambient_only:1;screened_residual_action:0")
    emit("lower_gcd_sequence=z219:(3,);z218:(2,);z217:(7,);z216:(8,18,24,36,72);z215:(1,5);alternating_free_factor_law:0")
    emit(f"v4_farey=v4_translations:{v4_farey_packet[0]};atom_swap:{v4_farey_packet[1]};actual_edge:{v4_farey_packet[2]};atom_mediant:{v4_farey_packet[3]};z220_cover:{v4_farey_packet[4]};divisor_cover_controls:{len(v4_farey_packet[5])}")
    emit(f"prime_two_state_boundary={prime_boundary_packet};Div3_states:2;Div4_states:3;unpointed_flip_not_join_preserving:1")
    emit(f"z219_sorted_position_support=census:{z219_position_packet[0]};edges:{z219_position_packet[1]};unweighted_graph_aut:{z219_position_packet[2]};weighted_graph_aut:{z219_position_packet[3]};deleted_position3_unweighted_aut_count:{len(z219_position_packet[4])};deleted_position3_weighted_aut:{z219_position_packet[5]};position3_pattern_mass:9;physical_supplier_labels_discarded:1")
    emit(f"supplier_hostile=body:{supplier_packet[0][0]};L:{supplier_packet[0][1]};separate:{supplier_packet[0][2]}/{supplier_packet[0][3]};coupled:{supplier_packet[0][4]}/{supplier_packet[0][5]}")
    emit(f"cross_level_arithmetic=body:{supplier_packet[1][0]};L:{supplier_packet[1][1]};ds:{supplier_packet[1][2]};q222:{supplier_packet[1][3]};q221:{supplier_packet[1][4]}")
    emit(f"legal_sidecar_control=common_carrier:6;S:{sidecar_packet[0]};R:{sidecar_packet[1]};S2:1;R3:1;noncommute:1")
    emit(f"external_action_hostile=Div30_B3;S:{external_action_packet[0]};R:{external_action_packet[1]};lattice_automorphisms:1;S2:1;R3:1;noncommute:1")
    emit("scope_internal=pointed_divisor_intervals_are_median_partial_cubes;internal_join_restorations_form_a_commuting_idempotent_band_with_trivial_group_image")
    emit("scope_external=external_lattice_automorphisms_can_exist_but_are_not_internal_restorations_or_physical_carrier_maps;unpointed_XOR_flips_destroy_pointed_join_semantics")
    emit("scope_supplier=q_and_sorted_denominator_positions_record_top_nonemptiness_and_truncated_valuation_height_but_forget_physical_supplier_labels_intersections_and_order;position_graph_automorphisms_are_not_carrier_symmetries")
    emit("scope_physical=no_PSL2_action_or_C2_C3_torsor_without_one_common_carrier_and_coefficient_phase_owner_ancestry_preservation;no_THM2056_Farey_flank_claim;no_new_LRC_decrement")
    emit(f"semantic_sha256={semantic}")
    emit("all_exact_controls=PASS")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    if args.output is not None:
        args.output.write_text(text, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
