#!/usr/bin/env python3
"""Exact primary for the c=34 AP-centred Hamming-six two-fibre obstruction.

Every labelled support and hereditary order word is visited.  Orders 1 and 2
are retained exactly as thick fibres of Z/34 -> Z/2.  Each order-17/order-34
provider is then allowed to choose its unit independently outside the exact
anchor union.  This is a sound upper relaxation, and its owner bound is at
most 29<34 on every scalar-surviving owner obligation.

Tournament Analysis uses owner obligations rather than runners.  It compares
the Z/2 and Z/17 gauges, with coordinate order as the tie Hamiltonian path; an
alternate three-vertex tournament compares quotient flags.  These completed
tournaments are telemetry only: the absolute Z/2 owner bound is the proof.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 34
ORDERS = (1, 2, 17, 34)
UNITS = {
    d: ((0,) if d == 1 else tuple(u for u in range(1, d) if gcd(u, d) == 1))
    for d in ORDERS
}
FULL = (1 << C) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, order: int, unit: int) -> int:
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def literal_crt_base(label: int, order: int, unit: int) -> int:
    answers = tuple(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )
    require(len(answers) == 1, "literal CRT uniqueness")
    return answers[0]


def local_mask(label: int, order: int, unit: int, owner: int) -> int:
    base = crt_base(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            mask |= 1 << sheet
    return mask


def analytic_cardinality(label: int, order: int, owner: int) -> int:
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    hits = sum(value % P == target for value in range(-order + 1, order + 1))
    return (C // order) * hits


def build_tables():
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    choices = {}
    rows = {}
    base_digest = sha256()
    mask_digest = sha256()
    for label in range(1, P):
        for order_index, order in enumerate(ORDERS):
            for unit in UNITS[order]:
                algebraic = crt_base(label, order, unit)
                literal = literal_crt_base(label, order, unit)
                require(algebraic == literal, "algebraic/literal CRT mismatch")
                base_digest.update(bytes((label, order, unit)))
                base_digest.update(algebraic.to_bytes(2, "little"))
            for owner in range(1, P):
                card = analytic_cardinality(label, order, owner)
                masks = tuple(sorted({
                    local_mask(label, order, unit, owner)
                    for unit in UNITS[order]
                }))
                require(all(mask.bit_count() == card for mask in masks),
                        "mask/cardinality mismatch")
                for mask in masks:
                    mask_digest.update(bytes((label, order, owner)))
                    mask_digest.update(mask.to_bytes(8, "little"))
                    for sheet in range(C):
                        require(
                            ((mask >> sheet) & 1)
                            == ((mask >> ((sheet + order) % C)) & 1),
                            "effective-order periodicity mismatch",
                        )
                    if order in (1, 2):
                        require(all(
                            ((mask >> sheet) & 1)
                            == ((mask >> ((sheet + 2) % C)) & 1)
                            for sheet in range(C)
                        ), "Z/2 anchor periodicity mismatch")
                cards[label, order_index, owner] = card
                choices[label, order, owner] = masks
        rows[label] = tuple(
            analytic_cardinality(label, order, 1) for order in ORDERS
        )
    for label in range(1, P):
        for owner in range(1, P):
            ratio = label * pow(owner, -1, P) % P
            require(
                tuple(int(cards[label, i, owner]) for i in range(len(ORDERS)))
                == rows[ratio],
                "ratio covariance mismatch",
            )
    return cards, choices, rows, base_digest.hexdigest(), mask_digest.hexdigest()


def hereditary(index_word) -> bool:
    word = tuple(ORDERS[i] for i in index_word)
    via_lcm = all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )
    via_carriers = (
        sum(order % 2 == 0 for order in word) >= 2
        and sum(order % 17 == 0 for order in word) >= 2
    )
    require(via_lcm == via_carriers, "hereditary grammar mismatch")
    return via_lcm


def build_grammar():
    words = []
    weighted = 0
    digest = sha256()
    for index_word in product(range(len(ORDERS)), repeat=6):
        if not hereditary(index_word):
            continue
        word = tuple(ORDERS[i] for i in index_word)
        weight = prod(len(UNITS[order]) for order in word)
        words.append(index_word)
        weighted += weight
        digest.update(bytes(index_word))
        digest.update(weight.to_bytes(8, "little"))
    require(len(words) == 57**2 == 3249, "grammar count")
    require(weighted == (2**6 - 1 - 6) * (17**6 - 1 - 6 * 16),
            "weighted grammar count")
    return np.asarray(words, dtype=np.int8), weighted, digest.hexdigest()


def scalar_census(words, cards):
    histogram = np.zeros(7, dtype=np.int64)
    support_histogram = Counter()
    survivors = []
    profiles = Counter()
    capacity_vectors = set()
    literal_survivors = 0
    digest = sha256()
    for support_tuple in combinations(range(1, P), 6):
        support = np.asarray(support_tuple, dtype=np.int8)
        capacities = np.zeros((len(words), 6), dtype=np.int16)
        for provider, label in enumerate(support):
            capacities += cards[label, words[:, provider, None], support[None, :]]
        feasible = (capacities >= C).sum(axis=1)
        histogram += np.bincount(feasible, minlength=7)
        indices = np.flatnonzero(feasible == 6)
        support_histogram[len(indices)] += 1
        for row_index in indices:
            word = tuple(ORDERS[int(i)] for i in words[row_index])
            capacity = tuple(int(v) for v in capacities[row_index])
            survivors.append((support_tuple, word, capacity))
            profile = tuple(word.count(order) for order in ORDERS)
            profiles[profile] += 1
            capacity_vectors.add(capacity)
            literal_survivors += prod(len(UNITS[order]) for order in word)
            digest.update(bytes(support_tuple))
            digest.update(bytes(word))
            digest.update(bytes(capacity))
    return {
        "histogram": histogram,
        "support_histogram": support_histogram,
        "survivors": tuple(survivors),
        "profiles": profiles,
        "capacity_vectors": capacity_vectors,
        "literal_survivors": literal_survivors,
        "digest": digest.hexdigest(),
    }


def relaxation_bound(support, word, owner, choices, anchors):
    bank = frozenset((0,))
    for label, order in zip(support, word):
        if order not in anchors:
            continue
        bank = frozenset(
            partial | option
            for partial in bank
            for option in choices[label, order, owner]
        )
    best = 0
    for anchor_union in bank:
        outside = FULL ^ anchor_union
        bound = anchor_union.bit_count()
        for label, order in zip(support, word):
            if order in anchors:
                continue
            bound += max(
                (mask & outside).bit_count()
                for mask in choices[label, order, owner]
            )
        best = max(best, bound)
    return best, len(bank)


def flag_census(survivors, choices):
    flags = {
        "Z2": frozenset((1, 2)),
        "Z17": frozenset((1, 17)),
        "mixed": frozenset((1, 2, 17)),
    }
    bounds = {name: Counter() for name in flags}
    banks = {name: Counter() for name in flags}
    live = {name: Counter() for name in flags}
    own_order_bound = Counter()
    gauge_histogram = Counter()
    flag_rank_histogram = Counter()
    digests = {name: sha256() for name in flags}
    for support, word, capacities in survivors:
        row_bounds = {name: [] for name in flags}
        row_banks = {name: [] for name in flags}
        for owner in support:
            for name, anchors in flags.items():
                bound, bank_size = relaxation_bound(
                    support, word, owner, choices, anchors
                )
                row_bounds[name].append(bound)
                row_banks[name].append(bank_size)
                bounds[name][bound] += 1
                banks[name][bank_size] += 1
                digests[name].update(bytes(support))
                digests[name].update(bytes(word))
                digests[name].update(bytes((owner, bound)))
                digests[name].update(bank_size.to_bytes(4, "little"))
        for name in flags:
            live[name][sum(value >= C for value in row_bounds[name])] += 1
        for owner_index in range(6):
            own_order_bound[(word[owner_index], row_bounds["Z2"][owner_index])] += 1
            names = tuple(flags)
            ranking = tuple(sorted(
                names,
                key=lambda name: (
                    row_bounds[name][owner_index],
                    row_banks[name][owner_index],
                    names.index(name),
                ),
            ))
            flag_rank_histogram[ranking] += 1

        def edges(name):
            result = set()
            ties = 0
            for left, right in combinations(range(6), 2):
                left_key = (
                    row_bounds[name][left] >= C,
                    row_bounds[name][left],
                    capacities[left],
                    row_banks[name][left],
                )
                right_key = (
                    row_bounds[name][right] >= C,
                    row_bounds[name][right],
                    capacities[right],
                    row_banks[name][right],
                )
                if left_key == right_key:
                    ties += 1
                    winner, loser = left, right
                elif left_key > right_key:
                    winner, loser = left, right
                else:
                    winner, loser = right, left
                result.add((winner, loser))
            scores = sorted(sum(a == vertex for a, _b in result) for vertex in range(6))
            require(scores == list(range(6)), "owner tournament not transitive")
            return result, ties

        edges2, ties2 = edges("Z2")
        edges17, ties17 = edges("Z17")
        flips = sum(
            ((left, right) in edges2) != ((left, right) in edges17)
            for left, right in combinations(range(6), 2)
        )
        gauge_histogram[(ties2, ties17, flips)] += 1
    return {
        "bounds": bounds,
        "banks": banks,
        "live": live,
        "own_order_bound": own_order_bound,
        "gauge_histogram": gauge_histogram,
        "flag_rank_histogram": flag_rank_histogram,
        "digests": {name: value.hexdigest() for name, value in digests.items()},
    }


def fmt(counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main() -> None:
    cards, choices, rows, base_digest, mask_digest = build_tables()
    words, weighted, grammar_digest = build_grammar()
    scalar = scalar_census(words, cards)
    flags = flag_census(scalar["survivors"], choices)

    require(tuple(int(value) for value in scalar["histogram"])
            == (33_112, 387_180, 1_283_004, 1_099_712, 180_156, 18_360, 552),
            "scalar histogram")
    require(scalar["support_histogram"]
            == Counter({0: 888, 15: 24, 16: 12}), "support histogram")
    require(len(scalar["survivors"]) == 552, "survivor count")
    require(len({support for support, _word, _cap in scalar["survivors"]}) == 36,
            "survivor support count")
    require(scalar["profiles"] == Counter({
        (0, 2, 2, 2): 216,
        (0, 2, 3, 1): 144,
        (0, 2, 1, 3): 144,
        (0, 2, 0, 4): 36,
        (0, 2, 4, 0): 12,
    }), "multiplicity profiles")
    require(len(scalar["capacity_vectors"]) == 356, "capacity vector count")
    require(scalar["literal_survivors"] == 36_175_872,
            "literal survivor count")
    require(flags["banks"]["Z2"] == Counter({1: 3_312}), "Z2 bank count")
    require(flags["bounds"]["Z2"] == Counter({
        25: 48, 26: 1_056, 27: 408, 28: 1_728, 29: 72,
    }), "Z2 bound histogram")
    require(flags["live"]["Z2"] == Counter({0: 552}), "Z2 terminality")
    require(max(flags["bounds"]["Z2"]) == 29, "Z2 ceiling")
    require(sum(flags["gauge_histogram"].values()) == 552,
            "owner tournament count")

    print("### c=34 AP-centred H6 exact two-fibre primary ###")
    print("orders:", ORDERS)
    print("unit counts:", tuple(len(UNITS[d]) for d in ORDERS))
    print("grammar: at least two 2-carriers and at least two 17-carriers")
    print("hereditary words:", len(words))
    print("labelled support/order contexts:", 924 * len(words))
    print("literal state words/support:", weighted)
    print("raw labelled state words:", 924 * weighted)
    print("grammar digest:", grammar_digest)
    print("CRT-base digest:", base_digest)
    print("literal-mask digest:", mask_digest)
    print("cardinality vectors by ratio; columns D=1,2,17,34:")
    for ratio in range(1, P):
        print(f"  r={ratio:2d}: {rows[ratio]}")
    print("scalar feasible-owner histogram:",
          " ".join(f"{i}:{int(v)}" for i, v in enumerate(scalar["histogram"])))
    print("scalar support-survivor histogram:", fmt(scalar["support_histogram"]))
    print("scalar survivors:", len(scalar["survivors"]))
    print("scalar survivor supports:",
          len({support for support, _word, _cap in scalar["survivors"]}))
    print("scalar multiplicity profiles:", fmt(scalar["profiles"]))
    print("scalar capacity vectors:", len(scalar["capacity_vectors"]))
    print("literal scalar-survivor state words:", scalar["literal_survivors"])
    print("scalar digest:", scalar["digest"])
    for name in ("Z2", "Z17", "mixed"):
        print(f"{name} bank sizes:", fmt(flags["banks"][name]))
        print(f"{name} bounds:", fmt(flags["bounds"][name]))
        print(f"{name} live owners/context:", fmt(flags["live"][name]))
        print(f"{name} digest:", flags["digests"][name])
    print("Z2 owner-order/bound:", fmt(flags["own_order_bound"]))
    print("owner tournament (Z2 ties,Z17 ties,gauge flips):",
          fmt(flags["gauge_histogram"]))
    print("owner tournaments: score (0,1,2,3,4,5), triangles 0, SCCs 1^6, HP 1")
    print("flag strength rankings:", fmt(flags["flag_rank_histogram"]))
    print("tournament vertices: owner obligations; alternate vertices: quotient flags")
    print("observable/switch: proof cost in Z2 versus Z17 gauge; coordinate tie path")
    print("preserved: owner-local absolute upper bound; destroyed: shared units/needle overlaps")
    print("Z2 ceiling: 29 < 34 on every scalar-surviving owner obligation")
    print("VERDICT: primitive proper AP-centred common-scale-34 H6 face EMPTY")


if __name__ == "__main__":
    main()
