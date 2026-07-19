#!/usr/bin/env python3
"""Exact primary for the common-scale-35 AP-centred Hamming-six face.

The program exhausts all labelled supports and all hereditary effective-order
words.  Scalar owner capacity leaves 216 rows.  At each of their 1,296 owner
obligations it retains every order-1/order-5 mask exactly and independently
maximizes the deviations of the order-7/order-35 masks outside that anchor.
This sound upper relaxation is at most 31, strictly below the 35 sheets.

Tournament Analysis uses owner obligations as vertices.  The observable is
the owner-local proof-cost key in the Z/5 or Z/7 quotient gauge; coordinate
order supplies the tie Hamiltonian path.  This tournament is telemetry only:
the absolute Z/5 ceiling is the proof-bearing predicate.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 35
ORDERS = (1, 5, 7, 35)
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
        sum(order % 5 == 0 for order in word) >= 2
        and sum(order % 7 == 0 for order in word) >= 2
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
    require(
        weighted == (5**6 - 1 - 6 * 4) * (7**6 - 1 - 6 * 6),
        "unit-weighted grammar count",
    )
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
            profiles[tuple(word.count(order) for order in ORDERS)] += 1
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
        "Z5": frozenset((1, 5)),
        "Z7": frozenset((1, 7)),
        "mixed": frozenset((1, 5, 7)),
    }
    bounds = {name: Counter() for name in flags}
    banks = {name: Counter() for name in flags}
    live = {name: Counter() for name in flags}
    own_order_bound = Counter()
    gauge_histogram = Counter()
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
            own_order_bound[(word[owner_index], row_bounds["Z5"][owner_index])] += 1

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
            scores = sorted(sum(a == vertex for a, _ in result) for vertex in range(6))
            require(scores == list(range(6)), "owner tournament not transitive")
            return result, ties

        edges5, ties5 = edges("Z5")
        edges7, ties7 = edges("Z7")
        flips = sum(
            ((left, right) in edges5) != ((left, right) in edges7)
            for left, right in combinations(range(6), 2)
        )
        gauge_histogram[(ties5, ties7, flips)] += 1
    return {
        "bounds": bounds,
        "banks": banks,
        "live": live,
        "own_order_bound": own_order_bound,
        "gauge_histogram": gauge_histogram,
        "digests": {name: digest.hexdigest() for name, digest in digests.items()},
    }


def fmt(counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main() -> None:
    cards, choices, rows, base_digest, mask_digest = build_tables()
    words, weighted, grammar_digest = build_grammar()
    scalar = scalar_census(words, cards)
    flags = flag_census(scalar["survivors"], choices)

    require(tuple(int(value) for value in scalar["histogram"]) == (
        8_522, 1_135_668, 1_246_950, 485_440, 114_720, 10_560, 216,
    ), "scalar histogram")
    require(scalar["support_histogram"] == Counter({0: 900, 9: 24}),
            "support histogram")
    require(len(scalar["survivors"]) == 216, "scalar survivor count")
    require(len({support for support, _word, _cap in scalar["survivors"]}) == 24,
            "scalar survivor support count")
    require(scalar["profiles"] == Counter({
        (0, 2, 2, 2): 96,
        (0, 2, 1, 3): 96,
        (0, 2, 0, 4): 24,
    }), "multiplicity profiles")
    require(len(scalar["capacity_vectors"]) == 130, "capacity vector count")
    require(scalar["literal_survivors"] == 286_654_464,
            "literal survivor count")
    require(flags["bounds"]["Z5"] == Counter({
        27: 576, 28: 336, 29: 240, 30: 96, 31: 48,
    }), "Z5 bound histogram")
    require(flags["live"]["Z5"] == Counter({0: 216}), "Z5 terminality")
    require(max(flags["bounds"]["Z5"]) == 31 < C, "Z5 strict ceiling")
    require(sum(flags["gauge_histogram"].values()) == 216,
            "owner tournament count")
    require(grammar_digest ==
            "c2ab055e4dceeafbb7b15aedd1429cde3dc44c6de41ec64ecc8011d2dbb9cce9",
            "grammar digest")
    require(base_digest ==
            "d397fc4e12b528261da7f6f5bccdc80ce1cb7bc720a467229153f27ab4d981f8",
            "CRT-base digest")
    require(mask_digest ==
            "9d2c2d2e7ac16f8bbb45dfa36dda1df82a428cffd6913532a0eefb171ebd2f15",
            "literal-mask digest")
    require(scalar["digest"] ==
            "db8aa0ce7e195ba1ae4024747073022ef5b4d51eb149cf54071a8293f3c9e467",
            "scalar digest")
    require(flags["digests"] == {
        "Z5": "25bc8abff0f769dd80d5c4e8302507586109c4da28e8c1d84c875637ef46235c",
        "Z7": "7acb48a927c55952e604b40fc0ff96ffecf7e686066dd169db7584b49baf0347",
        "mixed": "314d0ae90377f40af17ca87594633052db47a7d0dca8576202d5609867d3cdb6",
    }, "flag digests")

    print("### c=35 AP-centred H6 exact five-fibre primary ###")
    print("orders:", ORDERS)
    print("unit counts:", tuple(len(UNITS[d]) for d in ORDERS))
    print("grammar: at least two 5-carriers and at least two 7-carriers")
    print("hereditary words:", len(words))
    print("labelled support/order contexts:", 924 * len(words))
    print("literal state words/support:", weighted)
    print("raw labelled state words:", 924 * weighted)
    print("grammar digest:", grammar_digest)
    print("CRT-base digest:", base_digest)
    print("literal-mask digest:", mask_digest)
    print("cardinality vectors by ratio; columns D=1,5,7,35:")
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
    for name in ("Z5", "Z7", "mixed"):
        print(f"{name} bank sizes:", fmt(flags["banks"][name]))
        print(f"{name} bounds:", fmt(flags["bounds"][name]))
        print(f"{name} live owners/context:", fmt(flags["live"][name]))
        print(f"{name} digest:", flags["digests"][name])
    print("Z5 owner-order/bound:", fmt(flags["own_order_bound"]))
    print("owner tournament (Z5 ties,Z7 ties,gauge flips):",
          fmt(flags["gauge_histogram"]))
    print("owner tournaments: scores (0,1,2,3,4,5), triangles 0, SCCs 1^6, HP 1")
    print("tournament vertices: owner obligations; gauge: Z5 versus Z7 proof cost")
    print("preserved: owner-local absolute upper bound; destroyed: shared units/needle overlaps")
    print("Z5 ceiling: 31 < 35 on every scalar-surviving owner obligation")
    print("VERDICT: primitive proper AP-centred common-scale-35 H6 face EMPTY")


if __name__ == "__main__":
    main()
