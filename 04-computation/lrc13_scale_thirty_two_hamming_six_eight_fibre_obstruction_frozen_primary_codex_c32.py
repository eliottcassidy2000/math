#!/usr/bin/env python3
"""Frozen exact primary for the c=32 AP-centred Hamming-six Z/8 obstruction.

The proof-facing object is THM-994's anchor/nonanchor upper relaxation.  Every
effective-order 1, 2, 4, or 8 provider is retained exactly as a union of the
eight thick fibres of Z/32 -> Z/8.  Each order-16 or order-32 provider is then
allowed to choose its exact-order residue independently so as to maximize its
contribution outside the retained anchor union.  This only enlarges the
attainable literal union.  The resulting bound reaches 32 at no more than two
of the six owners in every labelled scalar row.

NumPy is used only for exact int16 batching of the 924*12281 scalar contexts.
CRT bases, literal masks, cyclic owner gauges, anchor banks, upper bounds,
hashes, multiplication orbits, and tournament telemetry use deterministic
Python integers and immutable tuples.  No orbit quotient discards a row.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 32
ORDERS = (1, 2, 4, 8, 16, 32)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL = (1 << C) - 1
MOD8_ANCHORS = frozenset((1, 2, 4, 8))

# Filled and enforced after the first deterministic freeze.
EXPECTED_BASE_DIGEST = "68b88f331468047800d59975eaabd0bdf02af4567bc70ffe8a4afb0984a12b74"
EXPECTED_MASK_DIGEST = "f715f217aae1be03edc31d51a4eb4380024e815c67cc73983bd321beff8c5b81"
EXPECTED_GRAMMAR_DIGEST = "ea035511ad805d5211005c01f3c230882cf5a5349a5a0a11830221bc7538740f"
EXPECTED_SCALAR_DIGEST = "aba8c8c9d589103d83c135f512bec5375ded7997adc166c77b9051cc2b2c106a"
EXPECTED_RELAXATION_DIGEST = "c8f577c67d078c3c0fc70f642ecf63a5ba7ab89e93056997002387bd53fda079"
EXPECTED_CONTEXT_DIGEST = "f08d8ddfeeba5d1d532e3e08a112f23210b0709a1dacf3575deef369a5ea4e43"


def require(condition: bool, message: str) -> None:
    """An optimization-stable assertion."""
    if not condition:
        raise RuntimeError(message)


def format_counter(counter: Counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def literal_crt_base(label: int, order: int, unit: int) -> int:
    answers = tuple(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )
    require(len(answers) == 1, "literal CRT representative is not unique")
    return answers[0]


def algebraic_crt_base(label: int, order: int, unit: int) -> int:
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def literal_mask(label: int, order: int, unit: int, owner: int) -> int:
    base = literal_crt_base(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            mask |= 1 << sheet
    return mask


def rotate32(mask: int, amount: int) -> int:
    amount %= C
    if amount == 0:
        return mask
    return ((mask << amount) | (mask >> (C - amount))) & FULL


def analytic_cardinality(ratio: int, order: int) -> int:
    target = order * ratio % P
    one_period = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * one_period


EXPECTED_CARDINALITY_ROWS = {
    1: (32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    2: (16, 0, 0, 0, 0, 16, 16, 0, 0, 0, 0, 0),
    4: (8, 0, 8, 8, 0, 8, 8, 0, 8, 8, 0, 0),
    8: (8, 4, 4, 8, 4, 4, 4, 4, 8, 4, 4, 4),
    16: (6, 4, 4, 6, 6, 4, 4, 6, 6, 4, 4, 4),
    32: (5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4),
}


def build_tables():
    """Build normalized options, then audit every labelled owner gauge."""
    normalized = {}
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    base_digest = sha256()
    mask_digest = sha256()
    labelled_masks = 0

    for ratio in range(1, P):
        for order_index, order in enumerate(ORDERS):
            options = tuple(sorted({
                literal_mask(ratio, order, unit, 1)
                for unit in UNITS[order]
            }))
            expected = analytic_cardinality(ratio, order)
            require(
                {mask.bit_count() for mask in options} == {expected},
                "normalized mask cardinality depends on exact-order residue",
            )
            normalized[ratio, order] = options
            cards[ratio, order_index, 1] = expected

    for label in range(1, P):
        for order_index, order in enumerate(ORDERS):
            for owner in range(1, P):
                ratio = label * pow(owner, -1, P) % P
                cards[label, order_index, owner] = analytic_cardinality(
                    ratio, order
                )
                shift = next(
                    candidate
                    for candidate in range(C)
                    if (pow(owner, -1, P) + P * candidate) % C == 1
                )
                actual_options = []
                for unit in UNITS[order]:
                    literal = literal_crt_base(label, order, unit)
                    algebraic = algebraic_crt_base(label, order, unit)
                    require(literal == algebraic, "literal/algebraic CRT mismatch")
                    if owner == 1:
                        base_digest.update(bytes((label, order, unit)))
                        base_digest.update(literal.to_bytes(2, "little"))

                    actual = literal_mask(label, order, unit, owner)
                    expected = rotate32(
                        literal_mask(ratio, order, unit, 1), shift
                    )
                    require(actual == expected, "cyclic owner-gauge mismatch")
                    require(
                        actual.bit_count() == cards[label, order_index, owner],
                        "literal mask/cardinality mismatch",
                    )
                    actual_options.append(actual)
                    labelled_masks += 1
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(actual.to_bytes(4, "little"))

                    for sheet in range(C):
                        require(
                            ((actual >> sheet) & 1)
                            == ((actual >> ((sheet + order) % C)) & 1),
                            "order-period law failed",
                        )
                        if order in MOD8_ANCHORS:
                            require(
                                ((actual >> sheet) & 1)
                                == ((actual >> ((sheet + 8) % C)) & 1),
                                "mod-eight anchor is not a thick-fibre union",
                            )

                require(
                    len(set(actual_options)) == len(normalized[ratio, order]),
                    "owner gauge changed option-bank size",
                )

    rows = {
        order: tuple(
            analytic_cardinality(ratio, order) for ratio in range(1, P)
        )
        for order in ORDERS
    }
    require(rows == EXPECTED_CARDINALITY_ROWS, "cardinality ledger changed")
    require(labelled_masks == 4_608, "labelled literal-mask count changed")
    base_hex = base_digest.hexdigest()
    mask_hex = mask_digest.hexdigest()
    if EXPECTED_BASE_DIGEST != "PENDING":
        require(base_hex == EXPECTED_BASE_DIGEST, "CRT-base digest changed")
    if EXPECTED_MASK_DIGEST != "PENDING":
        require(mask_hex == EXPECTED_MASK_DIGEST, "literal-mask digest changed")
    return {
        "normalized": normalized,
        "cards": cards,
        "rows": rows,
        "labelled_masks": labelled_masks,
        "base_digest": base_hex,
        "mask_digest": mask_hex,
    }


def hereditary(index_word: tuple[int, ...]) -> bool:
    order_word = tuple(ORDERS[index] for index in index_word)
    return all(
        lcm(*(order_word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )


def valuation_grammar(index_word: tuple[int, ...]) -> bool:
    return sum(ORDERS[index] == 32 for index in index_word) >= 2


def build_grammar():
    words = []
    weighted_states = 0
    digest = sha256()
    order32_histogram = Counter()
    for index_word in product(range(len(ORDERS)), repeat=6):
        via_lcm = hereditary(index_word)
        via_valuation = valuation_grammar(index_word)
        require(via_lcm == via_valuation, "lcm and valuation grammars disagree")
        if not via_lcm:
            continue
        words.append(index_word)
        order_word = tuple(ORDERS[index] for index in index_word)
        weight = prod(len(UNITS[order]) for order in order_word)
        weighted_states += weight
        order32_histogram[order_word.count(32)] += 1
        digest.update(bytes(index_word))
        digest.update(weight.to_bytes(8, "little"))

    require(
        len(words) == 6**6 - 5**6 - 6 * 5**5 == 12_281,
        "hereditary word count changed",
    )
    require(
        weighted_states == 32**6 - 7 * 16**6 == 956_301_312,
        "literal state-word count changed",
    )
    require(
        order32_histogram
        == Counter({2: 9_375, 3: 2_500, 4: 375, 5: 30, 6: 1}),
        "order-thirty-two grammar histogram changed",
    )
    digest_hex = digest.hexdigest()
    if EXPECTED_GRAMMAR_DIGEST != "PENDING":
        require(digest_hex == EXPECTED_GRAMMAR_DIGEST,
                "hereditary grammar digest changed")
    return {
        "words": np.asarray(words, dtype=np.int8),
        "weighted_states": weighted_states,
        "order32_histogram": order32_histogram,
        "digest": digest_hex,
    }


EXPECTED_MULTIPLICITIES = Counter({
    (0, 0, 1, 1, 2, 2): 672,
    (0, 2, 0, 1, 1, 2): 432,
    (0, 0, 1, 0, 3, 2): 384,
    (0, 0, 1, 2, 1, 2): 312,
    (0, 2, 0, 2, 0, 2): 216,
    (0, 2, 0, 0, 2, 2): 216,
    (0, 0, 1, 0, 1, 4): 168,
    (0, 3, 1, 0, 0, 2): 144,
    (0, 2, 0, 1, 0, 3): 144,
    (0, 2, 0, 0, 1, 3): 144,
    (0, 1, 3, 0, 0, 2): 96,
    (0, 2, 1, 0, 1, 2): 96,
    (0, 0, 1, 0, 0, 5): 72,
    (0, 0, 2, 1, 1, 2): 72,
    (0, 0, 1, 1, 0, 4): 48,
    (0, 0, 1, 1, 1, 3): 48,
    (0, 0, 1, 0, 2, 3): 48,
    (0, 0, 2, 2, 0, 2): 36,
    (0, 2, 0, 0, 0, 4): 36,
    (0, 0, 4, 0, 0, 2): 30,
    (0, 0, 3, 1, 0, 2): 12,
    (0, 0, 3, 0, 1, 2): 12,
    (0, 0, 1, 3, 0, 2): 12,
})


def scalar_census(words: np.ndarray, cards: np.ndarray):
    survivors = []
    feasible_owner_histogram = np.zeros(7, dtype=np.int64)
    support_histogram = Counter()
    multiplicities = Counter()
    capacity_vectors = set()
    minimum_slack = Counter()
    maximum_slack = Counter()
    tight_owner_histogram = Counter()
    digest = sha256()

    for support_tuple in combinations(range(1, P), 6):
        support = np.asarray(support_tuple, dtype=np.int8)
        capacities = np.zeros((len(words), 6), dtype=np.int16)
        for provider, label in enumerate(support):
            capacities += cards[
                label, words[:, provider, None], support[None, :]
            ]
        feasible_count = (capacities >= C).sum(axis=1)
        feasible_owner_histogram += np.bincount(feasible_count, minlength=7)
        survivor_indices = np.flatnonzero(feasible_count == 6)
        support_histogram[len(survivor_indices)] += 1
        for row_index in survivor_indices:
            order_word = tuple(
                ORDERS[int(index)] for index in words[row_index]
            )
            capacity = tuple(int(value) for value in capacities[row_index])
            survivors.append((support_tuple, order_word, capacity))
            multiplicities[tuple(
                order_word.count(order) for order in ORDERS
            )] += 1
            capacity_vectors.add(capacity)
            minimum_slack[min(capacity) - C] += 1
            maximum_slack[max(capacity) - C] += 1
            tight_owner_histogram[sum(value == C for value in capacity)] += 1
            digest.update(bytes(support_tuple))
            digest.update(bytes(order_word))
            digest.update(bytes(capacity))

    require(
        tuple(int(value) for value in feasible_owner_histogram)
        == (76_548, 2_800_212, 4_692_582, 2_946_408,
            743_040, 85_404, 3_450),
        "scalar feasible-owner histogram changed",
    )
    require(
        support_histogram
        == Counter({
            0: 640, 1: 96, 2: 48, 3: 30, 17: 24, 22: 24,
            33: 12, 34: 24, 38: 24, 54: 2,
        }),
        "support survivor histogram changed",
    )
    require(len(survivors) == 3_450, "scalar survivor count changed")
    require(
        len({support for support, _word, _capacity in survivors}) == 284,
        "scalar survivor support count changed",
    )
    require(all(1 not in word for _support, word, _capacity in survivors),
            "order-one scalar survivor appeared")
    require(multiplicities == EXPECTED_MULTIPLICITIES,
            "scalar multiplicity profiles changed")
    require(len(capacity_vectors) == 1_649,
            "scalar capacity-vector count changed")
    require(
        minimum_slack == Counter({0: 1_932, 1: 384, 2: 954, 3: 180}),
        "scalar minimum-slack histogram changed",
    )
    require(
        maximum_slack
        == Counter({1: 144, 2: 504, 3: 576, 4: 348, 5: 504,
                    6: 480, 7: 528, 9: 150, 10: 216}),
        "scalar maximum-slack histogram changed",
    )
    require(
        tight_owner_histogram
        == Counter({0: 1_518, 1: 24, 2: 960, 3: 756, 4: 192}),
        "scalar tight-owner histogram changed",
    )
    literal_survivor_words = sum(
        prod(len(UNITS[order]) for order in order_word)
        for _support, order_word, _capacity in survivors
    )
    require(literal_survivor_words == 621_084_672,
            "literal scalar-survivor unit-word count changed")
    digest_hex = digest.hexdigest()
    if EXPECTED_SCALAR_DIGEST != "PENDING":
        require(digest_hex == EXPECTED_SCALAR_DIGEST,
                "scalar survivor digest changed")
    return {
        "survivors": tuple(survivors),
        "feasible_owner_histogram": tuple(
            int(value) for value in feasible_owner_histogram
        ),
        "support_histogram": support_histogram,
        "multiplicities": multiplicities,
        "capacity_vectors": capacity_vectors,
        "minimum_slack": minimum_slack,
        "maximum_slack": maximum_slack,
        "tight_owner_histogram": tight_owner_histogram,
        "literal_survivor_words": literal_survivor_words,
        "digest": digest_hex,
    }


def normalized_owner_key(support, order_word, owner):
    owner_inverse = pow(owner, -1, P)
    return tuple(sorted(
        (label * owner_inverse % P, order)
        for label, order in zip(support, order_word)
    ))


def anchor_union_bank(key, masks):
    bank = frozenset((0,))
    for ratio, order in key:
        if order not in MOD8_ANCHORS:
            continue
        bank = frozenset(
            partial | option
            for partial in bank
            for option in masks[ratio, order]
        )
    return bank


def relaxation_bound(key, anchor_bank, masks):
    """Return U8 and the number of checked transverse option inequalities."""
    best = 0
    option_checks = 0
    for anchor_union in anchor_bank:
        bound = anchor_union.bit_count()
        for ratio, order in key:
            if order in MOD8_ANCHORS:
                continue
            outside_sizes = tuple(
                (option & ~anchor_union & FULL).bit_count()
                for option in masks[ratio, order]
            )
            maximum = max(outside_sizes)
            require(all(size <= maximum for size in outside_sizes),
                    "transverse maximum is not an upper bound")
            option_checks += len(outside_sizes)
            bound += maximum
        best = max(best, bound)
    return best, option_checks


def evaluate_quotient(unique_keys, masks):
    bounds = {}
    bank_sizes = {}
    total_anchor_unions = 0
    total_option_checks = 0
    digest = sha256()
    for key in unique_keys:
        bank = anchor_union_bank(key, masks)
        bound, option_checks = relaxation_bound(key, bank, masks)
        bounds[key] = bound
        bank_sizes[key] = len(bank)
        total_anchor_unions += len(bank)
        total_option_checks += option_checks
        for ratio, order in key:
            digest.update(bytes((ratio, order)))
        digest.update(bytes((bound,)))
        digest.update(len(bank).to_bytes(2, "little"))
        for anchor_union in sorted(bank):
            digest.update(anchor_union.to_bytes(4, "little"))
    digest_hex = digest.hexdigest()
    if EXPECTED_RELAXATION_DIGEST != "PENDING":
        require(digest_hex == EXPECTED_RELAXATION_DIGEST,
                "relaxation digest changed")
    return {
        "bounds": bounds,
        "bank_sizes": bank_sizes,
        "total_anchor_unions": total_anchor_unions,
        "total_option_checks": total_option_checks,
        "digest": digest_hex,
    }


def multiply_context(context, multiplier):
    return tuple(sorted(
        (multiplier * label % P, order) for label, order in context
    ))


def orbit_histogram(contexts):
    contexts = set(contexts)
    remaining = set(contexts)
    histogram = Counter()
    while remaining:
        context = min(remaining)
        orbit = {
            multiply_context(context, multiplier)
            for multiplier in range(1, P)
        }
        require(orbit <= contexts, "multiplication orbit escapes its class")
        remaining -= orbit
        histogram[len(orbit)] += 1
    return histogram


def tournament_fingerprint(owner_rows):
    """Complete the owner-obligation preorder along the coordinate tie path."""
    adjacency = [[False] * 6 for _ in range(6)]
    ties = 0
    flips = 0
    for left, right in combinations(range(6), 2):
        if owner_rows[left] == owner_rows[right]:
            ties += 1
            winner, loser = left, right
        elif owner_rows[left] > owner_rows[right]:
            winner, loser = left, right
        else:
            flips += 1
            winner, loser = right, left
        adjacency[winner][loser] = True

    scores = tuple(sorted(sum(row) for row in adjacency))
    triangles = sum(
        int(
            (adjacency[a][b] and adjacency[b][c] and adjacency[c][a])
            or (adjacency[a][c] and adjacency[c][b] and adjacency[b][a])
        )
        for a, b, c in combinations(range(6), 3)
    )
    paths = [[0] * 6 for _mask in range(1 << 6)]
    for last in range(6):
        paths[1 << last][last] = 1
    for mask in range(1, 1 << 6):
        for last in range(6):
            if not (mask >> last) & 1:
                continue
            previous_mask = mask ^ (1 << last)
            for previous in range(6):
                if ((previous_mask >> previous) & 1
                        and adjacency[previous][last]):
                    paths[mask][last] += paths[previous_mask][previous]
    require(scores == (0, 1, 2, 3, 4, 5), "tournament scores changed")
    require(triangles == 0, "owner tournament has a directed triangle")
    require(sum(paths[-1]) == 1, "owner tournament Hamiltonian paths changed")
    return ties, flips


def fibre_census(survivors, masks):
    key_rows = []
    context_rows = []
    for support, order_word, _capacities in survivors:
        keys = tuple(
            normalized_owner_key(support, order_word, owner)
            for owner in support
        )
        key_rows.append(keys)
        context_rows.append(tuple(sorted(zip(support, order_word))))

    key_frequency = Counter(key for keys in key_rows for key in keys)
    unique_keys = tuple(sorted(key_frequency))
    require(len(unique_keys) == 1_725, "normalized owner-key count changed")
    require(set(key_frequency.values()) == {12},
            "owner keys do not have uniform multiplicity twelve")
    quotient = evaluate_quotient(unique_keys, masks)
    bound_by_key = quotient["bounds"]
    bank_by_key = quotient["bank_sizes"]

    owner_bound = Counter()
    owner_bank = Counter()
    capacity_minus_bound = Counter()
    live_contexts = Counter()
    dead_margin = Counter()
    live_class_contexts = {0: set(), 1: set(), 2: set()}
    tournament_ties = Counter()
    tournament_flips = Counter()
    context_digest = sha256()

    for (support, order_word, capacities), keys, context in zip(
        survivors, key_rows, context_rows
    ):
        bounds = tuple(bound_by_key[key] for key in keys)
        bank_sizes = tuple(bank_by_key[key] for key in keys)
        live = sum(value >= C for value in bounds)
        require(live <= 2, "a scalar row has more than two live owners")
        live_contexts[live] += 1
        live_class_contexts[live].add(context)

        owner_rows = tuple(
            (
                int(bounds[index] >= C),
                bounds[index],
                capacities[index],
                bank_sizes[index],
            )
            for index in range(6)
        )
        ties, flips = tournament_fingerprint(owner_rows)
        tournament_ties[ties] += 1
        tournament_flips[flips] += 1

        for bound, bank_size, capacity in zip(bounds, bank_sizes, capacities):
            require(bound <= capacity, "relaxation exceeds scalar capacity")
            owner_bound[bound] += 1
            owner_bank[bank_size] += 1
            capacity_minus_bound[capacity - bound] += 1
            if bound < C:
                dead_margin[C - bound] += 1

        context_digest.update(bytes(support))
        context_digest.update(bytes(order_word))
        context_digest.update(bytes(capacities))
        context_digest.update(bytes(bounds))
        for bank_size in bank_sizes:
            context_digest.update(bank_size.to_bytes(2, "little"))

    require(
        owner_bound
        == Counter({
            20: 24, 21: 396, 22: 420, 23: 564, 24: 1_416,
            25: 2_352, 26: 3_252, 27: 3_780, 28: 3_708,
            29: 2_208, 30: 1_260, 31: 480, 32: 624,
            33: 192, 34: 24,
        }),
        "mod-eight owner-bound histogram changed",
    )
    require(
        owner_bank
        == Counter({
            1: 8_064, 2: 4_008, 3: 636, 4: 3_576, 6: 2_088,
            7: 240, 8: 192, 9: 36, 10: 684, 12: 984,
            14: 168, 26: 24,
        }),
        "mod-eight anchor-bank histogram changed",
    )
    require(
        live_contexts == Counter({0: 2_802, 1: 456, 2: 192}),
        "mod-eight live-owner/context histogram changed",
    )
    require(
        capacity_minus_bound
        == Counter({
            2: 1_104, 3: 480, 4: 1_560, 5: 2_400, 6: 3_504,
            7: 2_400, 8: 2_976, 9: 960, 10: 2_448,
            11: 336, 12: 1_548, 13: 216, 14: 360,
            15: 288, 16: 108, 20: 12,
        }),
        "scalar/relaxation loss histogram changed",
    )
    require(
        dead_margin
        == Counter({
            1: 480, 2: 1_260, 3: 2_208, 4: 3_708,
            5: 3_780, 6: 3_252, 7: 2_352, 8: 1_416,
            9: 564, 10: 420, 11: 396, 12: 24,
        }),
        "dead-owner deficit histogram changed",
    )
    require(
        tournament_ties
        == Counter({0: 1_536, 1: 888, 2: 732, 3: 204, 4: 48, 7: 42}),
        "tournament tie histogram changed",
    )
    require(
        tournament_flips
        == Counter({
            0: 11, 1: 18, 2: 53, 3: 117, 4: 253, 5: 372,
            6: 582, 7: 651, 8: 562, 9: 396, 10: 220,
            11: 127, 12: 59, 13: 19, 14: 8, 15: 2,
        }),
        "tournament flip histogram changed",
    )

    orbit_histograms = {
        live: orbit_histogram(contexts)
        for live, contexts in live_class_contexts.items()
    }
    require(
        orbit_histograms
        == {
            0: Counter({6: 3, 12: 232}),
            1: Counter({12: 38}),
            2: Counter({12: 16}),
        },
        "live-class multiplication orbits changed",
    )

    context_hex = context_digest.hexdigest()
    require(context_hex == EXPECTED_CONTEXT_DIGEST,
            "labelled context-bound digest changed")
    return {
        "unique_keys": unique_keys,
        "key_frequency": key_frequency,
        "owner_bound": owner_bound,
        "owner_bank": owner_bank,
        "capacity_minus_bound": capacity_minus_bound,
        "dead_margin": dead_margin,
        "live_contexts": live_contexts,
        "orbit_histograms": orbit_histograms,
        "tournament_ties": tournament_ties,
        "tournament_flips": tournament_flips,
        "total_anchor_unions": quotient["total_anchor_unions"],
        "total_option_checks": quotient["total_option_checks"],
        "relaxation_digest": quotient["digest"],
        "context_digest": context_hex,
    }


def main() -> None:
    tables = build_tables()
    grammar = build_grammar()
    scalar = scalar_census(grammar["words"], tables["cards"])
    fibre = fibre_census(scalar["survivors"], tables["normalized"])

    print("scale-thirty-two AP-centred Hamming-six Z/8 certificate")
    print("scope primitive proper common-scale H6 owner-local gate only")
    print("status FROZEN-PRIMARY THM-1096 independent-referee-pending")
    print("orders", ORDERS, "unit counts", tuple(len(UNITS[d]) for d in ORDERS))
    print("numpy exact-batching version", np.__version__)
    print("hereditary grammar at least two order-32 providers")
    print(
        "hereditary words", len(grammar["words"]),
        "state words/support", grammar["weighted_states"],
        "labelled support/order contexts", 924 * len(grammar["words"]),
        "raw labelled states", 924 * grammar["weighted_states"],
    )
    print("order32-count grammar histogram", format_counter(grammar["order32_histogram"]))
    for order in ORDERS:
        print(f"D{order} ratio-cardinalities", tables["rows"][order])
    print("labelled literal masks audited", tables["labelled_masks"])
    print(
        "scalar feasible-owner histogram",
        " ".join(
            f"{index}:{value}"
            for index, value in enumerate(scalar["feasible_owner_histogram"])
        ),
    )
    print("scalar supports-by-survivor-count", format_counter(scalar["support_histogram"]))
    print(
        "scalar survivors", len(scalar["survivors"]),
        "supports", len({row[0] for row in scalar["survivors"]}),
        "literal unit words", scalar["literal_survivor_words"],
        "multiplicity profiles", len(scalar["multiplicities"]),
        "capacity vectors", len(scalar["capacity_vectors"]),
    )
    print("scalar multiplicity histogram", format_counter(scalar["multiplicities"]))
    print("scalar minimum-slack histogram", format_counter(scalar["minimum_slack"]))
    print("scalar maximum-slack histogram", format_counter(scalar["maximum_slack"]))
    print("scalar tight-owner histogram", format_counter(scalar["tight_owner_histogram"]))
    print(
        "normalized owner keys", len(fibre["unique_keys"]),
        "uniform labelled multiplicity", min(fibre["key_frequency"].values()),
    )
    print("mod8 anchor orders", tuple(sorted(MOD8_ANCHORS)), "fibres 8x4")
    print("mod8 anchor-bank histogram", format_counter(fibre["owner_bank"]))
    print("mod8 owner-bound histogram", format_counter(fibre["owner_bound"]))
    print("mod8 live-owner/context histogram", format_counter(fibre["live_contexts"]))
    print("dead-owner deficit histogram", format_counter(fibre["dead_margin"]))
    print("scalar-minus-relaxation histogram", format_counter(fibre["capacity_minus_bound"]))
    print("anchor unions audited", fibre["total_anchor_unions"])
    print("transverse option inequalities audited", fibre["total_option_checks"])
    print(
        "proof implication literal-cover => U8>=32; every scalar row has at least four owners with U8<32"
    )
    print("live-class multiplication orbits", fibre["orbit_histograms"])
    print("tournament vertices owner obligations; pair observable lex(live,U8,scalar-capacity,anchor-bank-size); harder key wins; coordinate-order tie Hamiltonian path")
    print("tournament fingerprints all 3450 transitive: scores 0,1,2,3,4,5; cycles 0; singleton SCCs; Hamiltonian paths 1")
    print("tournament tie-edge histogram", format_counter(fibre["tournament_ties"]))
    print("tournament flip-edge histogram", format_counter(fibre["tournament_flips"]))
    print("Kakeya/toothpick carrier has eight thick four-point Z/8 fibres; order-16/32 providers are independently relaxed transverse needles")
    print("carrier preserves the absolute owner-cover implication and threshold; destroys within-fibre positions, transverse overlaps, shared exact units, and literal witnesses only in the upper-bound direction")
    print("alternate-carrier audit runners/providers, gaps, fixed sections/boundaries, wall crossings, isolated residues, cover arcs, Fourier modes, matroid circuits, Fano points, chi7 colours, and completed tournaments all lose the absolute threshold or higher mask overlap")
    print("SHA256 CRT-bases", tables["base_digest"])
    print("SHA256 literal-masks", tables["mask_digest"])
    print("SHA256 grammar", grammar["digest"])
    print("SHA256 scalar", scalar["digest"])
    print("SHA256 normalized-relaxation", fibre["relaxation_digest"])
    print("SHA256 labelled-context-bounds", fibre["context_digest"])
    print("frontier verdict scalar-empty no; Z/8-relaxation all-six empty yes")
    print("scope caveat closes only c=32 primitive proper AP-centred common-sheet H6; finite H5, non-AP/deep sheets, uniform n=12 sporadic emptiness, and LRC14 remain open")


if __name__ == "__main__":
    main()
