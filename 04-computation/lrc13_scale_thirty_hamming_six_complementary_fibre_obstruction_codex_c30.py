#!/usr/bin/env python3
"""Exact primary for the c=30 AP-centred Hamming-six obstruction.

The proof-facing object is a *pair* of sound upper relaxations.  The first
retains every provider whose effective order divides 6 as an exact union of
the six thick fibres of Z/30 -> Z/6.  The second does the same for the ten
thick fibres of Z/30 -> Z/10.  In either relaxation, each transverse provider
is allowed to choose its exact-order residue independently so as to maximise
its contribution outside the retained anchor union.  This only enlarges the
attainable union, so a literal 30-sheet cover must pass both bounds.

The full hereditary divisor grammar and all 924 labelled supports are scanned.
The mod-six relaxation alone leaves 120 all-six-owner scalar rows.  Every
owner of those 120 rows fails the mod-ten relaxation, and globally at most two
owners of any row pass both.  Hence no common unit word can cover at all six
owners.

NumPy is used only for exact int16 batching of the scalar capacity gate.
Literal CRT search, masks, anchor banks, upper relaxations, hashes, orbits, and
tournament telemetry use deterministic Python integers and immutable tuples.
Multiplication orbits and completed tournaments are telemetry only; no row is
discarded through either quotient.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 30
ORDERS = (1, 2, 3, 5, 6, 10, 15, 30)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL = (1 << C) - 1
MOD6_ANCHORS = frozenset(order for order in ORDERS if 6 % order == 0)
MOD10_ANCHORS = frozenset(order for order in ORDERS if 10 % order == 0)


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


def rotate30(mask: int, amount: int) -> int:
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


def build_tables():
    """Build normalized options, then audit every labelled owner gauge."""
    normalized = {}
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    base_digest = sha256()
    mask_digest = sha256()

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
                    expected = rotate30(
                        literal_mask(ratio, order, unit, 1), shift
                    )
                    require(actual == expected, "cyclic owner-gauge mismatch")
                    require(
                        actual.bit_count() == cards[label, order_index, owner],
                        "literal mask/cardinality mismatch",
                    )
                    actual_options.append(actual)
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(actual.to_bytes(4, "little"))

                    for modulus, anchors in (
                        (6, MOD6_ANCHORS), (10, MOD10_ANCHORS)
                    ):
                        if order not in anchors:
                            continue
                        for sheet in range(C):
                            require(
                                ((actual >> sheet) & 1)
                                == ((actual >> ((sheet + modulus) % C)) & 1),
                                "anchor mask is not a union of thick fibres",
                            )

                require(
                    len(set(actual_options)) == len(normalized[ratio, order]),
                    "owner gauge changed the option-bank size",
                )

    rows = {
        order: tuple(
            analytic_cardinality(ratio, order) for ratio in range(1, P)
        )
        for order in ORDERS
    }
    expected_rows = {
        1: (30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        2: (15, 0, 0, 0, 0, 15, 15, 0, 0, 0, 0, 0),
        3: (10, 0, 0, 10, 10, 0, 0, 10, 10, 0, 0, 0),
        5: (6, 6, 6, 0, 6, 6, 6, 6, 0, 6, 6, 0),
        6: (5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0),
        10: (6, 6, 6, 3, 3, 6, 6, 3, 3, 6, 6, 3),
        15: (6, 4, 4, 4, 4, 6, 6, 4, 4, 4, 4, 4),
        30: (5, 4, 5, 5, 4, 5, 5, 4, 5, 5, 4, 4),
    }
    require(rows == expected_rows, "ratio-cardinality table changed")
    base_hex = base_digest.hexdigest()
    mask_hex = mask_digest.hexdigest()
    require(
        base_hex
        == "c94c62de38c0238f686e34cfc77293359e8638418485e59e49bed509f992dc3e",
        "CRT-base digest changed",
    )
    require(
        mask_hex
        == "dcff92077f5e320927a98c7f5a86df40e9c6baed2075e7cad46542d40d3e8da9",
        "literal-mask digest changed",
    )
    return {
        "normalized": normalized,
        "cards": cards,
        "rows": rows,
        "base_digest": base_hex,
        "mask_digest": mask_hex,
    }


def hereditary(index_word: tuple[int, ...]) -> bool:
    order_word = tuple(ORDERS[index] for index in index_word)
    return all(
        lcm(*(order_word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )


def prime_flag_grammar(index_word: tuple[int, ...]) -> bool:
    order_word = tuple(ORDERS[index] for index in index_word)
    return all(
        sum(order % prime == 0 for order in order_word) >= 2
        for prime in (2, 3, 5)
    )


def build_grammar():
    words = []
    weighted_states = 0
    digest = sha256()
    prime_count_histogram = Counter()
    for index_word in product(range(len(ORDERS)), repeat=6):
        via_lcm = hereditary(index_word)
        via_flags = prime_flag_grammar(index_word)
        require(via_lcm == via_flags, "lcm and prime-flag grammars disagree")
        if not via_lcm:
            continue
        words.append(index_word)
        order_word = tuple(ORDERS[index] for index in index_word)
        weight = prod(len(UNITS[order]) for order in order_word)
        weighted_states += weight
        prime_count_histogram[tuple(
            sum(order % prime == 0 for order in order_word)
            for prime in (2, 3, 5)
        )] += 1
        digest.update(bytes(index_word))
        digest.update(weight.to_bytes(8, "little"))

    require(len(words) == 185_193, "hereditary word count")
    require(weighted_states == 636_667_200, "weighted state-word count")
    digest_hex = digest.hexdigest()
    require(
        digest_hex
        == "d5ec0f6ac4c8f74ed399fbce0e8b50f232789f82619f556f7606b2e41fb6d55f",
        "hereditary grammar digest changed",
    )
    return {
        "words": np.asarray(words, dtype=np.int8),
        "weighted_states": weighted_states,
        "prime_count_histogram": prime_count_histogram,
        "digest": digest_hex,
    }


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
        == (
            1_401_966,
            36_143_640,
            66_874_158,
            49_478_260,
            15_326_622,
            1_839_636,
            54_050,
        ),
        "scalar feasible-owner histogram changed",
    )
    expected_support_histogram = Counter({
        0: 152,
        2: 112,
        4: 24,
        8: 24,
        12: 120,
        14: 6,
        17: 96,
        22: 96,
        30: 24,
        54: 96,
        68: 24,
        122: 24,
        137: 24,
        208: 24,
        253: 24,
        294: 12,
        296: 12,
        545: 24,
        549: 6,
    })
    require(
        support_histogram == expected_support_histogram,
        "support survivor histogram changed",
    )
    require(len(survivors) == 54_050, "scalar survivor count changed")
    require(
        len({support for support, _word, _capacity in survivors}) == 772,
        "scalar support count changed",
    )
    literal_survivor_words = sum(
        prod(len(UNITS[order]) for order in order_word)
        for _support, order_word, _capacity in survivors
    )
    require(
        literal_survivor_words == 64_678_912,
        "literal scalar-survivor unit-word count changed",
    )
    require(len(multiplicities) == 244, "multiplicity-profile count changed")
    require(len(capacity_vectors) == 28_965, "capacity-vector count changed")
    require(
        minimum_slack
        == Counter({0: 39_152, 1: 9_672, 2: 4_248, 3: 744, 4: 198, 5: 36}),
        "minimum-slack histogram changed",
    )
    require(
        tight_owner_histogram
        == Counter({0: 14_898, 1: 14_376, 2: 14_604, 3: 7_368,
                    4: 2_592, 5: 84, 6: 128}),
        "tight-owner histogram changed",
    )
    digest_hex = digest.hexdigest()
    require(
        digest_hex
        == "46ac0f0015d7a1731470cf1444e91f4c41df4a4ba9eac4ce62555aaaadfbe73d",
        "scalar survivor digest changed",
    )
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


def anchor_union_bank(key, anchor_orders, masks):
    bank = frozenset((0,))
    for ratio, order in key:
        if order not in anchor_orders:
            continue
        bank = frozenset(
            partial | option
            for partial in bank
            for option in masks[ratio, order]
        )
    return bank


def relaxation_bound(key, anchor_orders, anchor_bank, masks):
    """Upper-bound every literal union compatible with this owner key."""
    best = 0
    for anchor_union in anchor_bank:
        bound = anchor_union.bit_count()
        for ratio, order in key:
            if order in anchor_orders:
                continue
            bound += max(
                (option & ~anchor_union).bit_count()
                for option in masks[ratio, order]
            )
        best = max(best, bound)
    return best


def evaluate_quotient(unique_keys, anchor_orders, masks):
    bounds = {}
    bank_sizes = {}
    for key in unique_keys:
        bank = anchor_union_bank(key, anchor_orders, masks)
        bank_sizes[key] = len(bank)
        bounds[key] = relaxation_bound(key, anchor_orders, bank, masks)
    return bounds, bank_sizes


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
    require(scores == (0, 1, 2, 3, 4, 5), "tournament score sequence")
    require(triangles == 0, "owner tournament has a directed triangle")
    return ties, flips


def complementary_fibre_census(survivors, masks):
    key_rows = []
    context_rows = []
    for support, order_word, capacities in survivors:
        keys = tuple(
            normalized_owner_key(support, order_word, owner)
            for owner in support
        )
        key_rows.append(keys)
        context_rows.append(tuple(sorted(zip(support, order_word))))

    key_frequency = Counter(key for keys in key_rows for key in keys)
    unique_keys = tuple(sorted(key_frequency))
    require(len(unique_keys) == 27_025, "normalized owner-key count changed")
    require(
        set(key_frequency.values()) == {12},
        "normalized owner keys do not have uniform multiplicity twelve",
    )

    bound6, bank6 = evaluate_quotient(
        unique_keys, MOD6_ANCHORS, masks
    )
    bound10, bank10 = evaluate_quotient(
        unique_keys, MOD10_ANCHORS, masks
    )

    owner_bound6 = Counter()
    owner_bound10 = Counter()
    owner_common_bound = Counter()
    owner_bank6 = Counter()
    owner_bank10 = Counter()
    live6_contexts = Counter()
    live10_contexts = Counter()
    common_live_contexts = Counter()
    residual_order_profiles = Counter()
    residual_u10_scores = Counter()
    residual_contexts = []
    tournament_ties = Counter()
    tournament_flips = Counter()
    relaxation_digest = sha256()
    residual_digest = sha256()

    for (support, order_word, capacities), keys, context in zip(
        survivors, key_rows, context_rows
    ):
        bounds6 = tuple(bound6[key] for key in keys)
        bounds10 = tuple(bound10[key] for key in keys)
        common = tuple(min(left, right) for left, right in zip(bounds6, bounds10))

        live6 = sum(value >= C for value in bounds6)
        live10 = sum(value >= C for value in bounds10)
        live_common = sum(value >= C for value in common)
        live6_contexts[live6] += 1
        live10_contexts[live10] += 1
        common_live_contexts[live_common] += 1

        if live6 == 6:
            residual_contexts.append(context)
            residual_order_profiles[tuple(
                order_word.count(order) for order in ORDERS
            )] += 1
            residual_u10_scores.update(bounds10)
            require(
                max(bounds10) < C,
                "a mod-six all-owner residual has a mod-ten-live owner",
            )
            residual_digest.update(bytes(support))
            residual_digest.update(bytes(order_word))
            residual_digest.update(bytes(bounds6))
            residual_digest.update(bytes(bounds10))

        owner_rows = tuple(
            (
                int(common[index] >= C),
                common[index],
                bounds6[index],
                bounds10[index],
                capacities[index],
            )
            for index in range(6)
        )
        ties, flips = tournament_fingerprint(owner_rows)
        tournament_ties[ties] += 1
        tournament_flips[flips] += 1

        for key, left, right, both in zip(keys, bounds6, bounds10, common):
            owner_bound6[left] += 1
            owner_bound10[right] += 1
            owner_common_bound[both] += 1
            owner_bank6[bank6[key]] += 1
            owner_bank10[bank10[key]] += 1

        relaxation_digest.update(bytes(support))
        relaxation_digest.update(bytes(order_word))
        relaxation_digest.update(bytes(bounds6))
        relaxation_digest.update(bytes(bounds10))
        relaxation_digest.update(bytes(common))

    require(
        live6_contexts == Counter({0: 45_110, 1: 7_536, 2: 1_284, 6: 120}),
        "mod-six live-owner histogram changed",
    )
    require(
        live10_contexts
        == Counter({
            0: 33_944,
            1: 10_344,
            2: 6_288,
            3: 2_232,
            4: 744,
            5: 132,
            6: 366,
        }),
        "mod-ten live-owner histogram changed",
    )
    require(
        common_live_contexts == Counter({0: 45_998, 1: 6_852, 2: 1_200}),
        "complementary live-owner histogram changed",
    )
    require(len(residual_contexts) == 120, "mod-six residual count changed")
    require(max(common_live_contexts) == 2, "an all-six common row survived")
    require(
        residual_order_profiles
        == Counter({
            (0, 0, 0, 1, 0, 3, 0, 2): 48,
            (0, 0, 0, 0, 0, 4, 0, 2): 48,
            (0, 0, 0, 0, 0, 2, 3, 1): 24,
        }),
        "mod-six residual order profiles changed",
    )
    residual_orbits = orbit_histogram(residual_contexts)
    all_orbits = orbit_histogram(context_rows)
    require(
        residual_u10_scores
        == Counter({23: 48, 24: 120, 25: 192, 26: 240, 27: 72, 28: 48}),
        "mod-six residual mod-ten scores changed",
    )
    require(residual_orbits == Counter({12: 10}), "residual orbit census changed")
    require(
        all_orbits == Counter({4: 2, 6: 9, 12: 4_499}),
        "scalar context orbit census changed",
    )
    require(
        owner_bank6
        == Counter({1: 114_588, 2: 123_024, 3: 42_996, 4: 35_760,
                    5: 5_244, 7: 2_688}),
        "mod-six anchor-bank histogram changed",
    )
    require(
        owner_bank10
        == Counter({
            1: 81_792, 4: 142_704, 10: 42_660, 14: 2_856,
            16: 25_752, 28: 6_072, 32: 384, 34: 4_848,
            38: 1_680, 40: 2_064, 50: 4_320, 52: 2_988,
            58: 2_880, 66: 144, 70: 528, 75: 108, 86: 984,
            94: 1_224, 98: 240, 100: 72,
        }),
        "mod-ten anchor-bank histogram changed",
    )
    require(
        owner_bound6
        == Counter({
            16: 144, 17: 144, 18: 384, 19: 840, 20: 7_632,
            21: 9_264, 22: 26_496, 23: 60_972, 24: 104_388,
            25: 31_368, 26: 42_480, 27: 22_464, 28: 6_600,
            29: 300, 30: 10_608, 31: 72, 32: 96, 34: 48,
        }),
        "mod-six owner-bound histogram changed",
    )
    require(
        owner_bound10
        == Counter({
            18: 288, 19: 1_656, 20: 6_720, 21: 6_888, 22: 17_976,
            23: 29_664, 24: 67_800, 25: 39_408, 26: 70_092,
            27: 20_400, 28: 21_324, 29: 6_636, 30: 17_328,
            31: 3_756, 32: 3_204, 33: 3_720, 34: 4_344,
            35: 1_812, 36: 624, 37: 168, 38: 336, 39: 36,
            42: 48, 43: 72,
        }),
        "mod-ten owner-bound histogram changed",
    )
    require(
        owner_common_bound
        == Counter({
            16: 144, 17: 144, 18: 384, 19: 1_944, 20: 8_928,
            21: 10_584, 22: 28_680, 23: 64_668, 24: 104_052,
            25: 30_000, 26: 38_640, 27: 21_768, 28: 5_112,
            30: 9_252,
        }),
        "common minimum-bound histogram changed",
    )
    require(
        tournament_ties
        == Counter({0: 17_184, 1: 23_232, 2: 8_088, 3: 2_886,
                    4: 1_608, 6: 660, 7: 264, 15: 128}),
        "tournament tie histogram changed",
    )
    require(
        tournament_flips
        == Counter({
            0: 283, 1: 425, 2: 1_102, 3: 2_349, 4: 4_208,
            5: 6_369, 6: 8_357, 7: 9_039, 8: 8_270, 9: 6_114,
            10: 3_827, 11: 2_149, 12: 997, 13: 430, 14: 120, 15: 11,
        }),
        "tournament flip histogram changed",
    )
    relaxation_hex = relaxation_digest.hexdigest()
    residual_hex = residual_digest.hexdigest()
    require(
        relaxation_hex
        == "22ec2ba4da34d3e669bb7564975b640f3cb4263fa0e75e55aa6353737e98da75",
        "complementary-relaxation digest changed",
    )
    require(
        residual_hex
        == "6adaeac9e329b35d1d5328c3490ca1c567d9e629ce2f574a2c91784e1cefda12",
        "mod-six residual digest changed",
    )

    return {
        "unique_keys": unique_keys,
        "key_frequency": key_frequency,
        "owner_bound6": owner_bound6,
        "owner_bound10": owner_bound10,
        "owner_common_bound": owner_common_bound,
        "owner_bank6": owner_bank6,
        "owner_bank10": owner_bank10,
        "live6_contexts": live6_contexts,
        "live10_contexts": live10_contexts,
        "common_live_contexts": common_live_contexts,
        "residual_order_profiles": residual_order_profiles,
        "residual_u10_scores": residual_u10_scores,
        "residual_orbits": residual_orbits,
        "all_orbits": all_orbits,
        "tournament_ties": tournament_ties,
        "tournament_flips": tournament_flips,
        "relaxation_digest": relaxation_hex,
        "residual_digest": residual_hex,
    }


def main() -> None:
    tables = build_tables()
    grammar = build_grammar()
    scalar = scalar_census(grammar["words"], tables["cards"])
    fibre = complementary_fibre_census(
        scalar["survivors"], tables["normalized"]
    )

    print("scale-thirty AP-centred Hamming-six complementary-fibre certificate")
    print("scope primitive proper common-scale H6 owner-local gate only")
    print("artifact FROZEN-PRIMARY THM-1090 theorem-status-in-canonical-page")
    print("orders", ORDERS, "unit counts", tuple(len(UNITS[d]) for d in ORDERS))
    print("numpy exact-batching version", np.__version__)
    print("hereditary grammar at least two providers for each prime flag 2,3,5")
    print(
        "hereditary words", len(grammar["words"]),
        "state words/support", grammar["weighted_states"],
        "labelled support/order contexts", 924 * len(grammar["words"]),
        "raw labelled states", 924 * grammar["weighted_states"],
    )
    print("prime-count triples", format_counter(grammar["prime_count_histogram"]))
    for order in ORDERS:
        print(f"D{order} ratio-cardinalities", tables["rows"][order])
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
    print("scalar minimum-slack histogram", format_counter(scalar["minimum_slack"]))
    print("scalar maximum-slack histogram", format_counter(scalar["maximum_slack"]))
    print("scalar tight-owner histogram", format_counter(scalar["tight_owner_histogram"]))
    print(
        "normalized owner keys", len(fibre["unique_keys"]),
        "uniform labelled multiplicity", min(fibre["key_frequency"].values()),
    )
    print("mod6 anchor orders", tuple(sorted(MOD6_ANCHORS)), "fibres 6x5")
    print("mod10 anchor orders", tuple(sorted(MOD10_ANCHORS)), "fibres 10x3")
    print("mod6 anchor-bank histogram", format_counter(fibre["owner_bank6"]))
    print("mod10 anchor-bank histogram", format_counter(fibre["owner_bank10"]))
    print("mod6 owner-bound histogram", format_counter(fibre["owner_bound6"]))
    print("mod10 owner-bound histogram", format_counter(fibre["owner_bound10"]))
    print("common min-bound histogram", format_counter(fibre["owner_common_bound"]))
    print("mod6 live-owner/context histogram", format_counter(fibre["live6_contexts"]))
    print("mod10 live-owner/context histogram", format_counter(fibre["live10_contexts"]))
    print("common live-owner/context histogram", format_counter(fibre["common_live_contexts"]))
    print(
        "mod6 all-six residual rows", sum(fibre["residual_order_profiles"].values()),
        "support/order multiplication orbits", format_counter(fibre["residual_orbits"]),
    )
    print("mod6 residual order profiles", format_counter(fibre["residual_order_profiles"]))
    print("mod6 residual mod10-score histogram", format_counter(fibre["residual_u10_scores"]))
    print("all scalar context multiplication orbits", format_counter(fibre["all_orbits"]))
    print(
        "proof implication literal-cover => U6>=30 and U10>=30; "
        "at most two owners satisfy both in every labelled scalar row"
    )
    print("tournament vertices owner obligations; pair observable lex(common-status,min(U6,U10),U6,U10,scalar-capacity); harder key wins; coordinate-order tie Hamiltonian path")
    print("tournament fingerprints all 54050 transitive: scores 0,1,2,3,4,5; cycles 0; singleton SCCs; Hamiltonian paths 1")
    print("tournament tie-edge histogram", format_counter(fibre["tournament_ties"]))
    print("tournament flip-edge histogram", format_counter(fibre["tournament_flips"]))
    print("Kakeya/toothpick carrier mod6 has six thick five-point fibres and mod10 has ten thick three-point fibres; nonanchors are independently relaxed transverse needles")
    print("carrier preserves absolute owner cover implication and the terminal threshold; destroys within-fibre positions, transverse overlaps, shared exact units, and literal witnesses, all in the upper-bound direction")
    print("paired-quotient insight mod6 forgets the 5-coordinate and mod10 forgets the 3-coordinate; their viability predicates, not a completed tournament, separate the scalar survivors")
    print("alternate-carrier audit providers/runners, gaps, fixed sections/boundaries, wall events, isolated residues, cover arcs, Fourier modes, matroid circuits, Fano points, chi7 colours, and completed tournaments each lose the absolute threshold or paired thick-fibre incidence")
    print("SHA256 CRT-bases", tables["base_digest"])
    print("SHA256 literal-masks", tables["mask_digest"])
    print("SHA256 grammar", grammar["digest"])
    print("SHA256 scalar", scalar["digest"])
    print("SHA256 complementary-relaxation", fibre["relaxation_digest"])
    print("SHA256 mod6-residual", fibre["residual_digest"])
    print("frontier verdict scalar-empty no; mod6-alone empty no; mod10-alone empty no; complementary mod6+mod10 all-six empty yes")
    print("scope caveat closes only c=30 primitive proper AP-centred common-sheet H6; finite H5, non-AP/deep sheets, uniform n=12 sporadic emptiness, and LRC14 remain open")


if __name__ == "__main__":
    main()
