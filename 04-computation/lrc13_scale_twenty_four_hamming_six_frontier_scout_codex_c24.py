#!/usr/bin/env python3
"""Exact scale-24 AP-centred Hamming-six owner-deficit certificate.

All arithmetic is integral.  NumPy is used only to batch the 100,543,212
support/order capacity rows; literal CRT masks and owner-local reachable unions
are independently reconstructed with Python integers and immutable sets.  No
support, unit word, or owner is quotiented out.  Multiplication orbits and the
completed owner tournament are telemetry only.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 24
ORDERS = (1, 2, 3, 4, 6, 8, 12, 24)
UNITS = {
    order: ((0,) if order == 1 else tuple(u for u in range(1, order + 1) if gcd(u, order) == 1))
    for order in ORDERS
}
FULL_MASK = (1 << C) - 1
MASK64 = (1 << 64) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base_algebraic(label: int, order: int, unit: int) -> int:
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def crt_base_literal(label: int, order: int, unit: int) -> int:
    for value in range(P * order):
        if value % P == order * label % P and value % order == unit % order:
            return value
    raise AssertionError("literal CRT search failed")


def local_mask(label: int, order: int, unit: int, owner: int) -> int:
    base = crt_base_algebraic(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    result = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            result |= 1 << sheet
    return result


def analytic_cardinality(label: int, order: int, owner: int) -> int:
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    period = sum(value % P == target for value in range(-order + 1, order + 1))
    return (C // order) * period


def build_tables():
    masks = {}
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    digest = sha256()
    for label in range(1, P):
        for order_index, order in enumerate(ORDERS):
            for owner in range(1, P):
                card = analytic_cardinality(label, order, owner)
                cards[label, order_index, owner] = card
                for unit in UNITS[order]:
                    require(
                        crt_base_algebraic(label, order, unit)
                        == crt_base_literal(label, order, unit),
                        "algebraic/literal CRT mismatch",
                    )
                    mask = local_mask(label, order, unit, owner)
                    require(mask.bit_count() == card, "mask/cardinality mismatch")
                    masks[label, order, unit, owner] = mask
                    digest.update(bytes((label, order, unit, owner)))
                    digest.update(mask.to_bytes(4, "little"))
    return masks, cards, digest.hexdigest()


def hereditary(index_word: tuple[int, ...]) -> bool:
    word = tuple(ORDERS[index] for index in index_word)
    return all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )


def hereditary_prime_power(index_word: tuple[int, ...]) -> bool:
    word = tuple(ORDERS[index] for index in index_word)
    return (
        sum(order % 8 == 0 for order in word) >= 2
        and sum(order % 3 == 0 for order in word) >= 2
    )


def build_order_words():
    words = []
    digest = sha256()
    for index_word in product(range(len(ORDERS)), repeat=6):
        via_lcm = hereditary(index_word)
        require(via_lcm == hereditary_prime_power(index_word), "hereditary audit")
        if via_lcm:
            words.append(index_word)
            digest.update(bytes(index_word))
    return np.asarray(words, dtype=np.int8), digest.hexdigest()


def scalar_bank(words: np.ndarray, cards: np.ndarray):
    bank = []
    support_histogram = Counter()
    feasible_owner_histogram = np.zeros(7, dtype=np.int64)
    multiplicities = Counter()
    capacities_seen = set()
    minimum_slack = Counter()
    maximum_slack = Counter()
    tight_owners = Counter()
    scalar_digest = sha256()

    for support_tuple in combinations(range(1, P), 6):
        support = np.asarray(support_tuple, dtype=np.int8)
        capacities = np.zeros((len(words), 6), dtype=np.int16)
        for provider, label in enumerate(support):
            capacities += cards[label, words[:, provider, None], support[None, :]]
        feasible_count = (capacities >= C).sum(axis=1)
        feasible_owner_histogram += np.bincount(feasible_count, minlength=7)
        survivor_indices = np.flatnonzero(feasible_count == 6)
        support_histogram[len(survivor_indices)] += 1
        for row_index in survivor_indices:
            index_word = tuple(int(x) for x in words[row_index])
            order_word = tuple(ORDERS[index] for index in index_word)
            capacity_word = tuple(int(x) for x in capacities[row_index])
            bank.append((support_tuple, order_word, capacity_word))
            multiplicities[tuple(order_word.count(order) for order in ORDERS)] += 1
            capacities_seen.add(capacity_word)
            minimum_slack[min(capacity_word) - C] += 1
            maximum_slack[max(capacity_word) - C] += 1
            tight_owners[sum(value == C for value in capacity_word)] += 1
            scalar_digest.update(bytes(support_tuple))
            scalar_digest.update(bytes(order_word))
            for value in capacity_word:
                scalar_digest.update(value.to_bytes(2, "little"))

    return {
        "bank": bank,
        "support_histogram": support_histogram,
        "feasible_owner_histogram": feasible_owner_histogram,
        "multiplicities": multiplicities,
        "capacities_seen": capacities_seen,
        "minimum_slack": minimum_slack,
        "maximum_slack": maximum_slack,
        "tight_owners": tight_owners,
        "digest": scalar_digest.hexdigest(),
    }


def mix64(value: int) -> int:
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9 & MASK64
    value = (value ^ (value >> 27)) * 0x94D049BB133111EB & MASK64
    return value ^ (value >> 31)


def owner_local(support, word, owner, masks):
    reachable = frozenset({0})
    for label, order in zip(support, word):
        choices = frozenset(masks[label, order, unit, owner] for unit in UNITS[order])
        reachable = frozenset(partial | choice for partial in reachable for choice in choices)
    maximum = 0
    mask_sum = 0
    mask_xor = 0
    for mask in reachable:
        maximum = max(maximum, mask.bit_count())
        mask_sum = (mask_sum + mask) & MASK64
        mask_xor ^= mix64(mask)
    feasible = FULL_MASK in reachable
    require(feasible == (maximum == C), "owner feasibility/maximum mismatch")
    return feasible, maximum, len(reachable), mask_sum, mask_xor


def tournament_summary(capacities, local_rows):
    adjacency = [[False] * 6 for _ in range(6)]
    ties = 0
    flips = 0
    for left in range(6):
        for right in range(left + 1, 6):
            left_key = (local_rows[left][0], local_rows[left][1], capacities[left])
            right_key = (local_rows[right][0], local_rows[right][1], capacities[right])
            if left_key == right_key:
                ties += 1
                winner, loser = left, right
            elif left_key > right_key:
                winner, loser = left, right
            else:
                winner, loser = right, left
                flips += 1
            adjacency[winner][loser] = True
    scores = tuple(sorted(sum(row) for row in adjacency))
    triangles = 0
    for a, b, c in combinations(range(6), 3):
        triangles += int(
            (adjacency[a][b] and adjacency[b][c] and adjacency[c][a])
            or (adjacency[a][c] and adjacency[c][b] and adjacency[b][a])
        )
    require(scores == (0, 1, 2, 3, 4, 5) and triangles == 0, "tournament audit")
    return ties, flips


def multiplication_orbits(supports):
    remaining = set(supports)
    histogram = Counter()
    while remaining:
        support = min(remaining)
        orbit = {
            tuple(sorted(multiplier * label % P for label in support))
            for multiplier in range(1, P)
        }
        require(orbit <= supports, "multiplication-invariance")
        remaining -= orbit
        histogram[len(orbit)] += 1
    return histogram


def fmt(counter: Counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main() -> None:
    masks, cards, mask_digest = build_tables()
    words, order_digest = build_order_words()
    state_words_per_support = sum(
        prod(len(UNITS[ORDERS[index]]) for index in word) for word in words
    )
    scalar = scalar_bank(words, cards)
    bank = scalar["bank"]
    support_set = {support for support, _, _ in bank}

    feasible_contexts = Counter()
    maximum_union = Counter()
    minimum_owner_maximum = Counter()
    owner_profiles = set()
    reachable_counts = Counter()
    reachable_total = 0
    reachable_maximum = 0
    feasible_rows = 0
    owner_digest = sha256()
    tie_histogram = Counter()
    flip_histogram = Counter()

    for support, word, capacities in bank:
        local_rows = []
        for owner in support:
            row = owner_local(support, word, owner, masks)
            feasible, maximum, count, mask_sum, mask_xor = row
            local_rows.append(row)
            feasible_rows += feasible
            maximum_union[maximum] += 1
            reachable_counts[count] += 1
            reachable_total += count
            reachable_maximum = max(reachable_maximum, count)
            owner_digest.update(bytes(support))
            owner_digest.update(bytes(word))
            owner_digest.update(bytes((owner, int(feasible), maximum)))
            owner_digest.update(count.to_bytes(4, "little"))
            owner_digest.update(mask_sum.to_bytes(8, "little"))
            owner_digest.update(mask_xor.to_bytes(8, "little"))
        feasible_count = sum(row[0] for row in local_rows)
        owner_vector = tuple(row[1] for row in local_rows)
        feasible_contexts[feasible_count] += 1
        minimum_owner_maximum[min(owner_vector)] += 1
        owner_profiles.add(owner_vector)
        ties, flips = tournament_summary(capacities, local_rows)
        tie_histogram[ties] += 1
        flip_histogram[flips] += 1

    require(len(words) == 108_813, "hereditary word count")
    require(state_words_per_support == 167_165_952, "literal state count")
    require(len(bank) == 66_984 and len(support_set) == 854, "scalar bank")
    require(len(scalar["multiplicities"]) == 202, "multiplicity profiles")
    require(
        feasible_contexts == Counter({0: 64_962, 1: 1_800, 2: 192, 4: 30}),
        "owner deficit",
    )
    require(
        maximum_union
        == Counter(
            {12: 72, 14: 2136, 15: 1644, 16: 15876, 17: 24420, 18: 76296,
             19: 94872, 20: 104592, 21: 53040, 22: 24948, 23: 1704, 24: 2304}
        ),
        "maximum-union histogram",
    )
    require(len(owner_profiles) == 20_302, "owner profiles")
    require(reachable_maximum == 7_728, "reachable maximum")

    orbit_histogram = multiplication_orbits(support_set)
    nonzero_support_histogram = Counter(scalar["support_histogram"])
    zero_supports = nonzero_support_histogram.pop(0)

    print("scale-twenty-four AP-centred Hamming-six owner-deficit certificate")
    print("divisor grammar 1,2,3,4,6,8,12,24; literal states 24")
    print(
        f"supports 924; hereditary order words {len(words)}; labelled order contexts "
        f"{924 * len(words)}"
    )
    print(
        f"state words/support {state_words_per_support}; raw labelled states "
        f"{924 * state_words_per_support}"
    )
    print(f"numpy exact-batching version {np.__version__}")
    print(f"mask SHA256 {mask_digest}; order SHA256 {order_digest}")
    print(
        f"scalar contexts {len(bank)} on {len(support_set)} supports; "
        f"multiplicity patterns {len(scalar['multiplicities'])}; "
        f"scalar-bank SHA256 {scalar['digest']}"
    )
    print("scalar capacity feasible-owner histogram", " ".join(
        f"{i}:{int(value)}" for i, value in enumerate(scalar["feasible_owner_histogram"])
    ))
    print("scalar-empty supports", zero_supports)
    print("contexts-per-support histogram", fmt(nonzero_support_histogram))
    print("multiplication orbit-size histogram", fmt(orbit_histogram), "(telemetry; no quotient)")
    print(f"capacity vectors {len(scalar['capacities_seen'])}")
    print("minimum scalar-slack histogram", fmt(scalar["minimum_slack"]))
    print("maximum scalar-slack histogram", fmt(scalar["maximum_slack"]))
    print("tight-owner/context histogram", fmt(scalar["tight_owners"]))
    print(f"owner-local rows {6 * len(bank)}; feasible rows {feasible_rows}")
    print("feasible-owner/context histogram", fmt(feasible_contexts))
    print("maximum reachable sheet-union histogram", fmt(maximum_union))
    print("minimum owner maximum/context histogram", fmt(minimum_owner_maximum))
    print(f"distinct owner max-union vectors {len(owner_profiles)}")
    print(
        f"reachable banks {6 * len(bank)}; total masks {reachable_total}; "
        f"distinct bank-size bins {len(reachable_counts)}; maximum bank {reachable_maximum}"
    )
    print(f"owner-summary SHA256 {owner_digest.hexdigest()}")
    print("owner-local all-six contexts", feasible_contexts.get(6, 0))
    print("tournament pair observable exact ordered (feasible,max-union,capacity) owner summaries; lexicographic switch and coordinate tie Hamiltonian path")
    print("tournament fingerprints all 66984 transitive: scores 0,1,2,3,4,5; cycles 0; SCCs 6; Hamiltonian paths 1")
    print("tournament tie-edge histogram", fmt(tie_histogram))
    print("tournament edge-flip histogram", fmt(flip_histogram))
    print("challenged vertices owner obligations preserve the terminal deficit through their feasibility/max-union vector and the scalar gate through capacities; the tournament loses thresholds and magnitudes, while provider, divisor, residue, isolated sheet, and wall-event vertices lose shared-unit incidence")
    print("frontier verdict scalar-empty no; owner-local all-six empty yes; next legal common scale 25")
    print("local D24 mask table at owner one (units 1,5,7,11,13,17,19,23; ratios 1..12 in hex)")
    for unit in UNITS[24]:
        print(f"  e={unit}: " + " ".join(
            f"{masks[label, 24, unit, 1]:x}" for label in range(1, P)
        ))


if __name__ == "__main__":
    main()
