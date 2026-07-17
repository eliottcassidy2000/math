#!/usr/bin/env python3
"""Independent exact referee for the common-scale-nine Hamming-six bank.

This implementation deliberately uses a different representation from the
C++ frontier certificate.  It builds masks from CRT congruences, searches the
473 hereditary order words as tuples of actual divisors, uses Python sets for
owner-local mask reachability, and packs the six owners only in the final
3,048,192-word replay.  Integer arithmetic is exact throughout.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd


P = 13
C = 9
LABELS = tuple(range(1, P))
FULL = (1 << C) - 1
FULL_PACKED = (1 << (6 * C)) - 1
STATES = (
    (1, 0),
    (3, 1),
    (3, 2),
    (9, 1),
    (9, 2),
    (9, 4),
    (9, 5),
    (9, 7),
    (9, 8),
)
STATE_INDICES = {
    divisor: tuple(i for i, state in enumerate(STATES) if state[0] == divisor)
    for divisor in (1, 3, 9)
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, divisor: int, unit: int) -> int:
    return next(
        value
        for value in range(P * divisor)
        if value % P == divisor * label % P and value % divisor == unit
    )


def sheet_mask(label: int, state: int, owner: int) -> int:
    divisor, unit = STATES[state]
    base = crt_base(label, divisor, unit)
    owner_inverse = pow(owner, -1, P)
    return sum(
        1 << sheet
        for sheet in range(C)
        if -divisor
        < centered(base * (owner_inverse + P * sheet), P * divisor)
        <= divisor
    )


MASK = {
    (label, state, owner): sheet_mask(label, state, owner)
    for label in LABELS
    for state in range(len(STATES))
    for owner in LABELS
}
MASK_DIGEST = sha256(
    b"".join(
        MASK[label, state, owner].to_bytes(2, "little")
        for label in LABELS
        for state in range(len(STATES))
        for owner in LABELS
    )
).hexdigest()


# Unit-independent mask cardinality is the exact scalar relaxation.
CARDINALITY = {}
for divisor in (1, 3, 9):
    for ratio in LABELS:
        values = {
            MASK[ratio, state, 1].bit_count()
            for state in STATE_INDICES[divisor]
        }
        require(len(values) == 1, "mask cardinality depends on a unit")
        CARDINALITY[divisor, ratio] = values.pop()


ORDER_WORDS = tuple(
    word
    for word in product((1, 3, 9), repeat=6)
    if word.count(9) >= 2
)
require(len(ORDER_WORDS) == 473, "hereditary order grammar mismatch")


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


def scalar_capacity(labels: tuple[int, ...], orders: tuple[int, ...]) -> bool:
    return all(
        sum(CARDINALITY[divisor, ratio(provider, owner)]
            for provider, divisor in zip(labels, orders))
        >= C
        for owner in labels
    )


def owner_locally_feasible(
    labels: tuple[int, ...], orders: tuple[int, ...], owner: int
) -> bool:
    reachable = {0}
    for provider, divisor in zip(labels, orders):
        provider_masks = {
            MASK[provider, state, owner]
            for state in STATE_INDICES[divisor]
        }
        reachable = {
            partial | provider_mask
            for partial in reachable
            for provider_mask in provider_masks
        }
    return FULL in reachable


def signed_doubling_cycle(labels: tuple[int, ...]) -> bool:
    edges = {
        (provider, owner)
        for provider in labels
        for owner in labels
        if provider != owner and owner * pow(provider, -1, P) % P in (2, 11)
    }
    if len(edges) != 6:
        return False
    successor = {}
    for vertex in labels:
        outgoing = tuple(owner for provider, owner in edges if provider == vertex)
        incoming = tuple(provider for provider, owner in edges if owner == vertex)
        if len(outgoing) != 1 or len(incoming) != 1:
            return False
        successor[vertex] = outgoing[0]
    orbit = [labels[0]]
    for _ in range(5):
        orbit.append(successor[orbit[-1]])
    return len(set(orbit)) == 6 and successor[orbit[-1]] == orbit[0]


def cycle_order(labels: tuple[int, ...]) -> tuple[int, ...]:
    require(signed_doubling_cycle(labels), "cycle order requested off the cycle bank")
    successor = {
        provider: owner
        for provider in labels
        for owner in labels
        if provider != owner and owner * pow(provider, -1, P) % P in (2, 11)
    }
    answer = [labels[0]]
    for _ in range(5):
        answer.append(successor[answer[-1]])
    return tuple(answer)


def mixed_orbit(labels: tuple[int, ...], orders: tuple[int, ...]) -> bool:
    d3 = {label for label, divisor in zip(labels, orders) if divisor == 3}
    d9 = {label for label, divisor in zip(labels, orders) if divisor == 9}
    return any(
        d3 == {anchor, 5 * anchor % P}
        and d9 == {2 * anchor % P, 11 * anchor % P,
                   3 * anchor % P, 10 * anchor % P}
        for anchor in d3
    )


supports = tuple(combinations(LABELS, 6))
scalar_bank = []
local_bank = []
scalar_pattern = Counter()
for labels in supports:
    for orders in ORDER_WORDS:
        if not scalar_capacity(labels, orders):
            continue
        scalar_bank.append((labels, orders))
        scalar_pattern[tuple(orders.count(divisor) for divisor in (1, 3, 9))] += 1
        if all(owner_locally_feasible(labels, orders, owner) for owner in labels):
            local_bank.append((labels, orders))

require(len(scalar_bank) == 1_186, "scalar-capacity bank mismatch")
require(len({labels for labels, _ in scalar_bank}) == 316,
        "scalar-capacity support count mismatch")
require(
    scalar_pattern
    == {
        (0, 0, 6): 82,
        (0, 2, 4): 474,
        (0, 3, 3): 132,
        (0, 4, 2): 330,
        (1, 3, 2): 168,
    },
    "scalar order-pattern histogram mismatch",
)
require(len(local_bank) == 76, "owner-local bank mismatch")
require(len({labels for labels, _ in local_bank}) == 76,
        "owner-local supports are not unique")


all_d9 = tuple(context for context in local_bank if context[1] == (9,) * 6)
mixed = tuple(context for context in local_bank if context[1] != (9,) * 6)
signed_supports = {labels for labels in supports if signed_doubling_cycle(labels)}
require(len(all_d9) == 64 and {labels for labels, _ in all_d9} == signed_supports,
        "all-D9 bank is not the signed-doubling bank")
require(len(mixed) == 12 and all(mixed_orbit(*context) for context in mixed),
        "mixed bank is not the claimed multiplicative orbit")

LOCAL_BANK_DIGEST = sha256(
    "\n".join(
        f"{','.join(map(str, labels))}|{','.join(map(str, orders))}"
        for labels, orders in local_bank
    ).encode()
).hexdigest()


def packed_provider_masks(labels: tuple[int, ...]):
    return tuple(
        tuple(
            sum(MASK[provider, state, owner] << (C * owner_index)
                for owner_index, owner in enumerate(labels))
            for state in range(len(STATES))
        )
        for provider in labels
    )


def satisfaction_subset(packed_cover: int) -> int:
    return sum(
        (1 << owner) if (packed_cover >> (C * owner)) & FULL == FULL else 0
        for owner in range(6)
    )


def fibre_profile(labels: tuple[int, ...], orders: tuple[int, ...]) -> Counter[int]:
    packed = packed_provider_masks(labels)
    profile = Counter()
    for states in product(*(STATE_INDICES[divisor] for divisor in orders)):
        cover = 0
        for provider in range(6):
            cover |= packed[provider][states[provider]]
        profile[satisfaction_subset(cover)] += 1
    return profile


global_words_tested = 0
global_survivors = 0
for labels, orders in local_bank:
    profile = fibre_profile(labels, orders)
    global_words_tested += sum(profile.values())
    global_survivors += profile[63]
    by_count = Counter()
    for subset, multiplicity in profile.items():
        by_count[subset.bit_count()] += multiplicity
    owner_sizes = tuple(
        sum(multiplicity for subset, multiplicity in profile.items()
            if subset & (1 << owner))
        for owner in range(6)
    )
    pair_sizes = {
        (i, j): sum(multiplicity for subset, multiplicity in profile.items()
                    if subset & (1 << i) and subset & (1 << j))
        for i, j in combinations(range(6), 2)
    }

    if orders == (9,) * 6:
        require(by_count == {0: 44_100, 1: 2_520, 2: 36},
                "all-D9 satisfaction profile mismatch")
        require(owner_sizes == (432,) * 6, "all-D9 owner size mismatch")
        cycle = cycle_order(labels)
        position = {label: i for i, label in enumerate(cycle)}
        for (i, j), intersection in pair_sizes.items():
            distance = abs(position[labels[i]] - position[labels[j]])
            distance = min(distance, 6 - distance)
            require(intersection == (12 if distance == 3 else 0),
                    "all-D9 pair nerve mismatch")
    else:
        require(by_count == {0: 2_928, 1: 1_776, 2: 336, 3: 144},
                "mixed satisfaction profile mismatch")
        require(owner_sizes == tuple(1_152 if divisor == 3 else 144
                                     for divisor in orders),
                "mixed owner size mismatch")
        for (i, j), intersection in pair_sizes.items():
            if orders[i] == orders[j] == 3:
                require(intersection == 192, "mixed D3-D3 intersection mismatch")
            elif orders[i] == 3 or orders[j] == 3:
                require(intersection in (24, 48), "mixed cross intersection mismatch")
            else:
                require(intersection == (144 if labels[i] + labels[j] == P else 0),
                        "mixed D9 sub-nerve mismatch")


require(global_words_tested == 3_048_192, "reduced fibre census mismatch")
require(global_survivors == 0, "a common-sheet state word survived")


print("scale-nine independent Python referee")
print(f"mask digest {MASK_DIGEST}")
print("hereditary order words 473; labelled contexts 437052")
print("hereditary state words 521964; labelled contexts 482294736")
print("scalar capacity 1186 contexts on 316 supports")
print("owner-local bank 76 contexts on 76 supports")
print("owner-local split all-D9 signed-C6:64 mixed-orbit:12")
print(f"owner-local bank digest {LOCAL_BANK_DIGEST}")
print("reduced literal state words 3048192")
print("global common-sheet survivors 0")
print("all-D9 nerve 3K2; owner 432; antipodal pair 12; profile 44100,2520,36")
print("mixed D9 sub-nerve 2K2; D9 owner 144; antipodal pair 144")
print("mixed profile 2928,1776,336,144")
