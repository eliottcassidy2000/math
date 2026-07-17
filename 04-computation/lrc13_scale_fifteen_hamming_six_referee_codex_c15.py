#!/usr/bin/env python3
"""Independent exact referee for common-scale-fifteen Hamming-six packets.

The full owner-local reachability sets are reconstructed with exact Python
integers.  The faithful terminal datum is the six-owner feasibility subset (or
the stronger vector of maximum reachable union sizes), not a runner graph.
For Tournament Analysis we orient owners by decreasing maximum union size and
break ties by label order.  This always gives a transitive completion and is
reported only to demonstrate what the ranking forgets: the absolute threshold
fifteen and the fact that at least two owners fail in every context.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 15
LABELS = tuple(range(1, P))
DIVISORS = (1, 3, 5, 15)
STATES = tuple(
    (divisor, unit)
    for divisor in DIVISORS
    for unit in range(divisor)
    if (divisor == 1 and unit == 0) or gcd(unit, divisor) == 1
)
STATE_INDICES = {
    divisor: tuple(
        index for index, state in enumerate(STATES) if state[0] == divisor
    )
    for divisor in DIVISORS
}
FULL = (1 << C) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def algebraic_crt_base(label: int, divisor: int, unit: int) -> int:
    if divisor == 1:
        answer = label
    else:
        correction = unit * pow(P, -1, divisor) % divisor
        answer = (divisor * label + P * correction) % (P * divisor)
    require(answer % P == divisor * label % P, "CRT label mismatch")
    require(answer % divisor == unit % divisor, "CRT unit mismatch")
    return answer


def sheet_mask(label: int, state_index: int, owner: int) -> int:
    divisor, unit = STATES[state_index]
    base = algebraic_crt_base(label, divisor, unit)
    owner_inverse = pow(owner, -1, P)
    return sum(
        1 << sheet
        for sheet in range(C)
        if -divisor
        < centered(base * (owner_inverse + P * sheet), P * divisor)
        <= divisor
    )


MASK = {
    (provider, state, owner): sheet_mask(provider, state, owner)
    for provider in LABELS
    for state in range(len(STATES))
    for owner in LABELS
}
MASK_DIGEST = sha256(
    b"".join(
        MASK[provider, state, owner].to_bytes(2, "little")
        for provider in LABELS
        for state in range(len(STATES))
        for owner in LABELS
    )
).hexdigest()


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


CARDINALITY = {}
for divisor in DIVISORS:
    for provider_owner_ratio in LABELS:
        sizes = {
            MASK[provider_owner_ratio, state, 1].bit_count()
            for state in STATE_INDICES[divisor]
        }
        require(len(sizes) == 1, "mask cardinality depends on unit")
        CARDINALITY[divisor, provider_owner_ratio] = sizes.pop()

for provider in LABELS:
    for owner in LABELS:
        for divisor in DIVISORS:
            literal_sizes = {
                MASK[provider, state, owner].bit_count()
                for state in STATE_INDICES[divisor]
            }
            require(
                literal_sizes
                == {CARDINALITY[divisor, ratio(provider, owner)]},
                "provider-owner ratio reduction failed",
            )


ORDER_WORDS = tuple(
    word
    for word in product(DIVISORS, repeat=6)
    if all(
        lcm(*(word[:omitted] + word[omitted + 1 :])) == C
        for omitted in range(6)
    )
)
STATE_WORDS_PER_SUPPORT = sum(
    prod(len(STATE_INDICES[divisor]) for divisor in word)
    for word in ORDER_WORDS
)
ORDER_DIGEST = sha256(
    "\n".join(",".join(map(str, word)) for word in ORDER_WORDS).encode()
).hexdigest()


def scalar_capacity(
    labels: tuple[int, ...], orders: tuple[int, ...]
) -> bool:
    return all(
        sum(
            CARDINALITY[divisor, ratio(provider, owner)]
            for provider, divisor in zip(labels, orders)
        )
        >= C
        for owner in labels
    )


def owner_reachable(
    labels: tuple[int, ...], orders: tuple[int, ...], owner: int
) -> set[int]:
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
    return reachable


def order_pattern(orders: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(orders.count(divisor) for divisor in DIVISORS)


SUPPORTS = tuple(combinations(LABELS, 6))
SCALAR_BANK = []
SCALAR_PATTERNS = Counter()
FEASIBLE_OWNER_COUNT = Counter()
FEASIBLE_SUBSET_SIZE = Counter()
OWNER_MAXIMUM = Counter()
TOURNAMENT_FLIPS = Counter()
MAXIMUM_VECTOR_DIGESTER = sha256()
for labels in SUPPORTS:
    for orders in ORDER_WORDS:
        if not scalar_capacity(labels, orders):
            continue
        SCALAR_BANK.append((labels, orders))
        SCALAR_PATTERNS[order_pattern(orders)] += 1
        reachable = tuple(
            owner_reachable(labels, orders, owner) for owner in labels
        )
        maxima = tuple(max(mask.bit_count() for mask in masks) for masks in reachable)
        feasible = tuple(FULL in masks for masks in reachable)
        feasible_count = sum(feasible)
        FEASIBLE_OWNER_COUNT[feasible_count] += 1
        FEASIBLE_SUBSET_SIZE[sum(1 << index for index, hit in enumerate(feasible) if hit)] += 1
        OWNER_MAXIMUM.update(maxima)
        MAXIMUM_VECTOR_DIGESTER.update(bytes(maxima))

        # Pair observable: larger maximum reachable union wins.  Gauge and tie
        # Hamiltonian order: increasing owner label.  The completion sorts by
        # (-maximum,label), so it is transitive; flips measure deviation from
        # the declared label path but play no role in the proof.
        flips = sum(
            maxima[first] < maxima[second]
            for first, second in combinations(range(6), 2)
        )
        TOURNAMENT_FLIPS[flips] += 1


SCALAR_SUPPORTS = {labels for labels, _ in SCALAR_BANK}


def multiply_support(
    labels: tuple[int, ...], multiplier: int
) -> tuple[int, ...]:
    return tuple(sorted(multiplier * label % P for label in labels))


remaining = set(SCALAR_SUPPORTS)
SUPPORT_ORBITS = []
while remaining:
    representative = min(remaining)
    orbit = {
        multiply_support(representative, multiplier)
        for multiplier in LABELS
    }
    require(orbit <= remaining, "multiplication orbit leaves scalar bank")
    SUPPORT_ORBITS.append((representative, len(orbit)))
    remaining -= orbit


EXPECTED_PATTERNS = {
    (0, 0, 2, 4): 12,
    (0, 1, 3, 2): 48,
    (0, 1, 4, 1): 84,
    (0, 2, 1, 3): 144,
    (0, 2, 2, 2): 276,
    (0, 2, 3, 1): 24,
    (0, 2, 4, 0): 6,
    (0, 3, 0, 3): 12,
    (0, 3, 1, 2): 36,
    (0, 3, 2, 1): 36,
    (0, 3, 3, 0): 12,
    (0, 4, 0, 2): 18,
    (0, 4, 1, 1): 792,
    (0, 4, 2, 0): 180,
    (1, 3, 1, 1): 336,
    (1, 3, 2, 0): 168,
}
require(len(STATES) == 15, "state alphabet mismatch")
require(len(SUPPORTS) == 924, "support census mismatch")
require(len(ORDER_WORDS) == 3_249, "hereditary order-word census mismatch")
require(STATE_WORDS_PER_SUPPORT == 11_169_600, "state-word census mismatch")
require(len(SCALAR_BANK) == 2_184, "scalar bank mismatch")
require(len(SCALAR_SUPPORTS) == 462, "scalar support census mismatch")
require(SCALAR_PATTERNS == EXPECTED_PATTERNS, "scalar pattern mismatch")
require(
    Counter(size for _, size in SUPPORT_ORBITS) == {6: 1, 12: 38},
    "scalar support multiplication-orbit census mismatch",
)
require(
    FEASIBLE_OWNER_COUNT == {0: 750, 1: 456, 2: 912, 4: 66},
    "locally feasible owner-count histogram mismatch",
)
require(
    OWNER_MAXIMUM == {11: 804, 12: 7_512, 13: 1_812, 14: 432, 15: 2_544},
    "owner maximum-union histogram mismatch",
)
require(max(FEASIBLE_OWNER_COUNT) == 4, "a context is feasible at too many owners")


SCALAR_BANK_DIGEST = sha256(
    "\n".join(
        f"{','.join(map(str, labels))}|{','.join(map(str, orders))}"
        for labels, orders in SCALAR_BANK
    ).encode()
).hexdigest()


print("scale-fifteen independent Python referee")
print("divisor grammar 1,3,5,15; literal states 15")
print(f"supports {len(SUPPORTS)}")
print(f"hereditary order words {len(ORDER_WORDS)}")
print(f"hereditary state words per support {STATE_WORDS_PER_SUPPORT}")
print(f"hereditary labelled state contexts {len(SUPPORTS) * STATE_WORDS_PER_SUPPORT}")
print(f"scalar contexts {len(SCALAR_BANK)} on {len(SCALAR_SUPPORTS)} supports across {len(SCALAR_PATTERNS)} patterns")
print(f"locally feasible owners per context {dict(sorted(FEASIBLE_OWNER_COUNT.items()))}")
print(f"maximum reachable sheet-union histogram {dict(sorted(OWNER_MAXIMUM.items()))}")
print(f"scalar support multiplication orbits {len(SUPPORT_ORBITS)} sizes {','.join(str(size) for _, size in SUPPORT_ORBITS)}")
print(f"tournament pair observable maximum-union comparison; owner-label gauge and tie paths within equality blocks; transitive scores 0,1,2,3,4,5; cycles 0; SCC 1^6; paths 1")
print(f"tournament edge-flip histogram {dict(sorted(TOURNAMENT_FLIPS.items()))}")
print("faithful terminal state is the six-owner feasibility subset, strengthened by the exact maximum-union vector")
print("the tournament ranks owners but forgets threshold fifteen and cannot certify that at least two owners fail")
print(f"mask_sha256 {MASK_DIGEST}")
print(f"order_sha256 {ORDER_DIGEST}")
print(f"scalar_bank_sha256 {SCALAR_BANK_DIGEST}")
print(f"maximum_vector_sha256 {MAXIMUM_VECTOR_DIGESTER.hexdigest()}")
