#!/usr/bin/env python3
"""Exact exploratory scout for common-scale-fourteen Hamming-six packets.

The script reconstructs the hereditary divisor grammar, scalar owner capacity,
owner-local mask reachability, and every surviving global unit fibre.  Runner
supports are useful algebraic coordinates but do not preserve simultaneous
coverage; owner obligations do preserve the terminal predicate only through
their full intersection nerve.  The reported pair graph and a deterministic
tournament completion are therefore telemetry, not the certificate itself.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd, lcm, prod


P = 13
C = 14
LABELS = tuple(range(1, P))
DIVISORS = tuple(divisor for divisor in range(1, C + 1) if C % divisor == 0)
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


def crt_base(label: int, divisor: int, unit: int) -> int:
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


CARDINALITY = {}
for divisor in DIVISORS:
    for provider_owner_ratio in LABELS:
        sizes = {
            MASK[provider_owner_ratio, state, 1].bit_count()
            for state in STATE_INDICES[divisor]
        }
        require(len(sizes) == 1, "mask cardinality depends on unit")
        CARDINALITY[divisor, provider_owner_ratio] = sizes.pop()


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


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


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
LOCAL_BANK = []
SCALAR_PATTERNS = Counter()
LOCAL_OWNER_COUNT = Counter()
OWNER_MAXIMUM = Counter()
for labels in SUPPORTS:
    for orders in ORDER_WORDS:
        if not scalar_capacity(labels, orders):
            continue
        SCALAR_BANK.append((labels, orders))
        SCALAR_PATTERNS[order_pattern(orders)] += 1
        owner_reachability = tuple(
            owner_reachable(labels, orders, owner) for owner in labels
        )
        owner_feasibility = tuple(FULL in reachable for reachable in owner_reachability)
        LOCAL_OWNER_COUNT[sum(owner_feasibility)] += 1
        OWNER_MAXIMUM.update(
            max(mask.bit_count() for mask in reachable)
            for reachable in owner_reachability
        )
        if all(owner_feasibility):
            LOCAL_BANK.append((labels, orders))


SCALAR_SUPPORTS = {labels for labels, _ in SCALAR_BANK}


def multiply_support(
    labels: tuple[int, ...], multiplier: int
) -> tuple[int, ...]:
    return tuple(sorted(multiplier * label % P for label in labels))


remaining_supports = set(SCALAR_SUPPORTS)
SUPPORT_ORBITS = []
while remaining_supports:
    representative = min(remaining_supports)
    orbit = {
        multiply_support(representative, multiplier)
        for multiplier in LABELS
    }
    require(orbit <= remaining_supports, "scalar support orbit leaves the bank")
    SUPPORT_ORBITS.append((representative, len(orbit)))
    remaining_supports -= orbit


require(len(STATES) == 14, "state alphabet mismatch")
require(len(SUPPORTS) == 924, "support census mismatch")
require(len(ORDER_WORDS) == 3_249, "hereditary order-word census mismatch")
require(STATE_WORDS_PER_SUPPORT == 6_703_884, "state-word census mismatch")
require(len(SCALAR_BANK) == 576, "scalar-capacity bank mismatch")
require(len(SCALAR_SUPPORTS) == 36, "scalar support census mismatch")
require(
    SCALAR_PATTERNS
    == {
        (0, 2, 0, 4): 36,
        (0, 2, 1, 3): 144,
        (0, 2, 2, 2): 216,
        (0, 2, 3, 1): 144,
        (0, 2, 4, 0): 36,
    },
    "scalar order-pattern histogram mismatch",
)
require(
    tuple(size for _, size in SUPPORT_ORBITS) == (12, 12, 12),
    "scalar support multiplication-orbit census mismatch",
)
require(LOCAL_OWNER_COUNT == {0: 576}, "an owner-local row unexpectedly survived")
require(
    OWNER_MAXIMUM == {10: 96, 11: 1_056, 12: 2_304},
    "owner-local maximum-union histogram mismatch",
)
require(not LOCAL_BANK, "owner-local bank is nonempty")


def packed_provider_masks(labels: tuple[int, ...]):
    return tuple(
        tuple(
            sum(
                MASK[provider, state, owner] << (C * owner_index)
                for owner_index, owner in enumerate(labels)
            )
            for state in range(len(STATES))
        )
        for provider in labels
    )


def satisfaction_subset(packed_cover: int) -> int:
    return sum(
        1 << owner
        for owner in range(6)
        if (packed_cover >> (C * owner)) & FULL == FULL
    )


def fibre_profile(
    labels: tuple[int, ...], orders: tuple[int, ...]
) -> Counter[int]:
    packed = packed_provider_masks(labels)
    profile = Counter()
    for states in product(*(STATE_INDICES[divisor] for divisor in orders)):
        cover = 0
        for provider in range(6):
            cover |= packed[provider][states[provider]]
        profile[satisfaction_subset(cover)] += 1
    return profile


GLOBAL_WORDS = 0
GLOBAL_SURVIVORS = 0
PROFILE_SPECIES = Counter()
PROFILE_DIGESTER = sha256()
PAIR_GRAPH_SPECIES = Counter()
TOURNAMENT_SPECIES = Counter()
for labels, orders in LOCAL_BANK:
    profile = fibre_profile(labels, orders)
    words = sum(profile.values())
    GLOBAL_WORDS += words
    GLOBAL_SURVIVORS += profile[63]
    PROFILE_DIGESTER.update(
        b"".join(profile[subset].to_bytes(8, "little") for subset in range(64))
    )
    by_count = tuple(
        sum(value for subset, value in profile.items() if subset.bit_count() == k)
        for k in range(7)
    )
    owner_sizes = tuple(
        sum(value for subset, value in profile.items() if subset & (1 << owner))
        for owner in range(6)
    )
    pair_edges = tuple(
        (i, j)
        for i, j in combinations(range(6), 2)
        if any(
            value and subset & (1 << i) and subset & (1 << j)
            for subset, value in profile.items()
        )
    )
    PROFILE_SPECIES[(by_count, owner_sizes)] += 1
    PAIR_GRAPH_SPECIES[pair_edges] += 1

    # Complete the symmetric pair observable along the label order.  This is
    # deliberately lossy telemetry; the full profile above is the proof data.
    tournament = set()
    for i, j in combinations(range(6), 2):
        if (i, j) in pair_edges:
            tournament.add((i, j) if labels[i] < labels[j] else (j, i))
        else:
            tournament.add((i, j))
    scores = tuple(sorted(sum((i, j) in tournament for j in range(6)) for i in range(6)))
    triangles = sum(
        all(sum((i, j) in tournament for j in triple if j != i) == 1 for i in triple)
        for triple in combinations(range(6), 3)
    )
    paths = sum(
        all((path[i], path[i + 1]) in tournament for i in range(5))
        for path in permutations(range(6))
    )
    TOURNAMENT_SPECIES[(scores, triangles, paths)] += 1


LOCAL_BANK_DIGEST = sha256(
    "\n".join(
        f"{','.join(map(str, labels))}|{','.join(map(str, orders))}"
        for labels, orders in LOCAL_BANK
    ).encode()
).hexdigest()


print("scale_fourteen_hamming_six_frontier_scout")
print(f"divisors={','.join(map(str, DIVISORS))}")
print(f"states={len(STATES)}")
print(f"supports={len(SUPPORTS)}")
print(f"hereditary_order_words={len(ORDER_WORDS)}")
print(f"state_words_per_support={STATE_WORDS_PER_SUPPORT}")
print(f"raw_labelled_contexts={len(SUPPORTS) * STATE_WORDS_PER_SUPPORT}")
print(f"scalar_contexts={len(SCALAR_BANK)}")
print(f"scalar_supports={len(SCALAR_SUPPORTS)}")
print(f"scalar_patterns={dict(sorted(SCALAR_PATTERNS.items()))}")
print(f"scalar_support_orbits={SUPPORT_ORBITS}")
print(f"locally_feasible_owner_count={dict(sorted(LOCAL_OWNER_COUNT.items()))}")
print(f"owner_maximum_union={dict(sorted(OWNER_MAXIMUM.items()))}")
print(f"owner_local_contexts={len(LOCAL_BANK)}")
print(f"owner_local_supports={len({labels for labels, _ in LOCAL_BANK})}")
print(f"owner_local_patterns={dict(sorted(Counter(order_pattern(o) for _, o in LOCAL_BANK).items()))}")
print(f"global_unit_words={GLOBAL_WORDS}")
print(f"global_survivors={GLOBAL_SURVIVORS}")
print(f"profile_species={dict(PROFILE_SPECIES)}")
print(f"pair_graph_species={dict(PAIR_GRAPH_SPECIES)}")
print(f"tournament_species={dict(TOURNAMENT_SPECIES)}")
print("tournament_analysis=not_applicable_after_owner_local_elimination; runner and owner completions erase the sheet-cardinality deficit")
print(f"mask_sha256={MASK_DIGEST}")
print(f"local_bank_sha256={LOCAL_BANK_DIGEST}")
print(f"profile_sha256={PROFILE_DIGESTER.hexdigest()}")
