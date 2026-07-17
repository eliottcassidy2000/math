#!/usr/bin/env python3
"""Independent exact referee for the common-scale-ten Hamming-six bank.

Unlike the C++ certificate, this implementation stores actual divisors in the
order words, constructs CRT masks as Python integers, uses set-valued union-mask
reachability owner by owner, and packs all six owners only for the final fibre
replay.  All arithmetic and all comparisons are exact.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from math import lcm, prod


P = 13
C = 10
LABELS = tuple(range(1, P))
DIVISORS = (1, 2, 5, 10)
FULL = (1 << C) - 1
STATES = (
    (1, 0),
    (2, 1),
    (5, 1),
    (5, 2),
    (5, 3),
    (5, 4),
    (10, 1),
    (10, 3),
    (10, 7),
    (10, 9),
)
STATE_INDICES = {
    divisor: tuple(i for i, state in enumerate(STATES) if state[0] == divisor)
    for divisor in DIVISORS
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


CARDINALITY = {}
for divisor in DIVISORS:
    for ratio in LABELS:
        values = {
            MASK[ratio, state, 1].bit_count()
            for state in STATE_INDICES[divisor]
        }
        require(len(values) == 1, "mask cardinality depends on the unit")
        CARDINALITY[divisor, ratio] = values.pop()


ORDER_WORDS = tuple(
    word
    for word in product(DIVISORS, repeat=6)
    if all(lcm(*(word[:omitted] + word[omitted + 1 :])) == C
           for omitted in range(6))
)
require(len(ORDER_WORDS) == 3_249, "hereditary order grammar mismatch")
HEREDITARY_STATE_WORDS = sum(
    prod(len(STATE_INDICES[divisor]) for divisor in word)
    for word in ORDER_WORDS
)
require(HEREDITARY_STATE_WORDS == 889_200,
        "hereditary state-word census mismatch")


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


def scalar_capacity(labels: tuple[int, ...], orders: tuple[int, ...]) -> bool:
    return all(
        sum(
            CARDINALITY[divisor, ratio(provider, owner)]
            for provider, divisor in zip(labels, orders)
        )
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


def order_pattern(orders: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(orders.count(divisor) for divisor in DIVISORS)


def sign_transversal(labels: tuple[int, ...]) -> bool:
    return {min(label, P - label) for label in labels} == set(range(1, 7))


SUPPORTS = tuple(combinations(LABELS, 6))
SCALAR_BANK = []
LOCAL_BANK = []
SCALAR_PATTERN = Counter()
for labels in SUPPORTS:
    for orders in ORDER_WORDS:
        if not scalar_capacity(labels, orders):
            continue
        SCALAR_BANK.append((labels, orders))
        SCALAR_PATTERN[order_pattern(orders)] += 1
        if all(owner_locally_feasible(labels, orders, owner) for owner in labels):
            LOCAL_BANK.append((labels, orders))

require(len(SCALAR_BANK) == 1_200, "scalar-capacity bank mismatch")
require(len({labels for labels, _ in SCALAR_BANK}) == 388,
        "scalar-capacity support count mismatch")
require(
    SCALAR_PATTERN
    == {
        (0, 0, 0, 6): 64,
        (0, 0, 3, 3): 120,
        (0, 2, 0, 4): 36,
        (0, 2, 2, 2): 48,
        (0, 2, 3, 1): 48,
        (0, 2, 4, 0): 12,
        (0, 3, 0, 3): 344,
        (0, 3, 1, 2): 336,
        (0, 3, 2, 1): 144,
        (0, 3, 3, 0): 48,
    },
    "scalar order-pattern histogram mismatch",
)
require(len(LOCAL_BANK) == 64, "owner-local bank mismatch")
require(all(orders == (10,) * 6 for _, orders in LOCAL_BANK),
        "owner-local bank contains a non-D10 word")
LOCAL_SUPPORTS = {labels for labels, _ in LOCAL_BANK}
SIGN_TRANSVERSALS = {labels for labels in SUPPORTS if sign_transversal(labels)}
require(LOCAL_SUPPORTS == SIGN_TRANSVERSALS and len(SIGN_TRANSVERSALS) == 64,
        "owner-local bank is not the sign-transversal bank")

LOCAL_BANK_DIGEST = sha256(
    "\n".join(
        f"{','.join(map(str, labels))}|{','.join(map(str, orders))}"
        for labels, orders in LOCAL_BANK
    ).encode()
).hexdigest()


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


def fibre_profile(labels: tuple[int, ...], orders: tuple[int, ...]) -> Counter[int]:
    packed = packed_provider_masks(labels)
    profile = Counter()
    for states in product(*(STATE_INDICES[divisor] for divisor in orders)):
        cover = 0
        for provider in range(6):
            cover |= packed[provider][states[provider]]
        profile[satisfaction_subset(cover)] += 1
    return profile


def zero_pair(first: int, second: int) -> bool:
    return second * pow(first, -1, P) % P in (2, 6, 7, 11)


GLOBAL_WORDS = 0
GLOBAL_SURVIVORS = 0
PROFILE_DIGESTER = sha256()
for labels, orders in LOCAL_BANK:
    profile = fibre_profile(labels, orders)
    GLOBAL_WORDS += sum(profile.values())
    GLOBAL_SURVIVORS += profile[63]
    PROFILE_DIGESTER.update(
        b"".join(profile[subset].to_bytes(2, "little") for subset in range(64))
    )
    by_count = Counter()
    for subset, multiplicity in profile.items():
        by_count[subset.bit_count()] += multiplicity
    require(by_count == {0: 3_940, 1: 120, 2: 36},
            "owner-satisfaction profile mismatch")
    owner_sizes = tuple(
        sum(value for subset, value in profile.items() if subset & (1 << owner))
        for owner in range(6)
    )
    require(owner_sizes == (32,) * 6, "owner-obligation size mismatch")
    for i, j in combinations(range(6), 2):
        intersection = sum(
            value
            for subset, value in profile.items()
            if subset & (1 << i) and subset & (1 << j)
        )
        require(
            intersection == (0 if zero_pair(labels[i], labels[j]) else 4),
            "owner-obligation pair nerve mismatch",
        )

require(GLOBAL_WORDS == 262_144 and GLOBAL_SURVIVORS == 0,
        "reduced literal replay mismatch")
PROFILE_DIGEST = PROFILE_DIGESTER.hexdigest()


PROJECTIVE_CYCLE = (1, 2, 4, 5, 3, 6)
PROJECTIVE_SUCCESSOR = {
    PROJECTIVE_CYCLE[i]: PROJECTIVE_CYCLE[(i + 1) % 6]
    for i in range(6)
}


def tournament_fingerprint(labels: tuple[int, ...]):
    projective = lambda label: min(label, P - label)
    sparse = {
        (first, second)
        for first in labels
        for second in labels
        if first != second
        and PROJECTIVE_SUCCESSOR[projective(first)] == projective(second)
    }
    require(len(sparse) == 6, "gauged sparse relation is not C6")
    tie_path = next(
        path
        for path in permutations(labels)
        if all(
            (path[i], path[i + 1]) not in sparse
            and (path[i + 1], path[i]) not in sparse
            for i in range(5)
        )
    )
    position = {label: i for i, label in enumerate(tie_path)}
    tournament = set()
    flips = 0
    for first, second in combinations(labels, 2):
        if (first, second) in sparse:
            edge = (first, second)
        elif (second, first) in sparse:
            edge = (second, first)
        else:
            edge = ((first, second) if position[first] < position[second]
                    else (second, first))
        tournament.add(edge)
        flips += edge in sparse and position[edge[0]] > position[edge[1]]
    scores = tuple(sorted(
        sum((vertex, other) in tournament for other in labels)
        for vertex in labels
    ))
    triangles = sum(
        all(sum((vertex, other) in tournament for other in triple) == 1
            for vertex in triple)
        for triple in combinations(labels, 3)
    )

    def reachable(source: int, target: int) -> bool:
        seen = {source}
        frontier = [source]
        while frontier:
            vertex = frontier.pop()
            for other in labels:
                if (vertex, other) in tournament and other not in seen:
                    seen.add(other)
                    frontier.append(other)
        return target in seen

    unused = set(labels)
    components = []
    while unused:
        seed = min(unused)
        component = {
            vertex for vertex in labels
            if reachable(seed, vertex) and reachable(vertex, seed)
        }
        components.append(len(component))
        unused -= component
    scc = tuple(sorted(components, reverse=True) + [0] * (6 - len(components)))
    paths = sum(
        all((path[i], path[i + 1]) in tournament for i in range(5))
        for path in permutations(labels)
    )
    return scores, triangles, scc, flips, paths


TOURNAMENTS = Counter(tournament_fingerprint(labels) for labels in LOCAL_SUPPORTS)
require(len(TOURNAMENTS) == 5, "tournament joint fingerprint mismatch")
FLIPS = Counter()
PATHS = Counter()
for fingerprint, multiplicity in TOURNAMENTS.items():
    scores, triangles, scc, flips, paths = fingerprint
    require(scores == (1, 2, 2, 3, 3, 4) and triangles == 6
            and scc == (6, 0, 0, 0, 0, 0),
            "unexpected tournament invariant")
    FLIPS[flips] += multiplicity
    PATHS[paths] += multiplicity
require(FLIPS == {2: 8, 3: 52, 4: 4}
        and PATHS == {29: 32, 31: 20, 37: 12},
        "tournament telemetry mismatch")


print("scale-ten independent Python referee")
print(f"mask digest {MASK_DIGEST}")
print("hereditary order words 3249; labelled contexts 3002076")
print("hereditary state words 889200; labelled contexts 821620800")
print("scalar capacity 1200 contexts on 388 supports")
print("owner-local bank 64 all-D10 contexts on 64 sign-transversal supports")
print(f"owner-local bank digest {LOCAL_BANK_DIGEST}")
print("reduced literal state words 262144; global common-sheet survivors 0")
print("obligation profile 3940,120,36; owners 32; nerve K6\\C6 with edge size 4")
print(f"profile digest {PROFILE_DIGEST}")
print("tournament joint fingerprints 5; flips 2:8,3:52,4:4; paths 29:32,31:20,37:12")
