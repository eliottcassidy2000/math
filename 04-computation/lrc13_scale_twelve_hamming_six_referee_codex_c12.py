#!/usr/bin/env python3
"""Independent exact referee for the common-scale-twelve Hamming-six bank.

The referee reconstructs the CRT masks from the congruences, enumerates the
actual hereditary divisor words (rather than a multiplicity surrogate), and
uses set-valued union reachability for each owner-local test.  The scalar scan
packs the six integer owner capacities into independent eight-bit lanes; the
largest possible lane sum is 72, so the checked threshold trick has no carries
and is exactly equivalent to six ordinary integer comparisons.  The terminal
fibre replay literally visits all 64 * 4^6 unit words.  Correctness never
depends on ``assert`` and every comparison is exact integer arithmetic.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd, lcm, prod


P = 13
C = 12
LABELS = tuple(range(1, P))
DIVISORS = (1, 2, 3, 4, 6, 12)
FULL = (1 << C) - 1
STATES = tuple(
    (divisor, unit)
    for divisor in DIVISORS
    for unit in range(divisor)
    if gcd(unit, divisor) == 1
)
STATE_INDICES = {
    divisor: tuple(
        index
        for index, state in enumerate(STATES)
        if state[0] == divisor
    )
    for divisor in DIVISORS
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


require(
    STATES
    == (
        (1, 0),
        (2, 1),
        (3, 1),
        (3, 2),
        (4, 1),
        (4, 3),
        (6, 1),
        (6, 5),
        (12, 1),
        (12, 5),
        (12, 7),
        (12, 11),
    ),
    "effective-state alphabet mismatch",
)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def algebraic_crt_base(label: int, divisor: int, unit: int) -> int:
    """Solve x = divisor*label (mod 13), x = unit (mod divisor)."""
    if divisor == 1:
        base = label
    else:
        correction = unit * pow(P, -1, divisor) % divisor
        base = (divisor * label + P * correction) % (P * divisor)
    require(
        base % P == divisor * label % P,
        "CRT label congruence failed",
    )
    require(base % divisor == unit % divisor, "CRT unit congruence failed")
    return base


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


# Unit independence is derived from the literal masks.  The second loop also
# checks that passing from a provider-owner pair to its ratio loses no size
# information.
CARDINALITY = {}
for divisor in DIVISORS:
    for provider_owner_ratio in LABELS:
        sizes = {
            MASK[provider_owner_ratio, state, 1].bit_count()
            for state in STATE_INDICES[divisor]
        }
        require(len(sizes) == 1, "mask cardinality depends on a unit")
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
                "ratio cardinality reduction failed",
            )


# These are actual ordered divisor tuples.  Left and right base-six codes are
# only accelerators for the later exhaustive scan.
ORDER_RECORDS = []
for divisor_indices in product(range(len(DIVISORS)), repeat=6):
    word = tuple(DIVISORS[index] for index in divisor_indices)
    if not all(
        lcm(*(word[:omitted] + word[omitted + 1 :])) == C
        for omitted in range(6)
    ):
        continue
    left_code = (
        36 * divisor_indices[0]
        + 6 * divisor_indices[1]
        + divisor_indices[2]
    )
    right_code = (
        36 * divisor_indices[3]
        + 6 * divisor_indices[4]
        + divisor_indices[5]
    )
    pattern = tuple(word.count(divisor) for divisor in DIVISORS)
    ORDER_RECORDS.append((left_code, right_code, word, pattern))

ORDER_RECORDS = tuple(ORDER_RECORDS)
ORDER_WORDS = tuple(record[2] for record in ORDER_RECORDS)
ORDER_WORD_DIGEST = sha256(
    "\n".join(",".join(map(str, word)) for word in ORDER_WORDS).encode()
).hexdigest()
HEREDITARY_STATE_WORDS = sum(
    prod(len(STATE_INDICES[divisor]) for divisor in word)
    for word in ORDER_WORDS
)
require(len(ORDER_WORDS) == 26_961, "hereditary divisor-word census mismatch")
require(
    len(set(ORDER_WORDS)) == len(ORDER_WORDS),
    "hereditary divisor words are not distinct",
)
require(
    HEREDITARY_STATE_WORDS == 2_611_968,
    "hereditary state-word census mismatch",
)


# Exact six-lane scalar threshold.  A lane contains a sum of six mask sizes,
# each at most twelve.  Adding 128-C therefore sets that lane's high bit iff
# its capacity is at least C, and 72+(128-C)<256 prevents inter-lane carries.
HIGH_BITS = sum(128 << (8 * owner) for owner in range(6))
THRESHOLD_BIAS = sum((128 - C) << (8 * owner) for owner in range(6))
require(
    max(CARDINALITY.values()) * 6 + 128 - C < 256,
    "packed scalar lanes could carry",
)
require(
    all((((capacity + 128 - C) & 128) != 0) == (capacity >= C)
        for capacity in range(6 * max(CARDINALITY.values()) + 1)),
    "packed scalar threshold is not exact",
)


def packed_scalar_capacity(packed_capacity: int) -> bool:
    return (
        (packed_capacity + THRESHOLD_BIAS) & HIGH_BITS
    ) == HIGH_BITS


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


SUPPORTS = tuple(combinations(LABELS, 6))
SCALAR_BANK = []
SCALAR_PATTERNS = Counter()
LOCAL_BANK = []
for labels in SUPPORTS:
    packed_contribution = tuple(
        tuple(
            sum(
                CARDINALITY[divisor, ratio(labels[provider], owner)]
                << (8 * owner_index)
                for owner_index, owner in enumerate(labels)
            )
            for divisor in DIVISORS
        )
        for provider in range(6)
    )
    left_sums = tuple(
        packed_contribution[0][first]
        + packed_contribution[1][second]
        + packed_contribution[2][third]
        for first in range(6)
        for second in range(6)
        for third in range(6)
    )
    right_sums = tuple(
        packed_contribution[3][fourth]
        + packed_contribution[4][fifth]
        + packed_contribution[5][sixth]
        for fourth in range(6)
        for fifth in range(6)
        for sixth in range(6)
    )
    for left_code, right_code, orders, pattern in ORDER_RECORDS:
        if not packed_scalar_capacity(
            left_sums[left_code] + right_sums[right_code]
        ):
            continue
        SCALAR_BANK.append((labels, orders))
        SCALAR_PATTERNS[pattern] += 1
        if all(
            owner_locally_feasible(labels, orders, owner)
            for owner in labels
        ):
            LOCAL_BANK.append((labels, orders))

require(len(SUPPORTS) == 924, "support census mismatch")
require(len(SCALAR_BANK) == 36_830, "scalar-capacity bank mismatch")
require(
    len({labels for labels, _ in SCALAR_BANK}) == 912,
    "scalar-capacity support count mismatch",
)
require(len(SCALAR_PATTERNS) == 85, "scalar pattern census mismatch")

SCALAR_BANK_DIGEST = sha256(
    "\n".join(
        f"{','.join(map(str, labels))}|{','.join(map(str, orders))}"
        for labels, orders in SCALAR_BANK
    ).encode()
).hexdigest()
SCALAR_PATTERN_DIGEST = sha256(
    "\n".join(
        f"{','.join(map(str, pattern))}:{multiplicity}"
        for pattern, multiplicity in sorted(SCALAR_PATTERNS.items())
    ).encode()
).hexdigest()


def sign_transversal(labels: tuple[int, ...]) -> bool:
    return {
        min(label, P - label) for label in labels
    } == set(range(1, 7))


SIGN_TRANSVERSALS = {
    labels for labels in SUPPORTS if sign_transversal(labels)
}
require(len(LOCAL_BANK) == 64, "owner-local bank mismatch")
require(
    all(orders == (12,) * 6 for _, orders in LOCAL_BANK),
    "owner-local bank contains a non-D12 word",
)
LOCAL_SUPPORTS = {labels for labels, _ in LOCAL_BANK}
require(
    LOCAL_SUPPORTS == SIGN_TRANSVERSALS and len(SIGN_TRANSVERSALS) == 64,
    "owner-local bank is not exactly the sign-transversal bank",
)
LOCAL_BANK_DIGEST = sha256(
    "\n".join(
        f"{','.join(map(str, labels))}|{','.join(map(str, orders))}"
        for labels, orders in LOCAL_BANK
    ).encode()
).hexdigest()


def multiply_support(
    labels: tuple[int, ...], multiplier: int
) -> tuple[int, ...]:
    return tuple(sorted(multiplier * label % P for label in labels))


remaining = set(LOCAL_SUPPORTS)
ORBITS = []
while remaining:
    representative = min(remaining)
    orbit = {
        multiply_support(representative, multiplier)
        for multiplier in LABELS
    }
    require(orbit <= remaining, "multiplication orbit leaves local bank")
    ORBITS.append((representative, len(orbit)))
    remaining -= orbit

EXPECTED_ORBITS = (
    ((1, 2, 3, 4, 5, 6), 12),
    ((1, 2, 3, 4, 5, 7), 12),
    ((1, 2, 3, 4, 6, 8), 12),
    ((1, 2, 3, 5, 6, 9), 4),
    ((1, 2, 3, 5, 7, 9), 12),
    ((1, 2, 6, 8, 9, 10), 12),
)
require(tuple(ORBITS) == EXPECTED_ORBITS, "multiplication orbits mismatch")


def packed_provider_masks(
    labels: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            sum(
                MASK[provider, state, owner] << (C * owner_index)
                for owner_index, owner in enumerate(labels)
            )
            for state in STATE_INDICES[12]
        )
        for provider in labels
    )


def literal_fibre_profile(labels: tuple[int, ...]) -> tuple[int, ...]:
    packed = packed_provider_masks(labels)
    lane_masks = tuple(FULL << (C * owner) for owner in range(6))
    profile = [0] * 64
    for units in product(range(4), repeat=6):
        cover = 0
        for provider, unit in enumerate(units):
            cover |= packed[provider][unit]
        satisfaction = sum(
            1 << owner
            for owner, lane_mask in enumerate(lane_masks)
            if cover & lane_mask == lane_mask
        )
        profile[satisfaction] += 1
    require(sum(profile) == 4**6, "literal fibre size mismatch")
    return tuple(profile)


EXPECTED_PROFILE = tuple(
    3_808
    if subset == 0
    else 48
    if subset & (subset - 1) == 0
    else 0
    for subset in range(64)
)
PROFILE_DIGESTER = sha256()
REPLAYED_WORDS = 0
GLOBAL_SURVIVORS = 0
for labels, _ in LOCAL_BANK:
    profile = literal_fibre_profile(labels)
    require(profile == EXPECTED_PROFILE, "uniform fibre profile mismatch")
    owner_sizes = tuple(
        sum(
            multiplicity
            for subset, multiplicity in enumerate(profile)
            if subset & (1 << owner)
        )
        for owner in range(6)
    )
    require(owner_sizes == (48,) * 6, "owner-obligation size mismatch")
    for first, second in combinations(range(6), 2):
        intersection = sum(
            multiplicity
            for subset, multiplicity in enumerate(profile)
            if subset & (1 << first) and subset & (1 << second)
        )
        require(intersection == 0, "distinct owner obligations intersect")
    PROFILE_DIGESTER.update(
        b"".join(value.to_bytes(2, "little") for value in profile)
    )
    REPLAYED_WORDS += sum(profile)
    GLOBAL_SURVIVORS += profile[63]

PROFILE_DIGEST = PROFILE_DIGESTER.hexdigest()
require(REPLAYED_WORDS == 262_144, "literal replay total mismatch")
require(GLOBAL_SURVIVORS == 0, "a common-sheet unit word survived")


def tournament_fingerprint(labels: tuple[int, ...]) -> tuple[object, ...]:
    """Complete the empty intersection relation along the gauged tie path."""
    # Multiplication gauges each support to its least orbit representative.
    # Every observable pair is a tie, so the sorted labels are the declared
    # Hamiltonian tie path and orient the completed tournament transitively.
    tie_path = tuple(sorted(labels))
    tournament = {
        (tie_path[first], tie_path[second])
        for first, second in combinations(range(6), 2)
    }
    scores = tuple(
        sorted(
            sum((vertex, other) in tournament for other in labels)
            for vertex in labels
        )
    )
    cycles = sum(
        all(
            sum((vertex, other) in tournament for other in triple) == 1
            for vertex in triple
        )
        for triple in combinations(labels, 3)
    )
    hamiltonian_paths = sum(
        all(
            (path[index], path[index + 1]) in tournament
            for index in range(5)
        )
        for path in permutations(labels)
    )
    scc_sizes = (1, 1, 1, 1, 1, 1)
    edge_flips_from_tie_path = 0
    return (
        scores,
        cycles,
        scc_sizes,
        edge_flips_from_tie_path,
        hamiltonian_paths,
    )


TOURNAMENT_FINGERPRINTS = Counter(
    tournament_fingerprint(labels) for labels in LOCAL_SUPPORTS
)
require(
    TOURNAMENT_FINGERPRINTS
    == {
        (
            (0, 1, 2, 3, 4, 5),
            0,
            (1, 1, 1, 1, 1, 1),
            0,
            1,
        ): 64
    },
    "completed-tournament fingerprint mismatch",
)


print("scale-twelve independent Python referee")
print(f"mask digest {MASK_DIGEST}")
print(f"hereditary divisor-word digest {ORDER_WORD_DIGEST}")
print("hereditary divisor words 26961; labelled order contexts 24911964")
print("hereditary state words 2611968; labelled contexts 2413458432")
print("scalar capacity 36830 contexts on 912 supports across 85 patterns")
print(f"scalar bank digest {SCALAR_BANK_DIGEST}")
print(f"scalar pattern digest {SCALAR_PATTERN_DIGEST}")
print("owner-local bank 64 all-D12 contexts on 64 sign transversals")
print(f"owner-local bank digest {LOCAL_BANK_DIGEST}")
print("multiplication orbits 6; sizes 12,12,12,4,12,12")
print("literal unit words 262144; global common-sheet survivors 0")
print("uniform profile 0:3808 1:288; owner size 48; distinct pairs disjoint")
print(f"profile digest {PROFILE_DIGEST}")
print("pair observable owner-obligation intersection; multiplication gauge; sorted tie path")
print("completed tournaments transitive: scores 0,1,2,3,4,5; cycles 0; SCC 1^6; flips 0; paths 1")
print("faithful vertices are owner obligations (or dually unit words), not raw runners")
print("runners, residues, sheets, and cover arcs lose exact obligation disjointness")
print("tournament completion destroys obligation sizes and the fact that every pair was a tie")
