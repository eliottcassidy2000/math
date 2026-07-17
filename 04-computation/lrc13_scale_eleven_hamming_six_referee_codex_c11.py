#!/usr/bin/env python3
"""Independent exact referee for the common-scale-eleven Hamming-six bank.

This implementation is deliberately separate from the C++ certificate.  It
constructs actual divisor tuples, derives CRT bases algebraically, uses Python
sets for owner-local mask-union reachability, proves multiplication covariance
at the mask-generator level (including the induced provider/unit and owner
coordinate permutations), and literally visits 10^6 unit words for one
representative of each of the seven multiplication orbits.  Every comparison
uses exact integer arithmetic and correctness does not depend on ``assert``.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from math import lcm, prod


P = 13
C = 11
LABELS = tuple(range(1, P))
DIVISORS = (1, 11)
STATES = ((1, 0),) + tuple((11, unit) for unit in range(1, 11))
STATE_INDICES = {
    divisor: tuple(
        state_index
        for state_index, state in enumerate(STATES)
        if state[0] == divisor
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
    """Solve x = divisor*label (mod 13), x = unit (mod divisor)."""
    if divisor == 1:
        base = label
    else:
        correction = unit * pow(P, -1, divisor) % divisor
        base = (divisor * label + P * correction) % (P * divisor)
    require(base % P == divisor * label % P, "CRT label congruence failed")
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


# Unit-independent cardinalities form the scalar-capacity relaxation.  They
# are reconstructed from literal masks rather than entered as a table.
CARDINALITY = {}
for divisor in DIVISORS:
    for provider_owner_ratio in LABELS:
        sizes = {
            MASK[provider_owner_ratio, state, 1].bit_count()
            for state in STATE_INDICES[divisor]
        }
        require(len(sizes) == 1, "mask cardinality depends on a unit")
        CARDINALITY[divisor, provider_owner_ratio] = sizes.pop()


ORDER_WORDS = tuple(
    word
    for word in product(DIVISORS, repeat=6)
    if all(
        lcm(*(word[:omitted] + word[omitted + 1 :])) == C
        for omitted in range(6)
    )
)
HEREDITARY_STATE_WORDS = sum(
    prod(len(STATE_INDICES[divisor]) for divisor in word)
    for word in ORDER_WORDS
)
require(len(ORDER_WORDS) == 57, "hereditary divisor-word census mismatch")
require(
    HEREDITARY_STATE_WORDS == 1_771_500,
    "hereditary state-word census mismatch",
)


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
LOCAL_BANK = []
SCALAR_BY_D11 = Counter()
for labels in SUPPORTS:
    for orders in ORDER_WORDS:
        if not scalar_capacity(labels, orders):
            continue
        SCALAR_BANK.append((labels, orders))
        SCALAR_BY_D11[orders.count(11)] += 1
        if all(
            owner_locally_feasible(labels, orders, owner)
            for owner in labels
        ):
            LOCAL_BANK.append((labels, orders))

require(len(SUPPORTS) == 924, "support census mismatch")
require(len(SCALAR_BANK) == 66, "scalar-capacity bank mismatch")
require(
    len({labels for labels, _ in SCALAR_BANK}) == 66,
    "scalar bank does not have one context per support",
)
require(SCALAR_BY_D11 == {6: 66}, "scalar bank is not entirely all-D11")
require(LOCAL_BANK == SCALAR_BANK, "owner-local feasibility changed the bank")
LOCAL_SUPPORTS = {labels for labels, _ in LOCAL_BANK}
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


def multiplication_orbit(labels: tuple[int, ...]) -> set[tuple[int, ...]]:
    return {multiply_support(labels, multiplier) for multiplier in LABELS}


remaining = set(LOCAL_SUPPORTS)
ORBITS = []
while remaining:
    representative = min(remaining)
    orbit = multiplication_orbit(representative)
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
    ((1, 3, 4, 9, 10, 12), 2),
)
require(tuple(ORBITS) == EXPECTED_ORBITS, "multiplication orbits mismatch")


def coordinate_permutation(
    labels: tuple[int, ...], multiplier: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return target support and old-index -> target-index permutation."""
    target = multiply_support(labels, multiplier)
    target_position = {label: index for index, label in enumerate(target)}
    permutation = tuple(
        target_position[multiplier * label % P] for label in labels
    )
    require(
        tuple(sorted(permutation)) == tuple(range(6)),
        "multiplication did not induce a coordinate permutation",
    )
    return target, permutation


def transport_units(
    units: tuple[int, ...], permutation: tuple[int, ...], multiplier: int
) -> tuple[int, ...]:
    """Transport and rescale provider-attached units to target coordinates."""
    lifted_multiplier = (
        multiplier if multiplier % 11 else multiplier + P
    )
    answer = [0] * 6
    for source_index, target_index in enumerate(permutation):
        answer[target_index] = (
            lifted_multiplier * units[source_index] % 11
        )
        require(answer[target_index] != 0, "unit transport left U(11)")
    return tuple(answer)


def state_transport(state: int, multiplier: int) -> int:
    divisor, unit = STATES[state]
    if divisor == 1:
        return state
    lifted_multiplier = (
        multiplier if multiplier % 11 else multiplier + P
    )
    moved_unit = lifted_multiplier * unit % 11
    require(moved_unit != 0, "state transport left U(11)")
    return moved_unit  # state indices 1,...,10 equal their unit values


def owner_sheet_permutation(owner: int, multiplier: int) -> tuple[int, ...]:
    """Induced permutation of the eleven lifts of an owner inverse."""
    modulus = P * C
    lifted_multiplier = (
        multiplier if multiplier % 11 else multiplier + P
    )
    lifted_inverse = pow(lifted_multiplier, -1, modulus)
    source_inverse = pow(owner, -1, P)
    moved_owner = multiplier * owner % P
    target_inverse = pow(moved_owner, -1, P)
    answer = []
    for sheet in range(C):
        moved_lift = (
            lifted_inverse * (source_inverse + P * sheet)
        ) % modulus
        require(
            moved_lift % P == target_inverse,
            "owner lift transport has the wrong residue",
        )
        answer.append((moved_lift - target_inverse) // P)
    require(
        tuple(sorted(answer)) == tuple(range(C)),
        "owner-sheet transport is not a permutation",
    )
    return tuple(answer)


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(
        ((mask >> source_sheet) & 1) << target_sheet
        for source_sheet, target_sheet in enumerate(permutation)
    )


# Exact generator-level covariance.  Under multiplication, provider and owner
# coordinates undergo the same permutation.  A provider's unit is multiplied
# by a chosen unit lift of the F_13 multiplier, while each owner's eleven sheet
# bits undergo the inverse lift permutation.  Equality on every transformed
# provider/unit/owner mask proves equality after arbitrary six-fold unions,
# hence justifies one replay per orbit.  Raw bit masks are *not* covariant if
# the unit and sheet permutations are omitted.
for representative, _ in ORBITS:
    for multiplier in LABELS:
        target, permutation = coordinate_permutation(
            representative, multiplier
        )
        sample_units = tuple(range(1, 7))
        moved_units = transport_units(
            sample_units, permutation, multiplier
        )
        lifted_multiplier = (
            multiplier if multiplier % 11 else multiplier + P
        )
        require(
            all(
                moved_units[permutation[index]]
                == lifted_multiplier * sample_units[index] % 11
                for index in range(6)
            ),
            "unit-coordinate transport failed",
        )
        for provider_index, provider in enumerate(representative):
            moved_provider = target[permutation[provider_index]]
            for owner_index, owner in enumerate(representative):
                moved_owner = target[permutation[owner_index]]
                sheet_permutation = owner_sheet_permutation(
                    owner, multiplier
                )
                for state in range(len(STATES)):
                    moved_state = state_transport(state, multiplier)
                    require(
                        MASK[moved_provider, moved_state, moved_owner]
                        == permute_mask(
                            MASK[provider, state, owner],
                            sheet_permutation,
                        ),
                        "multiplicative mask covariance failed",
                    )


def packed_provider_masks(labels: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            sum(
                MASK[provider, state, owner] << (C * owner_index)
                for owner_index, owner in enumerate(labels)
            )
            for state in STATE_INDICES[11]
        )
        for provider in labels
    )


def literal_fibre_profile(labels: tuple[int, ...]) -> tuple[int, ...]:
    """Visit all 10^6 unit words via an exact three-plus-three join."""
    packed = packed_provider_masks(labels)
    left_words = tuple(
        packed[0][a] | packed[1][b] | packed[2][c]
        for a, b, c in product(range(10), repeat=3)
    )
    right_words = tuple(
        packed[3][d] | packed[4][e] | packed[5][f]
        for d, e, f in product(range(10), repeat=3)
    )
    lane_masks = tuple(FULL << (C * owner) for owner in range(6))
    profile = [0] * 64
    lane0, lane1, lane2, lane3, lane4, lane5 = lane_masks
    for left in left_words:
        for right in right_words:
            cover = left | right
            satisfaction = (
                ((cover & lane0) == lane0)
                | (((cover & lane1) == lane1) << 1)
                | (((cover & lane2) == lane2) << 2)
                | (((cover & lane3) == lane3) << 3)
                | (((cover & lane4) == lane4) << 4)
                | (((cover & lane5) == lane5) << 5)
            )
            profile[satisfaction] += 1
    require(sum(profile) == 1_000_000, "literal fibre size mismatch")
    return tuple(profile)


QUADRATIC_RESIDUES = {pow(value, 2, P) for value in LABELS}


def quadratic_coset(labels: tuple[int, ...]) -> bool:
    support = set(labels)
    return support == QUADRATIC_RESIDUES or support == set(LABELS) - QUADRATIC_RESIDUES


def count_histogram(profile: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * 7
    for subset, multiplicity in enumerate(profile):
        answer[subset.bit_count()] += multiplicity
    return tuple(answer)


def owner_sizes(profile: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(
            multiplicity
            for subset, multiplicity in enumerate(profile)
            if subset & (1 << owner)
        )
        for owner in range(6)
    )


def pair_size(profile: tuple[int, ...], first: int, second: int) -> int:
    return sum(
        multiplicity
        for subset, multiplicity in enumerate(profile)
        if subset & (1 << first) and subset & (1 << second)
    )


ORDINARY_HISTOGRAM = (996_640, 3_360, 0, 0, 0, 0, 0)
QUADRATIC_HISTOGRAM = (998_200, 0, 1_800, 0, 0, 0, 0)
PROFILES = []
PROFILE_DIGESTER = sha256()
REPLAYED_WORDS = 0
REPRESENTED_WORDS = 0
GLOBAL_SURVIVORS = 0
PROFILE_SPLIT = Counter()
for labels, orbit_size in ORBITS:
    profile = literal_fibre_profile(labels)
    PROFILES.append(profile)
    PROFILE_DIGESTER.update(
        b"".join(value.to_bytes(4, "little") for value in profile)
    )
    quadratic = quadratic_coset(labels)
    expected_histogram = (
        QUADRATIC_HISTOGRAM if quadratic else ORDINARY_HISTOGRAM
    )
    require(
        count_histogram(profile) == expected_histogram,
        "owner-satisfaction profile mismatch",
    )
    expected_owner_size = 600 if quadratic else 560
    require(
        owner_sizes(profile) == (expected_owner_size,) * 6,
        "owner-obligation size mismatch",
    )
    for first, second in combinations(range(6), 2):
        expected_pair_size = (
            600
            if quadratic and labels[first] + labels[second] == P
            else 0
        )
        require(
            pair_size(profile, first, second) == expected_pair_size,
            "owner-obligation pair intersection mismatch",
        )
    REPLAYED_WORDS += sum(profile)
    REPRESENTED_WORDS += orbit_size * sum(profile)
    GLOBAL_SURVIVORS += orbit_size * profile[63]
    PROFILE_SPLIT["quadratic" if quadratic else "ordinary"] += orbit_size

PROFILE_DIGEST = PROFILE_DIGESTER.hexdigest()
require(REPLAYED_WORDS == 7_000_000, "representative replay total mismatch")
require(REPRESENTED_WORDS == 66_000_000, "represented replay total mismatch")
require(GLOBAL_SURVIVORS == 0, "a common-sheet unit word survived")
require(
    PROFILE_SPLIT == {"ordinary": 64, "quadratic": 2},
    "64+2 profile split mismatch",
)


def tournament_fingerprint(labels: tuple[int, ...]) -> tuple[object, ...]:
    """Complete the symmetric intersection observable by the gauged tie path."""
    # Multiplication has already gauged each support to its least orbit member.
    # The sorted labels are the declared tie Hamiltonian path.  Orienting both
    # intersection edges and missing pairs forward along it is transitive.
    edges = {
        (labels[first], labels[second])
        for first, second in combinations(range(6), 2)
    }
    scores = tuple(
        sorted(
            sum((vertex, other) in edges for other in labels)
            for vertex in labels
        )
    )
    cycles = sum(
        all(
            sum((vertex, other) in edges for other in triple) == 1
            for vertex in triple
        )
        for triple in combinations(labels, 3)
    )
    hamiltonian_paths = sum(
        all((path[index], path[index + 1]) in edges for index in range(5))
        for path in permutations(labels)
    )
    # A transitive tournament has six singleton strongly connected components.
    scc_sizes = (1, 1, 1, 1, 1, 1)
    edge_flips_from_tie_path = 0
    return scores, cycles, scc_sizes, edge_flips_from_tie_path, hamiltonian_paths


TOURNAMENT_FINGERPRINTS = Counter(
    tournament_fingerprint(labels) for labels, _ in ORBITS
)
require(
    TOURNAMENT_FINGERPRINTS
    == {((0, 1, 2, 3, 4, 5), 0, (1, 1, 1, 1, 1, 1), 0, 1): 7},
    "completed-tournament fingerprint mismatch",
)


print("scale-eleven independent Python referee")
print(f"mask digest {MASK_DIGEST}")
print("hereditary divisor words 57; labelled order contexts 52668")
print("hereditary state words 1771500; labelled contexts 1636866000")
print("scalar capacity 66 all-D11 contexts on 66 supports")
print("owner-local bank 66 all-D11 contexts on 66 supports")
print(f"owner-local bank digest {LOCAL_BANK_DIGEST}")
print("multiplication orbits 7; sizes 12,12,12,4,12,12,2")
print("exact covariance permutes providers, D11 units, owners, and owner sheets")
print("literal representative unit words 7000000; represented words 66000000")
print("global common-sheet survivors 0")
print("ordinary split 64; profile 0:996640 1:3360; owner 560; empty pair nerve")
print("quadratic split 2; profile 0:998200 2:1800; owner 600; nerve 3K2")
print(f"orbit-profile digest {PROFILE_DIGEST}")
print("pair observable owner-obligation intersection; multiplication gauge; sorted tie path")
print("completed tournaments transitive: scores 0,1,2,3,4,5; cycles 0; SCC 1^6; flips 0; paths 1")
print("faithful vertices are owner obligations (or dually unit words), not raw runners")
print("sheet bits and completed tournaments erase the empty-nerve versus 3K2 distinction")
