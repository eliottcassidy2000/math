#!/usr/bin/env python3
"""Exact c=7 AP-centred Hamming-six common-sheet obstruction (THM-962).

The raw leave-one-out-compatible bank has 108,673,488 label/state contexts.
This verifier does not iterate that bank literally.  It proves the exact local
mask grammar, uses sheet capacity to force all six effective orders to seven,
classifies all 924 six-label supports by a unit-independent row-product
character, and dispatches the only two exceptional supports by a square-sum
identity.  A direct 93,312-word literal replay independently checks the two
exceptional unit fibres.

Tournament telemetry is included deliberately as a lossy quotient.  Its
pairwise observable is the relative lambda character, its binary switch is a
centred discrete logarithm in F_7^*, and ties are completed along the
increasing-label Hamiltonian path.  The two hard supports collapse to the
fully tied, transitive completion, demonstrating why the six-fold row product
and power sums—not a bare tournament—are the faithful carrier at scale seven.
"""

from __future__ import annotations

from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations, product
from operator import mul


P = 13
SCALE = 7
LABELS = tuple(range(1, P))
UNITS = tuple(range(1, SCALE))
STATES = ((1, 0),) + tuple((SCALE, unit) for unit in UNITS)
FULL_SHEET_MASK = (1 << SCALE) - 1
SQUARES = frozenset(pow(unit, 2, SCALE) for unit in UNITS)

# lambda(a) = ceil(a/2)^(-1) in F_7.  The paired display is the exact
# owner-normalized local grammar derived below from literal CRT masks.
LAMBDA_BY_RATIO = (
    0,
    1, 1,
    4, 4,
    5, 5,
    2, 2,
    3, 3,
    6, 6,
)
LOG_BASE_3 = {1: 0, 3: 1, 2: 2, 6: 3, 4: 4, 5: 5}

QUADRATIC_CYCLE = (1, 4, 3, 12, 9, 10)  # powers of 4 modulo 13
QUADRATIC_RESIDUES = tuple(sorted(QUADRATIC_CYCLE))
NONSQUARE_COSET = tuple(sorted(2 * value % P for value in QUADRATIC_CYCLE))
EXCEPTIONAL_SUPPORTS = (QUADRATIC_RESIDUES, NONSQUARE_COSET)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, order: int, unit: int) -> int:
    return next(
        value
        for value in range(P * order)
        if value % P == order * label % P
        and value % order == unit % order
    )


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


def sheet_mask(label: int, state: int, owner: int) -> int:
    order, unit = STATES[state]
    base = crt_base(label, order, unit)
    inverse_owner = pow(owner, -1, P)
    return sum(
        1 << sheet
        for sheet in range(SCALE)
        if -order
        < centered(base * (inverse_owner + P * sheet), P * order)
        <= order
    )


MASK = {
    (label, state, owner): sheet_mask(label, state, owner)
    for label in LABELS
    for state in range(len(STATES))
    for owner in LABELS
}


# Exact local grammar.  A D=1 state fills all sheets at its own owner and none
# elsewhere.  At a D=7 owner, all six self masks share one fixed sheet; their
# extra sheets label F_7^*.  In that owner-dependent gauge an off-owner state
# (7,e) contributes the singleton lambda(provider/owner)*e.
OWNER_GAUGES: dict[int, tuple[int, dict[int, int]]] = {}
for owner in LABELS:
    require(MASK[owner, 0, owner] == FULL_SHEET_MASK, "D=1 self mask is not full")
    require(
        all(MASK[provider, 0, owner] == 0 for provider in LABELS if provider != owner),
        "off-owner D=1 state supplies a sheet",
    )

    self_masks = tuple(MASK[owner, unit, owner] for unit in UNITS)
    fixed_mask = reduce(int.__and__, self_masks)
    require(fixed_mask.bit_count() == 1, "D=7 self masks lack a unique fixed sheet")
    fixed_sheet = fixed_mask.bit_length() - 1
    symbol_of_sheet: dict[int, int] = {}
    for unit, mask in zip(UNITS, self_masks):
        require(mask.bit_count() == 2, "D=7 self mask does not have size two")
        extra_mask = mask ^ fixed_mask
        require(extra_mask.bit_count() == 1, "D=7 self extra is not a singleton")
        symbol_of_sheet[extra_mask.bit_length() - 1] = unit
    require(
        len(symbol_of_sheet) == 6 and fixed_sheet not in symbol_of_sheet,
        "self extras do not gauge the six nonfixed sheets",
    )

    for provider in LABELS:
        if provider == owner:
            continue
        character = LAMBDA_BY_RATIO[ratio(provider, owner)]
        for unit in UNITS:
            mask = MASK[provider, unit, owner]
            require(mask.bit_count() == 1, "off-owner D=7 mask is not a singleton")
            sheet = mask.bit_length() - 1
            require(sheet != fixed_sheet, "off-owner D=7 mask hits the fixed sheet")
            require(
                symbol_of_sheet[sheet] == character * unit % SCALE,
                "literal D=7 mask disagrees with the lambda grammar",
            )
    OWNER_GAUGES[owner] = (fixed_sheet, symbol_of_sheet)


literal_mask_payload = bytes(
    MASK[label, state, owner]
    for label in LABELS
    for state in range(len(STATES))
    for owner in LABELS
)


# THM-860's hereditary lcm law is equivalent, at prime scale seven, to at
# least two coordinates having effective order seven.  The exact strata are
# small enough to enumerate directly.  If k is their number, any D=7 owner
# sees at most 2+(k-1)=k+1 sheets, so every stratum k<=5 dies by capacity.
STATE_WORDS_BY_D7 = Counter()
for word in product(range(len(STATES)), repeat=6):
    d7_count = sum(state != 0 for state in word)
    if d7_count >= 2:
        STATE_WORDS_BY_D7[d7_count] += 1
require(
    STATE_WORDS_BY_D7
    == {2: 540, 3: 4_320, 4: 19_440, 5: 46_656, 6: 46_656},
    "unexpected leave-one-out state strata",
)


def owner_symbols(
    labels: tuple[int, ...], units: tuple[int, ...], owner: int
) -> tuple[int, ...]:
    return tuple(
        LAMBDA_BY_RATIO[ratio(provider, owner)] * unit % SCALE
        for provider, unit in zip(labels, units)
    )


def symbolic_owner_cover(
    labels: tuple[int, ...], units: tuple[int, ...], owner: int
) -> bool:
    return set(owner_symbols(labels, units, owner)) == set(UNITS)


def literal_owner_cover(
    labels: tuple[int, ...], units: tuple[int, ...], owner: int
) -> bool:
    union = 0
    for provider, unit in zip(labels, units):
        union |= MASK[provider, unit, owner]
    return union == FULL_SHEET_MASK


def row_product(labels: tuple[int, ...], owner: int) -> int:
    return reduce(
        mul,
        (LAMBDA_BY_RATIO[ratio(provider, owner)] for provider in labels),
        1,
    ) % SCALE


# If all six owner rows are permutations of F_7^*, then their products all
# equal -1.  Cancelling the common product of the six units forces the six
# row_product values to agree.  Exhaust the 924 supports in this quotient.
PRODUCT_MULTIPLICITY_MAX = Counter()
CONSTANT_PRODUCT_SUPPORTS = []
SUPPORT_PAYLOAD = []
for labels in combinations(LABELS, 6):
    signature = tuple(row_product(labels, owner) for owner in labels)
    multiplicities = Counter(signature)
    PRODUCT_MULTIPLICITY_MAX[max(multiplicities.values())] += 1
    if len(multiplicities) == 1:
        CONSTANT_PRODUCT_SUPPORTS.append(labels)
    SUPPORT_PAYLOAD.append(
        f"{','.join(map(str, labels))}|{','.join(map(str, signature))}"
    )
require(
    PRODUCT_MULTIPLICITY_MAX == {2: 438, 3: 328, 4: 132, 5: 24, 6: 2},
    "unexpected product-signature multiplicities",
)
require(
    tuple(CONSTANT_PRODUCT_SUPPORTS) == EXCEPTIONAL_SUPPORTS,
    "constant row product has supports beyond the two quadratic cosets",
)


def signed_doubling_cycle(labels: tuple[int, ...]) -> bool:
    edges = {
        (provider, owner)
        for provider in labels
        for owner in labels
        if provider != owner and ratio(owner, provider) in (2, 11)
    }
    return (
        len(edges) == 6
        and all(sum((provider, owner) in edges for provider in labels) == 1 for owner in labels)
        and all(sum((provider, owner) in edges for owner in labels) == 1 for provider in labels)
    )


SIGNED_DOUBLING_SUPPORTS = tuple(
    labels for labels in combinations(LABELS, 6) if signed_doubling_cycle(labels)
)
require(len(SIGNED_DOUBLING_SUPPORTS) == 64, "unexpected signed-doubling support count")
require(
    not set(SIGNED_DOUBLING_SUPPORTS) & set(CONSTANT_PRODUCT_SUPPORTS),
    "a recurring signed-doubling support survives the scale-seven product character",
)


# On either exceptional coset, order the labels by powers of 4.  The squared
# character row is [1,2,2,1,2,2].  Thus, with z_i=e_i^2 and S=sum z_i, the
# owner-i square-sum equation is 2S-z_i-z_(i+3)=0.  Summing all six equations
# gives 3S=0, hence z_(i+3)=-z_i.  This is impossible because -1 is a
# nonsquare modulo seven.  The 3^6 square-word census below is a tiny exact
# executable check of the same contradiction.
SQUARED_CHARACTER_ROW = tuple(
    pow(LAMBDA_BY_RATIO[value], 2, SCALE) for value in QUADRATIC_CYCLE
)
require(
    SQUARED_CHARACTER_ROW == (1, 2, 2, 1, 2, 2),
    "unexpected quadratic-coset squared character row",
)
require(
    SQUARES == {1, 2, 4}
    and SQUARES.isdisjoint({(-value) % SCALE for value in SQUARES}),
    "the mod-seven square/nonsquare contradiction failed",
)


def square_sum_equations(z: tuple[int, ...]) -> tuple[int, ...]:
    total = sum(z) % SCALE
    direct = tuple(
        sum(SQUARED_CHARACTER_ROW[(j - i) % 6] * z[j] for j in range(6))
        % SCALE
        for i in range(6)
    )
    reduced = tuple((2 * total - z[i] - z[(i + 3) % 6]) % SCALE for i in range(6))
    require(direct == reduced, "square-sum circulant does not reduce as claimed")
    require(sum(direct) % SCALE == 3 * total % SCALE, "summed square law mismatch")
    return direct


SQUARE_WORDS_TESTED = 0
for z in product(sorted(SQUARES), repeat=6):
    SQUARE_WORDS_TESTED += 1
    require(any(square_sum_equations(z)), "an exceptional square word survived")
require(SQUARE_WORDS_TESTED == 3**6, "unexpected exceptional square-word census")


# Independent literal replay of both exceptional supports.  The profiles are
# not needed by the proof; they ensure the owner gauge and square obstruction
# agree with the raw seven-sheet masks on every exceptional unit word.
EXCEPTIONAL_PROFILE_VARIANTS = Counter()
for labels in EXCEPTIONAL_SUPPORTS:
    profile = Counter()
    for units in product(UNITS, repeat=6):
        satisfied = 0
        for owner in labels:
            symbolic = symbolic_owner_cover(labels, units, owner)
            literal = literal_owner_cover(labels, units, owner)
            require(symbolic == literal, "symbolic/literal exceptional owner mismatch")
            satisfied += literal
        profile[satisfied] += 1
    require(
        profile == {0: 44_712, 2: 1_728, 4: 216},
        "unexpected exceptional owner-satisfaction profile",
    )
    EXCEPTIONAL_PROFILE_VARIANTS[tuple(sorted(profile.items()))] += 1


# Tournament telemetry.  For an unordered pair a<b use
#   kappa(a,b)=lambda(b/a)/lambda(a/b) in F_7^*.
# Write kappa=3^d.  Exponents d=1,2 point a->b; d=4,5 point b->a; the
# self-inverse values d=0,3 are ties.  Ties follow the increasing-label
# Hamiltonian path.  This completion is intentionally not used in the proof.
def tournament_edge(a: int, b: int) -> tuple[tuple[int, int], bool]:
    require(a < b, "tournament edge expects increasing labels")
    forward = LAMBDA_BY_RATIO[ratio(b, a)]
    backward = LAMBDA_BY_RATIO[ratio(a, b)]
    observable = forward * pow(backward, -1, SCALE) % SCALE
    exponent = LOG_BASE_3[observable]
    if exponent in (1, 2):
        return (a, b), False
    if exponent in (4, 5):
        return (b, a), False
    return (a, b), True


def directed_reachable(
    source: int, target: int, labels: tuple[int, ...], edges: set[tuple[int, int]]
) -> bool:
    if source == target:
        return True
    stack = [source]
    seen = {source}
    while stack:
        vertex = stack.pop()
        for neighbor in labels:
            if (vertex, neighbor) not in edges or neighbor in seen:
                continue
            if neighbor == target:
                return True
            seen.add(neighbor)
            stack.append(neighbor)
    return False


def scc_sizes(labels: tuple[int, ...], edges: set[tuple[int, int]]) -> tuple[int, ...]:
    unused = set(labels)
    sizes = []
    while unused:
        seed = min(unused)
        component = {
            vertex
            for vertex in labels
            if directed_reachable(seed, vertex, labels, edges)
            and directed_reachable(vertex, seed, labels, edges)
        }
        unused -= component
        sizes.append(len(component))
    return tuple(sorted(sizes, reverse=True))


def tournament_fingerprint(labels: tuple[int, ...]):
    edges: set[tuple[int, int]] = set()
    ties = 0
    backward_edges = 0
    for a, b in combinations(labels, 2):
        edge, tie = tournament_edge(a, b)
        edges.add(edge)
        ties += tie
        backward_edges += edge == (b, a)
    scores = tuple(
        sorted(sum((vertex, other) in edges for other in labels) for vertex in labels)
    )
    triangles = sum(
        all(sum((vertex, other) in edges for other in triple) == 1 for vertex in triple)
        for triple in combinations(labels, 3)
    )
    hamiltonian_paths = sum(
        all((path[i], path[i + 1]) in edges for i in range(5))
        for path in permutations(labels)
    )
    return (
        scores,
        triangles,
        scc_sizes(labels, edges),
        ties,
        backward_edges,
        hamiltonian_paths,
    )


TOURNAMENT_FINGERPRINTS = Counter(
    tournament_fingerprint(labels) for labels in combinations(LABELS, 6)
)
TOURNAMENT_SCORE_HIST = Counter()
TOURNAMENT_TRIANGLE_HIST = Counter()
TOURNAMENT_SCC_HIST = Counter()
TOURNAMENT_TIE_HIST = Counter()
TOURNAMENT_BACKWARD_HIST = Counter()
TOURNAMENT_HAMILTONIAN_HIST = Counter()
for fingerprint, multiplicity in TOURNAMENT_FINGERPRINTS.items():
    score, triangles, scc, ties, backward, paths = fingerprint
    TOURNAMENT_SCORE_HIST[score] += multiplicity
    TOURNAMENT_TRIANGLE_HIST[triangles] += multiplicity
    TOURNAMENT_SCC_HIST[scc] += multiplicity
    TOURNAMENT_TIE_HIST[ties] += multiplicity
    TOURNAMENT_BACKWARD_HIST[backward] += multiplicity
    TOURNAMENT_HAMILTONIAN_HIST[paths] += multiplicity

EXCEPTIONAL_TOURNAMENTS = tuple(
    tournament_fingerprint(labels) for labels in EXCEPTIONAL_SUPPORTS
)
require(
    len(TOURNAMENT_FINGERPRINTS) == 245,
    "unexpected scale-seven tournament fingerprint count",
)
require(
    len(set(EXCEPTIONAL_TOURNAMENTS)) == 1
    and EXCEPTIONAL_TOURNAMENTS[0]
    == ((0, 1, 2, 3, 4, 5), 0, (1, 1, 1, 1, 1, 1), 15, 0, 1),
    "quadratic cosets do not collapse to the fully tied transitive tournament",
)


STATE_CONTEXTS_BY_D7 = {
    count: 924 * words for count, words in sorted(STATE_WORDS_BY_D7.items())
}
MIXED_CONTEXTS_KILLED_BY_CAPACITY = sum(
    contexts for count, contexts in STATE_CONTEXTS_BY_D7.items() if count < 6
)
ALL_D7_CONTEXTS = STATE_CONTEXTS_BY_D7[6]
PRODUCT_CONTEXTS_KILLED = (924 - len(EXCEPTIONAL_SUPPORTS)) * 6**6
EXCEPTIONAL_CONTEXTS_KILLED = len(EXCEPTIONAL_SUPPORTS) * 6**6

require(
    sum(STATE_CONTEXTS_BY_D7.values()) == 108_673_488,
    "unexpected full scale-seven context count",
)
require(
    MIXED_CONTEXTS_KILLED_BY_CAPACITY == 65_563_344,
    "unexpected mixed-order capacity count",
)
require(
    ALL_D7_CONTEXTS == 43_110_144
    and PRODUCT_CONTEXTS_KILLED == 43_016_832
    and EXCEPTIONAL_CONTEXTS_KILLED == 93_312,
    "unexpected all-order-seven decomposition",
)


def sorted_counter(counter: Counter) -> dict:
    return dict(sorted(counter.items(), key=lambda item: item[0]))


print("THM-962 SCALE-SEVEN HAMMING-SIX PRODUCT/SQUARE OBSTRUCTION")
print(f"states={STATES}")
print(f"state_words_by_D7={sorted_counter(STATE_WORDS_BY_D7)}")
print(f"state_contexts_by_D7={STATE_CONTEXTS_BY_D7}")
print(f"all_admissible_contexts={sum(STATE_CONTEXTS_BY_D7.values())}")
print(f"mixed_contexts_killed_by_capacity={MIXED_CONTEXTS_KILLED_BY_CAPACITY}")
print(f"lambda_ratio_pairs={(1, 4, 5, 2, 3, 6)}")
print(f"literal_mask_payload_sha256={sha256(literal_mask_payload).hexdigest()}")
print(f"all_D7_supports={924}")
print(f"all_D7_unit_contexts={ALL_D7_CONTEXTS}")
print(f"product_signature_max_multiplicity={sorted_counter(PRODUCT_MULTIPLICITY_MAX)}")
print(f"constant_product_supports={tuple(CONSTANT_PRODUCT_SUPPORTS)}")
print(f"signed_doubling_supports_product_killed={len(SIGNED_DOUBLING_SUPPORTS)}")
print(f"product_contexts_killed={PRODUCT_CONTEXTS_KILLED}")
print(f"support_product_payload_sha256={sha256(chr(10).join(SUPPORT_PAYLOAD).encode()).hexdigest()}")
print(f"squared_character_row={SQUARED_CHARACTER_ROW}")
print(f"square_words_tested_per_coset={SQUARE_WORDS_TESTED}")
print(f"exceptional_contexts_literal_replayed={EXCEPTIONAL_CONTEXTS_KILLED}")
print(f"exceptional_profile_variants={dict(EXCEPTIONAL_PROFILE_VARIANTS)}")
print("full_common_sheet_contexts=0")
print("tournament_pair_observable=lambda(b/a)/lambda(a/b)")
print("tournament_switch=log_3_exponents_1,2_forward_4,5_backward_0,3_tied")
print("tournament_tie_Hamiltonian_path=increasing_labels")
print(f"tournament_joint_fingerprints={len(TOURNAMENT_FINGERPRINTS)}")
print(f"tournament_score_hist={sorted_counter(TOURNAMENT_SCORE_HIST)}")
print(f"tournament_directed_triangle_hist={sorted_counter(TOURNAMENT_TRIANGLE_HIST)}")
print(f"tournament_SCC_hist={sorted_counter(TOURNAMENT_SCC_HIST)}")
print(f"tournament_tie_edge_hist={sorted_counter(TOURNAMENT_TIE_HIST)}")
print(f"tournament_lex_backward_edge_hist={sorted_counter(TOURNAMENT_BACKWARD_HIST)}")
print(f"tournament_Hamilton_path_hist={sorted_counter(TOURNAMENT_HAMILTONIAN_HIST)}")
print(f"exceptional_tournament_fingerprint={EXCEPTIONAL_TOURNAMENTS[0]}")
print("PASS: every primitive proper AP-centred scale-seven H6 common-sheet context is impossible")
