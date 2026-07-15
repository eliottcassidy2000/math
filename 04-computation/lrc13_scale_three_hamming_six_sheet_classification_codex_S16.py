#!/usr/bin/env python3
"""Exact c=3 primitive AP-centred H6 common-sheet classification.

This verifier deliberately stops before the unbounded metric component
recursion.  It proves the finite sheet bank, its signed-provider/matching-code
description, the exact symmetry orbits, a metric-quotient counterexample, and
the exact first layer of the proposed THM-857 recursion.

Only integer arithmetic and fractions.Fraction are used.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path


P = 13
SCALE = 3
FULL_MASK = (1 << SCALE) - 1
UNITS = (1, 2)  # +1,-1 in F_3
MULTIPLIERS = tuple(range(1, P))
EDGE_SIGN = {4: -1, 5: -1, 8: 1, 9: 1}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, order: int, unit: int) -> int:
    return next(
        u
        for u in range(P * order)
        if u % P == order * label % P and u % order == unit % order
    )


def sheet_mask(label: int, order: int, unit: int, owner: int) -> int:
    u = crt_base(label, order, unit)
    inverse_owner = pow(owner, -1, P)
    answer = 0
    for sheet in range(SCALE):
        z = centered(u * (inverse_owner + P * sheet), P * order)
        if -order < z <= order:
            answer |= 1 << sheet
    return answer


MASK = {
    (label, order, unit, owner): sheet_mask(label, order, unit, owner)
    for label in range(1, P)
    for order in (1, 3)
    for unit in ((0,) if order == 1 else UNITS)
    for owner in range(1, P)
}


def direct_cover(order_three: tuple[int, ...], units: tuple[int, ...]) -> bool:
    """Test the nonautomatic D=3 owners with literal CRT sheet masks."""
    for owner in order_three:
        mask = reduce(
            int.__or__,
            (
                MASK[label, 3, unit, owner]
                for label, unit in zip(order_three, units)
            ),
            0,
        )
        if mask != FULL_MASK:
            return False
    return True


def nae_cover(order_three: tuple[int, ...], units: tuple[int, ...]) -> bool:
    """Equivalent signed-provider NAE test after the self sheet is removed."""
    signs = {
        label: (1 if unit == 1 else -1)
        for label, unit in zip(order_three, units)
    }
    for owner in order_three:
        literals = []
        for provider in order_three:
            if provider == owner:
                continue
            ratio = provider * pow(owner, -1, P) % P
            if ratio in EDGE_SIGN:
                literals.append(EDGE_SIGN[ratio] * signs[provider])
        if 1 not in literals or -1 not in literals:
            return False
    return True


# The D=1 owner covers all three sheets at itself and no sheet at another
# replacement owner.  Hence the common-sheet condition depends only on C,
# the set of D=3 labels, and its unit word.  Check the direct/symbolic identity
# before decorating C by the D=1 labels.
VALID_C_UNITS: dict[tuple[int, ...], tuple[tuple[int, ...], ...]] = {}
DIRECT_NAE_CHECKS = 0
for size in range(1, 7):
    for order_three in combinations(range(1, P), size):
        valid = []
        for units in product(UNITS, repeat=size):
            direct = direct_cover(order_three, units)
            symbolic = nae_cover(order_three, units)
            require(direct == symbolic, "direct mask / signed NAE mismatch")
            DIRECT_NAE_CHECKS += 1
            if direct:
                valid.append(units)
        if valid:
            VALID_C_UNITS[order_three] = tuple(valid)


def make_context(
    order_three: tuple[int, ...],
    order_one: tuple[int, ...],
    units: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    decorated = {label: (3, unit) for label, unit in zip(order_three, units)}
    decorated.update({label: (1, 0) for label in order_one})
    labels = tuple(sorted(decorated))
    return (
        labels,
        tuple(decorated[label][0] for label in labels),
        tuple(decorated[label][1] for label in labels),
    )


CONTEXTS = set()
PRESENTATIONS = set()
for order_three, unit_words in VALID_C_UNITS.items():
    outside = tuple(label for label in range(1, P) if label not in order_three)
    for order_one in combinations(outside, 6 - len(order_three)):
        for units in unit_words:
            row = make_context(order_three, order_one, units)
            CONTEXTS.add(row)
            PRESENTATIONS.add(row[:2])


def order_type(row) -> tuple[int, int]:
    return row[1].count(1), row[1].count(3)


PRESENTATION_COUNTS = Counter(order_type(row) for row in PRESENTATIONS)
CONTEXT_COUNTS = Counter(order_type(row) for row in CONTEXTS)
UNITS_PER_PRESENTATION = Counter(
    (order_type(presentation), sum(row[:2] == presentation for row in CONTEXTS))
    for presentation in PRESENTATIONS
)

require(DIRECT_NAE_CHECKS == 94_448, "unexpected direct/NAE check count")
require(len(PRESENTATIONS) == 212, "unexpected presentation count")
require(len(CONTEXTS) == 1_504, "unexpected context count")
require(
    PRESENTATION_COUNTS == {(2, 4): 84, (1, 5): 84, (0, 6): 44},
    "unexpected presentation split",
)
require(
    CONTEXT_COUNTS == {(2, 4): 336, (1, 5): 672, (0, 6): 496},
    "unexpected context split",
)
require(
    UNITS_PER_PRESENTATION
    == {
        ((2, 4), 4): 84,
        ((1, 5), 8): 84,
        ((0, 6), 8): 26,
        ((0, 6), 16): 18,
    },
    "unexpected unit-fibre distribution",
)


# Structural descriptions of the two mixed banks.
H = frozenset((1, 5, 8, 12))
H_COSETS = {frozenset(a * x % P for x in H) for a in MULTIPLIERS}
EXPECTED_K4 = set()
for quartet in H_COSETS:
    outside = tuple(label for label in range(1, P) if label not in quartet)
    for pair in combinations(outside, 2):
        labels = tuple(sorted(quartet | set(pair)))
        orders = tuple(3 if label in quartet else 1 for label in labels)
        EXPECTED_K4.add((labels, orders))

EXPECTED_K5 = set()
for quartet in H_COSETS:
    forward = {2 * x % P for x in quartet}
    for flag in forward:
        order_three = quartet | {flag}
        for order_one in range(1, P):
            if order_one in order_three:
                continue
            labels = tuple(sorted(order_three | {order_one}))
            orders = tuple(3 if label in order_three else 1 for label in labels)
            EXPECTED_K5.add((labels, orders))

require(len(EXPECTED_K4) == 84, "k=4 structural count mismatch")
require(len(EXPECTED_K5) == 84, "k=5 structural count mismatch")
require(
    EXPECTED_K4 == {row for row in PRESENTATIONS if order_type(row) == (2, 4)},
    "k=4 structural classification mismatch",
)
require(
    EXPECTED_K5 == {row for row in PRESENTATIONS if order_type(row) == (1, 5)},
    "k=5 structural classification mismatch",
)


MATCHING_CODES = {
    (1, 5, 8, 12): (((1, 12, 1), (5, 8, 1)), ()),
    (1, 2, 5, 8, 12): (((1, 12, 1), (5, 8, 1)), (2,)),
    (1, 2, 3, 4, 9, 10): (((1, 2, 1), (4, 9, 1), (3, 10, 1)), ()),
    (1, 2, 3, 5, 8, 12): (((5, 8, 1), (1, 12, 1)), (2, 3)),
    (1, 2, 5, 6, 8, 11): (((1, 6, -1), (5, 8, 1), (2, 11, 1)), ()),
    (1, 2, 5, 8, 11, 12): (((5, 8, 1), (1, 12, 1)), (2, 11)),
    (1, 3, 4, 9, 10, 12): (((4, 9, 1), (3, 10, 1), (1, 12, 1)), ()),
}

for order_three, (equations, _) in MATCHING_CODES.items():
    predicted = set()
    for units in product(UNITS, repeat=len(order_three)):
        signs = {
            label: (1 if unit == 1 else -1)
            for label, unit in zip(order_three, units)
        }
        if all(signs[a] * signs[b] == rhs for a, b, rhs in equations):
            predicted.add(units)
    require(
        predicted == set(VALID_C_UNITS[order_three]),
        f"matching-code mismatch at {order_three}",
    )


def multiply_context(row, multiplier: int):
    labels, orders, units = row
    triples = sorted(
        (multiplier * label % P, order, unit)
        for label, order, unit in zip(labels, orders, units)
    )
    return (
        tuple(x[0] for x in triples),
        tuple(x[1] for x in triples),
        tuple(x[2] for x in triples),
    )


def multiply_presentation(row, multiplier: int):
    labels, orders = row
    pairs = sorted(
        (multiplier * label % P, order) for label, order in zip(labels, orders)
    )
    return tuple(x[0] for x in pairs), tuple(x[1] for x in pairs)


def reflect_context(row):
    labels, orders, units = row
    return (
        labels,
        orders,
        tuple(3 - unit if order == 3 else 0 for order, unit in zip(orders, units)),
    )


def invert_context(row):
    labels, orders, units = row
    triples = sorted(
        (pow(label, -1, P), order, unit)
        for label, order, unit in zip(labels, orders, units)
    )
    return (
        tuple(x[0] for x in triples),
        tuple(x[1] for x in triples),
        tuple(x[2] for x in triples),
    )


def invert_presentation(row):
    labels, orders = row
    pairs = sorted(
        (pow(label, -1, P), order) for label, order in zip(labels, orders)
    )
    return tuple(x[0] for x in pairs), tuple(x[1] for x in pairs)


def orbit_partition(space, group, action):
    unseen = set(space)
    answer = []
    while unseen:
        seed = min(unseen)
        orbit = {action(seed, group_element) for group_element in group}
        require(orbit <= space, "action leaves the claimed space")
        representative = min(orbit)
        stabilizer = tuple(
            group_element
            for group_element in group
            if action(representative, group_element) == representative
        )
        answer.append((representative, len(orbit), stabilizer))
        unseen -= orbit
    return tuple(sorted(answer))


PRESENTATION_ORBITS = orbit_partition(
    PRESENTATIONS, MULTIPLIERS, multiply_presentation
)
CONTEXT_ORBITS = orbit_partition(CONTEXTS, MULTIPLIERS, multiply_context)
SHEET_GROUP = tuple((multiplier, reflection) for multiplier in MULTIPLIERS for reflection in (0, 1))


def sheet_action(row, group_element):
    multiplier, reflection = group_element
    image = multiply_context(row, multiplier)
    return reflect_context(image) if reflection else image


SHEET_ORBITS = orbit_partition(CONTEXTS, SHEET_GROUP, sheet_action)


# One explicit representative for each C_12 x C_2 sheet orbit.  The order of
# each word agrees with the displayed increasing D=3 label set.  These are
# normalized first by the D=3 set and then by its stabilizer and global flip;
# they need not be the globally lexicographically least context tuple.
CANONICAL_SHEET_REPS = []
for order_one in (
    (2, 3),
    (2, 4),
    (2, 6),
    (2, 7),
    (2, 9),
    (2, 11),
    (4, 6),
    (4, 9),
):
    for units in ((1, 1, 1, 1), (1, 2, 2, 1)):
        CANONICAL_SHEET_REPS.append(make_context((1, 5, 8, 12), order_one, units))
for order_one in ((3,), (4,), (6,), (7,), (9,), (10,), (11,)):
    for units in (
        (1, 1, 1, 1, 1),
        (1, 1, 2, 2, 1),
        (1, 2, 1, 1, 1),
        (1, 2, 2, 2, 1),
    ):
        CANONICAL_SHEET_REPS.append(
            make_context((1, 2, 5, 8, 12), order_one, units)
        )
for order_three, words in (
    (
        (1, 2, 3, 4, 9, 10),
        ((1, 1, 1, 1, 1, 1), (1, 1, 1, 2, 2, 1), (1, 1, 2, 1, 1, 2), (1, 1, 2, 2, 2, 2)),
    ),
    (
        (1, 2, 3, 5, 8, 12),
        (
            (1, 1, 1, 1, 1, 1),
            (1, 1, 1, 2, 2, 1),
            (1, 1, 2, 1, 1, 1),
            (1, 1, 2, 2, 2, 1),
            (1, 2, 1, 1, 1, 1),
            (1, 2, 1, 2, 2, 1),
            (1, 2, 2, 1, 1, 1),
            (1, 2, 2, 2, 2, 1),
        ),
    ),
    (
        (1, 2, 5, 6, 8, 11),
        ((1, 1, 1, 2, 1, 1), (1, 1, 2, 2, 2, 1), (1, 2, 1, 2, 1, 2), (1, 2, 2, 2, 2, 2)),
    ),
    (
        (1, 2, 5, 8, 11, 12),
        (
            (1, 1, 1, 1, 1, 1),
            (1, 1, 1, 1, 2, 1),
            (1, 1, 2, 2, 1, 1),
            (1, 1, 2, 2, 2, 1),
            (1, 2, 1, 1, 2, 1),
            (1, 2, 2, 2, 2, 1),
        ),
    ),
    (
        (1, 3, 4, 9, 10, 12),
        ((1, 1, 1, 1, 1, 1), (1, 1, 2, 2, 1, 1)),
    ),
):
    for units in words:
        CANONICAL_SHEET_REPS.append(make_context(order_three, (), units))

SHEET_ORBIT_ID = {}
for orbit_id, (representative, _, _) in enumerate(SHEET_ORBITS):
    for group_element in SHEET_GROUP:
        SHEET_ORBIT_ID[sheet_action(representative, group_element)] = orbit_id
require(len(SHEET_ORBIT_ID) == len(CONTEXTS), "sheet orbit map is incomplete")
require(len(CANONICAL_SHEET_REPS) == 68, "canonical sheet representative count mismatch")
require(
    all(row in CONTEXTS for row in CANONICAL_SHEET_REPS),
    "canonical sheet representative is not a valid context",
)
require(
    len({SHEET_ORBIT_ID[row] for row in CANONICAL_SHEET_REPS}) == 68,
    "canonical sheet representatives do not meet every orbit once",
)


def orbit_size_summary(orbits):
    return {
        kind: dict(
            sorted(
                Counter(
                    size for representative, size, _ in orbits if order_type(representative) == kind
                ).items()
            )
        )
        for kind in sorted({order_type(representative) for representative, _, _ in orbits})
    }


PRESENTATION_ORBIT_SUMMARY = orbit_size_summary(PRESENTATION_ORBITS)
CONTEXT_ORBIT_SUMMARY = orbit_size_summary(CONTEXT_ORBITS)
SHEET_ORBIT_SUMMARY = orbit_size_summary(SHEET_ORBITS)

require(len(PRESENTATION_ORBITS) == 20, "presentation orbit count mismatch")
require(len(CONTEXT_ORBITS) == 136, "context orbit count mismatch")
require(len(SHEET_ORBITS) == 68, "sheet orbit count mismatch")
require(
    PRESENTATION_ORBIT_SUMMARY
    == {
        (0, 6): {2: 1, 6: 1, 12: 3},
        (1, 5): {12: 7},
        (2, 4): {6: 2, 12: 6},
    },
    "presentation orbit-size mismatch",
)
require(
    CONTEXT_ORBIT_SUMMARY
    == {
        (0, 6): {2: 2, 6: 10, 12: 36},
        (1, 5): {12: 56},
        (2, 4): {6: 8, 12: 24},
    },
    "context orbit-size mismatch",
)
require(
    SHEET_ORBIT_SUMMARY
    == {
        (0, 6): {4: 1, 12: 5, 24: 18},
        (1, 5): {24: 28},
        (2, 4): {12: 4, 24: 12},
    },
    "sheet orbit-size mismatch",
)
require(
    all(not any(reflection for _, reflection in stabilizer) for _, _, stabilizer in SHEET_ORBITS),
    "unexpected multiplier-reflection stabilizer",
)

INVERTIBLE_PRESENTATIONS = {
    row for row in PRESENTATIONS if invert_presentation(row) in PRESENTATIONS
}
INVERTIBLE_CONTEXTS = {row for row in CONTEXTS if invert_context(row) in CONTEXTS}
require(len(INVERTIBLE_PRESENTATIONS) == 86, "inversion presentation count mismatch")
require(len(INVERTIBLE_CONTEXTS) == 352, "inversion context count mismatch")
require(
    Counter(order_type(row) for row in INVERTIBLE_PRESENTATIONS)
    == {(2, 4): 84, (0, 6): 2},
    "inversion presentation split mismatch",
)
require(
    Counter(order_type(row) for row in INVERTIBLE_CONTEXTS)
    == {(2, 4): 336, (0, 6): 16},
    "inversion context split mismatch",
)


def mult_set(labels: tuple[int, ...], multiplier: int) -> tuple[int, ...]:
    return tuple(sorted(multiplier * label % P for label in labels))


def order_three_set_orbits():
    answer = {}
    for size in (4, 5, 6):
        space = {labels for labels in VALID_C_UNITS if len(labels) == size}
        unseen = set(space)
        rows = []
        while unseen:
            seed = min(unseen)
            orbit = {mult_set(seed, multiplier) for multiplier in MULTIPLIERS}
            representative = min(orbit)
            stabilizer = tuple(
                multiplier
                for multiplier in MULTIPLIERS
                if mult_set(representative, multiplier) == representative
            )
            rows.append(
                (
                    representative,
                    len(orbit),
                    stabilizer,
                    len(VALID_C_UNITS[representative]),
                )
            )
            unseen -= orbit
        answer[size] = tuple(sorted(rows))
    return answer


C_SET_ORBITS = order_three_set_orbits()
require([len(C_SET_ORBITS[size]) for size in (4, 5, 6)] == [1, 1, 5], "C-set orbit mismatch")


def graph_fingerprint(labels: tuple[int, ...]):
    edges = {
        (provider, owner): EDGE_SIGN[provider * pow(owner, -1, P) % P]
        for provider in labels
        for owner in labels
        if provider != owner and provider * pow(owner, -1, P) % P in EDGE_SIGN
    }
    indegrees = Counter(owner for provider, owner in edges)
    outdegrees = Counter(provider for provider, owner in edges)
    pair_counts = Counter(
        int((a, b) in edges) + int((b, a) in edges)
        for a, b in combinations(labels, 2)
    )

    def reachable(start: int) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            here = stack.pop()
            for a, b in edges:
                if a == here and b not in seen:
                    seen.add(b)
                    stack.append(b)
        return seen

    reach = {label: reachable(label) for label in labels}
    remaining = set(labels)
    scc_sizes = []
    while remaining:
        seed = min(remaining)
        component = {
            label
            for label in remaining
            if label in reach[seed] and seed in reach[label]
        }
        scc_sizes.append(len(component))
        remaining -= component
    return (
        len(edges),
        tuple(sorted(Counter(edges.values()).items())),
        tuple(sorted(indegrees[label] for label in labels)),
        tuple(sorted(outdegrees[label] for label in labels)),
        tuple(pair_counts.get(number, 0) for number in (0, 1, 2)),
        tuple(sorted(scc_sizes, reverse=True)),
    )


def circle_norm(numerator: int, denominator: int) -> F:
    residue = numerator % denominator
    return F(min(residue, denominator - residue), denominator)


def exact_maximin(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    denominators = {2 * speed for speed in speeds}
    denominators |= {
        denominator
        for u, v in combinations(speeds, 2)
        for denominator in (u + v, abs(u - v))
        if denominator
    }
    best = F(0)
    witnesses = set()
    for denominator in sorted(denominators):
        for numerator in range(denominator):
            value = min(
                circle_norm(speed * numerator, denominator) for speed in speeds
            )
            witness = F(numerator, denominator)
            if value > best:
                best = value
                witnesses = {witness}
            elif value == best:
                witnesses.add(witness)
    return best, tuple(sorted(witnesses))


def least_proper_packet(row) -> tuple[int, ...]:
    labels, orders, units = row
    replacements = []
    for label, order, unit in zip(labels, orders, units):
        if order == 1:
            replacements.append(3 * label + 39)
        else:
            replacements.append(crt_base(label, 3, unit))
    core = [3 * label for label in range(1, P) if label not in labels]
    return tuple(sorted(core + replacements))


METRIC_CONTEXT = make_context(
    (1, 2, 3, 4, 9, 10), (), (1, 1, 1, 1, 1, 1)
)
METRIC_BASE = least_proper_packet(METRIC_CONTEXT)
METRIC_REFLECTION = least_proper_packet(reflect_context(METRIC_CONTEXT))
METRIC_TIMES_TWO = least_proper_packet(multiply_context(METRIC_CONTEXT, 2))
METRIC_VALUES = {
    "base": exact_maximin(METRIC_BASE),
    "flip": exact_maximin(METRIC_REFLECTION),
    "multiplier2": exact_maximin(METRIC_TIMES_TWO),
}
require(METRIC_VALUES["base"] == (F(1, 4), (F(7, 20), F(13, 20))), "base metric mismatch")
require(METRIC_VALUES["flip"] == (F(7, 26), (F(1, 52), F(51, 52))), "flip metric mismatch")
require(
    METRIC_VALUES["multiplier2"] == (F(14, 65), (F(32, 65), F(33, 65))),
    "multiplier metric mismatch",
)


def danger_intervals(speed: int) -> tuple[tuple[F, F], ...]:
    radius = F(1, P * speed)
    intervals = []
    for k in range(speed):
        lo = F(k, speed) - radius
        hi = F(k, speed) + radius
        if lo < 0:
            intervals.extend(((F(0), hi), (lo + 1, F(1))))
        elif hi > 1:
            intervals.extend(((lo, F(1)), (F(0), hi - 1)))
        else:
            intervals.append((lo, hi))
    return tuple(intervals)


DANGER_CACHE = {}


def longest_strict_safe_component(speeds: tuple[int, ...]) -> F:
    intervals = sorted(
        interval
        for speed in speeds
        for interval in DANGER_CACHE.setdefault(speed, danger_intervals(speed))
    )
    merged = []
    for lo, hi in intervals:
        if not merged or lo > merged[-1][1]:
            merged.append([lo, hi])
        elif hi > merged[-1][1]:
            merged[-1][1] = hi
    last = F(0)
    longest = F(0)
    for lo, hi in merged:
        longest = max(longest, lo - last)
        last = max(last, hi)
    return max(longest, F(1) - last)


ROOTS = {row[0] for row in CONTEXTS}
ROOTS_BY_TYPE = {
    kind: {row[0] for row in CONTEXTS if order_type(row) == kind}
    for kind in ((2, 4), (1, 5), (0, 6))
}
ROOT_LENGTH = {
    labels: longest_strict_safe_component(
        tuple(3 * label for label in range(1, P) if label not in labels)
    )
    for labels in ROOTS
}


def first_candidate_counts(row) -> tuple[int, int, F]:
    labels, orders, units = row
    length = ROOT_LENGTH[labels]
    cap = F(132, 13) / length
    cap_floor = cap.numerator // cap.denominator
    order_one = 0
    order_three = 0
    for label, order, unit in zip(labels, orders, units):
        if order == 1:
            base = 3 * label + 39
            if base <= cap:
                order_one += (cap_floor - base) // 39 + 1
        else:
            base = crt_base(label, 3, unit)
            if base <= cap:
                order_three += (cap_floor - base) // 39 + 1
    return order_one, order_three, cap


FIRST_ROWS = {row: first_candidate_counts(row) for row in CONTEXTS}
FIRST_BY_TYPE = defaultdict(Counter)
for row, (order_one, order_three, cap) in FIRST_ROWS.items():
    kind = order_type(row)
    FIRST_BY_TYPE[kind]["contexts"] += 1
    FIRST_BY_TYPE[kind]["order_one"] += order_one
    FIRST_BY_TYPE[kind]["order_three"] += order_three
    total = order_one + order_three
    if "minimum" not in FIRST_BY_TYPE[kind]:
        FIRST_BY_TYPE[kind]["minimum"] = total
        FIRST_BY_TYPE[kind]["maximum"] = total
    else:
        FIRST_BY_TYPE[kind]["minimum"] = min(
            FIRST_BY_TYPE[kind]["minimum"], total
        )
        FIRST_BY_TYPE[kind]["maximum"] = max(
            FIRST_BY_TYPE[kind]["maximum"], total
        )

FIRST_TOTAL = sum(a + b for a, b, cap in FIRST_ROWS.values())
FIRST_ORDER_ONE = sum(a for a, b, cap in FIRST_ROWS.values())
FIRST_ORDER_THREE = sum(b for a, b, cap in FIRST_ROWS.values())
ROOT_CAPS = [cap for a, b, cap in FIRST_ROWS.values()]

EXPECTED_FIRST = {
    (2, 4): {
        "contexts": 336,
        "order_one": 10_440,
        "order_three": 22_268,
        "minimum": 53,
        "maximum": 161,
    },
    (1, 5): {
        "contexts": 672,
        "order_one": 10_200,
        "order_three": 54_396,
        "minimum": 66,
        "maximum": 152,
    },
    (0, 6): {
        "contexts": 496,
        "order_one": 0,
        "order_three": 49_608,
        "minimum": 68,
        "maximum": 183,
    },
}
require(len(ROOTS) == 110, "distinct metric-root count mismatch")
require(
    {kind: len(roots) for kind, roots in ROOTS_BY_TYPE.items()}
    == {(2, 4): 84, (1, 5): 66, (0, 6): 44},
    "metric-root split mismatch",
)
require(
    {kind: dict(values) for kind, values in FIRST_BY_TYPE.items()} == EXPECTED_FIRST,
    "first-layer split mismatch",
)
require(FIRST_TOTAL == 146_912, "first-layer total mismatch")
require(FIRST_ORDER_ONE == 20_640, "first-layer D=1 mismatch")
require(FIRST_ORDER_THREE == 126_272, "first-layer D=3 mismatch")
require(min(ROOT_CAPS) == F(4_752, 13), "minimum root cap mismatch")
require(max(ROOT_CAPS) == F(1_188), "maximum root cap mismatch")

# These are planning extrapolations from THM-861's c=2 frozen census, not
# c=3 node counts.  The exact c=3 first layer is the numerator being scaled.
DEPTH_TWO_PROJECTION = F(FIRST_TOTAL * 641_866, 6_266)
FULL_TREE_PROJECTION = F(len(CONTEXTS) * 41_882_982, 64)
require(DEPTH_TWO_PROJECTION == F(47_148_908_896, 3_133), "depth-two projection mismatch")
require(FULL_TREE_PROJECTION == 984_250_077, "full-tree projection mismatch")


def presentation_line(row) -> str:
    return ",".join(f"{label}:{order}" for label, order in zip(*row))


def context_line(row) -> str:
    return ",".join(
        f"{label}:{order}:{unit}"
        for label, order, unit in zip(row[0], row[1], row[2])
    )


PRESENTATION_PAYLOAD = (
    "\n".join(presentation_line(row) for row in sorted(PRESENTATIONS)) + "\n"
).encode()
CONTEXT_PAYLOAD = (
    "\n".join(context_line(row) for row in sorted(CONTEXTS)) + "\n"
).encode()
FIRST_PAYLOAD = (
    "\n".join(
        f"{context_line(row)}:{a}:{b}:{cap.numerator}/{cap.denominator}"
        for row, (a, b, cap) in sorted(FIRST_ROWS.items())
    )
    + "\n"
).encode()


print("SCALE_THREE_HAMMING_SIX_SHEET_CLASSIFICATION")
print("arithmetic=integer+rational floating_point=none metric_recursion=not_run")
print(f"direct_NAE_states={DIRECT_NAE_CHECKS} mismatches=0")
print(f"all_order_one_presentations_excluded_nonprimitive={len(tuple(combinations(range(1, P), 6)))}")
print(f"presentations={len(PRESENTATIONS)} by_type={dict(sorted(PRESENTATION_COUNTS.items()))}")
print(f"contexts={len(CONTEXTS)} by_type={dict(sorted(CONTEXT_COUNTS.items()))}")
print(f"units_per_presentation={dict(sorted(UNITS_PER_PRESENTATION.items()))}")
print(f"presentation_payload_sha256={sha256(PRESENTATION_PAYLOAD).hexdigest()}")
print(f"context_payload_sha256={sha256(CONTEXT_PAYLOAD).hexdigest()}")
print(
    f"presentation_multiplier_orbits={len(PRESENTATION_ORBITS)} "
    f"size_by_type={PRESENTATION_ORBIT_SUMMARY}"
)
print(
    f"context_multiplier_orbits={len(CONTEXT_ORBITS)} "
    f"size_by_type={CONTEXT_ORBIT_SUMMARY}"
)
print(
    f"context_C12xC2_orbits={len(SHEET_ORBITS)} "
    f"size_by_type={SHEET_ORBIT_SUMMARY}"
)
print("global_unit_flip_fixed=0 multiplier_flip_stabilizers=0")
print("canonical_C12xC2_representatives=68 split=k4:16,k5:28,k6:24")
print(
    "canonical_scheme k4:C=1,5,8,12 B=23|24|26|27|29|2.11|46|49 "
    "words=1111|1221"
)
print(
    "canonical_scheme k5:C=1,2,5,8,12 B=3|4|6|7|9|10|11 "
    "words=11111|11221|12111|12221"
)
print("canonical_scheme k6:representatives_and_words_are_the_five_matching_code_rows")
print(
    f"inversion_presentations={len(INVERTIBLE_PRESENTATIONS)}/{len(PRESENTATIONS)} "
    f"by_type={dict(sorted(Counter(order_type(row) for row in INVERTIBLE_PRESENTATIONS).items()))} "
    f"fixed={sum(invert_presentation(row) == row for row in INVERTIBLE_PRESENTATIONS)}"
)
print(
    f"inversion_contexts={len(INVERTIBLE_CONTEXTS)}/{len(CONTEXTS)} "
    f"by_type={dict(sorted(Counter(order_type(row) for row in INVERTIBLE_CONTEXTS).items()))} "
    f"fixed={sum(invert_context(row) == row for row in INVERTIBLE_CONTEXTS)}"
)
print(
    "inversion_then_flip_fixed="
    + str(
        sum(
            reflect_context(invert_context(row)) == row
            for row in INVERTIBLE_CONTEXTS
        )
    )
)

for size in (4, 5, 6):
    for labels, orbit_size, stabilizer, solutions in C_SET_ORBITS[size]:
        print(
            f"C_orbit k={size} C={labels} orbit={orbit_size} "
            f"stabilizer={stabilizer} unit_solutions={solutions} "
            f"graph={graph_fingerprint(labels)}"
        )
for labels, (equations, free) in MATCHING_CODES.items():
    print(
        f"matching_code C={labels} equations={equations} free={free} "
        f"dimension={len(labels) - len(equations)}"
    )

print(
    f"metric_guardrail base_packet={METRIC_BASE} M={METRIC_VALUES['base'][0]} "
    f"witnesses={METRIC_VALUES['base'][1]}"
)
print(
    f"metric_guardrail flip_packet={METRIC_REFLECTION} M={METRIC_VALUES['flip'][0]} "
    f"witnesses={METRIC_VALUES['flip'][1]}"
)
print(
    f"metric_guardrail multiplier2_packet={METRIC_TIMES_TWO} "
    f"M={METRIC_VALUES['multiplier2'][0]} witnesses={METRIC_VALUES['multiplier2'][1]}"
)

print(
    f"distinct_metric_roots={len(ROOTS)} "
    f"by_type={dict((kind, len(ROOTS_BY_TYPE[kind])) for kind in ((2, 4), (1, 5), (0, 6)))}"
)
for kind in ((2, 4), (1, 5), (0, 6)):
    values = EXPECTED_FIRST[kind]
    print(
        f"first_layer type={kind} contexts={values['contexts']} "
        f"D1={values['order_one']} D3={values['order_three']} "
        f"total={values['order_one'] + values['order_three']} "
        f"per_context_min={values['minimum']} per_context_max={values['maximum']}"
    )
print(
    f"first_layer total={FIRST_TOTAL} D1={FIRST_ORDER_ONE} D3={FIRST_ORDER_THREE} "
    f"root_cap_range={min(ROOT_CAPS)}..{max(ROOT_CAPS)}"
)
print(f"first_layer_payload_sha256={sha256(FIRST_PAYLOAD).hexdigest()}")
print(
    f"planning_only_c2_depth2_projection={DEPTH_TWO_PROJECTION} "
    f"floor={DEPTH_TWO_PROJECTION.numerator // DEPTH_TWO_PROJECTION.denominator}"
)
print(f"planning_only_c2_full_tree_projection={FULL_TREE_PROJECTION}")
print("verdict=sheet_bank_classified_metric_bank_open")
print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")
