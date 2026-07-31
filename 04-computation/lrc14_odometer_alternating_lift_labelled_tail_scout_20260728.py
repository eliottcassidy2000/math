#!/usr/bin/env python3
"""Exact labelled scout for an alternating odometer-twisted LRC handoff.

Put ``R=13^6``, ``S=13^5``, ``k=S+2``, and

    T_-(x)={13x-k/R},       T_+(x)={13x+k/R}.

The points ``x_+=1/2+k/(14R)`` and ``x_-=1/2-k/(14R)`` form an exact
two-cycle.  Unlike the smaller ``k=14`` control, their inherited stored
clock edges are nonconstant and glue as ``4->3->4``.  This script inserts
the actual THM-2584/THM-2640 rails, present packets, predecessor carries,
half-digits, delayed guard sectors, private roots, and primitive-unit test.

There is a positive three-event (two-handoff) interval.  One fully typed
right-half path has states

    (shallow,owner,c,h,kappa,r)
      (4,3, 7,7,1, 2),
      (3,4, 5,9,1,11),
      (4,3,11,8,1,10).

The root handoff must retain the intermediate D-root.  It is

    2 --D-root 2 --(+9)--> 11,
   11 --D-root 6 --(+4)--> 10,

not the ill-typed shortcut which adds the odometer step directly to the
current private root.

The delayed coordinate supplies a sharp stopping mechanism.  For every
integer lift k,

    y={Rx}  ==>  y(T_k x)={13y};

the integer translation disappears.  In fact the union of the two inherited
guard-sector words is positive for three event states and exactly empty for
four *before choosing any clock*.  Every clock-cut delayed word is a subset,
so the four-event stopping law is uniform over the full intrinsic-clock
alphabet and every sequence of integer odometer lifts.  The sharper
4,3,4,3 ledger is retained to locate the positive witness.  Two independent
exact paths verify both statements: a Fraction cylinder recursion and a
fixed refined-integer-grid periodic intersection.  Thus the scout gives a
positive three-event packet, but no four-event (hence no seven-event)
continuation in the inherited delayed carrier.

There is also a much smaller analytic certificate hidden inside that exact
intersection.  The union of the two raw guard sectors cancels the guard cut
and is exactly

    D_(13^3) intersect D_14^c intersect D_27^c intersect ...
                    intersect D_(2*13^5)^c.

Already the first target tooth and the first unit safe factor suffice.  If
``I=[-1/14,1/14)`` and ``z={13^3 y}`` is represented in ``I``, then the
event-0,1,2 target conditions ``z,13z,13^2z in I`` have no wrap branch and
force ``z in I/13^2``.  Since ``14<13^2``, this implies ``14z in I``,
contradicting the event-3 safe factor at speed 14.  Thus the fourth-event
zero does not use the event-3 target tooth, the other four unit speeds, the
``2*13^5`` tooth, a guard sector, a clock label, or an odometer lift.  The
large exact replay below independently checks this minimal carrier and its
positive depths one through three.

Scope is finite exact support and unit typing for the stated inherited
packets.  It is not an endpoint-current transport theorem, row exclusion,
holonomy trivialization, scalar closure, or LRC(14).
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from itertools import product

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_predecessor_carry_private_root_atlas_thm2640 as private
import lrc14_successor_halfcell_carry_no_go_thm2623 as prior


P = 13
Q7 = 7
R = P**6
S = P**5
K = S + 2
T = cross.T
GRID_STEP = T // (14 * R)
CENTRE = T // 2
PLUS = CENTRE + K * GRID_STEP
MINUS = CENTRE - K * GRID_STEP
WITNESS_OFFSET = 1_504_594
DELTAS = ((-2 * K) % P, (2 * K) % P)

EXPECTED_TAIL = (
    (41_217, Fraction(524_901_093_439, 11_455_265_301_480)),
    (12_610, Fraction(797_217, 743_328_586)),
    (5_214, Fraction(28_677, 878_479_238)),
    (0, Fraction(0)),
)

EXPECTED_RAW_TAIL = (
    (47_484, Fraction(604_725_613_249, 11_455_265_301_480)),
    (16_244, Fraction(513_351, 371_664_293)),
    (6_776, Fraction(2_662, 62_748_517)),
    (0, Fraction(0)),
)

EXPECTED_MINIMAL_TAIL = (
    (1_886, Fraction(6, 49)),
    (1_606, Fraction(22_187, 2_798_978)),
    (1_452, Fraction(1_452, 2_599_051)),
    (0, Fraction(0)),
)

EXPECTED_RAW_SECTOR_LANGUAGE = {
    1: {
        (0,): (32_578, Fraction(207_428_797_703, 5_727_632_650_740)),
        (1,): (14_906, Fraction(2_707_797, 163_368_920)),
    },
    2: {
        (0, 0): (4_728, Fraction(447_838, 1_114_992_879)),
        (0, 1): (6_600, Fraction(18_975, 33_787_663)),
        (1, 0): (3_420, Fraction(46_301, 159_284_697)),
        (1, 1): (1_496, Fraction(4_301, 33_787_663)),
    },
    3: {
        (0, 0, 0): (1_232, Fraction(484, 62_748_517)),
        (0, 0, 1): (528, Fraction(1_452, 439_239_619)),
        (0, 1, 0): (2_420, Fraction(6_655, 439_239_619)),
        (0, 1, 1): (660, Fraction(1_815, 439_239_619)),
        (1, 0, 0): (880, Fraction(2_420, 439_239_619)),
        (1, 0, 1): (396, Fraction(1_089, 439_239_619)),
        (1, 1, 0): (528, Fraction(1_452, 439_239_619)),
        (1, 1, 1): (132, Fraction(363, 439_239_619)),
    },
    4: {},
}

EXPECTED_SECTOR_DEPTH3 = {
    (0, 0, 0): (682, Fraction(3_751, 878_479_238)),
    (0, 0, 1): (396, Fraction(1_089, 439_239_619)),
    (0, 1, 0): (2_024, Fraction(5_566, 439_239_619)),
    (0, 1, 1): (660, Fraction(1_815, 439_239_619)),
    (1, 0, 0): (550, Fraction(3_025, 878_479_238)),
    (1, 0, 1): (330, Fraction(1_815, 878_479_238)),
    (1, 1, 0): (440, Fraction(1_210, 439_239_619)),
    (1, 1, 1): (132, Fraction(363, 439_239_619)),
}

EXPECTED_ALTERNATING_SECTOR_LANGUAGE = {
    1: {
        (0,): (28_402, Fraction(180_836_059_223, 5_727_632_650_740)),
        (1,): (12_815, Fraction(48_885_587, 3_430_747_320)),
    },
    2: {
        (0, 0): (3_024, Fraction(573_073, 2_229_985_758)),
        (0, 1): (5_808, Fraction(16_698, 33_787_663)),
        (1, 0): (2_458, Fraction(233_020, 1_114_992_879)),
        (1, 1): (1_320, Fraction(3_795, 33_787_663)),
    },
    3: EXPECTED_SECTOR_DEPTH3,
    4: {},
}

EXPECTED_X_INTERVAL = (
    Fraction(4_286_847_956_371_849, 8_480_502_984_583_084),
    Fraction(4_286_847_956_371_861, 8_480_502_984_583_084),
)
EXPECTED_W_INTERVAL = (
    Fraction(42_841_277, 1_756_958_476),
    Fraction(42_841_289, 1_756_958_476),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge_intervals(intervals):
    """Merge sorted or unsorted half-open intervals."""
    out = []
    for left, right in sorted(intervals):
        require(left < right, "degenerate interval entered a support union")
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def interval_mass(intervals):
    return sum((right - left for left, right in intervals), Fraction(0))


def inside_unweighted(point, intervals):
    starts = [left for left, _ in intervals]
    index = bisect_right(starts, point) - 1
    return index >= 0 and intervals[index][0] <= point < intervals[index][1]


def containing_unweighted(point, intervals):
    starts = [left for left, _ in intervals]
    index = bisect_right(starts, point) - 1
    require(index >= 0 and intervals[index][0] <= point < intervals[index][1],
            "point left its asserted unweighted component")
    return intervals[index]


def containing_weighted(point, pieces):
    starts = [left for left, _, _ in pieces]
    index = bisect_right(starts, point) - 1
    require(index >= 0 and pieces[index][0] <= point < pieces[index][1],
            "point left its asserted weighted component")
    return pieces[index]


def prefix_intervals(prefix):
    starts, lengths, _ = prefix
    return [(left, left + length) for left, length in zip(starts, lengths)]


def fraction_pullback_intersection(current, target, factor):
    """Intersect current x-cylinders with ``{factor*x} in target``.

    The current cylinders are already small, so this local branch recursion
    never constructs all ``factor`` global preimages.  All endpoints remain
    exact Fractions.
    """
    starts = [left for left, _ in target]
    out = []
    for left, right in current:
        image_left = factor * left
        image_right = factor * right
        first_branch = image_left.numerator // image_left.denominator
        last_branch = image_right.numerator // image_right.denominator
        for branch in range(first_branch, last_branch + 1):
            local_left = max(Fraction(0), image_left - branch)
            local_right = min(Fraction(1), image_right - branch)
            if local_left >= local_right:
                continue
            index = max(0, bisect_right(starts, local_left) - 1)
            while index < len(target) and target[index][0] < local_right:
                overlap_left = max(local_left, target[index][0])
                overlap_right = min(local_right, target[index][1])
                if overlap_left < overlap_right:
                    out.append(((overlap_left + branch) / factor,
                                (overlap_right + branch) / factor))
                index += 1
    return merge_intervals(out)


def fraction_tail(sets, depth, sectors=None):
    """Exact source-y cylinders for a 4,3,4,3,... delayed itinerary."""
    require(depth >= 1, "tail depth must be positive")
    if sectors is None:
        current = sets[4]
    else:
        require(len(sectors) == depth, "sector word has wrong length")
        current = sets[sectors[0], 4]
    for index in range(1, depth):
        ell = 3 if index % 2 else 4
        target = sets[ell] if sectors is None else sets[sectors[index], ell]
        current = fraction_pullback_intersection(
            current, target, P**index
        )
        if not current:
            break
    return current


def periodic_integer_restrict(current, target, period):
    """Restrict integer-grid intervals by a periodic target union."""
    starts = [left for left, _ in target]
    out = []
    for left, right in current:
        first_period = left // period
        last_period = (right - 1) // period
        for block in range(first_period, last_period + 1):
            block_left = block * period
            local_left = max(0, left - block_left)
            local_right = min(period, right - block_left)
            if local_left >= local_right:
                continue
            index = max(0, bisect_right(starts, local_left) - 1)
            while index < len(target) and target[index][0] < local_right:
                overlap_left = max(local_left, target[index][0])
                overlap_right = min(local_right, target[index][1])
                if overlap_left < overlap_right:
                    out.append((block_left + overlap_left,
                                block_left + overlap_right))
                index += 1
    return merge_intervals(out)


def integer_itinerary(base_sets):
    """Independent fixed-grid replay for a list of base-grid target sets."""
    depth = len(base_sets)
    require(depth >= 1, "integer itinerary must be nonempty")
    scale0 = P ** (depth - 1)
    current = [(scale0 * left, scale0 * right)
               for left, right in base_sets[0]]
    for index in range(1, depth):
        scale = P ** (depth - 1 - index)
        target = [(scale * left, scale * right)
                  for left, right in base_sets[index]]
        current = periodic_integer_restrict(
            current, target, scale * T
        )
        if not current:
            break
    return current, scale0 * T


def integer_tail(base_sets, depth):
    """Independent fixed-grid replay for clocks ``4,3,4,3,...``."""
    return integer_itinerary([
        base_sets[4 if index % 2 == 0 else 3]
        for index in range(depth)
    ])


def fraction_repeated_tail(base_set, depth):
    """Fraction replay of one set under y -> {13y} for ``depth`` states."""
    current = base_set
    for index in range(1, depth):
        current = fraction_pullback_intersection(
            current, base_set, P**index
        )
        if not current:
            break
    return current


def fraction_sector_language(sector_sets, depth, clock_itinerary=None):
    """All positive binary sector words at one requested depth."""
    require(depth >= 1, "sector-language depth must be positive")
    if clock_itinerary is not None:
        require(len(clock_itinerary) == depth,
                "clock itinerary has wrong sector-language depth")
    language = {}
    for sectors in product(range(2), repeat=depth):
        if clock_itinerary is None:
            current = sector_sets[sectors[0]]
        else:
            current = sector_sets[sectors[0], clock_itinerary[0]]
        for index in range(1, depth):
            target = (sector_sets[sectors[index]]
                      if clock_itinerary is None
                      else sector_sets[sectors[index], clock_itinerary[index]])
            current = fraction_pullback_intersection(
                current, target, P**index
            )
            if not current:
                break
        if current:
            language[sectors] = (len(current), interval_mass(current))
    return language


def clock_label(point, iterate):
    phase = Fraction((P**iterate * point) % T, T)
    return int((Q7 * phase + Fraction(1, 2)) // 1) % Q7


def state_digits(point):
    quotient, future_coordinate = divmod(R * point, T)
    h = P * future_coordinate // T
    kappa = 2 * P * future_coordinate // T - 2 * h
    return quotient % P, h, kappa, future_coordinate


def comb_component(point, speed, low, high):
    require(0 <= low < high <= 182, "comb component must not wrap")
    unit = T // (182 * speed)
    step = 182 * unit
    base = low * unit
    branch = (point - base) // step
    left = base + branch * step
    right = left + (high - low) * unit
    require(left <= point < right, "point is not in requested comb component")
    return left, right


def root_memberships(module, point):
    phase = module.C3 * point % T
    out = []
    for edge in (0, 1):
        for root in range(1, P):
            low = 14 * root - 13 if edge == 0 else 14 * root
            high = low + 13
            if low * T <= 182 * phase < high * T:
                out.append((edge, root))
    return tuple(out)


def unit_vector(module, pair_prefixes, present, present_starts,
                rail, sector, edge, carry, h, kappa):
    root = (2 * carry + (2 * h + kappa) // P + (edge == 0)) % P
    require(root != 0, "zero root cannot be a private unit")
    values = []
    for ell5 in range(Q7):
        overlap = cross.old.intersect_weighted_union(
            rail[3], present[ell5, (-h) % P],
            present_starts[ell5, (-h) % P]
        )
        low = 14 * root - 13 if edge == 0 else 14 * root
        high = low + 13
        half = cross.old.intersect_weighted_comb(
            overlap, module.C3, 182, low, high
        )
        carry_values = private.delayed_carry_pair(
            half, pair_prefixes[sector][ell5][h], {}
        )
        values.append(carry_values[carry][kappa])
    values = tuple(values)
    unit = private.is_unit(values, root, cross.GLOBAL_CONTENT)
    raw_mod = tuple((value // cross.GLOBAL_CONTENT) % P for value in values)
    normalized = tuple(
        value * pow(root, -1, P) % P for value in raw_mod
    )
    reduced = tuple(
        (normalized[index] - normalized[-1]) % P
        for index in range(Q7 - 1)
    )
    determinant = cross.old.sat.multiplication_determinant_7(reduced)
    require(unit == bool(determinant), "unit predicate/determinant drift")
    return root, raw_mod, determinant


def main():
    require(T % (14 * R) == 0 and GRID_STEP > 0,
            "canonical grid stopped resolving the alternating lift")
    require(K % P == 2 and DELTAS == (9, 4),
            "odometer quotient steps changed")

    (module, _prefixes, _whole_prefixes, _digit_masses, rails,
     present, present_starts) = cross.build_carrier_data()
    pair_prefixes = private.build_pair_prefixes(module)

    # ------------------------------------------------------------------
    # Uniform delayed skew product and the sharp depth-three/depth-four cut
    # ------------------------------------------------------------------
    sector_words = prior.sector_words(module)
    raw_sector_fraction = {
        sector: [(Fraction(left, T), Fraction(right, T))
                 for left, right in word]
        for sector, word in enumerate(sector_words)
    }
    raw_union_base = merge_intervals(
        sector_words[0] + sector_words[1]
    )
    guard_free_base = module.make_comb(
        module.C2, 182, -13, 13
    )
    for unit_index in module.UNIT_IDX:
        guard_free_base = module.subtract_comb(
            guard_free_base, module.W[unit_index], 182, -13, 13
        )
    guard_free_base = module.subtract_comb(
        guard_free_base, module.C3, 182, -13, 13
    )
    require(raw_union_base == guard_free_base,
            "the two guard sectors stopped partitioning the guard-free word")
    raw_union_fraction = [
        (Fraction(left, T), Fraction(right, T))
        for left, right in raw_union_base
    ]
    require((len(sector_words[0]), len(sector_words[1]),
             len(raw_union_base)) == (32_578, 14_906, 47_484),
            "raw delayed sector-word census changed")

    # Minimal human-readable stopping certificate.  Only D_(13^3) and the
    # q_1=14 safe factor are retained.  On the centered half-open danger
    # interval I=[-1/14,1/14), the equality 13/14=1-1/14 makes every
    # apparent wrap branch land exactly on the excluded right endpoint.
    # Hence z,13z,13^2z in I imply z in I/13^2.  The terminal q_1 factor is
    # then dangerous because 14/13^2<1.  The interval replay is an
    # independent positive/hostile control for that two-factor proof.
    require(module.C2 == P**3 and module.W[module.UNIT_IDX[0]] == P + 1,
            "minimal target/unit scale identity changed")
    require(Fraction(P, P + 1) == 1 - Fraction(1, P + 1)
            and P + 1 < P**2,
            "minimal no-wrap/contraction inequality changed")
    minimal_base = module.make_comb(module.C2, 182, -13, 13)
    minimal_base = module.subtract_comb(
        minimal_base, module.W[module.UNIT_IDX[0]], 182, -13, 13
    )
    require(cross.old.intersect_sorted(raw_union_base, minimal_base)
            == raw_union_base,
            "raw delayed word escaped the minimal stopping carrier")
    minimal_fraction = [
        (Fraction(left, T), Fraction(right, T))
        for left, right in minimal_base
    ]
    minimal_tail_rows = []
    for depth in range(1, 5):
        support = fraction_repeated_tail(minimal_fraction, depth)
        mass = interval_mass(support)
        require((len(support), mass) == EXPECTED_MINIMAL_TAIL[depth - 1],
                f"minimal target/unit tail changed at depth {depth}")
        minimal_tail_rows.append((depth, len(support), mass))

    raw_tail_rows = []
    for depth in range(1, 5):
        fraction_support = fraction_repeated_tail(
            raw_union_fraction, depth
        )
        integer_support, refined_grid = integer_itinerary(
            [raw_union_base] * depth
        )
        fraction_as_integer = [
            (int(left * refined_grid), int(right * refined_grid))
            for left, right in fraction_support
        ]
        require(all(left.denominator == right.denominator == 1
                    for left, right in
                    ((left * refined_grid, right * refined_grid)
                     for left, right in fraction_support)),
                "raw Fraction cylinders missed the refined integer grid")
        require(fraction_as_integer == integer_support,
                "raw Fraction/integer tail replays disagree")
        mass = interval_mass(fraction_support)
        require((len(fraction_support), mass)
                == EXPECTED_RAW_TAIL[depth - 1],
                f"raw delayed tail changed at depth {depth}")
        raw_tail_rows.append((depth, len(fraction_support), mass))

    raw_sector_language = {
        depth: fraction_sector_language(raw_sector_fraction, depth)
        for depth in range(1, 5)
    }
    require(raw_sector_language == EXPECTED_RAW_SECTOR_LANGUAGE,
            "raw binary sector tail-language changed")

    sector_base = {}
    sector_fraction = {}
    for sector, word in enumerate(sector_words):
        for ell in (3, 4):
            intervals = module.subtract_comb(
                word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
            )
            require(cross.old.intersect_sorted(word, intervals) == intervals,
                    "clock-cut delayed word escaped its raw sector word")
            sector_base[sector, ell] = intervals
            sector_fraction[sector, ell] = [
                (Fraction(left, T), Fraction(right, T))
                for left, right in intervals
            ]
    require(
        {(ell, sector): len(sector_base[sector, ell])
         for ell in (3, 4) for sector in range(2)}
        == {(3, 0): 28_402, (3, 1): 12_815,
            (4, 0): 28_402, (4, 1): 12_815},
        "delayed sector interval census changed",
    )

    union_base = {
        ell: merge_intervals(
            sector_base[0, ell] + sector_base[1, ell]
        ) for ell in (3, 4)
    }
    union_fraction = {
        ell: [(Fraction(left, T), Fraction(right, T))
              for left, right in union_base[ell]]
        for ell in (3, 4)
    }
    require(all(len(union_base[ell]) == 41_217 for ell in (3, 4)),
            "sector union census changed")

    tail_rows = []
    for depth in range(1, 5):
        fraction_support = fraction_tail(union_fraction, depth)
        integer_support, refined_grid = integer_tail(union_base, depth)
        fraction_as_integer = [
            (int(left * refined_grid), int(right * refined_grid))
            for left, right in fraction_support
        ]
        require(all(left.denominator == right.denominator == 1
                    for left, right in
                    ((left * refined_grid, right * refined_grid)
                     for left, right in fraction_support)),
                "Fraction cylinders did not land on the refined grid")
        require(fraction_as_integer == integer_support,
                "independent Fraction/integer tail replays disagree")
        mass = interval_mass(fraction_support)
        require((len(fraction_support), mass) == EXPECTED_TAIL[depth - 1],
                f"alternating delayed tail changed at depth {depth}")
        tail_rows.append((depth, len(fraction_support), mass))

    alternating_sector_language = {
        depth: fraction_sector_language(
            sector_fraction, depth,
            tuple(4 if index % 2 == 0 else 3 for index in range(depth))
        ) for depth in range(1, 5)
    }
    require(alternating_sector_language
            == EXPECTED_ALTERNATING_SECTOR_LANGUAGE,
            "alternating-clock binary sector tail-language changed")

    # The raw union has no four-state orbit.  Since every clock-cut word is
    # a subset of one raw sector word, this proves the fourth-event zero for
    # every clock itinerary without enumerating 7^4 clock words.
    require(not raw_tail_rows[-1][1]
            and not raw_sector_language[4],
            "raw delayed carrier acquired a four-event orbit")

    # ------------------------------------------------------------------
    # Exact affine cycle, typed three-event witness, and unit census
    # ------------------------------------------------------------------
    lift_grid = K * T // R
    require((P * PLUS - lift_grid) % T == MINUS
            and (P * MINUS + lift_grid) % T == PLUS,
            "k=S+2 stopped giving the exact alternating two-cycle")
    require(
        ((clock_label(PLUS, 1), clock_label(PLUS, 2)),
         (clock_label(MINUS, 1), clock_label(MINUS, 2)))
        == ((4, 3), (3, 4)),
        "intrinsic stored clock edges changed",
    )

    points = (
        PLUS + WITNESS_OFFSET,
        MINUS + P * WITNESS_OFFSET,
        PLUS + P**2 * WITNESS_OFFSET,
    )
    require((P * points[0] - lift_grid) % T == points[1]
            and (P * points[1] + lift_grid) % T == points[2],
            "witness points do not follow the alternating affine handoff")

    configurations = (
        # owner, deep, shallow, carry, h, kappa, right root, D-right root
        (3, 0, 4, 7, 7, 1, 2, 2),
        (4, 12, 3, 5, 9, 1, 11, 6),
        (3, 0, 4, 11, 8, 1, 10, None),
    )
    expected_vectors = (
        (0, 0, 0, 0, 1, 6, 10),
        (0, 0, 8, 8, 0, 0, 0),
        (2, 0, 0, 0, 10, 4, 7),
    )
    expected_determinants = (2, 1, 1)
    expected_active_sources = (
        (1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12),
        (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12),
        (1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12),
    )
    expected_unit_sources = (
        (1, 3, 4, 5, 6, 8, 9, 10, 11, 12),
        (1, 2, 3, 4, 5, 7, 8, 9, 10, 12),
        (1, 3, 4, 5, 6, 8, 9, 10, 11, 12),
    )

    selected_rails = []
    unit_rows = []
    for event, (point, config) in enumerate(zip(points, configurations)):
        owner, deep, shallow, carry, h, kappa, right_root, d_root = config
        actual_carry, actual_h, actual_kappa, future_coordinate = (
            state_digits(point)
        )
        require((actual_carry, actual_h, actual_kappa)
                == (carry, h, kappa),
                f"event {event} physical carry/future labels changed")
        require(clock_label(point, 1) == shallow
                and clock_label(point, 2) == owner,
                f"event {event} stored clock labels changed")

        active = []
        for rail in rails:
            if rail[1] == owner and rail[2] == deep:
                try:
                    containing_weighted(point, rail[3])
                except RuntimeError:
                    continue
                active.append(rail)
        active_sources = tuple(rail[0] for rail in active)
        require(active_sources == expected_active_sources[event],
                f"event {event} active rail-source census changed")
        rail = next(row for row in active if row[0] == 1)
        selected_rails.append(rail)

        require(inside_unweighted(
            point, present[shallow, (-h) % P]
        ), f"event {event} left its present packet")
        require(root_memberships(module, point)
                == ((0, right_root + 1), (1, right_root)),
                f"event {event} private half-root overlap changed")
        require(inside_unweighted(
            future_coordinate,
            prefix_intervals(pair_prefixes[0][shallow][h][kappa])
        ), f"event {event} left delayed sector zero")
        require(not inside_unweighted(
            future_coordinate,
            prefix_intervals(pair_prefixes[1][shallow][h][kappa])
        ), f"event {event} unexpectedly entered delayed sector one")

        per_edge_sources = {}
        selected_vector = None
        selected_det = None
        for edge in (0, 1):
            sources = []
            for candidate in active:
                root, vector, determinant = unit_vector(
                    module, pair_prefixes, present, present_starts,
                    candidate, 0, edge, carry, h, kappa
                )
                expected_root = right_root + (edge == 0)
                require(root == expected_root,
                        f"event {event} private-root formula changed")
                if determinant:
                    sources.append(candidate[0])
                if candidate[0] == 1 and edge == 1:
                    selected_vector = vector
                    selected_det = determinant
            per_edge_sources[edge] = tuple(sources)
        require(per_edge_sources[0] == per_edge_sources[1]
                == expected_unit_sources[event],
                f"event {event} unit-source census changed")
        require(selected_vector == expected_vectors[event]
                and selected_det == expected_determinants[event],
                f"event {event} selected unit vector changed")
        unit_rows.append((active_sources, per_edge_sources[1],
                          selected_vector, selected_det))

        if d_root is not None:
            d_point = P * point % T
            require((1, d_root) in root_memberships(module, d_point),
                    f"event {event} intermediate D-root changed")
            require((d_root + DELTAS[event]) % P
                    == configurations[event + 1][6],
                    f"event {event} D-root/translation typing changed")

    require(configurations[0][4] - K % P == configurations[1][3]
            and (configurations[1][4] + K) % P
            == configurations[2][3],
            "predecessor carry covariance changed")

    # ------------------------------------------------------------------
    # Extract one exact positive interval supporting all selected labels.
    # The variable d is a displacement in the base T-grid coordinate from
    # WITNESS_OFFSET.  Every condition gives an exact half-open d-interval.
    # ------------------------------------------------------------------
    constraints = []

    def add_x_constraint(name, event, left, right, slope=None):
        if slope is None:
            slope = P**event
        point = points[event]
        constraints.append((Fraction(left - point, slope),
                            Fraction(right - point, slope),
                            name, event))

    for event, (point, config, rail) in enumerate(
            zip(points, configurations, selected_rails)):
        owner, deep, shallow, carry, h, kappa, right_root, d_root = config
        left, right, _weight = containing_weighted(point, rail[3])
        add_x_constraint("rail", event, left, right)

        left, right = containing_unweighted(
            point, present[shallow, (-h) % P]
        )
        add_x_constraint("present", event, left, right)

        predecessor = R * point // T
        add_x_constraint(
            "carry", event,
            Fraction(predecessor * T, R),
            Fraction((predecessor + 1) * T, R),
        )

        left, right = comb_component(
            point, module.C3, 14 * right_root, 14 * right_root + 13
        )
        add_x_constraint("current_root", event, left, right)

        future_coordinate = R * point % T
        delayed_component = containing_unweighted(
            future_coordinate,
            prefix_intervals(pair_prefixes[0][shallow][h][kappa])
        )
        constraints.append((
            Fraction(delayed_component[0] - future_coordinate,
                     R * P**event),
            Fraction(delayed_component[1] - future_coordinate,
                     R * P**event),
            "delayed_Q", event,
        ))

        if d_root is not None:
            d_point = P * point % T
            left, right = comb_component(
                d_point, module.C3, 14 * d_root, 14 * d_root + 13
            )
            constraints.append((
                Fraction(left - d_point, P**(event + 1)),
                Fraction(right - d_point, P**(event + 1)),
                "D_root", event,
            ))

    lower = max(row[0] for row in constraints)
    upper = min(row[1] for row in constraints)
    require(lower < upper, "selected labelled three-event interval is empty")
    require(
        [(name, event) for lo, hi, name, event in constraints if lo == lower]
        == [("delayed_Q", 2)]
        and [(name, event) for lo, hi, name, event in constraints if hi == upper]
        == [("delayed_Q", 2)],
        "positive interval stopped being limited by the third delayed atom",
    )

    w_interval = (
        Fraction(WITNESS_OFFSET) + lower,
        Fraction(WITNESS_OFFSET) + upper,
    )
    w_interval = tuple(value / (T // R) for value in w_interval)
    x_interval = (
        (Fraction(PLUS + WITNESS_OFFSET) + lower) / T,
        (Fraction(PLUS + WITNESS_OFFSET) + upper) / T,
    )
    require(w_interval == EXPECTED_W_INTERVAL,
            "normalized positive interval changed")
    require(x_interval == EXPECTED_X_INTERVAL,
            "physical positive interval changed")
    require(x_interval[1] - x_interval[0]
            == Fraction(3, 2_120_125_746_145_771),
            "positive interval width changed")

    independent_source_paths = 2 * 10**3
    same_source_paths = 2 * len(set(expected_unit_sources[0]).intersection(
        expected_unit_sources[1], expected_unit_sources[2]
    ))
    require((independent_source_paths, same_source_paths) == (2_000, 16),
            "unit-labelled path census changed")

    print("LRC14 odometer alternating-lift labelled tail scout")
    print("status=VERIFIED-EXACT SCOUT; inherited rail/present/carry/delayed/private-root/unit support")
    print(f"p={P}; R={R}; S={S}; k=S+2={K}; lifts={(-K, K)}; quotient_steps={DELTAS}")
    print(f"cycle_points=({Fraction(PLUS,T)},{Fraction(MINUS,T)}); stored_edges=((4,3),(3,4)); interfaces_glue=True")
    print("skew_product=y={R*x}; y(T_k x)={13*y} for every integer k; translation_is_invisible=True")
    print("tail_columns=(event_states,components,exact_mass)")
    print("minimal_D_13^3_and_safe_14_tail=")
    for row in minimal_tail_rows:
        print(row)
    print("minimal_analytic_certificate=z,13z,13^2z_in_I_implies_z_in_I/13^2; 14<13^2 forces_14z_in_I; event3_safe14_contradiction")
    print("raw_guard_sector_union_tail_before_clock_cut=")
    for row in raw_tail_rows:
        print(row)
    print("raw_binary_sector_language_columns=(depth,sector_word,components,exact_mass)")
    print(tuple(
        (depth, sectors, values[0], values[1])
        for depth in range(1, 5)
        for sectors, values in sorted(raw_sector_language[depth].items())
    ))
    print("alternating_clock_4_3_tail=")
    for row in tail_rows:
        print(row)
    print("alternating_binary_sector_language_columns=(depth,sector_word,components,exact_mass)")
    print(tuple(
        (depth, sectors, values[0], values[1])
        for depth in range(1, 5)
        for sectors, values in sorted(
            alternating_sector_language[depth].items()
        )
    ))
    print("depth4_raw_sector_words=0/16; depth4_full_clock_alphabet=EMPTY_by_subset; all_event_depths_ge_4=EMPTY; seven_event_tail=EMPTY")
    print("tail_automaton=all 2^d binary sector words survive for d=1,2,3; none survive for d=4; refined component transition has no recurrent SCC")
    print(f"witness_offset={WITNESS_OFFSET}; source_w_interval=[{w_interval[0]},{w_interval[1]}); source_x_interval=[{x_interval[0]},{x_interval[1]})")
    print(f"positive_interval_width={x_interval[1]-x_interval[0]}; limiting_constraint=event2_delayed_Q")
    print("event_columns=(shallow,owner,c,h,kappa,right_root,D_right_root,active_rails,unit_rails,right_unit_vector,determinant)")
    for config, row in zip(configurations, unit_rows):
        owner, deep, shallow, carry, h, kappa, root, d_root = config
        active_sources, unit_sources, vector, determinant = row
        print((shallow, owner, carry, h, kappa, root, d_root,
               len(active_sources), len(unit_sources), vector, determinant))
    print("carry_handoffs=(7 --h7-2--> 5, 5 --h9+2--> 11)")
    print("root_handoffs=(r2 --D:r2,+9--> r11, r11 --D:r6,+4--> r10)")
    print(f"unit_labelled_three_event_paths_independent_sources={independent_source_paths}; same_source={same_source_paths}")
    print("mechanism=three-event affine escape is real, but the k-independent pre-clock delayed y-tail is nilpotent at the fourth event")
    print("scope=finite exact support/unit scout; not endpoint transport, row exclusion, holonomy, scalar closure, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
