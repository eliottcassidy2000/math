#!/usr/bin/env python3
"""Exact periodic-tail scout for variable affine odometer handoffs.

This companion studies

    T_k(x) = {13*x + k/R},                 R=13^6,

with a variable integer carry k at each edge.  It deliberately does not
repeat the finite-depth two-cycle fibre-product scan.  Its object is the
exact periodic skew product before integration.

Write, on the R-fold circle,

    R*x_i = N_i + y_i,        N_i in Z/RZ, 0 <= y_i < 1,
    a_i   = floor(13*y_i).

Then every integer carry disappears from the terminal tail and acts only in
the odometer fibre:

    y_(i+1) = {13*y_i},
    N_(i+1) = 13*N_i + k_i + a_i                 (mod R).

Consequently an n-periodic tail has

    y_0 = A/(13^n-1),

and an exact period n iff no proper d<n satisfies

    (13^d-1)A = 0                              (mod 13^n-1).

After n steps the fibre closure is

    (13^n-1)N_0
      + sum_i 13^(n-1-i)(k_i+a_i) = 0          (mod R).

Since gcd(13^n-1,R)=1, every periodic tail and carry word has one and only
one odometer lift N_0.  Variable carries steer the coarse digit, never y.

There is a sharp reason not to keep iterating the attractive symmetric
two-cycle.  If

    x_+ = 1/2+u,  x_- = 1/2-u,
    x_+ --T_(-K)--> x_- --T_K--> x_+,

then (13+1)u=K/R, so

    y_+ = {1/2+K/14},       y_- = {1/2-K/14}.

Thus 14*y_+ and 14*y_- are integers for every integer K: symmetric
two-cycles are forced onto the speed-14 resonance by the factor 13+1.

The finite scout uses the smallest variable lawful alphabet around the
THM-2657 scale S=13^5,

    {-S-2,-S-1,S+1,S+2},

and exhausts every word of length 3,...,7.  Exact period-four and period-six
nonresonant intrinsic-clock cycles exist; odd lengths do not.  One explicit
period-four cycle lies strictly inside a common-source rail and the forced
present packet at every state.  Nevertheless every exact period-four and
period-six cycle in this universe misses the inherited delayed guard-safe
word already at its first factor, the target-a danger D_(c2).  This is a
finite exact carrier boundary, not a universal affine-handoff no-go.

No Bockstein unit, semantic endpoint, full transition mass, row exclusion,
or LRC(14) conclusion is asserted.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd

import lrc14_cross_time_target_future_diagonal_thm2616 as cross


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13
Q7 = 7
R = P**6
S = P**5
CENTRE = Fraction(1, 2)
CENTRAL_TOOTH = (Fraction(13, 28), Fraction(15, 28))
LIFT_ALPHABET = (-S - 2, -S - 1, S + 1, S + 2)
PERIODS = tuple(range(3, 8))
TARGET_A_SPEED = P**3

require(R == 4826809 and S == 371293, "odometer scale drift")
require(LIFT_ALPHABET == (-371295, -371294, 371294, 371295),
        "near-carry alphabet drift")
require(all(value % P for value in LIFT_ALPHABET),
        "near-carry alphabet contains a non-lawful quotient-zero lift")


def circle(value):
    return value % 1


def handoff(value, lift):
    return circle(P * value + Fraction(lift, R))


def clock(value):
    return int((Q7 * circle(value) + Fraction(1, 2)) // 1) % Q7


def shallow(value):
    return clock(P * value)


def owner(value):
    return clock(P * P * value)


def tail(value):
    return circle(R * value)


def coarse_digit(value):
    return int((R * value) // 1) % R


def predecessor_carry(value):
    return coarse_digit(value) % P


def future_digit(value):
    return int(P * tail(value))


def circle_distance(value):
    residue = circle(value)
    return min(residue, 1 - residue)


def exact_tail_period(value, upper_bound):
    iterate = value
    for period in range(1, upper_bound + 1):
        iterate = circle(P * iterate)
        if iterate == value:
            return period
    raise RuntimeError("tail did not return at the asserted word length")


def congruence_tail_period(value, upper_bound):
    """Independent proper-period test on A/(13^n-1)."""
    denominator = P**upper_bound - 1
    numerator = value * denominator
    require(numerator.denominator == 1,
            "tail is not represented on its asserted periodic lattice")
    residue = numerator.numerator % denominator
    for period in range(1, upper_bound + 1):
        if ((P**period - 1) * residue) % denominator == 0:
            return period
    raise RuntimeError("tail congruence did not close")


def periodic_orbit(lifts):
    """Return the unique central-branch periodic orbit for one carry word."""
    length = len(lifts)
    denominator = P**length - 1
    weighted_carry = sum(
        P**(length - 1 - index) * lift
        for index, lift in enumerate(lifts)
    )
    q0 = -Fraction(weighted_carry, denominator)
    q = q0
    states = []
    for lift in lifts:
        value = CENTRE + q / R
        states.append(value)
        next_q = P * q + lift
        require(handoff(value, lift) == CENTRE + next_q / R,
                "central recurrence disagrees with the circle handoff")
        q = next_q
    require(q == q0 and handoff(states[-1], lifts[-1]) == states[0],
            "periodic central lift failed to close")
    return tuple(states)


def audit_skew_product(states, lifts):
    """Check the tail/fibre recurrence and the unique N0 closure formula."""
    length = len(states)
    denominator = P**length - 1
    tails = tuple(tail(value) for value in states)
    digits = tuple(int(P * value) for value in tails)
    coarse = tuple(coarse_digit(value) for value in states)

    for index in range(length):
        following = (index + 1) % length
        require(tails[following] == circle(P * tails[index]),
                "integer carry leaked into the terminal tail")
        require(
            coarse[following]
            == (P * coarse[index] + lifts[index] + digits[index]) % R,
            "odometer fibre recurrence changed",
        )

    numerator = tails[0] * denominator
    require(numerator.denominator == 1,
            "periodic tail left the 1/(13^n-1) lattice")
    tail_numerator = numerator.numerator % denominator
    fibre_sum = sum(
        P**(length - 1 - index) * (lifts[index] + digits[index])
        for index in range(length)
    )
    require(((denominator * coarse[0] + fibre_sum) % R) == 0,
            "n-step odometer closure congruence failed")
    require(gcd(denominator, R) == 1,
            "tail multiplier unexpectedly shares an odometer factor")
    unique_coarse = -fibre_sum * pow(denominator, -1, R) % R
    require(unique_coarse == coarse[0],
            "periodic tail/carry word lost its unique odometer lift")
    return tails, digits, coarse, tail_numerator


def intrinsic_clock_cycle(states):
    edges = tuple((shallow(value), owner(value)) for value in states)
    nonconstant = all(left != right for left, right in edges)
    interfaces = all(
        edges[index][1] == edges[(index + 1) % len(edges)][0]
        for index in range(len(edges))
    )
    return edges, nonconstant and interfaces


def speed14_resonant(tails):
    flags = tuple(circle(14 * value) == 0 for value in tails)
    require(len(set(flags)) == 1,
            "speed-14 resonance was not tail-orbit invariant")
    return flags[0]


# The symbolic proof is in the module docstring.  These exact controls include
# both resonant and nontrivially moving tail presentations.
for symmetric_carry in (1, 2, 14, S + 1, S + 2):
    amplitude = Fraction(symmetric_carry, 14 * R)
    plus = CENTRE + amplitude
    minus = CENTRE - amplitude
    require(handoff(plus, -symmetric_carry) == minus
            and handoff(minus, symmetric_carry) == plus,
            "symmetric two-cycle identity changed")
    require(circle(14 * tail(plus)) == circle(14 * tail(minus)) == 0,
            "symmetric two-cycle escaped speed-14 resonance")


EXPECTED_CENSUS = {
    3: (64, 0, Counter()),
    4: (256, 32, Counter({4: 24, 2: 6, 1: 2})),
    5: (1024, 0, Counter()),
    6: (4096, 128, Counter({6: 120, 2: 6, 1: 2})),
    7: (16384, 0, Counter()),
}

CLOCK_CYCLES = {}
CENSUS_ROWS = {}
for length in PERIODS:
    total = 0
    central = 0
    compatible = []
    period_histogram = Counter()
    resonant = 0
    for lifts in product(LIFT_ALPHABET, repeat=length):
        total += 1
        states = periodic_orbit(lifts)
        tails, _digits, _coarse, _numerator = audit_skew_product(states, lifts)
        if all(CENTRAL_TOOTH[0] <= value < CENTRAL_TOOTH[1]
               for value in states):
            central += 1
        else:
            continue
        edges, clock_compatible = intrinsic_clock_cycle(states)
        alternating_signs = all(
            (lifts[index] < 0) != (lifts[(index + 1) % length] < 0)
            for index in range(length)
        )
        require(clock_compatible == alternating_signs,
                "near-carry intrinsic clocks stopped detecting sign alternation")
        if not clock_compatible:
            continue
        period = exact_tail_period(tails[0], length)
        require(period == congruence_tail_period(tails[0], length),
                "orbit and modular proper-period tests disagree")
        period_histogram[period] += 1
        is_resonant = speed14_resonant(tails)
        resonant += int(is_resonant)
        compatible.append((lifts, states, tails, edges, period, is_resonant))

    expected_total, expected_compatible, expected_periods = EXPECTED_CENSUS[length]
    require(total == central == expected_total,
            f"length-{length} near-carry orbit left the central tooth")
    require(len(compatible) == expected_compatible,
            f"length-{length} clock-compatible census changed")
    require(period_histogram == expected_periods,
            f"length-{length} proper-period audit changed")
    if length in (4, 6):
        require(resonant == 4,
                f"length-{length} speed-14 hostile count changed")
        require(all(not row[5] for row in compatible if row[4] == length),
                f"exact length-{length} tail became speed-14 resonant")
    else:
        require(not compatible,
                "odd near-carry word acquired a clock-compatible cycle")
    CENSUS_ROWS[length] = (total, len(compatible), period_histogram, resonant)
    CLOCK_CYCLES[length] = tuple(compatible)


# ---------------------------------------------------------------------------
# Current rail/present filters and the exact delayed-word hostile
# ---------------------------------------------------------------------------
(MODULE, PREFIXES, _WHOLE_PREFIXES, _DIGIT_MASSES, RAILS,
 PRESENT, _PRESENT_STARTS) = cross.build_carrier_data()
T = cross.T
require(MODULE.C1 == P and MODULE.C2 == TARGET_A_SPEED,
        "current word scale drift")


def interval_margin(intervals, value):
    coordinate = value * T
    margins = [
        min(coordinate - left, right - coordinate)
        for left, right in intervals
        if left < coordinate < right
    ]
    return max(margins) if margins else None


def weighted_margin(pieces, value):
    coordinate = value * T
    margins = [
        min(coordinate - left, right - coordinate)
        for left, right, weight in pieces
        if weight and left < coordinate < right
    ]
    return max(margins) if margins else None


def source_one_rail(value):
    wanted_owner = owner(value)
    matches = []
    for source, rail_clock, rail_digit, pieces in RAILS:
        if source != 1 or rail_clock != wanted_owner:
            continue
        margin = weighted_margin(pieces, value)
        if margin is not None:
            matches.append((rail_digit, margin))
    require(len(matches) <= 1,
            "point entered multiple source-one rails in one owner cell")
    return matches[0] if matches else None


def present_margin(value):
    h = future_digit(value)
    return interval_margin(PRESENT[shallow(value), (-h) % P], value)


def delayed_word_hit(value):
    ell = shallow(value)
    h = future_digit(value)
    starts, lengths, _prefix = PREFIXES[ell][h]
    coordinate = tail(value) * T
    index = bisect_right(starts, coordinate) - 1
    if index < 0:
        return False
    return starts[index] < coordinate < starts[index] + lengths[index]


def target_a_margin(value):
    return circle_distance(TARGET_A_SPEED * tail(value)) - Fraction(1, 14)


EXPECTED_FILTER_HISTOGRAMS = {
    4: Counter({(4, 3, 0): 16, (4, 4, 0): 4, (4, 2, 0): 4}),
    6: Counter({(6, 5, 0): 48, (6, 4, 0): 48,
                (6, 6, 0): 12, (6, 3, 0): 12}),
}
EXPECTED_MIN_TARGET_MARGINS = {
    4: Fraction(1433, 4080),
    6: Fraction(847565, 2413404),
}

FILTER_ROWS = {}
for length in (4, 6):
    histogram = Counter()
    target_margins = []
    exact_cycles = [row for row in CLOCK_CYCLES[length] if row[4] == length]
    for _lifts, states, _tails, _edges, _period, _resonant in exact_cycles:
        rail_support = tuple(source_one_rail(value) is not None
                             for value in states)
        present_support = tuple(present_margin(value) is not None
                                for value in states)
        word_support = tuple(delayed_word_hit(value) for value in states)
        margins = tuple(target_a_margin(value) for value in states)
        require(min(margins) > 0,
                "delayed-word hostile is not outside target-a danger")
        require(not any(word_support),
                "exact near-carry cycle unexpectedly entered the delayed word")
        target_margins.extend(margins)
        histogram[
            (sum(rail_support), sum(present_support), sum(word_support))
        ] += 1
    require(histogram == EXPECTED_FILTER_HISTOGRAMS[length],
            f"length-{length} rail/present/word filter histogram changed")
    require(min(target_margins) == EXPECTED_MIN_TARGET_MARGINS[length],
            f"length-{length} target-a hostile margin changed")
    FILTER_ROWS[length] = (histogram, min(target_margins))


# A period-four witness chosen to pass the rail and present filters at every
# state.  The last carry differs from the first positive carry, which breaks
# the symmetric speed-14 resonance and raises the exact tail period to four.
WITNESS_LIFTS = (-S - 2, S + 1, -S - 2, S + 2)
WITNESS_STATES = periodic_orbit(WITNESS_LIFTS)
WITNESS_TAILS, WITNESS_FUTURE_DIGITS, WITNESS_COARSE, _ = audit_skew_product(
    WITNESS_STATES, WITNESS_LIFTS
)
WITNESS_EDGES, WITNESS_CLOCK_OK = intrinsic_clock_cycle(WITNESS_STATES)
WITNESS_CARRIES = tuple(value % P for value in WITNESS_COARSE)
WITNESS_ROOT_STEPS = tuple(2 * lift % P for lift in WITNESS_LIFTS)
WITNESS_RAILS = tuple(source_one_rail(value) for value in WITNESS_STATES)
WITNESS_PRESENT_MARGINS = tuple(present_margin(value)
                                for value in WITNESS_STATES)
WITNESS_TARGET_MARGINS = tuple(target_a_margin(value)
                               for value in WITNESS_STATES)

require(WITNESS_STATES == (
    Fraction(69684274489, 137853665040),
    Fraction(68169392917, 137853665040),
    Fraction(69684274321, 137853665040),
    Fraction(68169390733, 137853665040),
), "period-four x witness changed")
require(WITNESS_TAILS == (
    Fraction(16489, 28560), Fraction(14437, 28560),
    Fraction(16321, 28560), Fraction(12253, 28560),
), "period-four tail witness changed")
require(WITNESS_CLOCK_OK and WITNESS_EDGES
        == ((4, 3), (3, 4), (4, 3), (3, 4)),
        "period-four intrinsic clock witness changed")
require(WITNESS_CARRIES == (7, 5, 7, 5)
        and WITNESS_FUTURE_DIGITS == (7, 6, 7, 5),
        "period-four carry/future digits changed")
require(WITNESS_ROOT_STEPS == (9, 2, 9, 4),
        "period-four quotient root steps changed")
require(tuple(label for label, _margin in WITNESS_RAILS) == (0, 12, 0, 12),
        "period-four common-source rail labels changed")
require(all(margin is not None and margin > 0
            for margin in WITNESS_PRESENT_MARGINS),
        "period-four witness left a forced present packet")
require(min(margin for _label, margin in WITNESS_RAILS) / T
        == Fraction(58260191, 137853665040),
        "period-four rail margin changed")
require(min(WITNESS_PRESENT_MARGINS) / T
        == Fraction(15143, 137853665040),
        "period-four present margin changed")
require(WITNESS_TARGET_MARGINS == (
    Fraction(1459, 4080), Fraction(1433, 4080),
    Fraction(12083, 28560), Fraction(1457, 4080),
), "period-four delayed target-a hostile changed")
require(exact_tail_period(WITNESS_TAILS[0], 4) == 4
        and not speed14_resonant(WITNESS_TAILS),
        "period-four witness lost exact nonresonant tail period")


def canonical_counter(counter):
    return tuple(sorted(counter.items()))


def main():
    print("LRC14 affine-odometer periodic-tail scout")
    print("status=VERIFIED-EXACT SCOUT; tail/clock/rail-present level")
    print(f"p={P}; R={R}; S={S}; lift_alphabet={LIFT_ALPHABET}")
    print("tail_skew=y_next={13y}; "
          "N_next=13N+k+floor(13y) mod R; unique_periodic_fibre=True")
    print("symmetric_two_cycle=speed14_resonant_FORCED_BY_13_PLUS_1")
    for length in PERIODS:
        total, compatible, periods, resonant = CENSUS_ROWS[length]
        print(
            f"period_search_n={length}:words={total}:central={total}:"
            f"clock_compatible={compatible}:"
            f"minimal_period_hist={canonical_counter(periods)}:"
            f"speed14_resonant={resonant}"
        )
    for length in (4, 6):
        histogram, margin = FILTER_ROWS[length]
        print(
            f"current_filters_n={length}:"
            f"(rail_hits,present_hits,word_hits)_hist="
            f"{canonical_counter(histogram)}:"
            f"min_target_a_failure_margin={margin}"
        )
    print(f"period4_witness_lifts={WITNESS_LIFTS}")
    print(f"period4_witness_x={WITNESS_STATES}")
    print(f"period4_witness_y={WITNESS_TAILS}")
    print(f"period4_witness_edges={WITNESS_EDGES}; "
          f"carries={WITNESS_CARRIES}; future_digits={WITNESS_FUTURE_DIGITS}; "
          f"root_steps={WITNESS_ROOT_STEPS}")
    print("period4_witness_common_source=1; "
          f"rail_digits={tuple(label for label, _margin in WITNESS_RAILS)}; "
          "rail_present_all_strict=True")
    print("period4_witness_min_margins="
          f"rail:{min(margin for _label, margin in WITNESS_RAILS) / T};"
          f"present:{min(WITNESS_PRESENT_MARGINS) / T}")
    print("delayed_word_hostile=every exact n=4,6 state fails first word "
          "factor D_c2; witness margins="
          f"{WITNESS_TARGET_MARGINS}")
    print("verdict=PASS: variable carries create exact nonresonant periods 4/6 "
          "and evade scalar-dilation clock nilpotence, but the inherited "
          "delayed word kills this finite near-carry alphabet")
    print("SCOPE: no full packet transition, unit transport, endpoint gluing, "
          "row exclusion, or LRC14")


if __name__ == "__main__":
    main()
